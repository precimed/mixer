#!/usr/bin/env python
'''
(c) 2018-2019 Oleksandr Frei, Alexey A. Shadrin
MiXeR software: Univariate and Bivariate Causal Mixture for GWAS
'''

import sys
import os
import argparse
import numpy as np
from numpy import ma
import scipy.optimize
import logging
import json
import scipy.stats
from scipy.interpolate import interp1d

from .libbgmg import LibBgmg
from .utils import UnivariateParams
from .utils import BivariateParams
from .utils import _log_exp_converter
from .utils import _logit_logistic_converter
from .utils import _arctanh_tanh_converter
from .utils import UnivariateParametrization_natural_axis		                     # diffevo, neldermead
from .utils import UnivariateParametrization_constPI_constSIG2BETA                   # inflation
from .utils import UnivariateParametrization_constPI		                         # infinitesimal
from .utils import UnivariateParametrization                                         # uncertainty
from .utils import BivariateParametrization_constUNIVARIATE_constRHOBETA_constPI     # inflation
from .utils import BivariateParametrization_constSIG2BETA_constSIG2ZERO_infPI_maxRG  # infinitesimal
from .utils import BivariateParametrization_constUNIVARIATE_natural_axis             # diffevo, neldermead
from .utils import BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO     # brute1, brent1
from .utils import BivariateParametrization_constUNIVARIATE                          # uncertainty
from .utils import _hessian_robust
from .utils import _max_rg
from .utils import _calculate_univariate_uncertainty
from .utils import _calculate_bivariate_uncertainty


__version__ = '1.0.0'
MASTHEAD = "***********************************************************************\n"
MASTHEAD += "* mixer.py: Univariate and Bivariate Causal Mixture for GWAS\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (c) 2018-2019 Oleksandr Frei, Alexey A. Shadrin\n"
MASTHEAD += "* Norwegian Centre for Mental Disorders Research / University of Oslo\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "***********************************************************************\n"

_cost_calculator_sampling = 0
_cost_calculator_gaussian = 1
_cost_calculator_convolve = 2

_cost_calculator = {
    'gaussian': _cost_calculator_gaussian, 
    'sampling': _cost_calculator_sampling,
    'convolve': _cost_calculator_convolve,
}

def enhance_optimize_result(r, cost_n, cost_df=None, cost_fast=None):
    # optimize_result - an instance of scipy.optimize.OptimizeResult
    # const_n - number of genetic variants effectively contributing to the cost (sum of weights)
    r['cost_n'] = cost_n
    r['cost_df'] = len(r.x) if (cost_df is None) else cost_df
    r['cost'] = r.fun
    r['BIC'] = np.log(r['cost_n']) * r['cost_df'] + 2 * r['cost']
    r['AIC'] =                   2 * r['cost_df'] + 2 * r['cost']
    if cost_fast is not None: r['cost_fast'] = cost_fast

def fix_and_validate_args(args):
    if 'fit_sequence' in args:
        if not args.fit_sequence:
            if not args.trait2_file: # univariate fit
                args.fit_sequence = ['diffevo-fast', 'neldermead']
            else:                    # bivariate fit
                args.fit_sequence = ['diffevo-fast', 'neldermead-fast', 'brute1', 'brent1']
        if args.trait2_file and (args.fit_sequence[0] != 'load'):
            if (not args.trait2_params_file) or (not args.trait1_params_file):
                raise ValueError('--trait1-params-file and --trait2-params-file are required for bivariate analysis (i.e. with --trait2-file argument), unless --fit-sequence starts with "load" step' )

    arg_dict = vars(args)
    chr2use_arg = arg_dict["chr2use"]
    chr2use = []
    for a in chr2use_arg.split(","):
        if "-" in a:
            start, end = [int(x) for x in a.split("-")]
            chr2use += [str(x) for x in range(start, end+1)]
        else:
            chr2use.append(a.strip())
    if np.any([not x.isdigit() for x in chr2use]): raise ValueError('Chromosome labels must be integer')
    arg_dict["chr2use"] = [int(x) for x in chr2use]

def convert_args_to_libbgmg_options(args, num_snp):
    libbgmg_options = {
        'r2min': args.r2min, 'kmax': args.kmax, 
        'max_causals': args.max_causals if (args.max_causals > 1) else (args.max_causals * num_snp),
        'num_components': 1 if (not args.trait2_file) else 3,
        'cache_tag_r2sum': args.cache_tag_r2sum, 'threads': args.threads, 'seed': args.seed,
        'cubature_rel_error': args.cubature_rel_error, 'cubature_max_evals':args.cubature_max_evals
        # 'z1max': args.z1max, 'z2max': args.z2max, 
    }
    return [(k, v) for k, v in libbgmg_options.items() if v is not None ]

# https://stackoverflow.com/questions/27433316/how-to-get-argparse-to-read-arguments-from-a-file-with-an-option-rather-than-pre
class LoadFromFile (argparse.Action):
    def __call__ (self, parser, namespace, values, option_string=None):
        with values as f:
            contents = f.read()

        data = parser.parse_args(contents.split(), namespace=namespace)
        for k, v in vars(data).items():
            if v and k != option_string.lstrip('-'):
                setattr(namespace, k, v)

def parser_fit_add_arguments(args, func, parser):
    parser.add_argument("--bim-file", type=str, default=None, help="Plink bim file. "
        "Defines the reference set of SNPs used for the analysis. "
        "Marker names must not have duplicated entries. "
        "May contain simbol '@', which will be replaced with the actual chromosome label. ")
    parser.add_argument("--frq-file", type=str, default=None, help="Plink frq file (alleles frequencies). "
        "May contain simbol '@', similarly to --bim-file argument. ")
    parser.add_argument("--plink-ld-bin", type=str, default=None, help="File with linkage disequilibrium information, "
        "converted from plink format as described in the README.md file. "
        "May contain simbol '@', similarly to --bim-file argument. ")
    parser.add_argument("--plink-ld-bin0", type=str, default=None, help="File with linkage disequilibrium information in an old format (deprecated)")
    parser.add_argument("--chr2use", type=str, default="1-22", help="Chromosome ids to use "
         "(e.g. 1,2,3 or 1-4,12,16-20). Chromosome must be labeled by integer, i.e. X and Y are not acceptable. ")
    parser.add_argument("--trait1-file", type=str, default=None, help="GWAS summary statistics for the first trait. ")
    parser.add_argument("--trait2-file", type=str, default="", help="GWAS summary statistics for the first trait. "
        "Specifying this argument triggers cross-trait analysis.")
    parser.add_argument('--fit-sequence', type=str, default=[], nargs='+',
        choices=['load', 'inflation', 'infinitesimal', 'diffevo', 'diffevo-fast', 'neldermead', 'neldermead-fast', 'brute1', 'brute1-fast', 'brent1', 'brent1-fast'],
        help="Specify fit sequence: "
             "'load' reads previosly fitted parameters from a file (--load-params-file); "
             "'diffevo' performs a single iteration of differential evolution, which is the right way to setup an initial approximation; "
             "'neldermead' applies Nelder-Mead downhill simplex search; "             
             "'brute1' applies only to bivariate optimization; it performs brute-force for one-dimentional optimization, searching optimal pi12 value constrained on genetic correlation (rg) and intercept (rho0); "
             "'brent1' is similar to brute1, but it uses brent method (inverse parabolic interpolation); "
             "'inflation' fits sig2zero (univariate) and rho_zero (bivariate), using fast cost function; this is quite special optimization step, typically useful for adjusting inflation parameters to another reference; "
             "'infinitesimal' fits a model with pi1=1 (univariate) or pi12=1 (bivariate) constrains, using fast cost function; this is quite special optimization step, typically used internally for AIC/BIC computation; " 
             "Note that bivariate fit is always constrained on univariate parameters, except for 'inflation' fit which adjust rho_zero and sig2_zero. "
             "The '...-fast' optimizations use fast cost function. "
             "Note that univariate optimization uses 'convolve' cost calculator, bivariate optimization uses 'sampling' cost calculator. "
             "Typical univariate sequence: 'diffevo-fast neldermead'"
             "Typical bivariate sequence: 'diffevo neldermead brute1 brent1'")
    parser.add_argument('--preliminary', default=False, action="store_true",
        help="perform an additional run using fast model with 'diffevo-fast nelderead-fast' to generate preliminary data. "
        "After preliminary run fit sequence is applied from scratch using full model.")

    parser.add_argument('--extract', type=str, default="", help="File with variants to include in the fit procedure")
    parser.add_argument('--exclude', type=str, default="", help="File with variants to exclude from the fit procedure")

    parser.add_argument('--randprune-n', type=int, default=64, help="Number of random pruning iterations")
    parser.add_argument('--randprune-r2', type=float, default=0.1, help="Threshold for random pruning")
    parser.add_argument('--kmax', type=int, default=20000, help="Number of sampling iterations")
    parser.add_argument('--seed', type=int, default=123, help="Random seed")

    parser.add_argument('--cache-tag-r2sum', default=False, action="store_true", help="enable tag-r2sum caching")
    parser.add_argument('--max-causals', type=float, default=0.03, help="upper limit for the total number of causal variants in the reference; a number between 0 and 1 represents a fraction of the total number SNPs in the reference")
    parser.add_argument('--r2min', type=float, default=0.05, help="r2 values below this threshold will contribute via infinitesimal model")
    parser.add_argument('--ci-alpha', type=float, default=None, help="significance level for the confidence interval estimation")
    parser.add_argument('--ci-samples', type=int, default=10000, help="number of samples in uncertainty estimation")
    parser.add_argument('--threads', type=int, default=None, help="specify how many threads to use (concurrency). None will default to the total number of CPU cores. ")
    parser.add_argument('--tol-x', type=float, default=1e-2, help="tolerance for the stop criteria in fminsearch optimization. ")
    parser.add_argument('--tol-func', type=float, default=1e-2, help="tolerance for the stop criteria in fminsearch optimization. ")
    parser.add_argument('--cubature-rel-error', type=float, default=1e-5, help="relative error for cubature stop criteria (applies to 'convolve' cost calculator). ")
    parser.add_argument('--cubature-max-evals', type=float, default=1000, help="max evaluations for cubature stop criteria (applies to 'convolve' cost calculator). "
        "Bivariate cubature require in the order of 10^4 evaluations and thus is much slower than sampling, therefore it is not exposed via mixer.py command-line interface. ")

    # NB! zmax won't apply in univariate fit because it uses convolution cost function. This must be fixed.
    # parser.add_argument('--z1max', type=float, default=None, help="right-censoring threshold for the first trait. ")
    # parser.add_argument('--z2max', type=float, default=None, help="right-censoring threshold for the second trait. ")

    parser.add_argument('--load-params-file', type=str, default=None, help="initial params for the optimization. ")
    parser.add_argument('--trait1-params-file', type=str, default=None, help="univariate params for the first trait (for the cross-trait analysis only). ")
    parser.add_argument('--trait2-params-file', type=str, default=None, help="univariate params for the second trait (for the cross-trait analysis only). ")

    parser.add_argument('--downsample-factor', default=50, type=int, help="Applies to --power-curve. "
        "'--downsample-factor N' imply that only 1 out of N available z-score values will be used in calculations.")
    parser.add_argument('--power-curve', default=False, action="store_true", help="generate power curves")
    parser.add_argument('--qq-plots', default=False, action="store_true", help="generate qq plot curves")    

    parser.set_defaults(func=func)

def parser_ld_add_arguments(args, func, parser):
    parser.add_argument("--plink-ld", type=str, default=None, help="Path to plink .ld.gz file to convert into BGMG binary format. ")
    parser.add_argument("--bim-file", type=str, default=None, help="Plink bim file. "
        "Defines the reference set of SNPs used for the analysis. "
        "Marker names must not have duplicated entries. "
        "May contain simbol '@', which will be replaced with the actual chromosome label. ")
    parser.add_argument("--chr2use", type=str, default="1-22", help="Chromosome ids to use "
         "(e.g. 1,2,3 or 1-4,12,16-20). Chromosome must be labeled by integer, i.e. X and Y are not acceptable. ")
    parser.set_defaults(func=func)

def parse_args(args):
    parser = argparse.ArgumentParser(description="MiXeR: Univariate and Bivariate Causal Mixture for GWAS.")

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--argsfile', type=open, action=LoadFromFile, default=None, help="file with additional command-line arguments")
    parent_parser.add_argument("--out", type=str, default="mixer", help="prefix for the output files")
    parent_parser.add_argument("--lib", type=str, default="libbgmg.so", help="path to libbgmg.so plugin")
    parent_parser.add_argument("--log", type=str, default=None, help="file to output log, defaults to <out>.log")
    
    subparsers = parser.add_subparsers()

    parser_fit_add_arguments(args=args, func=execute_fit_parser, parser=subparsers.add_parser("fit", parents=[parent_parser], help='fit MiXeR model'))
    parser_ld_add_arguments(args=args, func=execute_ld_parser, parser=subparsers.add_parser("ld", parents=[parent_parser], help='prepare files with linkage disequilibrium information'))

    return parser.parse_args(args)

def log_header(args, subparser_name, lib):
    defaults = vars(parse_args([subparser_name]))
    opts = vars(args)
    non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
    header = MASTHEAD
    header += "Call: \n"
    header += './mixer.py {} \\\n'.format(subparser_name)
    options = ['\t--'+x.replace('_','-')+' '+str(opts[x]).replace('\t', '\\t')+' \\' for x in non_defaults]
    header += '\n'.join(options).replace('True','').replace('False','')
    header = header[0:-1]+'\n'
    lib.log_message(header)

def load_univariate_params_file(fname):
    data = json.loads(open(fname).read())
    if 'analysis' not in data or data['analysis'] != 'univariate': raise ValueError('Invalid univariate results file')
    p = data['params']
    return UnivariateParams(pi=p['pi'], sig2_beta=p['sig2_beta'], sig2_zero=p['sig2_zero'])

def load_bivariate_params_file(fname):
    data = json.loads(open(fname).read())
    if 'analysis' not in data or data['analysis'] != 'bivariate': raise ValueError('Invalid bivariate results file')
    p = data['params']
    return (BivariateParams(pi=p['pi'], sig2_beta=p['sig2_beta'], sig2_zero=p['sig2_zero'], rho_beta=p['rho_beta'], rho_zero=p['rho_zero']),
            UnivariateParams(pi=p['pi'][0]+p['pi'][2], sig2_beta=p['sig2_beta'][0], sig2_zero=p['sig2_zero'][0]),
            UnivariateParams(pi=p['pi'][1]+p['pi'][2], sig2_beta=p['sig2_beta'][1], sig2_zero=p['sig2_zero'][1]))

def apply_univariate_fit_sequence(args, libbgmg, fit_sequence, init_params=None, trait=1):
    params=init_params; optimize_result_sequence=[]
    for fit_type in fit_sequence:
        libbgmg.log_message("fit_type=={}...".format(fit_type))

        if fit_type == 'load':
            params = load_univariate_params_file(args.load_params_file)
            optimize_result = {}
            libbgmg.log_message("fit_type==load: Done, {}".format(params))

        elif (fit_type == 'diffevo') or fit_type == ('diffevo-fast'):
            libbgmg.set_option('cost_calculator', _cost_calculator_convolve if (fit_type == 'diffevo') else _cost_calculator_gaussian)
            parametrization = UnivariateParametrization_natural_axis(lib=libbgmg, trait=trait)
            bounds_left = parametrization.params_to_vec(UnivariateParams(pi=5e-5, sig2_beta=5e-6, sig2_zero=0.9))
            bounds_right = parametrization.params_to_vec(UnivariateParams(pi=5e-1, sig2_beta=5e-2, sig2_zero=2.5))
            bounds4opt = [(l, r) for l, r in zip(bounds_left, bounds_right)]
            optimize_result = scipy.optimize.differential_evolution(lambda x: parametrization.calc_cost(x), bounds4opt,
                tol=0.01, mutation=(0.5, 1), recombination=0.7, atol=0, updating='immediate', polish=False, workers=1)  #, **global_opt_options)
            params = parametrization.vec_to_params(optimize_result.x)

        elif (fit_type == 'neldermead') or (fit_type == 'neldermead-fast'):
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "neldermead" fit'))
            libbgmg.set_option('cost_calculator', _cost_calculator_convolve if (fit_type == 'neldermead') else _cost_calculator_gaussian)
            parametrization = UnivariateParametrization_natural_axis(lib=libbgmg, trait=trait)
            optimize_result = scipy.optimize.minimize(lambda x: parametrization.calc_cost(x), parametrization.params_to_vec(params),
                method='Nelder-Mead', options={'maxiter':1200, 'fatol':1e-7, 'xatol':1e-4, 'adaptive':True})
            params = parametrization.vec_to_params(optimize_result.x)

        elif fit_type == 'inflation':
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "inflation" fit'))
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            optimizer = lambda func, x0: scipy.optimize.minimize(func, x0, method='Nelder-Mead')
            params, optimize_result = UnivariateParametrization_constPI_constSIG2BETA(
                init_sig2_zero=params._sig2_zero, const_params=params, lib=libbgmg, trait=trait).fit(optimizer)

        elif fit_type == 'infinitesimal':
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "infinitesimal" fit'))
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            optimizer = lambda func, x0: scipy.optimize.minimize(func, x0, method='Nelder-Mead')
            params, optimize_result = UnivariateParametrization_constPI(
                const_pi=1.0, init_sig2_zero=params._sig2_zero,
                init_sig2_beta=params._sig2_beta * params._pi,
                lib=libbgmg, trait=trait).fit(optimizer)
            # optimize_result['cost_df'] = 2    # happens automatically?

        else:
            libbgmg.log_message("Unable to apply {} in univariate fit".format(fit_type))

        libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
        if optimize_result:
            enhance_optimize_result(optimize_result, cost_n=np.sum(libbgmg.weights), cost_fast=params.cost(libbgmg, trait))
        optimize_result['params']=params.as_dict()   # params after optimization
        optimize_result_sequence.append((fit_type, optimize_result))
        libbgmg.log_message("fit_type=={} done ({}, {})".format(fit_type, params, optimize_result))

    if params == None: raise(RuntimeError('Empty --fit-sequence'))
    return params, optimize_result_sequence

def apply_bivariate_fit_sequence(args, libbgmg):
    if args.trait1_params_file or args.trait2_params_file:
        libbgmg.log_message("Loading univariate constrains for bivariate analysis...")
        params1 = load_univariate_params_file(args.trait1_params_file)
        params2 = load_univariate_params_file(args.trait2_params_file)
        libbgmg.log_message("trait1: {}".format(params1))
        libbgmg.log_message("trait2: {}".format(params2))

    params = None; optimize_result_sequence = []
    for fit_type in args.fit_sequence:
        libbgmg.log_message("fit_type=={}...".format(fit_type))

        if fit_type == 'load':
            params, params1, params2 = load_bivariate_params_file(args.load_params_file)
            if args.trait1_params_file or args.trait2_params_file:
                libbgmg.log_message("overwrite params previously loaded from --trait1-params-file, --trait2-params-file")
            libbgmg.log_message("fit_type==load... trait1 params: {}".format(params1))            
            libbgmg.log_message("fit_type==load... trait2 params: {}".format(params2))    
            libbgmg.log_message("fit_type=={} done ({})".format(fit_type, params))
            continue

        elif fit_type == 'inflation':
            if params == None: raise(RuntimeError('params == None, unable to proceed with "{}" fit'.format(fit_type)))
            params1, ors1 = apply_univariate_fit_sequence(args, libbgmg, ['inflation'], init_params=params1, trait=1)
            optimize_result_sequence.extend(ors1)
            params2, ors2 = apply_univariate_fit_sequence(args, libbgmg, ['inflation'], init_params=params2, trait=2)
            optimize_result_sequence.extend(ors2)
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            optimizer = lambda func, x0: scipy.optimize.minimize(func, x0, method='Nelder-Mead')
            params, optimize_result = BivariateParametrization_constUNIVARIATE_constRHOBETA_constPI(
                const_params1=params1, const_params2=params2,
                const_pi12=params._pi[2], const_rho_beta=params._rho_beta, init_rho_zero=params._rho_zero,
                lib=libbgmg).fit(optimizer)

        elif fit_type == 'infinitesimal':
            if params == None: raise(RuntimeError('params == None, unable to proceed with "{}" fit'.format(fit_type)))
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            optimizer = lambda func, x0: scipy.optimize.minimize(func, x0, method='Nelder-Mead')
            params, optimize_result = BivariateParametrization_constSIG2BETA_constSIG2ZERO_infPI_maxRG(
                const_sig2_beta=[params1._sig2_beta, params2._sig2_beta],
                const_sig2_zero=[params1._sig2_zero, params2._sig2_zero],
                max_rg = 1,
                init_rho_beta = params._rho_beta,
                init_rho_zero = params._rho_zero,
                lib=libbgmg).fit(optimizer)

        elif (fit_type == 'diffevo') or fit_type == ('diffevo-fast'):
            libbgmg.set_option('cost_calculator', _cost_calculator_sampling if (fit_type == 'diffevo') else _cost_calculator_gaussian)
            parametrization = BivariateParametrization_constUNIVARIATE_natural_axis(lib=libbgmg, const_params1=params1, const_params2=params2)
            max_pi12 = min(params1._pi, params2._pi)
            bounds_left = parametrization.params_to_vec(BivariateParams(params1=params1, params2=params2, pi12=0.05 * max_pi12, rho_beta=-0.95, rho_zero=-0.95))
            bounds_right = parametrization.params_to_vec(BivariateParams(params1=params1, params2=params2, pi12=0.95 * max_pi12, rho_beta=0.95, rho_zero=0.95))
            bounds4opt = [(l, r) for l, r in zip(bounds_left, bounds_right)]
            optimize_result = scipy.optimize.differential_evolution(lambda x: parametrization.calc_cost(x), bounds4opt,
                tol=0.01, mutation=(0.5, 1), recombination=0.7, atol=0, updating='immediate', polish=False, workers=1)  #, **global_opt_options)
            params = parametrization.vec_to_params(optimize_result.x)

        elif (fit_type == 'neldermead') or (fit_type == 'neldermead-fast'):
            if params == None: raise(RuntimeError('params == None, unable to proceed with "{}" fit'.format(fit_type)))
            libbgmg.set_option('cost_calculator', _cost_calculator_sampling if (fit_type == 'neldermead') else _cost_calculator_gaussian)
            parametrization = BivariateParametrization_constUNIVARIATE_natural_axis(lib=libbgmg, const_params1=params1, const_params2=params2)
            optimize_result = scipy.optimize.minimize(lambda x: parametrization.calc_cost(x), parametrization.params_to_vec(params),
                method='Nelder-Mead', options={'maxiter':1200, 'fatol':1e-7, 'xatol':1e-4, 'adaptive':True})
            params = parametrization.vec_to_params(optimize_result.x)

        elif (fit_type == 'brute1') or (fit_type == 'brute1-fast'):
            # brute1 optimization intentionally forgets previously fitted value of params._pi[2], to avoid being stuck in a local minimum.
            if params == None: raise(RuntimeError('params == None, unable to proceed with "{}" fit'.format(fit_type)))
            libbgmg.set_option('cost_calculator', _cost_calculator_sampling if (fit_type == 'brute1') else _cost_calculator_gaussian)
            parametrization = BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO(
                lib=libbgmg, const_params1=params1, const_params2=params2, const_rg=params._rg(), const_rho_zero=params._rho_zero)
            min_pi12 = parametrization._min_pi12; max_pi12 = parametrization._max_pi12
            params._pi[2] = 0.99 * min_pi12 + 0.01 * max_pi12; bounds_left = parametrization.params_to_vec(params)
            params._pi[2] = 0.01 * min_pi12 + 0.99 * max_pi12; bounds_right = parametrization.params_to_vec(params)
            bounds4opt = [(l, r) for l, r in zip(bounds_left, bounds_right)]
            x0, fval, grid, Jout = scipy.optimize.brute(lambda x: parametrization.calc_cost(x), bounds4opt, Ns=20, full_output=True)
            optimize_result = scipy.optimize.OptimizeResult(fun=fval, nit=1, nfev=len(grid),
                            success=True, message="Optimization terminated successfully.", x=x0)
            optimize_result['grid'] = grid
            optimize_result['Jout'] = Jout
            params = parametrization.vec_to_params(optimize_result.x)

        elif (fit_type == 'brent1') or (fit_type == 'brent1-fast'):
            if params == None: raise(RuntimeError('params == None, unable to proceed with "{}" fit'.format(fit_type)))
            libbgmg.set_option('cost_calculator', _cost_calculator_sampling if (fit_type == 'brent1') else _cost_calculator_gaussian)
            parametrization = BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO(
                lib=libbgmg, const_params1=params1, const_params2=params2, const_rg=params._rg(), const_rho_zero=params._rho_zero)
            min_pi12 = parametrization._min_pi12; max_pi12 = parametrization._max_pi12
            bracket_middle = parametrization.params_to_vec(params)
            params._pi[2] = min_pi12; bounds_left = parametrization.params_to_vec(params)
            params._pi[2] = max_pi12; bounds_right = parametrization.params_to_vec(params)
            bracket4opt = (bounds_left[0], bracket_middle[0], bounds_right[0])
            bracketfunc = (parametrization.calc_cost(bounds_left[0]), parametrization.calc_cost(bracket_middle[0]),parametrization.calc_cost(bounds_right[0]))
            libbgmg.log_message("bracket4opt = {}, {}".format(bracket4opt,bracketfunc))
            try:
                xmin, fval, iter, funcalls = scipy.optimize.brent(lambda x: parametrization.calc_cost(x), brack=bracket4opt, full_output=True)
                success=True; message = "Optimization terminated successfully."
            except ValueError as ve:
                xmin = bracket_middle[0]; fval = parametrization.calc_cost(bracket_middle[0])
                iter = 1; funcalls = 4 # three byscipy.optimize.brent, plus one in the line above
                success=False; message=str(ve)
            optimize_result = scipy.optimize.OptimizeResult(fun=fval, nit=iter, nfev=funcalls, x=[xmin], success=success, message=message)
            params = parametrization.vec_to_params(optimize_result.x)

        else:
            libbgmg.log_message("Unable to apply {} in bivariate fit".format(fit_type))

        libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
        # cost_df=9 --- nine free parameters (incl. univariate)
        enhance_optimize_result(optimize_result, cost_df=9, cost_n=np.sum(libbgmg.weights), cost_fast=params.cost(libbgmg))
        optimize_result['params']=params.as_dict()   # params after optimization
        optimize_result_sequence.append((fit_type, optimize_result))
        libbgmg.log_message("fit_type=={} done ({}, {})".format(fit_type, params, optimize_result))

    if params == None: raise(RuntimeError('Empty --fit-sequence'))
    return (params, params1, params2, optimize_result_sequence)

def calc_power_curve(libbgmg, params, trait_index, downsample):
    power_nvec = np.power(10, np.arange(3, 8, 0.1))

    original_weights = libbgmg.weights
    if not np.all(np.isfinite(original_weights)): raise(RuntimeError('undefined weights not supported'))
    if not np.all(np.isfinite(libbgmg.get_zvec(trait_index))): raise(RuntimeError('undefined weights not supported'))
    model_weights = libbgmg.weights

    mask = np.zeros((len(model_weights), ), dtype=bool)
    mask[range(0, len(model_weights), downsample)] = 1
    model_weights[~mask] = 0
    model_weights = model_weights/np.sum(model_weights)

    libbgmg.weights = model_weights     # temporary set downsampled weights
    power_svec = libbgmg.calc_univariate_power(trait_index, params._pi, params._sig2_zero, params._sig2_beta, 5.45, power_nvec)
    libbgmg.weights = original_weights  # restore original weights

    return power_nvec, power_svec

def calc_qq_plot(libbgmg, params, trait_index, downsample, mask=None):
    # mask can subset SNPs that are going into QQ curve, for example LDxMAF bin.
    if mask is None:
        mask = np.ones((libbgmg.num_tag, ), dtype=bool)

    # Empirical (data) QQ plot
    # Step 0. calculate weights for all data poitns
    # Step 1. convert zvec to -log10(pvec)
    # Step 2. sort pval from large (1.0) to small (0.0), and take empirical cdf of this distribution
    # Step 3. interpolate (data_x, data_y) curve, as we don't need 10M data points on QQ plots
    data_weights = libbgmg.weights                                      # step 0 
    data_weights = data_weights[mask] / np.sum(data_weights[mask])      # step 0
    zvec = libbgmg.get_zvec(trait_index)
    data_y = -np.log10(2*scipy.stats.norm.cdf(-np.abs(zvec[mask])))     # step 1
    si = np.argsort(data_y); data_y = data_y[si]                        # step 2
    data_x=-np.log10(np.flip(np.cumsum(np.flip(data_weights[si]))))     # step 2
    hv_z = np.linspace(0, np.min([np.max(np.abs(zvec)), 38.0]), 1000)   # step 3
    hv_logp = -np.log10(2*scipy.stats.norm.cdf(-hv_z))                  # step 3
    data_idx = np.not_equal(data_y, np.concatenate((data_y[1:], [np.Inf])))
    data_logpvec = interp1d(data_y[data_idx], data_x[data_idx],
                            bounds_error=False, fill_value=np.nan)(hv_logp) 
    #plt.plot(data_logpvec, hv_logp, '.')
    #plt.plot(data_x, data_y, '-')

    # Estimated (model) QQ plots
    model_weights = libbgmg.weights
    mask_indices = np.nonzero(mask)[0]
    model_mask = np.zeros((len(model_weights), ), dtype=bool)
    model_mask[mask_indices[range(0, len(mask_indices), downsample)]] = 1
    model_weights[~model_mask] = 0
    model_weights = model_weights/np.sum(model_weights)

    zgrid = np.arange(0, 38.0, 0.05, np.float32)

    original_weights = libbgmg.weights
    libbgmg.weights = model_weights
    pdf = libbgmg.calc_univariate_pdf(trait_index, params._pi, params._sig2_zero, params._sig2_beta, zgrid)
    libbgmg.weights = original_weights

    pdf = pdf / np.sum(model_weights)
    zgrid = np.concatenate((np.flip(-zgrid[1:]), zgrid))  # extend [0, 38] to [-38, 38]
    pdf = np.concatenate((np.flip(pdf[1:]), pdf))
    model_cdf = np.cumsum(pdf)  * (zgrid[1] - zgrid[0])
    model_cdf = 0.5 * (np.concatenate(([0.0], model_cdf[:-1])) + np.concatenate((model_cdf[:-1], [1.0])))
    model_logpvec = -np.log10(2*interp1d(-zgrid[zgrid<=0], model_cdf[zgrid<=0])(hv_z))
    #plt.plot(model_logpvec, hv_logp, '-')
    #plt.plot(data_logpvec, hv_logp, '-')
    #plt.plot(hv_logp, hv_logp, '--k')

    return hv_logp, data_logpvec, model_logpvec


# helper function to debug non-json searizable types...
def print_types(results, libbgmg):
    if isinstance(results, dict):
        for k, v in results.items():
            libbgmg.log_message('{}: {}'.format(k, type(v)))
            print_types(v, libbgmg)

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if callable(obj):
            return str(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.float32):
            return np.float64(obj)
        return json.JSONEncoder.default(self, obj)

def execute_ld_parser(args):
    libbgmg = LibBgmg(args.lib)
    fix_and_validate_args(args)

    libbgmg.init(args.bim_file, "", args.chr2use, "", "", "", "")
    libbgmg.convert_plink_ld(args.plink_ld, args.out + '.bin')
    libbgmg.log_message('Done')

def execute_fit_parser(args):
    libbgmg = LibBgmg(args.lib)
    fix_and_validate_args(args)

    # for univariate optimization, if fit steps involve "diffevo" or "neldermead", set "use_complete_tag_indices" (to enable convolution cost calculator)
    if not args.trait2_file:
        for fit_type in args.fit_sequence:
            if fit_type in ['diffevo', 'neldermead']:
                libbgmg.set_option('use_complete_tag_indices', 1)
                break

    libbgmg.init(args.bim_file, args.frq_file, args.chr2use, args.trait1_file, args.trait2_file,
        args.exclude, args.extract)

    for opt, val in convert_args_to_libbgmg_options(args, libbgmg.num_snp):
        libbgmg.set_option(opt, val)

    if args.plink_ld_bin0 is not None:
        libbgmg.set_option('ld_format_version', 0)
        args.plink_ld_bin = args.plink_ld_bin0
        args.plink_ld_bin0 = None

    for chr_label in args.chr2use: 
        libbgmg.set_ld_r2_coo_from_file(chr_label, args.plink_ld_bin.replace('@', str(chr_label)))
        libbgmg.set_ld_r2_csr(chr_label)

    libbgmg.set_weights_randprune(args.randprune_n, args.randprune_r2)
    libbgmg.set_option('diag', 0)

    totalhet = float(2.0 * np.dot(libbgmg.mafvec, 1.0 - libbgmg.mafvec))

    fit_sequence = args.fit_sequence.copy()
    for repeat in range(2):
        if (repeat == 0) and (not args.preliminary): continue # skip generating preliminary data
        args.fit_sequence = ['diffevo-fast', 'neldermead-fast'] if (repeat==0) else fit_sequence
        libbgmg.log_message('Apply fit sequence: {}{}...'.format(args.fit_sequence, ' (--preliminary)' if (repeat == 0) else '' ))

        results = {}
        results['options'] = vars(args).copy()
        results['options']['totalhet'] = totalhet
        results['options']['num_snp'] = float(libbgmg.num_snp)
        results['options']['num_tag'] = float(libbgmg.num_tag)
        results['options']['sum_weights'] = float(np.sum(libbgmg.weights))
        results['options']['trait1_nval'] = float(np.nanmedian(libbgmg.get_nvec(trait=1)))

        if not args.trait2_file:
            results['analysis'] = 'univariate'

            params, optimize_result = apply_univariate_fit_sequence(args, libbgmg, args.fit_sequence)
            results['params'] = params.as_dict()
            results['optimize'] = optimize_result

            libbgmg.log_message('Calculate AIC/BIC w.r.t. infinitesimal model (fast cost function)...')
            params_inft, optimize_result_inft = apply_univariate_fit_sequence(args, libbgmg, ['infinitesimal'], init_params=params)
            results['inft_params'] = params_inft.as_dict()
            results['inft_optimize'] = optimize_result_inft

            if args.ci_alpha and np.isfinite(args.ci_alpha):
                libbgmg.log_message("Uncertainty estimation...")
                libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
                results['ci'], _ = _calculate_univariate_uncertainty(UnivariateParametrization(params, libbgmg, trait=1), args.ci_alpha, totalhet, libbgmg.num_snp, args.ci_samples)
                for k, v in results['ci'].items():
                    libbgmg.log_message('{}: point_estimate={:.3g}, mean={:.3g}, median={:.3g}, std={:.3g}, ci=[{:.3g}, {:.3g}]'.format(k, v['point_estimate'], v['mean'], v['median'], v['std'], v['lower'], v['upper']))
                libbgmg.log_message("Uncertainty estimation done.")

            if args.power_curve:
                trait_index = 1
                power_nvec, power_svec = calc_power_curve(libbgmg, params, trait_index, args.downsample_factor)
                results['power'] = {'nvec': power_nvec, 'svec': power_svec}
            if args.qq_plots:
                trait_index = 1
                mask = np.ones((libbgmg.num_tag, ), dtype=bool)
                hv_logp, data_logpvec, model_logpvec = calc_qq_plot(libbgmg, params, trait_index, args.downsample_factor, mask)
                results['qqplot'] = {'hv_logp': hv_logp, 'data_logpvec': data_logpvec, 'model_logpvec': model_logpvec,
                                     'n_snps': int(np.sum(mask)), 'sum_data_weights': float(np.sum(libbgmg.weights[mask])),
                                     'title' : 'maf \\in [{:.3g},{:.3g}); L \\in [{:.3g},{:.3g})'.format(-np.inf,np.inf,-np.inf,np.inf)}
                mafvec = libbgmg.mafvec[libbgmg.defvec]
                tldvec = libbgmg.ld_tag_r2_sum
                maf_bins = np.concatenate(([-np.inf], np.quantile(mafvec, [1/3, 2/3]), [np.inf]))
                tld_bins = np.concatenate(([-np.inf], np.quantile(tldvec, [1/3, 2/3]), [np.inf]))
                results['qqplot_bins'] = []
                for i in range(0, 3):
                    for j in range(0, 3):
                        mask = ((mafvec>=maf_bins[i]) & (mafvec<maf_bins[i+1]) & (tldvec >= tld_bins[j]) &  (tldvec < tld_bins[j+1]))
                        hv_logp, data_logpvec, model_logpvec = calc_qq_plot(libbgmg, params, trait_index, args.downsample_factor, mask)
                        results['qqplot_bins'].append(
                            {'hv_logp': hv_logp, 'data_logpvec': data_logpvec, 'model_logpvec': model_logpvec,
                             'n_snps': int(np.sum(mask)), 'sum_data_weights': float(np.sum(libbgmg.weights[mask])),
                             'title' : 'maf \\in [{:.3g},{:.3g}); L \\in [{:.3g},{:.3g})'.format(maf_bins[i], maf_bins[i+1], tld_bins[j], tld_bins[j+1])}
                        )
        else:
            results['analysis'] = 'bivariate'
            results['options']['trait2_nval'] = float(np.nanmedian(libbgmg.get_nvec(trait=2)))

            params, params1, params2, optimize_result = apply_bivariate_fit_sequence(args, libbgmg)
            results['params'] = params.as_dict()
            results['optimize'] = optimize_result
            if args.ci_alpha and np.isfinite(args.ci_alpha):
                libbgmg.log_message("Uncertainty estimation...")
                libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
                _, ci_sample1 = _calculate_univariate_uncertainty(UnivariateParametrization(params1, libbgmg, trait=1), args.ci_alpha, totalhet, libbgmg.num_snp, args.ci_samples)
                _, ci_sample2 = _calculate_univariate_uncertainty(UnivariateParametrization(params2, libbgmg, trait=2), args.ci_alpha, totalhet, libbgmg.num_snp, args.ci_samples)
                results['ci'], _ = _calculate_bivariate_uncertainty(BivariateParametrization_constUNIVARIATE(
                    const_params1=params1, const_params2=params2, init_pi12=params._pi[2],
                    init_rho_beta=params._rho_beta, init_rho_zero=params._rho_zero,
                    lib=libbgmg), [ci_sample1, ci_sample2], args.ci_alpha, totalhet, libbgmg.num_snp, args.ci_samples)
                for k, v in results['ci'].items():
                    libbgmg.log_message('{}: point_estimate={:.3g}, mean={:.3g}, median={:.3g}, std={:.3g}, ci=[{:.3g}, {:.3g}]'.format(k, v['point_estimate'], v['mean'], v['median'], v['std'], v['lower'], v['upper']))
                libbgmg.log_message("Uncertainty estimation done.")

        with open(args.out + ('.preliminary' if (repeat == 0) else '')+ '.json', 'w') as outfile:  
            json.dump(results, outfile, cls=NumpyEncoder)

    libbgmg.set_option('diag', 0)
