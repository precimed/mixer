#!/usr/bin/env python
'''
(c) 2016-2020 Oleksandr Frei, Alexey A. Shadrin, Dominic Holland
MiXeR software: Univariate and Bivariate Causal Mixture for GWAS
'''

import sys
import os
import argparse
import numpy as np
import pandas as pd
from numpy import ma
import scipy.optimize
import logging
import json
import scipy.stats
import random
from scipy.interpolate import interp1d
import time

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
from .utils import _calculate_univariate_uncertainty_funcs
from .utils import _calculate_bivariate_uncertainty
from .utils import _calculate_bivariate_uncertainty_funcs

from .utils import calc_qq_data
from .utils import calc_qq_model
from .utils import calc_qq_plot
from .utils import calc_power_curve
from .utils import NumpyEncoder

__version__ = '1.2.0'
MASTHEAD = "***********************************************************************\n"
MASTHEAD += "* mixer.py: Univariate and Bivariate Causal Mixture for GWAS\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (c) 2016-2020 Oleksandr Frei, Alexey A. Shadrin, Dominic Holland\n"
MASTHEAD += "* Norwegian Centre for Mental Disorders Research / University of Oslo\n"
MASTHEAD += "* Center for Multimodal Imaging and Genetics / UCSD\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "***********************************************************************\n"

_cost_calculator_sampling = 0
_cost_calculator_gaussian = 1
_cost_calculator_convolve = 2
_cost_calculator_smplfast = 3

_cost_calculator = {
    'gaussian': _cost_calculator_gaussian, 
    'sampling': _cost_calculator_sampling,
    'convolve': _cost_calculator_convolve,
    'smplfast': _cost_calculator_smplfast,
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

def check_input_file(args, argname, chri=None):
    arg_dict = vars(args)
    if argname in arg_dict:
        if not arg_dict[argname]: raise ValueError("Missing required argument --{}".format(argname))
        argval = arg_dict[argname]
        if chri: argval=argval.replace('@', str(chri))
        if not os.path.isfile(argval): raise ValueError("Input file --{a} does not exist: {f}".format(a=argname, f=argval))

def fix_and_validate_args(args):
    check_input_file(args, 'trait1_file')
    check_input_file(args, 'trait2_file')
    check_input_file(args, 'trait1_params_file')
    check_input_file(args, 'trait2_params_file')
    check_input_file(args, 'load-params-file')

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

    for chri in arg_dict["chr2use"]:
        check_input_file(args, 'bim-file', chri)
        check_input_file(args, 'ld-file', chri)

def convert_args_to_libbgmg_options(args, num_snp):
    libbgmg_options = {
        'r2min': args.r2min if ('r2min' in args) else None,
        'kmax': args.kmax[0] if ('kmax' in args) else None, 
        'threads': args.threads[0] if ('threads' in args) else None,
        'seed': args.seed if ('seed' in args) else None,
        'cubature_rel_error': args.cubature_rel_error if ('cubature_rel_error' in args) else None,
        'cubature_max_evals': args.cubature_max_evals if ('cubature_max_evals' in args) else None,
        'z1max': args.z1max if ('z1max' in args) else None,
        'z2max': args.z2max if ('z2max' in args) else None, 
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

def parser_add_common_arguments(parser, num_traits):
    parser.add_argument("--bim-file", type=str, default=None, help="Plink bim file. "
        "Defines the reference set of SNPs used for the analysis. "
        "Marker names must not have duplicated entries. "
        "May contain simbol '@', which will be replaced with the actual chromosome label. ")
    parser.add_argument("--ld-file", type=str, default=None, help="File with linkage disequilibrium information, "
        "generated via 'mixer.py ld' command. "
        "May contain simbol '@', similarly to --bim-file argument. ")
    parser.add_argument("--chr2use", type=str, default="1-22", help="Chromosome ids to use "
         "(e.g. 1,2,3 or 1-4,12,16-20). Chromosome must be labeled by integer, i.e. X and Y are not acceptable. ")
    parser.add_argument("--trait1-file", type=str, default=None, help="GWAS summary statistics for the first trait. ")
    if num_traits==2: parser.add_argument("--trait2-file", type=str, default="", help="GWAS summary statistics for the second trait. ")
    parser.add_argument('--z1max', type=float, default=None, help="right-censoring threshold for the first trait. ")
    if num_traits==2: parser.add_argument('--z2max', type=float, default=None, help="right-censoring threshold for the second trait. ")
    parser.add_argument('--extract', type=str, default="", help="File with variants to include in the analysis")
    parser.add_argument('--exclude', type=str, default="", help="File with variants to exclude from the analysis")
    parser.add_argument('--randprune-n', type=int, default=64, help="Number of random pruning iterations")
    parser.add_argument('--randprune-r2', type=float, default=0.1, help="Threshold for random pruning")
    parser.add_argument('--seed', type=int, default=123, help="Random seed")
    parser.add_argument('--r2min', type=float, default=0.0, help="r2 values below this threshold will contribute via infinitesimal model")
    parser.add_argument('--threads', type=int, default=[None], nargs='+', help="specify how many threads to use (concurrency). None will default to the total number of CPU cores. ")
    parser.add_argument('--cubature-rel-error', type=float, default=1e-5, help="relative error for cubature stop criteria (applies to 'convolve' cost calculator). ")
    parser.add_argument('--cubature-max-evals', type=float, default=1000, help="max evaluations for cubature stop criteria (applies to 'convolve' cost calculator). "
        "Bivariate cubature require in the order of 10^4 evaluations and thus is much slower than sampling, therefore it is not exposed via mixer.py command-line interface. ")

def parser_fit_or_test_add_arguments(args, func, parser, do_fit, num_traits):
    parser_add_common_arguments(parser, num_traits)

    if do_fit and (num_traits==1):
        parser.add_argument('--analysis', type=str, default='fit1', help=argparse.SUPPRESS, choices=['fit1', 'fit2', 'test1', 'test2'])
        parser.add_argument('--fit-sequence', type=str, default=['diffevo-fast', 'neldermead'], nargs='+', help=argparse.SUPPRESS, choices=['diffevo', 'diffevo-fast', 'neldermead', 'neldermead-fast'])
    elif do_fit and (num_traits==2):
        parser.add_argument('--analysis', type=str, default='fit2', help=argparse.SUPPRESS, choices=['fit1', 'fit2', 'test1', 'test2'])
        parser.add_argument('--fit-sequence', type=str, default=['diffevo-fast', 'neldermead-fast', 'brute1', 'brent1'], nargs='+', help=argparse.SUPPRESS, choices=['diffevo-fast', 'neldermead-fast', 'brute1', 'brute1-fast', 'brent1', 'brent1-fast'])
        parser.add_argument('--trait1-params-file', type=str, default=None, help="univariate params for the first trait (for the cross-trait analysis only). ")
        parser.add_argument('--trait2-params-file', type=str, default=None, help="univariate params for the second trait (for the cross-trait analysis only). ")
    elif (not do_fit) and (num_traits==1):
        parser.add_argument('--analysis', type=str, default='test1', help=argparse.SUPPRESS, choices=['fit1', 'fit2', 'test1', 'test2'])
        parser.add_argument('--fit-sequence', type=str, default=['load', 'inflation'], nargs='+', help=argparse.SUPPRESS, choices=['load', 'inflation'])
    elif (not do_fit) and (num_traits==2):
        parser.add_argument('--analysis', type=str, default='test2', help=argparse.SUPPRESS, choices=['fit1', 'fit2', 'test1', 'test2'])
        parser.add_argument('--fit-sequence', type=str, default=['load', 'inflation'], nargs='+', help=argparse.SUPPRESS, choices=['load', 'inflation'])
    else:
        raise ValueError('internal error: invalid combination of do_fit and num_traits')

    if do_fit:
        parser.add_argument('--kmax', type=int, default=[20000], nargs='+', help="Number of sampling iterations")
    else:
        parser.add_argument('--kmax', type=int, default=[100], nargs='+', help="Number of sampling iterations")
        parser.add_argument('--load-params-file', type=str, default=None, help="params of the fitted model (for 'mixer.py test1' and 'mixer.py test2' runs). ")

    # all arguments below marked with argparse.SUPRESS option are internal and not recommended for a general use.
    # Valid options for --fit-sequence (remember to combine them in a sequence that makes sence):
    #       'load' reads previosly fitted parameters from a file (--load-params-file);
    #       'diffevo' performs a single iteration of differential evolution, which is the right way to setup an initial approximation;
    #       'neldermead' applies Nelder-Mead downhill simplex search;
    #       'brute1' applies only to bivariate optimization; it performs brute-force for one-dimentional optimization, searching optimal pi12 value constrained on genetic correlation (rg) and intercept (rho0);
    #       'brent1' is similar to brute1, but it uses brent method (inverse parabolic interpolation);
    #       'inflation' fits sig2zero (univariate) and rho_zero (bivariate), using fast cost function; this is quite special optimization step, typically useful for adjusting inflation parameters to another reference;
    #       'infinitesimal' fits a model with pi1=1 (univariate) or pi12=1 (bivariate) constrains, using fast cost function; this is quite special optimization step, typically used internally for AIC/BIC computation;
    # Note that bivariate fit is always constrained on univariate parameters, except for 'inflation' fit which adjust rho_zero and sig2_zero.
    # The '...-fast' optimizations use fast cost function.
    # Note that univariate optimization uses 'convolve' cost calculator, bivariate optimization uses 'sampling' cost calculator.
    # Typical univariate sequence: 'diffevo-fast neldermead'
    # Typical bivariate sequence: 'diffevo neldermead brute1 brent1'

    parser.add_argument('--diffevo-fast-repeats', type=int, default=20, help=argparse.SUPPRESS)      # repeat --diffevo-fast step this many times and choose the best run
    parser.add_argument('--downsample-factor', default=50, type=int, help=argparse.SUPPRESS)         # Applies to --power-curve and --qq-plots, --downsample-factor N' imply that only 1 out of N available z-score values will be used in model calculations.

    parser.add_argument('--power-curve', default=(not do_fit), action="store_true", help=argparse.SUPPRESS) # generate power curves
    parser.add_argument('--qq-plots', default=(not do_fit), action="store_true", help=argparse.SUPPRESS)    # generate qq plot curves

    # Replace with The Sandwich Estimator ? http://www.stat.umn.edu/geyer/5601/notes/sand.pdf and re-implement CI for bivariate analysis?
    if num_traits==1:
        parser.add_argument('--ci-alpha', type=float, default=None, help=argparse.SUPPRESS)              # significance level for the confidence interval estimation
        parser.add_argument('--ci-samples', type=int, default=10000, help=argparse.SUPPRESS)             # number of samples in uncertainty estimation
        parser.add_argument('--ci-power-samples', type=int, default=100, help=argparse.SUPPRESS)         # number of samples in power curves uncertainty estimation

    parser.set_defaults(func=func)

def parser_ld_add_arguments(args, func, parser):
    parser.add_argument("--bfile", type=str, default=None, help="Path to plink bfile. ")
    parser.add_argument('--r2min', type=float, default=0.05, help="r2 values above this threshold will be stored in sparse LD format")
    parser.add_argument('--ldscore-r2min', type=float, default=0.001, help="r2 values above this threshold (and below --r2min) will be stored as LD scores that contribute to the cost function via an infinitesimal model")
    parser.add_argument('--ld-window-kb', type=float, default=0, help="limit window similar to --ld-window-kb in 'plink r2'; 0 will disable this constraint")
    parser.add_argument('--ld-window', type=int, default=0, help="limit window similar to --ld-window in 'plink r2'; 0 will disable this constraint")
    parser.set_defaults(func=func)

def parser_snps_add_arguments(args, func, parser):
    parser.add_argument("--bim-file", type=str, default=None, help="Plink bim file. "
        "Defines the reference set of SNPs used for the analysis. "
        "Marker names must not have duplicated entries. "
        "May contain simbol '@', which will be replaced with the actual chromosome label. ")
    parser.add_argument("--ld-file", type=str, default=None, help="File with linkage disequilibrium information, "
        "generated via 'mixer.py ld' command. "
        "May contain simbol '@', similarly to --bim-file argument. ")
    parser.add_argument("--chr2use", type=str, default="1-22", help="Chromosome ids to use "
         "(e.g. 1,2,3 or 1-4,12,16-20). Chromosome must be labeled by integer, i.e. X and Y are not acceptable. ")
    parser.add_argument('--r2', type=float, default=0.8, help="r2 threshold for random prunning")
    parser.add_argument('--maf', type=float, default=0.05, help="maf threshold")
    parser.add_argument('--subset', type=int, default=2000000, help="number of SNPs to randomly select")
    parser.add_argument('--seed', type=int, default=123, help="Random seed")
    parser.set_defaults(func=func)

def parser_perf_add_arguments(args, func, parser):
    parser_add_common_arguments(parser, num_traits=2)
    parser.add_argument('--kmax', type=int, default=[20000, 2000, 200], nargs='+', help="Number of sampling iterations")
    parser.set_defaults(func=func)

def parse_args(args):
    parser = argparse.ArgumentParser(description="MiXeR: Univariate and Bivariate Causal Mixture for GWAS.")

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--argsfile', type=open, action=LoadFromFile, default=None, help="file with additional command-line arguments, e.g. those that are across all of your mixer.py runs (--lib, --bim-file and --ld-file)")
    parent_parser.add_argument("--out", type=str, default="mixer", help="prefix for the output files")
    parent_parser.add_argument("--lib", type=str, default="libbgmg.so", help="path to libbgmg.so plugin")
    parent_parser.add_argument("--log", type=str, default=None, help="file to output log, defaults to <out>.log")
    
    subparsers = parser.add_subparsers()

    parser_fit_or_test_add_arguments(args=args, func=execute_fit1_or_test1_parser, parser=subparsers.add_parser("fit1", parents=[parent_parser], help='fit univariate MiXeR model'), do_fit=True, num_traits=1)
    parser_fit_or_test_add_arguments(args=args, func=execute_fit1_or_test1_parser, parser=subparsers.add_parser("test1", parents=[parent_parser], help='test univariate MiXeR model'), do_fit=False, num_traits=1)
    parser_fit_or_test_add_arguments(args=args, func=execute_fit2_or_test2_parser, parser=subparsers.add_parser("fit2", parents=[parent_parser], help='fit bivariate MiXeR model'), do_fit=True, num_traits=2)
    parser_fit_or_test_add_arguments(args=args, func=execute_fit2_or_test2_parser, parser=subparsers.add_parser("test2", parents=[parent_parser], help='test bivariate MiXeR model'), do_fit=False, num_traits=2)

    parser_ld_add_arguments(args=args, func=execute_ld_parser, parser=subparsers.add_parser("ld", parents=[parent_parser], help='prepare files with linkage disequilibrium information'))
    parser_perf_add_arguments(args=args, func=execute_perf_parser, parser=subparsers.add_parser("perf", parents=[parent_parser], help='run performance evaluation of the MiXeR'))
    parser_snps_add_arguments(args=args, func=execute_snps_parser, parser=subparsers.add_parser("snps", parents=[parent_parser], help='generate random sets of SNPs'))

    return parser.parse_args(args)

def log_header(args, subparser_name, lib):
    defaults = vars(parse_args([subparser_name]))
    opts = vars(args)
    non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
    header = MASTHEAD
    header += "Call: \n"
    header += './mixer.py {} \\\n'.format(subparser_name)
    options = ['\t--'+x.replace('_','-')+' '+(' '.join([str(y) for y in opts[x]]) if isinstance(opts[x], list) else str(opts[x])).replace('\t', '\\t')+' \\' for x in non_defaults]
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
    return BivariateParams(pi=p['pi'], sig2_beta=p['sig2_beta'], sig2_zero=p['sig2_zero'], rho_beta=p['rho_beta'], rho_zero=p['rho_zero'])

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
            repeats = (1 if (fit_type == 'diffevo') else args.diffevo_fast_repeats)
            for repeat in range(repeats):
                optimize_result_tmp = scipy.optimize.differential_evolution(lambda x: parametrization.calc_cost(x), bounds4opt,
                    tol=0.01, mutation=(0.5, 1), recombination=0.7, seed=(args.seed + repeat), atol=0, updating='immediate', polish=False, workers=1)  #, **global_opt_options)
                params_tmp = parametrization.vec_to_params(optimize_result_tmp.x)
                libbgmg.log_message("--diffevo-fast-repeat={}: {}".format(repeat, params_tmp))
                if (repeat == 0) or (optimize_result_tmp.fun < optimize_result.fun):
                    optimize_result, params = optimize_result_tmp, params_tmp

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
    if args.analysis == 'fit2':
        libbgmg.log_message("Loading univariate constrains for bivariate analysis...")
        params1 = load_univariate_params_file(args.trait1_params_file)
        params2 = load_univariate_params_file(args.trait2_params_file)
        libbgmg.log_message("trait1: {}".format(params1))
        libbgmg.log_message("trait2: {}".format(params2))

    params = None; optimize_result_sequence = []
    for fit_type in args.fit_sequence:
        libbgmg.log_message("fit_type=={}...".format(fit_type))

        if fit_type == 'load':
            params = load_bivariate_params_file(args.load_params_file)
            params1 = params._params1()
            params2 = params._params2()
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
            repeats = (1 if (fit_type == 'diffevo') else args.diffevo_fast_repeats)
            for repeat in range(repeats):
                optimize_result_tmp = scipy.optimize.differential_evolution(lambda x: parametrization.calc_cost(x), bounds4opt,
                    tol=0.01, mutation=(0.5, 1), recombination=0.7, seed=(args.seed + repeat), atol=0, updating='immediate', polish=False, workers=1)  #, **global_opt_options)
                params_tmp = parametrization.vec_to_params(optimize_result_tmp.x)
                libbgmg.log_message("--diffevo-fast-repeat={}: {}".format(repeat, params_tmp))
                if (repeat == 0) or (optimize_result_tmp.fun < optimize_result.fun):
                    optimize_result, params = optimize_result_tmp, params_tmp

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

def calc_bivariate_pdf(libbgmg, params, downsample):
    original_weights = libbgmg.weights
    if not np.all(np.isfinite(original_weights)): raise(RuntimeError('undefined weights not supported'))
    if not np.all(np.isfinite(libbgmg.zvec1)): raise(RuntimeError('undefined weights not supported'))
    if not np.all(np.isfinite(libbgmg.zvec2)): raise(RuntimeError('undefined weights not supported'))

    model_weights = libbgmg.weights
    mask = np.zeros((len(model_weights), ), dtype=bool)
    mask[range(0, len(model_weights), downsample)] = 1
    model_weights[~mask] = 0
    model_weights = model_weights/np.sum(model_weights)

    zgrid = np.arange(-25, 25.00001, 0.05)
    [zgrid1, zgrid2] = np.meshgrid(zgrid, zgrid)
    zgrid1=zgrid1[zgrid>=0, :]; zgrid2=zgrid2[zgrid>=0, :]

    libbgmg.weights = model_weights     # temporary set downsampled weights
    pdf = params.pdf(libbgmg, zgrid)
    pdf = pdf/np.sum(model_weights)
    libbgmg.weights = original_weights  # restore original weights

    return zgrid, pdf

def calc_bivariate_qq(libbgmg, zgrid, pdf):
    zgrid_fine = np.arange(-25, 25.00001, 0.005)   # project to a finer grid
    pdf_fine=scipy.interpolate.interp2d(zgrid, zgrid, pdf)(zgrid_fine, zgrid_fine)
    pthresh_vec = [1, 0.1, 0.01, 0.001]
    zthresh_vec = -scipy.stats.norm.ppf(np.array(pthresh_vec)/2)

    zvec1=libbgmg.zvec1
    zvec2=libbgmg.zvec2

    # Regular grid (vertical axis of the QQ plots)
    hv_z = np.linspace(0, np.min([np.max(np.abs(np.concatenate((zvec1, zvec2)))), 38.0]), 1000)
    hv_logp = -np.log10(2*scipy.stats.norm.cdf(-hv_z))

    result = []
    for zthresh, pthresh in zip(zthresh_vec, pthresh_vec):
        mask = abs(zvec2)>=zthresh
        data_logpvec = calc_qq_data(zvec1[mask], libbgmg.weights[mask], hv_logp)

        pd_cond = np.sum(pdf_fine[abs(zgrid_fine) >= zthresh, :], axis=0)
        pd_cond = pd_cond / np.sum(pd_cond) / (zgrid_fine[1]-zgrid_fine[0])
        model_logpvec = calc_qq_model(zgrid_fine, pd_cond, hv_z)

        title = 'T1|T2|{}'.format(pthresh)
        result.append({'hv_logp': hv_logp, 'data_logpvec': data_logpvec, 'model_logpvec': model_logpvec,
                       'n_snps': int(np.sum(mask)), 'sum_data_weights': float(np.sum(libbgmg.weights[mask])), 'title' : title})

    for zthresh in zthresh_vec:
        mask = abs(zvec1)>=zthresh
        data_logpvec = calc_qq_data(zvec2[mask], libbgmg.weights[mask], hv_logp)

        pd_cond = np.sum(pdf_fine[:, abs(zgrid_fine) >= zthresh], axis=1)
        pd_cond = pd_cond / np.sum(pd_cond) / (zgrid_fine[1]-zgrid_fine[0])
        model_logpvec = calc_qq_model(zgrid_fine, pd_cond, hv_z)

        title = 'T2|T1|{}'.format(pthresh)
        result.append({'hv_logp': hv_logp, 'data_logpvec': data_logpvec, 'model_logpvec': model_logpvec,
                       'n_snps': int(np.sum(mask)), 'sum_data_weights': float(np.sum(libbgmg.weights[mask])), 'title' : title})
    return result

# helper function to debug non-json searizable types...
def print_types(results, libbgmg):
    if isinstance(results, dict):
        for k, v in results.items():
            libbgmg.log_message('{}: {}'.format(k, type(v)))
            print_types(v, libbgmg)

def execute_ld_parser(args):
    libbgmg = LibBgmg(args.lib)
    libbgmg.calc_ld_matrix(args.bfile, args.out, args.r2min, args.ldscore_r2min, args.ld_window, args.ld_window_kb)
    libbgmg.log_message('Done')

def initialize_mixer_plugin(args):
    libbgmg = LibBgmg(args.lib)
    libbgmg.init(args.bim_file, "", args.chr2use,
                 args.trait1_file if ('trait1_file' in args) else "",
                 args.trait2_file if ('trait2_file' in args) else "",
                 args.exclude if ('exclude' in args) else "",
                 args.extract if ('extract' in args) else "")

    for opt, val in convert_args_to_libbgmg_options(args, libbgmg.num_snp):
        libbgmg.set_option(opt, val)

    for chr_label in args.chr2use: 
        libbgmg.set_ld_r2_coo_from_file(chr_label, args.ld_file.replace('@', str(chr_label)))
        libbgmg.set_ld_r2_csr(chr_label)

    if ('randprune_n' in args) and ('randprune_r2' in args):
        libbgmg.set_weights_randprune(args.randprune_n, args.randprune_r2, exclude="", extract="")
    
    libbgmg.set_option('diag', 0)
    return libbgmg

def execute_perf_parser(args):
    fix_and_validate_args(args)
    libbgmg = initialize_mixer_plugin(args)
    params1 = UnivariateParams(pi=3e-3, sig2_beta=1e-5, sig2_zero=1.05)
    params12 = BivariateParams(params1=params1, params2=params1, pi12=2e-3, rho_beta=0.5, rho_zero=0.3)

    perf_data = []
    for threads in args.threads:
        for kmax in args.kmax:
            for costcalc in _cost_calculator:
                libbgmg.set_option('threads', threads)
                libbgmg.set_option('kmax', kmax)
                libbgmg.set_option('cost_calculator', _cost_calculator[costcalc])

                if (kmax!=np.max(args.kmax)) and (costcalc in ['gaussian', 'convolve']):
                    continue  # gaussian and convolve are calculated only for one value of kmax, e.g. for the largest one
                if (costcalc == 'convolve') and (args.trait2_file):
                    continue  # skip convolve for bivariate analysis

                start = time.time()
                cost = (params12.cost(libbgmg) if args.trait2_file else params1.cost(libbgmg, trait=1))
                end = time.time()
                perf_data.append((threads, kmax, costcalc, end-start, cost))
                libbgmg.log_message('threads={}, kmax={}, costcalt={} took {} seconds'.format(threads, kmax, costcalc, end-start))

    pd.DataFrame(perf_data, columns=['threads', 'kmax', 'costcalc', 'time_sec', 'cost']).to_csv(args.out + '.csv', sep='\t', index=False)
    libbgmg.log_message('Done')

def execute_snps_parser(args):
    fix_and_validate_args(args)
    if args.seed is not None: np.random.seed(args.seed)
    libbgmg = initialize_mixer_plugin(args)
    mafvec = np.minimum(libbgmg.mafvec, 1-libbgmg.mafvec)

    libbgmg.log_message('Load {}...'.format(args.bim_file))
    ref=pd.concat([pd.read_csv(args.bim_file.replace('@', str(chr_label)), sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split()) for chr_label in args.chr2use])
    libbgmg.log_message('{} SNPs in total'.format(len(ref)))

    # step0 - generate random values for clumping (this is a way to implement random pruning)
    buf = np.random.rand(libbgmg.num_tag, 1)
    
    # step1 - filter SNPs below MAF threshold
    buf[mafvec<args.maf] = np.nan
    libbgmg.log_message('{} SNPs pass --maf {} threshold'.format(np.sum(np.isfinite(buf)), args.maf))
    
    # step2 - select a random subset of SNPs
    indices = list(np.where(np.isfinite(buf))[0])
    sample = random.sample(indices, min(args.subset, len(indices)))
    mask = np.ones(buf.shape); mask[sample] = 0
    buf[mask == 1] = np.nan
    libbgmg.log_message('{} SNPs randomly selected (--subset)'.format(np.sum(np.isfinite(buf))))

    # step3 - further prune SNPs at certain r2 threshold
    buf_pruned = libbgmg.perform_ld_clump(args.r2, buf.flatten()) 
    libbgmg.log_message('{} SNPs pass random prunning at --r2 {} threshold'.format(np.sum(np.isfinite(buf_pruned)), args.r2))

    # step4 - save results
    ref.SNP[np.isfinite(buf_pruned)].to_csv(args.out, index=False, header=False)
    libbgmg.log_message('Result saved to {}'.format(args.out))
    
    libbgmg.log_message('Done')

def init_results_struct(libbgmg, args):
    results = {}
    results['options'] = vars(args).copy()
    results['options']['totalhet'] = float(2.0 * np.dot(libbgmg.mafvec, 1.0 - libbgmg.mafvec))
    results['options']['num_snp'] = float(libbgmg.num_snp)
    results['options']['num_tag'] = float(libbgmg.num_tag)
    results['options']['sum_weights'] = float(np.sum(libbgmg.weights))
    results['options']['trait1_nval'] = float(np.nanmedian(libbgmg.get_nvec(trait=1)))
    results['ci'] = {}
    return results

def execute_fit1_or_test1_parser(args):
    fix_and_validate_args(args)
    libbgmg = initialize_mixer_plugin(args)
    results = init_results_struct(libbgmg, args)
    results['analysis'] = 'univariate'
    totalhet = results['options']['totalhet']

    libbgmg.log_message('--fit-sequence: {}...'.format(args.fit_sequence))

    params, optimize_result = apply_univariate_fit_sequence(args, libbgmg, args.fit_sequence)
    results['params'] = params.as_dict()
    results['optimize'] = optimize_result

    libbgmg.log_message('Calculate AIC/BIC w.r.t. infinitesimal model (fast cost function)...')
    params_inft, optimize_result_inft = apply_univariate_fit_sequence(args, libbgmg, ['infinitesimal'], init_params=params)
    results['inft_params'] = params_inft.as_dict()
    results['inft_optimize'] = optimize_result_inft

    funcs, _ = _calculate_univariate_uncertainty_funcs(None, totalhet, libbgmg.num_snp)
    for func_name, func in funcs: results['ci'][func_name] = {'point_estimate': func(params)}

    ci_sample = None
    if args.ci_alpha and np.isfinite(args.ci_alpha):
        libbgmg.log_message("Uncertainty estimation...")
        libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
        results['ci'], ci_sample = _calculate_univariate_uncertainty(UnivariateParametrization(params, libbgmg, trait=1), args.ci_alpha, totalhet, libbgmg.num_snp, args.ci_samples)
        for k, v in results['ci'].items():
            libbgmg.log_message('{}: point_estimate={:.3g}, mean={:.3g}, median={:.3g}, std={:.3g}, ci=[{:.3g}, {:.3g}]'.format(k, v['point_estimate'], v['mean'], v['median'], v['std'], v['lower'], v['upper']))
        libbgmg.log_message("Uncertainty estimation done.")

    if args.power_curve:
        trait_index = 1
        results['power'] = calc_power_curve(libbgmg, params, trait_index, args.downsample_factor)
        if ci_sample is not None:
            power_ci = []
            for ci_index, ci_params in enumerate(ci_sample[:args.ci_power_samples]):
                libbgmg.log_message("Power curves uncertainty, {} of {}".format(ci_index, args.ci_power_samples))
                power_ci.append(calc_power_curve(libbgmg, ci_params, trait_index, args.downsample_factor))
            results['power_ci'] = power_ci

    if args.qq_plots:
        trait_index = 1
        mask = np.ones((libbgmg.num_tag, ), dtype=bool)
        results['qqplot'] = calc_qq_plot(libbgmg, params, trait_index, args.downsample_factor, mask,
            title='maf \\in [{:.3g},{:.3g}); L \\in [{:.3g},{:.3g})'.format(-np.inf,np.inf,-np.inf,np.inf))

        mafvec = libbgmg.mafvec[libbgmg.defvec]
        tldvec = libbgmg.ld_sum_r2[libbgmg.defvec]
        maf_bins = np.concatenate(([-np.inf], np.quantile(mafvec, [1/3, 2/3]), [np.inf]))
        tld_bins = np.concatenate(([-np.inf], np.quantile(tldvec, [1/3, 2/3]), [np.inf]))
        results['qqplot_bins'] = []
        for i in range(0, 3):
            for j in range(0, 3):
                mask = ((mafvec>=maf_bins[i]) & (mafvec<maf_bins[i+1]) & (tldvec >= tld_bins[j]) &  (tldvec < tld_bins[j+1]))
                results['qqplot_bins'].append(calc_qq_plot(libbgmg, params, trait_index, args.downsample_factor, mask,
                    title='maf \\in [{:.3g},{:.3g}); L \\in [{:.3g},{:.3g})'.format(maf_bins[i], maf_bins[i+1], tld_bins[j], tld_bins[j+1])))

    with open(args.out + '.json', 'w') as outfile:
        json.dump(results, outfile, cls=NumpyEncoder)

    libbgmg.set_option('diag', 0)
    libbgmg.log_message('Done')

def execute_fit2_or_test2_parser(args):
    fix_and_validate_args(args)
    libbgmg = initialize_mixer_plugin(args)
    results = init_results_struct(libbgmg, args)
    results['analysis'] = 'bivariate'
    results['options']['trait2_nval'] = float(np.nanmedian(libbgmg.get_nvec(trait=2)))
    totalhet = results['options']['totalhet']

    libbgmg.log_message('--fit-sequence: {}...'.format(args.fit_sequence))
    params, _, _, optimize_result = apply_bivariate_fit_sequence(args, libbgmg)
    results['params'] = params.as_dict()
    results['optimize'] = optimize_result

    funcs, _ = _calculate_bivariate_uncertainty_funcs(None, totalhet, libbgmg.num_snp)
    for func_name, func in funcs: results['ci'][func_name] = {'point_estimate': func(params)}

    if args.qq_plots:
        zgrid, pdf = calc_bivariate_pdf(libbgmg, params, args.downsample_factor)
        results['pdf_zgrid'] = zgrid
        results['pdf'] = pdf
        results['qqplot'] = calc_bivariate_qq(libbgmg, zgrid, pdf)

    with open(args.out + '.json', 'w') as outfile:
        json.dump(results, outfile, cls=NumpyEncoder)

    libbgmg.set_option('diag', 0)
    libbgmg.log_message('Done')
