#!/usr/bin/env python
'''
(c) 2018 Oleksandr Frei
Univariate and Bivariate Causal Mixture for GWAS
'''

import sys
import os
import argparse
import pandas as pd
import numpy as np
import scipy.optimize
import logging
import json

from libbgmg import LibBgmg
from utils import MixerOptimizeResult
from utils import UnivariateParams
from utils import BivariateParams
from utils import _log_exp_converter
from utils import _logit_logistic_converter
from utils import _arctanh_tanh_converter
from utils import UnivariateParametrization_constPI
from utils import UnivariateParametrization_constH2_constSIG2ZERO
from utils import UnivariateParametrization_constH2_constSIG2ZERO_boundedPI
from utils import UnivariateParametrization_constPI_constSIG2BETA
from utils import UnivariateParametrization
from utils import BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO
from utils import BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO_boundedPI
from utils import BivariateParametrization_constSIG2BETA_constSIG2ZERO_infPI_maxRG
from utils import BivariateParametrization_constUNIVARIATE
from utils import BivariateParametrization_constUNIVARIATE_natural_axis
from utils import BivariateParametrization_constUNIVARIATE_constRHOBETA_constPI
from utils import _hessian_robust
from utils import _max_rg
from utils import _calculate_univariate_uncertainty
from utils import _calculate_bivariate_uncertainty

__version__ = '1.0.0'
MASTHEAD = "***********************************************************************\n"
MASTHEAD += "* mixer.py: Univariate and Bivariate Causal Mixture for GWAS\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (C) 2018 Oleksandr Frei\n"
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

def fix_and_validate_args(args):
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
        'z1max': args.z1max, 'z2max': args.z2max, 'cubature_rel_error': args.cubature_rel_error, 'cubature_max_evals':args.cubature_max_evals
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

def parse_args(args):
    parser = argparse.ArgumentParser(description="MiXeR: Univariate and Bivariate Causal Mixture for GWAS.")

    parser.add_argument("--out", type=str, default="mixer", help="prefix for the output files")
    parser.add_argument("--lib", type=str, default="libbgmg.so", help="path to libbgmg.so plugin")
    parser.add_argument("--log", type=str, default=None, help="file to output log, defaults to <out>.log")
    parser.add_argument('--argsfile', type=open, action=LoadFromFile, default=None, help="file with additional command-line arguments")

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
         "(e.g. 1,2,3 or 1-4,12,16-20 or 19-22,X,Y). "
         "Chromosomes with non-integer ids should be indicated separately. ")
    parser.add_argument("--trait1-file", type=str, default=None, help="GWAS summary statistics for the first trait. ")
    parser.add_argument("--trait2-file", type=str, default="", help="GWAS summary statistics for the first trait. "
        "Specifying this argument triggers cross-trait analysis.")
    parser.add_argument('--full-cost-calc', type=str, default='sampling', choices=['sampling', 'convolve', 'gaussian'],
        help="calculator for the full cost function; note that 'convolve' is only recommended for the univariate fit.")
    parser.add_argument('--fit-sequence', type=str, default=['init', 'constrained'], nargs='+',
        choices=['load', 'init', 'inflation', 'infinitesimal', 'constrained', 'full'],
        help="Specify fit sequence: "
             "'load' reads previosly fitted parameters from a file; "
             "'init' initialize fit using 'fast model'; "
             "'inflation' fits sig2zero (univariate) and rho_zero (bivariate), typically useful when adjusting parameters to another reference; "
             "'infinitesimal' fits a model with pi1=1 (univariate) or pi12=1 (bivariate) constrains; " 
             "'constrained' fits a model with constrained h2 (univariate) or rg (bivariate); "
             "'full' jointly fits the three parameters: sig2_zero, sig2_beta, pi (univariate) or rho_zero, rho_beta, pi12 (bivariate); " 
             "Note that bivariate fit is always constrained on univariate parameters, except for 'inflation' fit which adjust rho_zero and sig2_zero. "
             "All fit sequences must start with 'init' or 'load'. "
             "Typical sequences: 'init', 'init constrained', 'init full', 'load inflation'")
    parser.add_argument('--preliminary', default=False, action="store_true",
        help="perform an additional run using fast model to generate preliminary data. "
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
    parser.add_argument('--ci-alpha', type=float, default=0.05, help="significance level for the confidence interval estimation")
    parser.add_argument('--ci-samples', type=int, default=10000, help="number of samples in uncertainty estimation")
    parser.add_argument('--threads', type=int, default=None, help="specify how many threads to use (concurrency). None will default to the total number of CPU cores. ")
    parser.add_argument('--tol-x', type=float, default=1e-2, help="tolerance for the stop criteria in fminsearch optimization. ")
    parser.add_argument('--tol-func', type=float, default=1e-2, help="tolerance for the stop criteria in fminsearch optimization. ")
    parser.add_argument('--cubature-rel-error', type=float, default=1e-6, help="relative error for cubature stop criteria (only relevant for --full-cost-calc convolve). ")
    parser.add_argument('--cubature-max-evals', type=float, default=100, help="max evaluations for cubature stop criteria (only relevant for --full-cost-calc convolve). Bivariate cubature require in the order of 10^4 evaluations and thus is much slower than sampling. ")

    parser.add_argument('--z1max', type=float, default=None, help="right-censoring threshold for the first trait. ")
    parser.add_argument('--z2max', type=float, default=None, help="right-censoring threshold for the second trait. ")

    parser.add_argument('--load-params-file', type=str, default=None, help="initial params for the optimization. ")
    parser.add_argument('--trait1-params-file', type=str, default=None, help="univariate params for the first trait (for the cross-trait analysis only). ")
    parser.add_argument('--trait2-params-file', type=str, default=None, help="univariate params for the second trait (for the cross-trait analysis only). ")

    return parser.parse_args(args)

def log_header(args, lib):
    defaults = vars(parse_args([]))
    opts = vars(args)
    non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
    header = MASTHEAD
    header += "Call: \n"
    header += './mixer.py \\\n'
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

def apply_univariate_fit_sequence(args, libbgmg, optimizer, scalar_optimizer, fit_sequence, init_params=None, trait=1):
    # 'load', 'init', 'inflation', 'infinitesimal', 'constrained', 'full'
    params=init_params; optimization_result=None
    for fit_type in fit_sequence:
        if fit_type == 'load':
            libbgmg.log_message("fit_type==load: {}".format(args.load_params_file))
            params = load_univariate_params_file(args.load_params_file)
            libbgmg.log_message("fit_type==load: Done, {}".format(params))
        elif fit_type == 'init':
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            libbgmg.log_message("fit_type==init: UnivariateParametrization_constPI.fit(const_pi=1.0) with 'fast model'...")
            params, details = UnivariateParametrization_constPI(
                const_pi=1.0, init_sig2_zero=np.ma.var(np.ma.masked_invalid(libbgmg.get_zvec(trait))),
                init_sig2_beta=1.0/np.nanmean(libbgmg.get_nvec(trait)),
                lib=libbgmg, trait=trait).fit(optimizer)
            libbgmg.log_message("fit_type==init: intermediate {}".format(params))
            libbgmg.log_message("fit_type==init: UnivariateParametrization_constH2_constSIG2ZERO.fit() with 'fast model'...")
            params, details = UnivariateParametrization_constH2_constSIG2ZERO(init_pi=0.01, const_params=params,
                lib=libbgmg, trait=trait).fit(optimizer)
            libbgmg.log_message("fit_type==init: intermediate {}".format(params))
            libbgmg.log_message("fit_type==init: UnivariateParametrization.fit() with 'fast model'...")
            params, details = UnivariateParametrization(init_params=params, lib=libbgmg, trait=trait).fit(optimizer)
            optimization_result = MixerOptimizeResult(details, cost_n=np.sum(libbgmg.weights))
            libbgmg.log_message("fit_type==init: Done, {}, {}".format(params, optimization_result))
        elif fit_type == 'constrained':
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "constrained" fit'))
            libbgmg.set_option('cost_calculator', _cost_calculator[args.full_cost_calc])
            libbgmg.log_message("fit_type==constrained, UnivariateParametrization_constH2_constSIG2ZERO.fit(), 'full model...'")
            params, details = UnivariateParametrization_constH2_constSIG2ZERO(init_pi=params._pi, const_params=params,
                lib=libbgmg, trait=trait).fit(optimizer)
            optimization_result = MixerOptimizeResult(details, cost_n=np.sum(libbgmg.weights))
            libbgmg.log_message("fit_type==constrained: Done, {}, {}".format(params, optimization_result))
        elif fit_type == 'full':
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "full" fit'))
            libbgmg.set_option('cost_calculator', _cost_calculator[args.full_cost_calc])
            libbgmg.log_message("fit_type==full: UnivariateParametrization.fit(), 'full model'...")
            params, details = UnivariateParametrization(init_params=params, lib=libbgmg, trait=trait).fit(optimizer)
            optimization_result = MixerOptimizeResult(details, cost_n=np.sum(libbgmg.weights))
            libbgmg.log_message("fit_type==full: Done, {}, {}".format(params, optimization_result))
        elif fit_type == 'inflation':
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "inflation" fit'))
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            libbgmg.log_message("fit_type==inflation: UnivariateParametrization_constPI_constSIG2BETA.fit(), 'fast model'...")
            params, details = UnivariateParametrization_constPI_constSIG2BETA(
                init_sig2_zero=params._sig2_zero, const_params=params, lib=libbgmg, trait=trait).fit(optimizer)
            optimization_result = MixerOptimizeResult(details, cost_n=np.sum(libbgmg.weights))
            libbgmg.log_message("fit_type==inflation: Done, {}, {}".format(params, optimization_result))
        elif fit_type == 'infinitesimal':
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "infinitesimal" fit'))
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            libbgmg.log_message("fit_type==infinitesimal: UnivariateParametrization_constPI.fit(const_pi=1.0) with 'fast model'...")
            params, details = UnivariateParametrization_constPI(
                const_pi=1.0, init_sig2_zero=params._sig2_zero,
                init_sig2_beta=params._sig2_beta * params._pi,
                lib=libbgmg, trait=trait).fit(optimizer)
            optimization_result = MixerOptimizeResult(details, cost_n=np.sum(libbgmg.weights))
            libbgmg.log_message("fit_type==infinitesimal: Done, {}, {}".format(params, optimization_result))
        else:
            libbgmg.log_message("Unable to apply {} in univariate fit".format(fit_type))

    if params == None: raise(RuntimeError('Empty --fit-sequence'))
    
    libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
    optimization_result._cost_fast = params.cost(libbgmg, trait)
    return params, optimization_result

def apply_bivariate_fit_sequence(args, libbgmg, optimizer, scalar_optimizer):
    # 'load', 'init', 'inflation', 'infinitesimal', 'constrained', 'full'

    if args.trait1_params_file or args.trait2_params_file:
        libbgmg.log_message("Loading univariate constrains for bivariate analysis...")
        params1 = load_univariate_params_file(args.trait1_params_file)
        params2 = load_univariate_params_file(args.trait2_params_file)
        libbgmg.log_message("trait1: {}".format(params1))
        libbgmg.log_message("trait2: {}".format(params2))

    if 'init.full' in '.'.join(args.fit_sequence):
        libbgmg.log_message("WARNING: 'full' fit may behave weird after 'init' in bivariate optimization; "
            "instead, we recommend to apply the following sequence: 'init'->'constrained'->'full'")

    params = None; optimization_result = None
    for fit_type in args.fit_sequence:
        if fit_type == 'load':
            libbgmg.log_message("fit_type==load: {}".format(args.load_params_file))
            params, params1, params2 = load_bivariate_params_file(args.load_params_file)
            libbgmg.log_message("fit_type==load: Done,   {}".format(params))
            libbgmg.log_message("fit_type==load: Trait1: {}".format(params1))            
            libbgmg.log_message("fit_type==load: Trait2: {}".format(params2))            
        elif fit_type == 'init':
            if (params1==None) or (params2==None): raise(RuntimeError('params1==None or params2==None, unable to proceed apply "init" fit'))
            libbgmg.log_message("fit_type==init: BivariateParametrization_constUNIVARIATE.fit(), 'fast model'...")
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            zcorr = np.ma.corrcoef(np.ma.masked_invalid(libbgmg.zvec1), np.ma.masked_invalid(libbgmg.zvec2))[0, 1]
            params, details = BivariateParametrization_constUNIVARIATE(
                const_params1=params1, const_params2=params2,
                init_pi12=min(params1._pi, params2._pi)*0.1,
                init_rho_beta=zcorr, init_rho_zero=zcorr,
                lib=libbgmg).fit(optimizer)
            optimization_result = MixerOptimizeResult(details, cost_n=np.sum(libbgmg.weights))
            libbgmg.log_message("fit_type==init: Done, {}, {}".format(params, optimization_result))
            
        elif fit_type == 'constrained':
            # Constrained optimization intgentionally does params._pi[2] to initialize,
            # because "fast cost" function used in 'init' sometimes leads to biased estimates,
            # which give to bad initial guess for the "full cost".
            # Hense we use scalar_optimizer (golden search on bounded interval) which does not require x0.

            if (params1==None) or (params2==None): raise(RuntimeError('params1==None or params2==None, unable to proceed apply "init" fit'))
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "constrained" fit'))
            max_pi12 = np.min([params1._pi, params2._pi])
            min_pi12 = abs(params._rg()) * np.sqrt(params1._pi * params2._pi)
            init_pi12 = np.sqrt(np.max([min_pi12, np.min([max_pi12/2.0, 1e-6])]) * max_pi12)
            libbgmg.log_message("fit_type==constrained: BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO.fit(init_pi12={}), 'full model'...".format(init_pi12))
            libbgmg.set_option('cost_calculator', _cost_calculator[args.full_cost_calc])
            params, details = BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO(
                const_params1=params1, const_params2=params2,
                const_rg=params._rg(), const_rho_zero=params._rho_zero,
                init_pi12=init_pi12,
                lib=libbgmg).fit(optimizer)
            optimization_result = MixerOptimizeResult(details, cost_n=np.sum(libbgmg.weights))
            libbgmg.log_message("fit_type==constrained: Done {}".format(params, optimization_result))

        elif fit_type == 'full':
            if (params1==None) or (params2==None): raise(RuntimeError('params1==None or params2==None, unable to proceed apply "init" fit'))
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "full" fit'))
            libbgmg.log_message("fit_type==full: BivariateParametrization_constUNIVARIATE.fit(), 'full model'...")
            libbgmg.set_option('cost_calculator', _cost_calculator[args.full_cost_calc])
            params, details = BivariateParametrization_constUNIVARIATE(
                const_params1=params1, const_params2=params2,
                init_pi12=params._pi[2], init_rho_beta=params._rho_beta, init_rho_zero=params._rho_zero,
                lib=libbgmg).fit(optimizer)
            optimization_result = MixerOptimizeResult(details, cost_n=np.sum(libbgmg.weights))
            libbgmg.log_message("fit_type==full: Done, {}, {}".format(params, optimization_result))

        elif fit_type == 'inflation':
            if (params1==None) or (params2==None): raise(RuntimeError('params1==None or params2==None, unable to proceed apply "init" fit'))
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "inflation" fit'))
            params1, _ = apply_univariate_fit_sequence(args, libbgmg, optimizer, scalar_optimizer, ['inflation'], init_params=params1, trait=1)
            params2, _ = apply_univariate_fit_sequence(args, libbgmg, optimizer, scalar_optimizer, ['inflation'], init_params=params2, trait=2)
            libbgmg.log_message("fit_type==inflation: BivariateParametrization_constUNIVARIATE_constRHOBETA_constPI.fit(), 'fast model'...")
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            params, details = BivariateParametrization_constUNIVARIATE_constRHOBETA_constPI(
                const_params1=params1, const_params2=params2,
                const_pi12=params._pi[2], const_rho_beta=params._rho_beta, init_rho_zero=params._rho_zero,
                lib=libbgmg).fit(optimizer)
            optimization_result = MixerOptimizeResult(details, cost_n=np.sum(libbgmg.weights))
            libbgmg.log_message("fit_type==inflation: Done, {}, {}".format(params, optimization_result))

        elif fit_type == 'infinitesimal':
            if (params1==None) or (params2==None): raise(RuntimeError('params1==None or params2==None, unable to proceed apply "init" fit'))
            if params == None: raise(RuntimeError('params == None, unable to proceed apply "infinitesimal" fit'))
            libbgmg.log_message("fit_type==infinitesimal: BivariateParametrization_constSIG2BETA_constSIG2ZERO_infPI_maxRG.fit(), 'fast model'...")
            libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
            params, details = BivariateParametrization_constSIG2BETA_constSIG2ZERO_infPI_maxRG(
                const_sig2_beta=[params1._sig2_beta, params2._sig2_beta],
                const_sig2_zero=[params1._sig2_zero, params2._sig2_zero],
                max_rg = 1,
                init_rho_beta = params._rho_beta,
                init_rho_zero = params._rho_zero,
                lib=libbgmg).fit(optimizer)
            optimization_result = MixerOptimizeResult(details, cost_n=np.sum(libbgmg.weights))
            libbgmg.log_message("fit_type==infinitesimal: Done, {}, {}".format(params, optimization_result))

        else:
            libbgmg.log_message("Unable to apply {} in bivariate fit".format(fit_type))

    if params == None: raise(RuntimeError('Empty --fit-sequence'))
    libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
    optimization_result._cost_fast = params.cost(libbgmg)
    return (params, params1, params2, optimization_result)

# helper function to debug non-json searizable types...
def print_types(results, libbgmg):
    if isinstance(results, dict):
        for k, v in results.items():
            libbgmg.log_message('{}: {}'.format(k, type(v)))
            print_types(v, libbgmg)

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    args = parse_args(sys.argv[1:])
    libbgmg = LibBgmg(args.lib)
    libbgmg.init_log(args.log if args.log else args.out + '.log')
    log_header(args, libbgmg)
    libbgmg.dispose()

    fix_and_validate_args(args)

    libbgmg.set_option('use_complete_tag_indices', 1)
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

    # TBD: specify stop criteria
    nelder_optimizer = lambda func, x0: scipy.optimize.minimize(func, x0, method='Nelder-Mead')
    #scalar_optimizer = lambda func, xLeft, xRight: scipy.optimize.minimize_scalar(func,  method='bounded', bounds=[xLeft, xRight])
    scalar_optimizer = lambda func, xLeft, xRight: scipy.optimize.minimize_scalar(func,  method='Brent', bracket=[xLeft, xRight], options={'maxiter':50})
   
    totalhet = float(2.0 * np.dot(libbgmg.mafvec, 1.0 - libbgmg.mafvec))

    full_cost_calc = args.full_cost_calc
    fit_sequence = args.fit_sequence.copy()
    for repeat in range(2):
        if (repeat == 0) and (not args.preliminary): continue
        args.full_cost_calc = 'gaussian' if (repeat == 0) else full_cost_calc
        if (repeat == 0) and (args.fit_sequence[-1] != 'full'): args.fit_sequence.append('full')
        else: args.fit_sequence = fit_sequence
        libbgmg.log_message('Apply fit sequence with {} full cost calculator{}...'.format(args.full_cost_calc, ' (preliminary data)' if (repeat == 0) else '' ))

        results = {}
        results['options'] = vars(args).copy()
        results['options']['totalhet'] = totalhet
        results['options']['num_snp'] = float(libbgmg.num_snp)
        results['options']['num_tag'] = float(libbgmg.num_tag)
        results['options']['sum_weights'] = float(np.sum(libbgmg.weights))
        results['options']['trait1_nval'] = float(np.nanmedian(libbgmg.get_nvec(trait=1)))

        if not args.trait2_file:
            results['analysis'] = 'univariate'

            params, optimization_result = apply_univariate_fit_sequence(args, libbgmg, nelder_optimizer, scalar_optimizer, args.fit_sequence)
            optimization_result._cost_df = 3       # three free parameters
            results['params'] = params.as_dict()
            results['optimization'] = optimization_result.as_dict()

            libbgmg.log_message('Calculate AIC/BIC w.r.t. infinitesimal model (fast cost function)...')
            params_inft, optimization_result_inft = apply_univariate_fit_sequence(args, libbgmg, nelder_optimizer, scalar_optimizer, ['infinitesimal'], init_params=params)
            optimization_result_inft._cost_df = 2  # two free parameters        
            results['inft_params'] = params_inft.as_dict()
            results['inft_optimization'] = optimization_result_inft.as_dict()

            if np.isfinite(args.ci_alpha):
                libbgmg.log_message("Uncertainty estimation...")
                libbgmg.set_option('cost_calculator', _cost_calculator_gaussian)
                results['ci'], _ = _calculate_univariate_uncertainty(UnivariateParametrization(params, libbgmg, trait=1), args.ci_alpha, totalhet, libbgmg.num_snp, args.ci_samples)
                for k, v in results['ci'].items():
                    libbgmg.log_message('{}: point_estimate={:.3g}, mean={:.3g}, median={:.3g}, std={:.3g}, ci=[{:.3g}, {:.3g}]'.format(k, v['point_estimate'], v['mean'], v['median'], v['std'], v['lower'], v['upper']))
                libbgmg.log_message("Uncertainty estimation done.")
        else:
            results['analysis'] = 'bivariate'
            results['options']['trait2_nval'] = float(np.nanmedian(libbgmg.get_nvec(trait=2)))

            params, params1, params2, optimization_result = apply_bivariate_fit_sequence(args, libbgmg, nelder_optimizer, scalar_optimizer)
            optimization_result._cost_df = 9  # nine free parameters (incl. univariate)
            results['params'] = params.as_dict()
            if np.isfinite(args.ci_alpha):
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
            json.dump(results, outfile)
