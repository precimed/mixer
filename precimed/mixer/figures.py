#!/usr/bin/env python
'''
(c) 2016-2020 Oleksandr Frei, Alexey A. Shadrin, Dominic Holland
MiXeR software: Univariate and Bivariate Causal Mixture for GWAS
'''

import argparse
import json
import os
import itertools
import glob

import pandas as pd
import numpy as np
from numpy import ma

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors

from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib_venn import venn2

from scipy.interpolate import interp1d
from scipy.stats import multivariate_normal

from .utils import BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO
from .utils import _calculate_bivariate_uncertainty_funcs
from .utils import BivariateParams

__version__ = '1.2.0'
MASTHEAD = "***********************************************************************\n"
MASTHEAD += "* mixer_figures.py: Visualization tools for MiXeR\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (c) 2016-2020 Oleksandr Frei, Alexey A. Shadrin, Dominic Holland\n"
MASTHEAD += "* Norwegian Centre for Mental Disorders Research / University of Oslo\n"
MASTHEAD += "* Center for Multimodal Imaging and Genetics / UCSD\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "***********************************************************************\n"

def make_qq_plot(qq, ci=True, ylim=7.3, xlim=7.3):
    hv_logp = np.array(qq['hv_logp']).astype(float)
    data_logpvec = np.array(qq['data_logpvec']).astype(float)
    model_logpvec = np.array(qq['model_logpvec']).astype(float)
    ylim_data = max(hv_logp[np.isfinite(data_logpvec)])
    model_logpvec[hv_logp > ylim_data]=np.nan
    if ci:
        with np.errstate(all='ignore'):
            q = 10**-data_logpvec; dq= 1.96*np.sqrt(q*(1-q)/qq['sum_data_weights'])
            y1=hv_logp
            x1=ma.filled(-ma.log10(q+dq), np.nan)  #left CI bound
            x2=ma.filled(-ma.log10(q-dq), np.nan)  #right CI bound
            if True:
                y2 = np.empty(hv_logp.shape); y2[:]=np.nan
                y2[x2<np.nanmax(x1)]=interp1d(x1, y1)(x2[x2<np.nanmax(x1)])                #upper CI bound
                y2[np.isnan(y2)]=ylim_data  
                plt.fill_between(x2, y1, y2, color=(0.1843, 0.3098, 0.3098), alpha=0.25)
            else:
                plt.plot(x1,hv_logp,x2,hv_logp)
    hData = plt.plot(data_logpvec, hv_logp)
    hModel = plt.plot(model_logpvec, hv_logp)
    hNull = plt.plot(hv_logp, hv_logp, 'k--')
    plt.ylim(0, ylim); plt.xlim(0, xlim)

def make_venn_plot(data, flip=False, factor='K', traits=['Trait1', 'Trait2'], colors=[0, 1], max_size=None, formatter=None):
    cm = plt.cm.get_cmap('tab10')

    if factor=='K': scale_factor=1000
    elif factor=='': scale_factor=1
    else: raise(ValueError('Unknow factor: {}'.format(factor)))

    n1 = data['ci']['nc1@p9']['point_estimate']/scale_factor; n1_se = None
    n2 = data['ci']['nc2@p9']['point_estimate']/scale_factor; n2_se = None
    n12 = data['ci']['nc12@p9']['point_estimate']/scale_factor; n12_se = None
    rg = data['ci']['rg']['point_estimate']

    if max_size is None: max_size = n1+n2+n12
    if flip: n1, n2 = n2, n1; n1_se, n2_se = n2_se, n1_se
    f = lambda x: x if x < 7 else x+1

    v = venn2(subsets = (n1, n2, n12), normalize_to=(n1+n2+n12)/max_size, set_labels = ("", ""))
    v.get_patch_by_id('100').set_color(cm.colors[f(colors[0])])
    v.get_patch_by_id('010').set_color(cm.colors[f(colors[1])])
    v.get_patch_by_id('110').set_color(cm.colors[7])   
    if formatter==None:
        if (n1_se is not None) and (n2_se is not None) and (n12_se is not None):
            formatter1 = '{:.2f}\n({:.2f})' if ((n1+n12+n2) < 1) else '{:.1f}\n({:.1f})'
        else:
            formatter1 = '{:.2f}' if ((n1+n12+n2) < 1) else '{:.1f}'
        formatter = [formatter1, formatter1, formatter1]
    v.get_label_by_id('100').set_text(formatter[0].format(n1, n1_se))
    v.get_label_by_id('010').set_text(formatter[1].format(n2, n2_se))
    v.get_label_by_id('110').set_text(formatter[2].format(n12, n12_se))

    plt.xlim([-0.75, 0.75]), plt.ylim([-0.7, 0.6])
    newline=''
    plt.title(traits[0] +' & ' + newline + traits[1], y=-0.18)

    clr = plt.cm.get_cmap('seismic')((rg+1)/2)
    plt.gca().add_patch(patches.Rectangle(((-abs(0.7*rg) if (rg < 0) else 0) , -0.7), abs(0.7 * rg), 0.15, fill=True, clip_on=False, color=clr))
    plt.gca().add_patch(patches.Rectangle((-0.70, -0.7), 1.4, 0.15, fill=False, clip_on=False))
    plt.gca().add_patch(patches.Rectangle((0, -0.7), 0, 0.15, fill=False, clip_on=False, linewidth=3))
    plt.gca().text(-0.35 if (rg>0) else 0.35, -0.7+0.15/2, '$r_g$={:.2f}'.format(rg), fontsize=11, horizontalalignment='center',         verticalalignment='center')

def make_strat_qq_plots(data, flip=False, traits=['Trait1', 'Trait2'], do_legend=True):
    cm = plt.cm.get_cmap('tab10')
    for i in np.array(range(0, 4)) + (4 if flip else 0):
        hData = plt.plot(data['qqplot'][i]['data_logpvec'], data['qqplot'][i]['hv_logp'], color=cm.colors[i % 4], linestyle='solid')
    hNull = plt.plot(data['qqplot'][i]['hv_logp'], data['qqplot'][i]['hv_logp'], 'k--')
    if do_legend: plt.legend(['All SNPs'] + '$P{0}\leq0.1$ $P{0}\leq0.01$ $P{0}\leq0.001$'.format('_{' + traits[1]+'}').split(), loc='lower right', fontsize=14, borderpad=0.2, frameon=False, borderaxespad=0.2, labelspacing=0.1)
    for i in np.array(range(0, 4)) + (4 if flip else 0):
        hModel = plt.plot(data['qqplot'][i]['model_logpvec'], data['qqplot'][i]['hv_logp'], color=cm.colors[i % 4], linestyle='dashed')
    plt.ylim(0, 7.3); plt.xlim(0, 7.3); 
    plt.title('{} | {}'.format(traits[0], traits[1]))

def merge_z_vs_z(df1, df2):
    import pandas as pd

    _N_CHR = 22
    # complementary bases
    COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # bases
    BASES = COMPLEMENT.keys()
    # true iff strand ambiguous
    STRAND_AMBIGUOUS = {''.join(x): x[0] == COMPLEMENT[x[1]]
                        for x in itertools.product(BASES, BASES)
                        if x[0] != x[1]}
    # SNPS we want to keep (pairs of alleles)
    VALID_SNPS = {x for x in map(lambda y: ''.join(y), itertools.product(BASES, BASES))
                  if x[0] != x[1] and not STRAND_AMBIGUOUS[x]}
    # T iff SNP 1 has the same alleles as SNP 2 (allowing for strand or ref allele flip).
    MATCH_ALLELES = {x for x in map(lambda y: ''.join(y), itertools.product(VALID_SNPS, VALID_SNPS))
                     # strand and ref match
                     if ((x[0] == x[2]) and (x[1] == x[3])) or
                     # ref match, strand flip
                     ((x[0] == COMPLEMENT[x[2]]) and (x[1] == COMPLEMENT[x[3]])) or
                     # ref flip, strand match
                     ((x[0] == x[3]) and (x[1] == x[2])) or
                     ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))}  # strand and ref flip
    # T iff SNP 1 has the same alleles as SNP 2 w/ ref allele flip.
    FLIP_ALLELES = {''.join(x):
                    ((x[0] == x[3]) and (x[1] == x[2])) or  # strand match
                    # strand flip
                    ((x[0] == COMPLEMENT[x[3]]) and (x[1] == COMPLEMENT[x[2]]))
                    for x in MATCH_ALLELES} 

    df1 = df1[['SNP', 'A1', 'A2', 'Z']].rename(columns={'Z': 'Z1'}).copy()
    df2 = df2[['SNP', 'A1', 'A2', 'Z']].rename(columns={'Z': 'Z2', 'A1': 'A1x', 'A2': 'A2x'}).copy()
    df = pd.merge(df1, df2, how='inner', on='SNP')
    df = df.dropna(how='any')
    alleles = df.A1 + df.A2 + df.A1x + df.A2x
    df = df[alleles.apply(lambda y: y in MATCH_ALLELES)]
    alleles = df.A1 + df.A2 + df.A1x + df.A2x
    flip_status = alleles.apply(lambda y: FLIP_ALLELES[y])
    df.Z2 *= (-1) ** alleles.apply(lambda y: FLIP_ALLELES[y])
    df = df.drop(['A1', 'A1x', 'A2', 'A2x'], axis=1)
    return df

def plot_z_vs_z_data(df, flip=False, traits=['Trait1', 'Trait2'], plot_limits=15, bins=100):
    '''
        # input can be generated as follows:
        import pandas as pd
        df1 = pd.read_csv(fname1, delim_whitespace=True, usecols=['SNP', 'A1', 'A2', 'Z'])
        df2 = pd.read_csv(fname2, delim_whitespace=True, usecols=['SNP', 'A1', 'A2', 'Z'])
        df = precimed.mixer.figures.merge_z_vs_z(df1, df2)
    '''
    plot_extent = [-plot_limits, plot_limits, -plot_limits, plot_limits]
    z1name, z2name = ('Z2', 'Z1') if flip else ('Z1', 'Z2')
    z, _, _ = np.histogram2d(df[z2name], df[z1name], bins=bins, range=[[-plot_limits, plot_limits], [-plot_limits, plot_limits]])
    im=plt.imshow(np.maximum(1,z),interpolation='none', origin='lower', cmap='hot', norm=matplotlib.colors.LogNorm(), vmin=1, vmax=1e4,extent=plot_extent)
    plt.xlabel('$z_{'+traits[0]+'}$')
    plt.ylabel('$z_{'+traits[1]+'}$')
    plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))

def plot_predicted_zscore(data, num_snps, flip=False, traits=['Trait1', 'Trait2'], plot_limits=15, bins=100):
    if ('pdf' not in data):
        print('Skip generating pdf plots, data not available. Did you include --qq-plots in your "python mixer.py fit" command?')
        return
    pdf = np.array(data['pdf'])
    zgrid = np.array(data['pdf_zgrid'])
    data_limits = max(zgrid)
    data_extent = [-data_limits, data_limits, -data_limits, data_limits]
    
    # zDATA histogram is 100x100 grid, while zBGMG histogram is 1000x1000 grid.
    # Therefore counts in zBGMG are 10 times smaller, and we need to adjust for it.
    nbins_to_pdfsize_scale = len(zgrid) / bins
    plot_step = zgrid[1] - zgrid[0]
    zBGMG = (nbins_to_pdfsize_scale*nbins_to_pdfsize_scale) * pdf * plot_step * plot_step * num_snps
    
    im=plt.imshow(np.maximum(1, zBGMG.T if flip else zBGMG),
                  interpolation='none', origin='lower', cmap='hot', norm=matplotlib.colors.LogNorm(),
                  vmin=1, vmax=1e4, extent=data_extent)
    
    plt.axis([-plot_limits, plot_limits, -plot_limits, plot_limits])
    plt.xlabel('$\\hat z_{'+traits[0]+'}$')
    plt.ylabel('$\\hat z_{'+traits[1]+'}$')
    plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))

def plot_causal_density(data, flip=False, traits=['Trait1', 'Trait2'], vmax=1e3, plot_limits=0.025):
    params = data['params']
    sb1, sb2 = params['sig2_beta'][0], params['sig2_beta'][1]
    pi1, pi2, pi12 = tuple(params['pi'])
    rho = max(min(params['rho_beta'], 0.98), -0.98)
    cov = rho * np.sqrt(sb1*sb2)
    factor = 1; sb_null = 1e-7
    rv1 = multivariate_normal([0, 0], [[sb1, 0], [0, sb_null]])
    rv2 = multivariate_normal([0, 0], [[sb_null, 0], [0, sb2]])
    rv12 = multivariate_normal([0, 0], [[sb1, cov], [cov, sb2]])
    grid_step=plot_limits/50
    x, y = np.mgrid[-plot_limits:plot_limits:grid_step, -plot_limits:plot_limits:grid_step]
    pos = np.empty(x.shape + (2,))
    pos[:, :, 0] = x; pos[:, :, 1] = y
    z=factor*1e7*grid_step*grid_step*(pi1*rv1.pdf(pos)+pi2*rv2.pdf(pos)+pi12*rv12.pdf(pos))
    plot_extent = [-plot_limits, plot_limits, -plot_limits, plot_limits]
    im=plt.imshow(np.maximum(1,z if flip else z.T),interpolation='none', origin='lower', cmap='magma', norm=matplotlib.colors.LogNorm(), vmin=1, vmax=vmax,extent=plot_extent)
    plt.xlabel('$\\beta_{'+traits[0]+'}$')
    plt.ylabel('$\\beta_{'+traits[1]+'}$')
    plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))

def extract_brute1_results(data):
    brute1_results = []
    if 'optimize' in data:
        for x in data['optimize']:
            if x[0]=='brute1':
                brute1_results=x[1]
                break
    return brute1_results

def extract_likelihood_function(data):
    funcs, stats = _calculate_bivariate_uncertainty_funcs(alpha=0.05, totalhet=data['options']['totalhet'], num_snps=data['options']['num_snp'])
    p=data['params']; params=BivariateParams(pi=p['pi'],sig2_beta=p['sig2_beta'], rho_beta=p['rho_beta'],sig2_zero=p['sig2_zero'],rho_zero=p['rho_zero'])
    parametrization = BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO(lib=None, const_params1=params._params1(), const_params2=params._params2(), const_rg=params._rg(), const_rho_zero=params._rho_zero)
    brute1_results = extract_brute1_results(data)
    if not brute1_results:
        return [], []

    return [dict(funcs)['nc12@p9'](parametrization.vec_to_params([x])) for x in brute1_results['grid']], brute1_results['Jout']

def plot_likelihood(data):
    like_x, like_y = extract_likelihood_function(data)
    if (not like_x) or (not like_y):
        print('--json argument does not contain brute1 optimization results, skip likelihood plot generation')
        return
    plt.plot(np.array(like_x)/1000, like_y - np.min(like_y))
    #plt.title('cost={}'.format(data['optimize'][-1][1]['fun']))
    plt.title('-log(L) + const')
    #plt.xlabel('n12')
    #plt.ylabel('-log(L)')

def make_power_plot(data_vec, colors=None, traits=None, power_thresh=None):
    if colors is None: colors = list(range(0, len(data_vec)))
    if traits is None: traits = ['trait{}'.format(i) for i in range(1, len(data_vec) + 1)]
    leg_labels = []
    current_n = []
    current_s = []
    future_n = []
    
    cm = plt.cm.get_cmap('tab10')
    for data, color, trait in zip(data_vec, colors, traits):
        ax=plt.plot(np.log10(data['power']['nvec']), data['power']['svec'], color=cm.colors[color % 10], linestyle='solid' )
        current_n.append(data['options']['trait1_nval'])
        cs = interp1d(np.log10(data['power']['nvec']),data['power']['svec'])(np.log10(data['options']['trait1_nval']))
        current_s.append(cs)

        if power_thresh is not None:
            future_n_val = np.power(10, float(interp1d(data['power']['svec'], np.log10(data['power']['nvec']))(power_thresh)))
            future_n.append(future_n_val)
        else:
            future_n.append(None)

        display_n = lambda x: '{}'.format(int(float('{:0.1e}'.format(x))))
        display_auto = lambda x: '{}K'.format(display_n(x/1000)) if (x < 1e6) else '{:0.1f}M'.format(x/1e6)

        if 'power_ci' in data:
            if power_thresh is not None:
                future_n_ci = [np.power(10, float(interp1d(data_power['svec'], np.log10(data_power['nvec']))(power_thresh))) for data_power in data['power_ci'] if data_power]
                leg_labels.append('{} {} ({})'.format(trait, display_auto(future_n_val), display_auto(np.std(future_n_ci))))
            else:
                current_s_ci = [float(interp1d(np.log10(data_power['nvec']),data_power['svec'])(np.log10(data['options']['trait1_nval']))) for data_power in data['power_ci'] if data_power]
                leg_labels.append('{} {:.1f}% ({:.1f}%)'.format(trait, 100 * cs, 100*np.std(current_s_ci)))
        else:
            if power_thresh is not None:
                leg_labels.append('{} ({})'.format(trait, display_auto(future_n_val)))
            else:
                leg_labels.append('{} ({:.1f}%)'.format(trait, 100 * cs))

    if power_thresh:
        plt.hlines([float(power_thresh)], 4, 8, linestyles='--')

    for cn, cs, fn, trait in zip(current_n, current_s, future_n, traits):
        plt.plot([np.log10(cn)], [cs], '*k', markersize=8)
        if power_thresh is not None:
            plt.plot([np.log10(fn)], [float(power_thresh)], '.k--', markersize=8)

    leg_labels.append('Current N')
    if power_thresh is not None:
        leg_labels.append('N at {}%'.format(int(100*float(power_thresh))))
    
    plt.legend(leg_labels, loc='lower right',frameon=False, numpoints=1)
    plt.xlabel('Sample size')
    plt.ylabel('Estimated percent variance explained\nby genome-wide significant SNPs')
    plt.xlim([4, 8])
    plt.ylim([0, 1])
    plt.locator_params(axis='x', nbins=5)
    plt.gca().set_xticklabels(labels=['10K', '100K', '1M', '10M', '100M'])

# https://stackoverflow.com/questions/27433316/how-to-get-argparse-to-read-arguments-from-a-file-with-an-option-rather-than-pre
class LoadFromFile (argparse.Action):
    def __call__ (self, parser, namespace, values, option_string=None):
        with values as f:
            contents = f.read()

        data = parser.parse_args(contents.split(), namespace=namespace)
        for k, v in vars(data).items():
            if v and k != option_string.lstrip('-'):
                setattr(namespace, k, v)

def parser_one_add_arguments(args, func, parser):
    parser.add_argument('--json', type=str, default=[""], nargs='+', help="json file from a univariate analysis. This argument does support wildcards (*) or a list with multiple space-separated arguments to process more than one .json file. This allows to generate a combined .csv table across many traits.")
    parser.add_argument('--trait1', type=str, default=[], nargs='+', help="name of the first trait")
    parser.add_argument('--power-thresh', type=str, default=None, help="threshold for power analysis, e.g. 0.9 or 0.5, to estimate corresponding N")
    parser.add_argument('--power-figsize', type=float, nargs='+', default=[], help="figure size for power plots")
    parser.set_defaults(func=func)

def parser_two_add_arguments(args, func, parser):
    parser.add_argument('--json', type=str, default=[], nargs='*', help="json file from a bivariate analysis, i.e. either 'mixer.py fit2' or 'mixer.py test2' step. This argument does support wildcards (*) to process multiple .json files (this allows to generate a combined .csv table across many cross-trait combinations, but it doesn't generate figures; to generate figures, use --json on a single file, or alternatively use --json-fit and --json-test. ")
    parser.add_argument('--json-fit', type=str, default="", help="json file from a bivariate analysis with 'mixer.py fit2' step. This argument does NOT support wildcards. Using --json-fit in conjunction with --json-test produces figures that contain both log-likelihood plots (based on fit2 results) and QQ plot (based on test2 results). When use --json-fit and --json-test, there is no need to specify --json argument. ")
    parser.add_argument('--json-test', type=str, default="", help="json file from a bivariate analysis with 'mixer.py test2' step). This argument does NOT support wildcards.")
    parser.add_argument('--trait1', type=str, default="trait1", help="name of the first trait")    
    parser.add_argument('--trait2', type=str, default="trait2", help="name of the second trait")    
    parser.add_argument('--trait1-file', type=str, default=None, help="summary statistics file for the first trait (optional parameter; use it only if you need to generate bivariate z-vs-z density plots)")    
    parser.add_argument('--trait2-file', type=str, default=None, help="summary statistics file for the second trait (optional parameter; see comment for --trait1-file")    
    parser.add_argument('--trait1-color', type=int, default=0, choices=list(range(9)), help="color for the venn diagram (first trait); 0-8, encoded as tab10 color palette (https://matplotlib.org/3.1.1/tutorials/colors/colormaps.html) excluding grey code which is reserved the polygenic overlap")    
    parser.add_argument('--trait2-color', type=int, default=1, choices=list(range(9)), help="color for the venn diagram (second trait)")    
    parser.add_argument('--flip', default=False, action="store_true", help="flip venn diagram and stratified QQ plots. Note that this arguments does not apply to --trait1 and --trait2 arguments, not to --trait1-color and --trait2-color.")
    parser.set_defaults(func=func)

def parse_args(args):
    parser = argparse.ArgumentParser(description="MiXeR visualization tools.")

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--argsfile', type=open, action=LoadFromFile, default=None, help="file with additional command-line arguments")
    parent_parser.add_argument("--out", type=str, default="mixer", help="prefix for the output files")
    parent_parser.add_argument('--ext', type=str, default=['png'], nargs='+', choices=['png', 'svg'], help="output extentions")
    parent_parser.add_argument('--zmax', type=float, default=10, help="limit for z-vs-z density plots")

    subparsers = parser.add_subparsers()
    parser_one_add_arguments(args=args, func=execute_one_parser, parser=subparsers.add_parser("one", parents=[parent_parser], help='produce figures for univariate analysis'))
    parser_two_add_arguments(args=args, func=execute_two_parser, parser=subparsers.add_parser("two", parents=[parent_parser], help='produce figures for cross-trait analysis'))

    return parser.parse_args(args)

def execute_two_parser(args):
    df_data = {}

    if len(args.json) == 0: args.json = [x for x in [args.json_fit, args.json_test] if x]

    files = glob.glob(args.json[0]) if (len(args.json) == 1) else args.json
    if len(files) == 0: raise(ValueError('no files detected, check --json {}'.format(args.json)))
    print('generate {}.csv from {} json files...'.format(args.out, len(files)))

    for fname in files:
        keys = 'dice pi1 pi2 pi12 nc1@p9 nc2@p9 nc12@p9 rho_zero rho_beta rg'.split()
        try:
            data = json.loads(open(fname).read())
            if 'dice' not in data['ci']:
                data['ci']['dice'] = {'point_estimate' : 2 * data['ci']['pi12']['point_estimate'] / (data['ci']['pi1u']['point_estimate'] + data['ci']['pi2u']['point_estimate'])}

            trait1 = os.path.basename(data['options']['trait1_file']).replace('.sumstats.gz', '')
            trait2 = os.path.basename(data['options']['trait2_file']).replace('.sumstats.gz', '')
            for k in keys:   # test that all keys are available
                val = data['ci'][k]['point_estimate']
        except:
            print('error reading from {}, skip'.format(fname))
            continue

        insert_key_to_dictionary_as_list(df_data, 'fname', fname)
        insert_key_to_dictionary_as_list(df_data, 'trait1', trait1)
        insert_key_to_dictionary_as_list(df_data, 'trait2', trait2)
        for k in keys: insert_key_to_dictionary_as_list(df_data, k, data['ci'][k]['point_estimate'])

        rho_beta = data['ci']['rho_beta']['point_estimate']
        fraction_concordant_within_shared = 2 * multivariate_normal([0, 0], [[1, rho_beta], [rho_beta, 1]]).cdf([0, 0])
        insert_key_to_dictionary_as_list(df_data, 'fraction_concordant_within_shared', fraction_concordant_within_shared)

        brute1_results = extract_brute1_results(data)
        if brute1_results:
            min_overlap = brute1_results['Jout'][0]
            max_overlap = brute1_results['Jout'][-1]
            best_cost = data['optimize'][-1][1]['fun']
            cost_n = data['options']['sum_weights']
            df_diff = -1  # fitting polygenic overlap require 1 extra parameter
            insert_key_to_dictionary_as_list(df_data, 'best_vs_min_AIC',              2 * df_diff + 2 * (min_overlap - best_cost))
            insert_key_to_dictionary_as_list(df_data, 'best_vs_min_BIC', np.log(cost_n) * df_diff + 2 * (min_overlap - best_cost))
            insert_key_to_dictionary_as_list(df_data, 'best_vs_max_AIC',              2 * df_diff + 2 * (max_overlap - best_cost))
            insert_key_to_dictionary_as_list(df_data, 'best_vs_max_BIC', np.log(cost_n) * df_diff + 2 * (max_overlap - best_cost))
        else:
            insert_key_to_dictionary_as_list(df_data, 'best_vs_min_AIC', None)
            insert_key_to_dictionary_as_list(df_data, 'best_vs_min_BIC', None)
            insert_key_to_dictionary_as_list(df_data, 'best_vs_max_AIC', None)
            insert_key_to_dictionary_as_list(df_data, 'best_vs_max_BIC', None)

    pd.DataFrame(df_data).to_csv(args.out+'.csv', index=False, sep='\t')
    print('Done.')

    if args.json_fit and args.json_test:
        data_fit = json.loads(open(args.json_fit).read())
        data_test = json.loads(open(args.json_test).read())
    elif len(files) == 1:
        data_fit = json.loads(open(files[0]).read())
        data_test = data_fit
    else:
        print('--json argument lists multiple files or is a wild-card (contains *), skip figures generation')
        return

    if 'qqplot' not in data_test:
        print('Skip generating stratified QQ plots, data not available.')

    if args.trait1_file and args.trait2_file:
        df1 = pd.read_table(args.trait1_file, delim_whitespace=True, usecols=['SNP', 'A1', 'A2', 'Z'])
        df2 = pd.read_table(args.trait2_file, delim_whitespace=True, usecols=['SNP', 'A1', 'A2', 'Z'])
        df = merge_z_vs_z(df1, df2)

        plt.figure(figsize=[12, 6])
        plt.subplot(2,4,1); make_venn_plot(data_fit, flip=args.flip, traits=[args.trait1, args.trait2], colors=[args.trait1_color, args.trait2_color])
        if 'qqplot' in data_test:
            plt.subplot(2,4,2); make_strat_qq_plots(data_test, flip=args.flip, traits=[args.trait1, args.trait2], do_legend=False)
            plt.subplot(2,4,3); make_strat_qq_plots(data_test, flip=(not args.flip), traits=[args.trait2, args.trait1], do_legend=True)
        plt.subplot(2,4,4); plot_likelihood(data_fit)
        plt.subplot(2,4,5); plot_causal_density(data_test, flip=args.flip, traits=[args.trait1, args.trait2])
        plt.subplot(2,4,6); plot_z_vs_z_data(df, flip=args.flip, plot_limits=args.zmax, traits=[args.trait1, args.trait2])
        plt.subplot(2,4,7); plot_predicted_zscore(data_test, len(df), flip=args.flip, plot_limits=args.zmax, traits=[args.trait1, args.trait2])
    else:
        plt.figure(figsize=[12, 3])
        plt.subplot(1,4,1); make_venn_plot(data_fit, flip=args.flip, traits=[args.trait1, args.trait2], colors=[args.trait1_color, args.trait2_color])
        if 'qqplot' in data_test:
            plt.subplot(1,4,2); make_strat_qq_plots(data_test, flip=args.flip, traits=[args.trait1, args.trait2], do_legend=False)
            plt.subplot(1,4,3); make_strat_qq_plots(data_test, flip=(not args.flip), traits=[args.trait2, args.trait1], do_legend=True)
        plt.subplot(1,4,4); plot_likelihood(data_fit)
        
    for ext in args.ext:
        plt.savefig(args.out + '.' + ext, bbox_inches='tight')
        print('Generated ' + args.out + '.' + ext)

def insert_key_to_dictionary_as_list(df_data, key, value):
    if key not in df_data:
        df_data[key] = []
    df_data[key].append(value)

def execute_one_parser(args):
    df_data = {}
    files = glob.glob(args.json[0]) if (len(args.json) == 1) else args.json
    if len(files) == 0: raise(ValueError('no files detected, check --json {}'.format(args.json)))
    print('generate {}.csv from {} json files...'.format(args.out, len(files)))
    for fname in files:
        keys = 'pi sig2_beta sig2_zero h2 nc@p9'.split()
        try:
            data = json.loads(open(fname).read())
            for k in keys:   # test that all keys are available
                val = data['ci'][k]['point_estimate']
                val = data['inft_optimize'][-1][1]['AIC'] 
                val = data['optimize'][-1][1]['AIC']
                val = data['inft_optimize'][-1][1]['BIC'] 
                val = data['optimize'][-1][1]['BIC']
        except:
            print('error reading from {}, skip'.format(fname))
            continue

        insert_key_to_dictionary_as_list(df_data, 'fname', fname)
        for k in keys:
            insert_key_to_dictionary_as_list(df_data, k, data['ci'][k]['point_estimate'])

        aic_diff = data['inft_optimize'][-1][1]['AIC'] - data['optimize'][-1][1]['AIC']
        bic_diff = data['inft_optimize'][-1][1]['BIC'] - data['optimize'][-1][1]['BIC']
        insert_key_to_dictionary_as_list(df_data, 'AIC', aic_diff)
        insert_key_to_dictionary_as_list(df_data, 'BIC', bic_diff)
        
    pd.DataFrame(df_data).to_csv(args.out+'.csv', index=False, sep='\t')
    print('Done.')

    data_list = []
    traits_list = []
    for fname, traitname in zip(files, files if (len(args.trait1) == 0) else args.trait1):
        try:
            data = json.loads(open(fname).read())
            if 'power' not in data: continue
            data_list.append(data)
            traits_list.append(traitname)
        except:
            print('error reading from {}, skip'.format(fname))
            continue

    if len(data_list) > 0:
        plt.figure(figsize=(tuple(args.power_figsize) if (len(args.power_figsize) > 0) else None))
        make_power_plot(data_list, traits=traits_list, power_thresh=args.power_thresh)
        for ext in args.ext:
            plt.savefig(args.out + '.power.' + ext, bbox_inches='tight')
            print('Generated ' + args.out + '.power.' + ext)
    else:
        print('Skip generating power plots, data not available. Did you include --power-curve in your "python mixer.py fit" command?')


    if len(files) > 1:
        print('--json argument is a wild-card (contains *), skip figures generation')        
        return

    data = json.loads(open(args.json[0]).read())
    if 'qqplot' in data:
        plt.figure()
        make_qq_plot(data['qqplot'], ci=True)
        for ext in args.ext:
            plt.savefig(args.out + '.qq.' + ext, bbox_inches='tight')
            print('Generated ' + args.out + '.qq.' + ext)
    else:
        print('Skip generating QQ plots, data not available. Did you include --qq-plots in your "python mixer.py fit" command?')

    if 'qqplot_bins' in data:
        plt.figure(figsize=[12, 12])
        for i in range(0, 3):
            for j in range(0, 3):
                plt.subplot(3,3,i*3+j+1)
                make_qq_plot(data['qqplot_bins'][i*3+j])
                plt.title(data['qqplot_bins'][i*3+j]['title'].replace(';', '\n'))
        for ext in args.ext:
            plt.savefig(args.out + '.qqbin.' + ext, bbox_inches='tight')
            print('Generated ' + args.out + '.qqbin.' + ext)            
    else:
        print('Skip generating qq plots (3x3 bins of MAF and LD score), data not available. Did you include --qq-plots in your "python mixer.py fit" command?')
  
