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
from matplotlib_venn import venn2_circles

from scipy.interpolate import interp1d
from scipy.stats import multivariate_normal

from .utils import BivariateParametrization_constUNIVARIATE_constRG_constRHOZERO
from .utils import _calculate_bivariate_uncertainty_funcs
from .utils import BivariateParams
from .utils import NumpyEncoder

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

def make_venn_plot(data, flip=False, factor='K', traits=['Trait1', 'Trait2'], colors=[0, 1], max_size=None, formatter=None, statistic=["point_estimate"], plot_rg=True):
    cm = plt.cm.get_cmap('tab10')

    if factor=='K': scale_factor=1000
    elif factor=='': scale_factor=1
    else: raise(ValueError('Unknow factor: {}'.format(factor)))

    with_std = (len(statistic) > 1)
    n1 = data['ci']['nc1@p9'][statistic[0]]/scale_factor; n1_se = data['ci']['nc1@p9'][statistic[1]]/scale_factor if with_std else None
    n2 = data['ci']['nc2@p9'][statistic[0]]/scale_factor; n2_se = data['ci']['nc2@p9'][statistic[1]]/scale_factor if with_std else None
    n12 = data['ci']['nc12@p9'][statistic[0]]/scale_factor; n12_se = data['ci']['nc12@p9'][statistic[1]]/scale_factor if with_std else None
    rg = data['ci']['rg'][statistic[0]]

    if max_size is None: max_size = n1+n2+n12
    if flip: n1, n2 = n2, n1; n1_se, n2_se = n2_se, n1_se
    f = lambda x: x if x < 7 else x+1

    v = venn2(subsets = (n1, n2, n12), normalize_to=(n1+n2+n12)/max_size, set_labels = ("", ""))
    v.get_patch_by_id('100').set_color(cm.colors[f(colors[0])])
    v.get_patch_by_id('010').set_color(cm.colors[f(colors[1])])
    v.get_patch_by_id('110').set_color(cm.colors[7])   
    c=venn2_circles(subsets = (n1, n2, n12), normalize_to=(n1+n2+n12)/max_size, linewidth=1.5, color="white")
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
    plt.title(traits[0] +' & ' + newline + traits[1], y=(-0.18 if plot_rg else 0))

    if plot_rg:
        clr = plt.cm.get_cmap('seismic')((rg+1)/2)
        plt.gca().add_patch(patches.Rectangle(((-abs(0.7*rg) if (rg < 0) else 0) , -0.7), abs(0.7 * rg), 0.15, fill=True, clip_on=False, color=clr))
        plt.gca().add_patch(patches.Rectangle((-0.70, -0.7), 1.4, 0.15, fill=False, clip_on=False))
        plt.gca().add_patch(patches.Rectangle((0, -0.7), 0, 0.15, fill=False, clip_on=False, linewidth=3))
        plt.gca().text(-0.35 if (rg>0) else 0.35, -0.7+0.13/2, '$r_g$={:.2f}'.format(rg), fontsize=11, horizontalalignment='center',         verticalalignment='center')

def make_strat_qq_plots(data, flip=False, traits=['Trait1', 'Trait2'], do_legend=True):
    cm = plt.cm.get_cmap('tab10')
    for i in np.array(range(0, 4)) + (4 if flip else 0):
        hData = plt.plot(data['qqplot'][i]['data_logpvec'], data['qqplot'][i]['hv_logp'], color=cm.colors[i % 4], linestyle='solid')
    hNull = plt.plot(data['qqplot'][i]['hv_logp'], data['qqplot'][i]['hv_logp'], 'k--')
    if do_legend: plt.legend(['All SNPs'] + '$P{0}\leq0.1$ $P{0}\leq0.01$ $P{0}\leq0.001$'.format('_{' + traits[1]+'}').split(), loc='lower right', fontsize=10, borderpad=0.2, frameon=False, borderaxespad=0.2, labelspacing=0.1)
    for i in np.array(range(0, 4)) + (4 if flip else 0):
        hModel = plt.plot(data['qqplot'][i]['model_logpvec'], data['qqplot'][i]['hv_logp'], color=cm.colors[i % 4], linestyle='dashed')
    plt.ylim(0, 7.3); plt.xlim(0, 7.3); 
    plt.title('{} | {}'.format(traits[0], traits[1]), fontsize=12)
    plt.xticks(np.arange(0, 7.3, step=1),fontsize=10)
    plt.yticks(np.arange(0, 7.3, step=1),fontsize=10)
    plt.xlabel('Expected -$log_{{10}}(p_{{{}}})$'.format(traits[0]), fontsize=12)
    plt.ylabel('Observed -$log_{{10}}(p_{{{}}})$'.format(traits[0]), fontsize=12)

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
    plt.xlabel('$z_{'+traits[0]+'}$', fontsize=12)
    plt.ylabel('$z_{'+traits[1]+'}$', fontsize=12, labelpad=-0.1)
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
    plt.xlabel('$\\hat z_{'+traits[0]+'}$', fontsize=12)
    plt.ylabel('$\\hat z_{'+traits[1]+'}$', fontsize=12, labelpad=-0.1)
    plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))

def plot_causal_density(data, flip=False, traits=['Trait1', 'Trait2'], vmax=1e3, plot_limits=0.025, statistic=['point_estimate']):
    if statistic[0] not in ['point_estimate', 'mean', 'median']:
        print('Unable to make plot_causal_density() for statistic=={}'.format(statistic[0]))
    sb1 = data['ci']['sig2_beta_T1'][statistic[0]]
    sb2 = data['ci']['sig2_beta_T2'][statistic[0]]
    pi1 = data['ci']['pi1'][statistic[0]]
    pi2 = data['ci']['pi2'][statistic[0]]
    pi12 = data['ci']['pi12'][statistic[0]]
    rho = max(min(data['ci']['rho_beta'][statistic[0]], 0.98), -0.98)
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
    plt.xlabel('$\\beta_{'+traits[0]+'}$', fontsize=12)
    plt.ylabel('$\\beta_{'+traits[1]+'}$', fontsize=12, labelpad=0.1)
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
    if 'likelihood' in data:
        cm = plt.cm.get_cmap('tab10')
        like_x = np.array([like_x for like_x, like_y in data['likelihood']])
        like_y = np.array([like_y for like_x, like_y in data['likelihood']])
        like_x_min = np.min(like_x)
        like_x_max = np.max(like_x)
        plot_x = np.arange(like_x_min, like_x_max, (like_x_max-like_x_min) / 100)
        plot_y = np.zeros((like_y.shape[0], len(plot_x)))
        for i in range(like_x.shape[0]):
            plot_y[i, :] = interp1d(like_x[i, :], like_y[i, :], bounds_error=False)(plot_x)
            plot_y[i, :] = plot_y[i, :] - np.nanmin(plot_y[i, :])
        plot_y_def = np.sum(np.isfinite(plot_y), 0)
        plot_y = np.nanmean(plot_y, 0)
        plot_y[plot_y_def < 3] = np.nan
        plt.plot(plot_x / 1000, plot_y)
        for like_x, like_y in data['likelihood']:
            plt.plot(np.array(like_x)/1000, like_y - np.min(like_y), linestyle='dotted', color=cm.colors[0], alpha=0.3)
    else:        
        like_x, like_y = extract_likelihood_function(data)
        if (not like_x) or (not like_y):
            print('--json argument does not contain brute1 optimization results, skip likelihood plot generation')
            return
        plt.plot(np.array(like_x)/1000, like_y - np.min(like_y))
        plt.title('-log(L) + const')
    plt.xlabel('Shared variant number [k]',fontsize=12)
    plt.ylabel('-log(L) + const',fontsize=12)
    plt.title('Log-likelihood',fontsize=13)

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

        print('HAS POWER?', 'power_ci' in data)
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
    plt.xlabel('Sample size (N)', fontsize=11)
    plt.ylabel('Estimated variance (%) explained\nby genome-wide significant SNPs', fontsize=11)
    plt.xlim([4, 8])
    plt.ylim([-0.017, 1])
    plt.locator_params(axis='x', nbins=5)
    plt.gca().set_xticklabels(labels=['10K', '100K', '1M', '10M', '100M'])
    plt.yticks(np.arange(0, 1.01, step=0.2))
    plt.axes().set_yticklabels(labels=['0', '20', '40', '60', '80', '100'])

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

def parser_combine_add_arguments(args, func, parser):
    parser.add_argument('--json', type=str, default=None, help="Path to json files from mixer runs. Must not contain wildcards. Must contain '@' sign, indicating the location of the repeat index.")
    parser.add_argument('--rep2use', type=str, default='1-20', help="Repeat indices to use, e.g. 1,2,3 or 1-4,12,16-20")
    parser.set_defaults(func=func)

def parse_args(args):
    parser = argparse.ArgumentParser(description="MiXeR visualization tools.")

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--argsfile', type=open, action=LoadFromFile, default=None, help="file with additional command-line arguments")
    parent_parser.add_argument("--out", type=str, default="mixer", help="prefix for the output files")
    parent_parser.add_argument('--ext', type=str, default=['png'], nargs='+', choices=['png', 'svg'], help="output extentions")
    parent_parser.add_argument('--zmax', type=float, default=10, help="limit for z-vs-z density plots")
    parent_parser.add_argument('--statistic', type=str, nargs='+', default=["point_estimate"], choices=["point_estimate", "mean", "median", "std", "min", "max"], help="Which statistic to show in the tables and on the Venn diagrams. Can have multiple values. In the case of venn diagram, the first value (typically 'point_estimate' or 'mean') indicate the size of the venn diagram; the second value (optional, typically 'std') allow to include error bars on the Venn diagramm.")

    subparsers = parser.add_subparsers()
    parser_one_add_arguments(args=args, func=execute_one_parser, parser=subparsers.add_parser("one", parents=[parent_parser], help='produce figures for univariate analysis'))
    parser_two_add_arguments(args=args, func=execute_two_parser, parser=subparsers.add_parser("two", parents=[parent_parser], help='produce figures for cross-trait analysis'))
    parser_combine_add_arguments(args=args, func=execute_combine_parser, parser=subparsers.add_parser("combine", parents=[parent_parser], help='combine .json files MiXeR runs (e.g. with different --extract setting)'))

    return parser.parse_args(args)

def execute_two_parser(args):
    df_data = {}

    if len(args.json) == 0: args.json = [x for x in [args.json_fit, args.json_test] if x]

    files = glob.glob(args.json[0]) if (len(args.json) == 1) else args.json
    if len(files) == 0: raise(ValueError('no files detected, check --json {}'.format(args.json)))
    print('generate {}.csv from {} json files...'.format(args.out, len(files)))

    for fname in files:
        keys = 'dice pi1 pi2 pi12 nc1@p9 nc2@p9 nc12@p9 rho_zero rho_beta rg fraction_concordant_within_shared'.split()
        try:
            data = json.loads(open(fname).read())
            if 'dice' not in data['ci']:
                data['ci']['dice'] = {'point_estimate' : 2 * data['ci']['pi12']['point_estimate'] / (data['ci']['pi1u']['point_estimate'] + data['ci']['pi2u']['point_estimate'])}
            if 'fraction_concordant_within_shared' not in data['ci']:
                rho_beta = data['ci']['rho_beta']['point_estimate']
                data['ci']['fraction_concordant_within_shared'] = {'point_estimate' : 2 * multivariate_normal([0, 0], [[1, rho_beta], [rho_beta, 1]]).cdf([0, 0])}

            trait1 = os.path.basename(data['options']['trait1_file']).replace('.sumstats.gz', '')
            trait2 = os.path.basename(data['options']['trait2_file']).replace('.sumstats.gz', '')
            for k in keys:   # test that all keys are available
                for stat in args.statistic:
                    val = data['ci'][k][stat]
        except:
            print('error reading from {}, skip'.format(fname))
            continue

        insert_key_to_dictionary_as_list(df_data, 'fname', fname)
        insert_key_to_dictionary_as_list(df_data, 'trait1', trait1)
        insert_key_to_dictionary_as_list(df_data, 'trait2', trait2)
        for k in keys:
            for stat in args.statistic:
                insert_key_to_dictionary_as_list(df_data, k if (stat=="point_estimate") else "{} ({})".format(k, stat), data['ci'][k][stat])

        brute1_results = extract_brute1_results(data)
        mskeys = 'best_vs_min_AIC best_vs_min_BIC best_vs_max_AIC best_vs_max_BIC'.split()
        if 'modelselection' in data:
            for mskey in mskeys: insert_key_to_dictionary_as_list(df_data, mskey, data['modelselection'][mskey])
        elif brute1_results:
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
            for mskey in mskeys: insert_key_to_dictionary_as_list(df_data, mskey, None)

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
        plt.figure(figsize=[12, 5.5])
        plt.subplot(2,4,1); make_venn_plot(data_fit, flip=args.flip, traits=[args.trait1, args.trait2], colors=[args.trait1_color, args.trait2_color], statistic=args.statistic)
        if 'qqplot' in data_test:
            plt.subplot(2,4,2); make_strat_qq_plots(data_test, flip=args.flip, traits=[args.trait1, args.trait2], do_legend=True)
            plt.subplot(2,4,3); make_strat_qq_plots(data_test, flip=(not args.flip), traits=[args.trait2, args.trait1], do_legend=True)
        plt.subplot(2,4,4); plot_likelihood(data_fit)
        plt.subplot(2,4,6); plot_causal_density(data_test, flip=args.flip, traits=[args.trait1, args.trait2], statistic=args.statistic)
        df1 = pd.read_table(args.trait1_file, delim_whitespace=True, usecols=['SNP', 'A1', 'A2', 'Z'])
        df2 = pd.read_table(args.trait2_file, delim_whitespace=True, usecols=['SNP', 'A1', 'A2', 'Z'])
        df = merge_z_vs_z(df1, df2)
        plt.subplot(2,4,7); plot_z_vs_z_data(df, flip=args.flip, plot_limits=args.zmax, traits=[args.trait1, args.trait2])
        plt.subplot(2,4,8); plot_predicted_zscore(data_test, len(df), flip=args.flip, plot_limits=args.zmax, traits=[args.trait1, args.trait2])
        plt.tight_layout(pad=1.5)
    else:
        plt.figure(figsize=[12, 3.5])
        plt.subplot(1,4,1); make_venn_plot(data_fit, flip=args.flip, traits=[args.trait1, args.trait2], colors=[args.trait1_color, args.trait2_color], statistic=args.statistic)
        if 'qqplot' in data_test:
            plt.subplot(1,4,2); make_strat_qq_plots(data_test, flip=args.flip, traits=[args.trait1, args.trait2], do_legend=True)
            plt.subplot(1,4,3); make_strat_qq_plots(data_test, flip=(not args.flip), traits=[args.trait2, args.trait1], do_legend=True)
        plt.subplot(1,4,4); plot_likelihood(data_fit)
        plt.tight_layout(pad=1.5)
        
    for ext in args.ext:
        plt.savefig(args.out + '.' + ext, bbox_inches='tight')
        print('Generated ' + args.out + '.' + ext)

def insert_key_to_dictionary_as_list(df_data, key, value):
    if key not in df_data:
        df_data[key] = []
    df_data[key].append(value)

def parse_val2use(val2use_arg):
    val2use = []
    for a in val2use_arg.split(","):
        if "-" in a:
            start, end = [int(x) for x in a.split("-")]
            val2use += [str(x) for x in range(start, end+1)]
        else:
            val2use.append(a.strip())
    if np.any([not x.isdigit() for x in val2use]): raise ValueError('Value labels must be integer: {}'.format(val2use_arg))
    return val2use

'''
python ~/github/mixer/precimed/mixer_figures.py combine --json PGC_SCZ_2014_EUR.fit.rep@.json  --out combined/PGC_SCZ_2014_EUR.fit
python ~/github/mixer/precimed/mixer_figures.py one --json combined/PGC_SCZ_2014_EUR.fit.json  --out combined/PGC_SCZ_2014_EUR.fit

python ~/github/mixer/precimed/mixer_figures.py combine --json PGC_SCZ_2014_EUR_vs_PGC_BIP_2016.fit.rep@.json  --out combined/PGC_SCZ_2014_EUR_vs_PGC_BIP_2016.fit
python ~/github/mixer/precimed/mixer_figures.py combine --json PGC_SCZ_2014_EUR_vs_PGC_BIP_2016.test.rep@.json  --out combined/PGC_SCZ_2014_EUR_vs_PGC_BIP_2016.test
python ~/github/mixer/precimed/mixer_figures.py two --json-fit combined/PGC_SCZ_2014_EUR_vs_PGC_BIP_2016.fit.json --json-test combined/PGC_SCZ_2014_EUR_vs_PGC_BIP_2016.test.json --out combined/PGC_SCZ_2014_EUR_vs_PGC_BIP_2016 --statistic mean std
'''

def execute_combine_parser(args):
    args.rep2use = parse_val2use(args.rep2use)
    failed_indices = []
    data_vec = []
    for rep in args.rep2use:
        try:
            data = json.loads(open(args.json.replace('@', str(rep)), 'r').read())
            if data:
                data_vec.append(data)
        except:
            failed_indices.append(rep)
    if failed_indices: print('WARNING: {}: results for {} runs are missing (rep {})'.format(args.out, len(failed_indices), ' '.join(failed_indices)))
    
    results = {'ci':{}, 'options':{}}

    for key in ['totalhet', 'num_snp']:
        values = [data['options'][key] for data in data_vec]
        results['options'][key] = np.mean(values)   

    for key in ['trait1_file', 'trait2_file']:
        values = [data['options'][key] for data in data_vec if (key in data['options'])]
        if not values: continue
        if len(set(values)) > 1: raise(ValueError('Input files have distinct value in "{}" field: {}'.format(key, ' '.join(set(values)))))
        results['options'][key] = values[0]

    for key in ['trait1_nval', 'trait2_nval']:
        values = [data['options'][key] for data in data_vec if (key in data['options'])]
        if not values: continue
        results['options'][key] = np.mean(values)

    values = [data['analysis'] for data in data_vec]
    if len(set(values)) > 1: raise(ValueError('Input files have distinct value in "analysis" field: {}'.format(' '.join(set(values)))))
    results['analysis'] = values[0]

    univariate_keys = ['pi', 'nc', 'nc@p9', 'sig2_beta', 'sig2_zero', 'h2']
    bivariate_keys = ['sig2_zero_T1', 'sig2_zero_T2', 'sig2_beta_T1', 'sig2_beta_T2', 'h2_T1', 'h2_T2', 'rho_zero', 'rho_beta', 'rg', 'pi1', 'pi2', 'pi12', 'pi1u', 'pi2u', 'dice', 'nc1', 'nc2', 'nc12', 'nc1u', 'nc2u', 'nc1@p9', 'nc2@p9', 'nc12@p9', 'nc1u@p9', 'nc2u@p9', 'totalpi', 'totalnc', 'totalnc@p9', 'pi1_over_totalpi', 'pi2_over_totalpi', 'pi12_over_totalpi', 'pi1_over_pi1u', 'pi2_over_pi2u', 'pi12_over_pi1u', 'pi12_over_pi2u', 'pi1u_over_pi2u', 'pi2u_over_pi1u']
    for key in (univariate_keys + bivariate_keys):
        values = [data['ci'][key]['point_estimate'] for data in data_vec if (key in data['ci'])]
        if values: results['ci'][key] = {'mean': np.mean(values), 'median':np.median(values), 'std': np.std(values), 'min': np.min(values), 'max': np.max(values)}
        if values and (key=='rho_beta'):
            values = [2 * multivariate_normal([0, 0], [[1, rho_beta], [rho_beta, 1]]).cdf([0, 0]) for rho_beta in values]
            results['ci']['fraction_concordant_within_shared'] = {'mean': np.mean(values), 'median':np.median(values), 'std': np.std(values), 'min': np.min(values), 'max': np.max(values)}

    if results['analysis'] == 'bivariate':
        for data in data_vec:
            like_x, like_y = extract_likelihood_function(data)
            if (not like_x) or (not like_y): continue
            if 'likelihood' not in results: results['likelihood'] = []
            results['likelihood'].append((like_x, like_y))

    combine_qqplots = lambda qqplots: {
        'hv_logp':np.mean(np.array([np.array(qq['hv_logp']).astype(float) for qq in qqplots]), 0),
        'data_logpvec':np.mean(np.array([np.array(qq['data_logpvec']).astype(float) for qq in qqplots]), 0),
        'model_logpvec':np.mean(np.array([np.array(qq['model_logpvec']).astype(float) for qq in qqplots]), 0),
        'sum_data_weights':np.mean(np.array([np.array(qq['sum_data_weights']).astype(float) for qq in qqplots]), 0)
    }

    qqplots = [data['qqplot'] for data in data_vec if ('qqplot' in data)]
    if qqplots:
        if results['analysis'] == 'univariate':
            results['qqplot'] = combine_qqplots(qqplots)
        elif results['analysis'] == 'bivariate':
            num_plots = len(qqplots[0])
            results['qqplot'] = []
            for index in range(num_plots):
                results['qqplot'].append(combine_qqplots([qq[index] for qq in qqplots]))

    if results['analysis'] == 'bivariate':
        best_vs_min_AIC = [];         best_vs_min_BIC = []
        best_vs_max_AIC = [];         best_vs_max_BIC = []
        for data in data_vec:
            brute1_results = extract_brute1_results(data)
            if not brute1_results: continue
            min_overlap = brute1_results['Jout'][0]
            max_overlap = brute1_results['Jout'][-1]
            best_cost = data['optimize'][-1][1]['fun']
            cost_n = data['options']['sum_weights']
            df_diff = -1  # fitting polygenic overlap require 1 extra parameter
            best_vs_min_AIC.append(              2 * df_diff + 2 * (min_overlap - best_cost))
            best_vs_min_BIC.append( np.log(cost_n) * df_diff + 2 * (min_overlap - best_cost))
            best_vs_max_AIC.append(              2 * df_diff + 2 * (max_overlap - best_cost))
            best_vs_max_BIC.append( np.log(cost_n) * df_diff + 2 * (max_overlap - best_cost))
        results['modelselection'] = {
            'best_vs_min_AIC' : np.mean(best_vs_min_AIC) if best_vs_min_AIC else None, 
            'best_vs_min_BIC' : np.mean(best_vs_min_BIC) if best_vs_min_BIC else None, 
            'best_vs_max_AIC' : np.mean(best_vs_max_AIC) if best_vs_max_AIC else None, 
            'best_vs_max_BIC' : np.mean(best_vs_max_BIC) if best_vs_max_BIC else None
        }

    if results['analysis'] == 'univariate':
        AIC = []; BIC = []
        for data in data_vec:        
            try:
                aic_diff = data['inft_optimize'][-1][1]['AIC'] - data['optimize'][-1][1]['AIC']
                bic_diff = data['inft_optimize'][-1][1]['BIC'] - data['optimize'][-1][1]['BIC']
            except:
                continue
            AIC.append(aic_diff)
            BIC.append(bic_diff)
        results['modelselection'] = {
            'mixture_vs_inft_AIC' : np.mean(AIC) if AIC else None, 
            'mixture_vs_inft_BIC' : np.mean(BIC) if BIC else None, 
        }

    data_pdf_vec = [np.array(data['pdf']) for data in data_vec if ('pdf' in data)]
    if data_pdf_vec: results['pdf'] = np.mean(np.array(data_pdf_vec), 0)
    data_pdf_zgrid_vec = [np.array(data['pdf_zgrid']) for data in data_vec if ('pdf_zgrid' in data)]
    if data_pdf_zgrid_vec: results['pdf_zgrid'] = np.mean(np.array(data_pdf_zgrid_vec), 0)

    data_power = [data['power'] for data in data_vec if ('power' in data)]
    if data_power:
        results['power_ci'] = data_power
        results['power'] = {}
        results['power']['svec'] = np.mean(np.array([np.array(power['svec']).astype(float) for power in data_power]), 0)
        results['power']['nvec'] = np.mean(np.array([np.array(power['nvec']).astype(float) for power in data_power]), 0)

    with open(args.out + '.json', 'w') as outfile:
        json.dump(results, outfile, cls=NumpyEncoder)

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
                for stat in args.statistic:
                    val = data['ci'][k][stat]
        except:
            print('error reading from {}, skip'.format(fname))
            continue

        insert_key_to_dictionary_as_list(df_data, 'fname', fname)
        for k in keys:
            for stat in args.statistic:
                insert_key_to_dictionary_as_list(df_data, k if (stat=="point_estimate") else "{} ({})".format(k, stat), data['ci'][k][stat])

        aic_diff = None; bic_diff = None
        if 'modelselection' in data:
            aic_diff = data['modelselection']['mixture_vs_inft_AIC']
            bic_diff = data['modelselection']['mixture_vs_inft_BIC']
        else:
            try:
                aic_diff = data['inft_optimize'][-1][1]['AIC'] - data['optimize'][-1][1]['AIC']
                bic_diff = data['inft_optimize'][-1][1]['BIC'] - data['optimize'][-1][1]['BIC']
            except:
                pass
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
        pd.concat([pd.DataFrame({'trait':[trait for i in data['power']['nvec']], 'nvec':data['power']['nvec'], 'svec':data['power']['svec']}) for data, trait in zip(data_list, traits_list)]).to_csv(args.out + '.power.csv', sep='\t', index=False)            
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
  
