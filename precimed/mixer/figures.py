#!/usr/bin/env python
'''
(c) 2018-2019 Oleksandr Frei, Alexey A. Shadrin
MiXeR software: Univariate and Bivariate Causal Mixture for GWAS
'''

import argparse
import json
import os
import itertools

import numpy as np
from numpy import ma

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors

from mpl_toolkits.axes_grid1 import make_axes_locatable

from matplotlib_venn import venn2

from scipy.interpolate import interp1d
from scipy.stats import multivariate_normal

__version__ = '1.0.0'
MASTHEAD = "***********************************************************************\n"
MASTHEAD += "* mixer_figures.py: Visualization tools for MiXeR\n"
MASTHEAD += "* Version {V}\n".format(V=__version__)
MASTHEAD += "* (c) 2018-2019 Oleksandr Frei, Alexey A. Shadrin\n"
MASTHEAD += "* Norwegian Centre for Mental Disorders Research / University of Oslo\n"
MASTHEAD += "* GNU General Public License v3\n"
MASTHEAD += "***********************************************************************\n"

def make_qq_plot(qq, ci=True, ylim=7.3, xlim=7.3):
    hv_logp = np.array(qq['hv_logp']).astype(float)
    data_logpvec = np.array(qq['data_logpvec']).astype(float)
    model_logpvec = np.array(qq['model_logpvec']).astype(float)
    ylim_data = max(hv_logp[np.isfinite(data_logpvec)])
    model_logpvec[hv_logp > ylim_data]=np.nan
    if ci: 
        q = 10**-data_logpvec; dq= 1.96*np.sqrt(q*(1-q)/qq['sum_data_weights']);
        y1=hv_logp
        x1=ma.filled(-ma.log10(q+dq), np.nan)  #left CI bound
        x2=ma.filled(-ma.log10(q-dq), np.nan)  #right CI bound
        if True:
            y2 = np.empty(hv_logp.shape); y2[:]=np.nan;
            y2[x2<np.nanmax(x1)]=interp1d(x1, y1)(x2[x2<np.nanmax(x1)])                #upper CI bound
            y2[np.isnan(y2)]=ylim_data  
            plt.fill_between(x2, y1, y2, color=(0.1843, 0.3098, 0.3098), alpha=0.25)
        else:
            plt.plot(x1,hv_logp,x2,hv_logp)
    hData = plt.plot(data_logpvec, hv_logp)
    hModel = plt.plot(model_logpvec, hv_logp)
    hNull = plt.plot(hv_logp, hv_logp, 'k--')
    plt.ylim(0, ylim); plt.xlim(0, xlim)

def make_venn_plot(data, flip=False, factor='K', traits=['Trait1', 'Trait2'], colors=[0, 1], max_size=None):
    cm = plt.cm.get_cmap('tab10')

    if factor=='K': scale_factor=1000
    elif factor=='': scale_factor=1
    else: raise(ValueError('Unknow factor: {}'.format(factor)))
        
    n1 = data['ci']['nc1@p9']['point_estimate']/scale_factor; n1_se = data['ci']['nc1@p9']['std']/scale_factor
    n2 = data['ci']['nc2@p9']['point_estimate']/scale_factor; n2_se = data['ci']['nc2@p9']['std']/scale_factor
    n12 = data['ci']['nc12@p9']['point_estimate']/scale_factor; n12_se = data['ci']['nc12@p9']['std']/scale_factor
    rg = data['ci']['rg']['point_estimate']

    if max_size is None: max_size = n1+n2+n12
    if flip: n1, n2 = n2, n1; n1_se, n2_se = n2_se, n1_se
    f = lambda x: x if x < 7 else x+1

    v = venn2(subsets = (n1, n2, n12), normalize_to=(n1+n2+n12)/max_size, set_labels = ("", ""))
    v.get_patch_by_id('100').set_color(cm.colors[f(colors[0])])
    v.get_patch_by_id('010').set_color(cm.colors[f(colors[1])])
    v.get_patch_by_id('110').set_color(cm.colors[7])   
    formatter = '{:.2f}\n({:.2f})' if ((n1+n12+n2) < 1) else '{:.1f}\n({:.1f})' 
    v.get_label_by_id('100').set_text(formatter.format(n1, n1_se))
    v.get_label_by_id('010').set_text(formatter.format(n2, n2_se))
    v.get_label_by_id('110').set_text(formatter.format(n12, n12_se))

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

def plot_z_vs_z_data(df, traits=['Trait1', 'Trait2'], plot_limits=15, bins=100):
    '''
        # input can be generated as follows:
        import pandas as pd
        df1 = pd.read_csv(fname1, delim_whitespace=True, usecols=['SNP', 'A1', 'A2', 'Z'])
        df2 = pd.read_csv(fname2, delim_whitespace=True, usecols=['SNP', 'A1', 'A2', 'Z'])
        df = precimed.mixer.figures.merge_z_vs_z(df1, df2)
    '''
    plot_extent = [-plot_limits, plot_limits, -plot_limits, plot_limits]
    z, _, _ = np.histogram2d(df['Z2'], df['Z1'], bins=bins, range=[[-plot_limits, plot_limits], [-plot_limits, plot_limits]])
    im=plt.imshow(np.maximum(1,z),interpolation='none', origin='lower', cmap='hot', norm=matplotlib.colors.LogNorm(), vmin=1, vmax=1e4,extent=plot_extent)
    plt.xlabel('$z_{'+traits[0]+'}$')
    plt.ylabel('$z_{'+traits[1]+'}$')
    plt.colorbar(im, cax=make_axes_locatable(plt.gca()).append_axes("right", size="5%", pad=0.05))

def plot_predicted_zscore(data, num_snps, flip=False, traits=['Trait1', 'Trait2'], plot_limits=15, bins=100):
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

def make_power_plot(data_vec, colors=None, traits=None):
    if colors is None: colors = list(range(0, len(data_vec)))
    if traits is None: traits = ['trait{}'.format(i) for i in range(1, len(data_vec) + 1)]
    leg_labels = []
    current_n = []
    current_s = []
    
    cm = plt.cm.get_cmap('tab10')
    for data, color, trait in zip(data_vec, colors, traits):
        ax=plt.plot(np.log10(data['power']['nvec']), data['power']['svec'], color=cm.colors[color % 10], linestyle='solid' )
        current_n.append(data['options']['trait1_nval'])
        cs = interp1d(np.log10(data['power']['nvec']),data['power']['svec'])(np.log10(data['options']['trait1_nval']))
        current_s.append(cs)
        leg_labels.append('{0} ({1:.1f}%)'.format(trait, 100 * cs))

    for cn, cs, trait in zip(current_n, current_s, traits):
        plt.plot([np.log10(cn)], [cs], '*k', markersize=8)
    leg_labels.append('Current N')
    
    plt.legend(leg_labels, loc='upper left',frameon=False, numpoints=1)
    plt.xlabel('Sample size')
    plt.ylabel('Estimated percent variance explained\nby genome-wide significant SNPs')
    plt.xlim([4, 8])
    plt.ylim([0, 1])
    plt.locator_params(axis='x', nbins=5)
    plt.axes().set_xticklabels(labels=['10K', '100K', '1M', '10M', '100M'])

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
    parser.add_argument('--json', type=str, default="", help="json file from univariate analysis")    
    parser.add_argument('--trait1', type=str, default="Trait1", help="name of the first trait")    
    parser.set_defaults(func=func)

def parser_two_add_arguments(args, func, parser):
    parser.add_argument('--json', type=str, default="", help="json file from cross-trait analysis")    
    parser.add_argument('--trait1', type=str, default="Trait1", help="name of the first trait")    
    parser.add_argument('--trait2', type=str, default="Trait2", help="name of the second trait")    
    parser.set_defaults(func=func)

def parse_args(args):
    parser = argparse.ArgumentParser(description="MiXeR visualization tools.")

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--argsfile', type=open, action=LoadFromFile, default=None, help="file with additional command-line arguments")
    parent_parser.add_argument("--out", type=str, default="mixer", help="prefix for the output files")
    parent_parser.add_argument('--ext', type=str, default=['png'], nargs='+', choices=['png', 'svg'], help="output extentions")

    subparsers = parser.add_subparsers()
    parser_one_add_arguments(args=args, func=execute_one_parser, parser=subparsers.add_parser("one", parents=[parent_parser], help='produce figures for univariate analysis'))
    parser_two_add_arguments(args=args, func=execute_two_parser, parser=subparsers.add_parser("two", parents=[parent_parser], help='produce figures for cross-trait analysis'))

    return parser.parse_args(args)

def execute_two_parser(args):
    data = json.loads(open(args.json).read())
    plt.figure()
    plt.figure(figsize=[12, 3])
    plt.subplot(1,3,1)
    make_venn_plot(data, flip=False, traits=['SCZ', 'BIP'])
    plt.subplot(1,3,2)
    make_strat_qq_plots(data, flip=False, traits=['SCZ', 'BIP'], do_legend=False)
    plt.subplot(1,3,3)
    make_strat_qq_plots(data, flip=True, traits=['BIP', 'SCZ'], do_legend=True)
    for ext in args.ext:
        plt.savefig(args.out + '.' + ext, bbox_inches='tight')

def execute_one_parser(args):
    data = json.loads(open(args.json).read())
    plt.figure()
    make_qq_plot(data['qqplot'], ci=True)
    for ext in args.ext:
        plt.savefig(args.out + '.' + ext, bbox_inches='tight')

