#!/usr/bin/env python
'''
(c) 2018-2019 Oleksandr Frei, Alexey A. Shadrin
MiXeR software: Univariate and Bivariate Causal Mixture for GWAS
'''

import argparse
import numpy as np
from numpy import ma
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import matplotlib.patches as patches
from scipy.interpolate import interp1d
import json
import os

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
    plt.gca().add_patch(patches.Rectangle(((-abs(0.7*row.rg) if (rg < 0) else 0) , -0.7), abs(0.7 * rg), 0.15, fill=True, clip_on=False, color=clr))
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

