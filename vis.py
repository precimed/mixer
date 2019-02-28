import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt

import os
from scipy.interpolate import interp1d
import scipy.stats as sstats
import itertools
import matplotlib as mpl
import math
import glob
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import multivariate_normal
from matplotlib.colors import LogNorm
from numpy import ma
from matplotlib_venn import venn2

NCKoef = 0.226  # see below for "calculate NCKoef"; this koef gives proportion of causal variants that explain 90% of heritability. 
                # it is specific to BGMG with single gaussian, with MAF specific model
DO_GCTA = False

#figures_folder = r'H:\Dropbox\analysis\2018_09_24_BGMG_manuscript\figures_and_tables_7pheno.model=full.run3b'
drive_letter = 'C';  # 'H'

figures_folder = '/home/oleksanf/vmshare/analysis/2019_02_11_MiXeR_display_items/figs'
if DO_GCTA: figures_folder = drive_letter + r':\Users\oleksanf\Dropbox\analysis\2018_09_24_BGMG_manuscript\figures_and_tables_laptop_2018_10_31_gcta'

tables_folder = figures_folder;
if not os.path.exists(tables_folder):  os.makedirs(tables_folder)
    
data_root = '/home/oleksanf/vmshare/data'    

format_string_UGMG_fit  = data_root + r"/LDSR/BGMG_result/{}"      +r".model=full.r2min=p05.randprune=n64p05.kmax=20000.run1.fit.json"
format_string_UGMG_test = data_root + r"/LDSR/BGMG_result/{}"      +r".model=full.r2min=p05.randprune=n64p05.kmax=20000.run1.test.json"
format_string_BGMG_fit  = data_root + r"/LDSR/BGMG_result/{}_vs_{}"+r".model=full.r2min=p05.randprune=n64p05.kmax=20000.run3bRGconst.fit.short.json"
format_string_BGMG_test = data_root + r"/LDSR/BGMG_result/{}_vs_{}"+r".model=full.r2min=p05.randprune=n64p05.kmax=20000.run3bRGconst.test.json"

format_string_UGMG_fit_censored  = data_root + r"/LDSR/BGMG_result/{}"      +r".model=full.r2min=p05.randprune=n64p05.kmax=20000.run2.fit.json"
format_string_UGMG_test_censored = data_root + r"/LDSR/BGMG_result/{}"      +r".model=full.r2min=p05.randprune=n64p05.kmax=20000.run2.test.json"
format_string_BGMG_fit_censored  = data_root + r"/LDSR/BGMG_result/{}_vs_{}"+r".model=full.r2min=p05.randprune=n64p05.kmax=20000.run2RGconst.fit.short.json"
format_string_BGMG_test_censored = data_root + r"/LDSR/BGMG_result/{}_vs_{}"+r".model=full.r2min=p05.randprune=n64p05.kmax=20000.run2RGconst.test.json"

format_string_LDSR_phase3_h2   = data_root + r"/MMIL/SUMSTAT/ANALYSIS/{}.h2.DEPRECATED.log"
format_string_LDSR_phase3_rg   = data_root +r"/MMIL/SUMSTAT/ANALYSIS/ldsr_rg.csv"
format_string_LDSR_h2   = data_root + r"/MMIL/SUMSTAT/ANALYSIS/{}.h2.log"
format_string_LDSR_rg   = data_root +r"/MMIL/SUMSTAT/ANALYSIS/ldsr_rg.DEPRECATED.csv"

folder_simu_bgmg_paper_examples = data_root + r'/run_simu_bgmg_paper_examples/final'

format_string_sumstats_noMHC = data_root + r"/MMIL/SUMSTAT/TMP/nomhc/{}.sumstats.gz"
format_string_sumstats_ldsr = data_root + r"/MMIL/SUMSTAT/TMP/ldsr/{}.sumstats.gz"

# TBD: there are some plots that hardcode the path - those must be configurable here...
globpat_SIMU_BGMG_11pifrac       = data_root + r'/SIMU_BGMG_11pifrac' + r'/*run3.bgmg.fit.short.json'
globpat_SIMU_BGMG_11pifrac_wave2 = data_root + r'/SIMU_BGMG_11pifrac_wave2' + r'/*run3.bgmg.fit.short.json'
globpat_SIMU_BGMG_spow2          = data_root + r'/SIMU_BGMG_spow2' + r'/*run3.bgmg.fit.short.json'

folder_SIMU_BGMG_GCTA = drive_letter + r":\work\SIMU_BGMG_11pifrac_gctaSigma"
folder_SIMU_BGMG_annotenrich = drive_letter + r':\work\SIMU_BGMG_annotenrich'
folder_SIMU_BGMG_subref = drive_letter + r':\work\SIMU_BGMG_subref'

nsnps_HAPGEN = 11015833;
nsnps_LDSR = 9997231;

if True:
    traits4_ordered = ['PGC_SCZ_2014_EUR', 'PGC_BIP_2016', 'SSGAC_EDU_2018_no23andMe', 'GIANT_HEIGHT_2018_UKB']
    traits4 = {'GIANT_HEIGHT_2018_UKB':'Height', 'PGC_BIP_2016':'Bipolar Disorder', 'PGC_SCZ_2014_EUR':'Schizophrenia', 'SSGAC_EDU_2018_no23andMe': 'Educational attainment'}
    traits4_short = {'GIANT_HEIGHT_2018_UKB':'HEIGHT', 'PGC_BIP_2016':'BIP', 'PGC_SCZ_2014_EUR':'SCZ', 'SSGAC_EDU_2018_no23andMe': 'EDU'}
    trait4_index_map = {trait:index for index, trait in enumerate([trait for trait in traits4_ordered])}

if True:
    traits_ordered = ['PGC_SCZ_2014_EUR', 'PGC_BIP_2016','PGC_MDD_2018_no23andMe',
                      'PGC_ASD_2017_iPSYCH','PGC_ADHD_2017_EUR','SSGAC_EDU_2018_no23andMe',
                      'GIANT_HEIGHT_2018_UKB']
    traits_short = {'PGC_BIP_2016':'BIP','PGC_SCZ_2014_EUR':'SCZ','PGC_MDD_2018_no23andMe':'MDD',
                      'PGC_ASD_2017_iPSYCH':'ASD','PGC_ADHD_2017_EUR':'ADHD','SSGAC_EDU_2018_no23andMe':'EDU',
                      'GIANT_HEIGHT_2018_UKB':'HEIGHT'}
    traits = {'GIANT_HEIGHT_2018_UKB':'Height',
              'PGC_BIP_2016':'Bipolar Disorder',
              'PGC_SCZ_2014_EUR':'Schizophrenia',
              'SSGAC_EDU_2018_no23andMe':'Educational attainment',
              'PGC_MDD_2018_no23andMe':'Major Depressive Disorder',
              'PGC_ASD_2017_iPSYCH':'Autism Spectrum Disorder',
              'PGC_ADHD_2017_EUR':'ADHD'}  # Attention deficit hyperactivity disorder
    traitPS_index_map = {trait:index for index, trait in enumerate([trait for trait in traits_ordered])}

traits_immuno = [('OKADA_RA_2014_EUR','Rheumatoid Arthritis','RA'),
                     ('IIBDGC_IBD_2017','Inflam. Bowel Disease','IBD'),
                     ('IIBDGC_CD_2017','Crohn Disease','CD'),
                     ('IIBDGC_UC_2017','Ulcerative Colitis','UC')]

traits_antro = [('EGG_BIRTHWEIGHT_2016','Birth Weight','BirthWeight'),
                ('GIANT_WHR_2015_EUR','Waist Hip Ratio','WHR'),
                ('GIANT_BMI_2015_EUR','Body Mass Index','BMI'),
                ('GIANT_HEIGHT_2018_UKB','Height','HEIGHT')]

traits4_psych = [('PGC_SCZ_2014_EUR','Schizophrenia','SCZ'),
                 ('PGC_BIP_2016','Bipolar Disorder','BIP'),
                 ('PGC_MDD_2018_no23andMe', 'Major Depressive Disorder', 'MDD'),
                 ('PGC_ASD_2017_iPSYCH', 'Autism Spectrum Disorder', 'ASD'),
                 ('PGC_ADHD_2017_EUR', 'ADHD', 'ADHD'),
                 ('SSGAC_EDU_2018_no23andMe','Educational attainment','EDU'),
                 ('GIANT_HEIGHT_2018_UKB','Height','HEIGHT')]

def create_traits_lists(full_list):
    traits_ordered = [idx for idx, full, short in full_list] 
    traits = {idx:full for idx, full, short in full_list}
    traits_short = {idx:short for idx, full, short in full_list}
    trait_index_map = {trait:index for index, trait in enumerate([trait for trait in traits_ordered])}
    return traits_ordered, traits, traits_short, trait_index_map

traits14_ordered, traits14, traits14_short, trait14_index_map = create_traits_lists(traits4_psych + [(x,y,z) for x,y,z in traits_antro if z != 'HEIGHT'] + traits_immuno)
traitsIM_ordered, traitsIM, traitsIM_short, traitIM_index_map = create_traits_lists(traits_immuno)
traitsAM_ordered, traitsAM, traitsAM_short, traitAM_index_map = create_traits_lists(traits_antro)

plt.rcParams.update({'mathtext.default': 'regular', 'font.size': 20 })

# Figures (#enabled means I made these plots on the new Ubuntu laptop)
DO_QQ_MODEL_DATA_vs_NULL = False      # enabled
DO_QQ_MODEL_vs_DATA = False
DO_QQ_BINS_MODEL_DATA_vs_NULL = False # enabled
DO_POWER_PLOT = False                 # enabled
DO_ADJ_TRAJECTORY = False
DO_ADJ_TRAJECTORY_BGMG = False
DO_VENN_DIAGRAMS = False              # enabled, works only after DO_BGMG_TABLE
DO_VENN_DIAGRAMS_SUPPL = False        # enabled 
DO_STRATIFIED_QQ=False                # enabled

DO_SIMU_QQ=False
DO_SIMU_QQ_BINS=False
DO_SIMU_STRATIFIED_QQ=False

DO_SIMU_UGMG=False                 # enabled  
DO_SIMU_BGMG=False                 # enabled (this does both "spow=0 & h2=0.1, 0.4, 0.7" and "h2=0.4 & spow=-0.25, -0.5, -0.75")
DO_SIMU_UGMG_ANNOTENRICH=False

# Tables
DO_UGMG_TABLE = False              # enabled
DO_BGMG_TABLE = False              # enabled
DO_GWAS_DATA_TABLE = False;
DO_SIMU_UGMG_TABLE = True
DO_SIMU_BGMG_TABLE = True
DO_SIMU_BGMG_ANNOTENRICH_TABLE=False
DO_SIMU_UGMG_SUBREF = False
table_index = 0; tables_writer = None

def concat_se(df, zeros_as_nan=True):
    df = df.copy()
    se_cols = [col[:-3] for col in df if col.endswith('_se')]
    for se_col in se_cols:
        df[se_col + '_se'] = ['{} ({})'.format(val, 'nan' if ((float(val_se)==0) and zeros_as_nan) else val_se) for val, val_se in zip(df[se_col].values, df[se_col + '_se'].values)]
        del df[se_col]
    df.columns = [(x.replace('_se', ' (se)') if x.endswith('_se') else x) for x in df.columns]
    
    return df

def concat_se_mean_std(df):
    df = df.copy()
    se_cols = [col[:-8] for col in df if col.endswith('_se_mean')]
    for se_col in se_cols:
        df[se_col] = ['{} ({} / {})'.format(val, 'nan' if (float(val_se_mean)==0) else val_se_mean,  'nan' if (float(val_std)==0) else val_std) for val, val_se_mean, val_std  in zip(df[se_col + '_mean'].values, df[se_col + '_se_mean'], df[se_col + '_std'].values)]
        del df[se_col+'_mean']
        del df[se_col+'_std']
        del df[se_col+'_se_mean']
    df.columns = [(x.replace('_se', ' (se/std)') if x.endswith('_se') else x) for x in df.columns]
    return df

def savefig(folder, filename, exts=['svg', 'png']):
    for ext in exts:
        plt.savefig(os.path.join(folder, ext, filename + '.' + ext), bbox_inches='tight')
        plt.savefig(os.path.join(folder, ext, filename + '.' + ext), bbox_inches='tight')
    
def json_loads(fname):
    data = json.loads(open(fname).read())
    if 'result' not in data: data={'result': data}  
    return data
    
def make_qq_plot(qq, ci=True):
    hv_logp = np.array(qq['hv_logp']).astype(float)
    data_logpvec = np.array(qq['data_logpvec']).astype(float)
    model_logpvec = np.array(qq['model_logpvec']).astype(float)
    ylim_data = max(hv_logp[np.isfinite(data_logpvec)])
    model_logpvec[hv_logp > ylim_data]=np.nan
    if ci: 
        q = 10**-data_logpvec; dq= 1.96*np.sqrt(q*(1-q)/qq['qq_options']['sum_data_weights']);
        y1=hv_logp
        x1=ma.filled(-ma.log10(q+dq), np.nan)  #left CI bound
        x2=ma.filled(-ma.log10(q-dq), np.nan)  #right CI bound
        if True:
            y2 = np.empty(hv_logp.shape); y2[:]=np.nan;
            y2[x2<max(x1)]=interp1d(x1, y1)(x2[x2<max(x1)])                #upper CI bound
            y2[np.isnan(y2)]=ylim_data  
            plt.fill_between(x2, y1, y2, color=(0.1843, 0.3098, 0.3098), alpha=0.25)
        else:
            plt.plot(x1,hv_logp,x2,hv_logp)
    hData = plt.plot(data_logpvec, hv_logp)
    hModel = plt.plot(model_logpvec, hv_logp)
    hNull = plt.plot(hv_logp, hv_logp, 'k--')

def plot_simu_bgmg_pi12(df_plot, do_title=True):
    h2_index=df_plot.true_h2.iloc[0]
    spow_index=df_plot.spow.iloc[0]
    pi1u=df_plot.true_pi1u.iloc[0]
    pi2u=df_plot.true_pi2u.iloc[0]
    scale_coef, scale_text = (1e4, '10^{-4}') if pi2u==3e-4 else (1e3, '10^{-3}')
    bars1 = df_plot['pi12_mean'].values*scale_coef
    bars2 = df_plot['true_pi12'].values*scale_coef
    yer1 = df_plot['pi12_se_mean'].values*scale_coef
    yer2 = df_plot['pi12_std'].values*scale_coef

    barWidth = 0.3
    r1 = np.arange(len(bars1)); r2 = [x + barWidth for x in r1] # # The x position of bars

    plt.bar(r1, bars1, width = barWidth, color = 'blue', edgecolor = 'black', yerr=yer2, capsize=7, label='estimated')
    plt.bar(r2, bars2, width = barWidth, color = 'cyan', edgecolor = 'black',            capsize=7, label='true')

    plt.ylabel('estimated $\pi_{12}\;(\\times ' + scale_text + '$)' )  # @10k
    plt.xlabel(     'true $\pi_{12}\;(\\times ' + scale_text + '$)' )
    if do_title: plt.title('${}=${}, $\pi^u_1=${:.0e}, $\pi^u_2=${:.0e}'.format('S' if (spow_index != 0) else h2, spow_index if (spow_index != 0) else h2_index, pi1u, pi2u).replace('e-0', 'e-'))
    plt.ylim((-0.0*pi2u*scale_coef, 1.05*pi2u*scale_coef))
    plt.gca().set_xticks([r + 0.5*barWidth for r in range(len(bars1))]);
    plt.gca().set_xticklabels(['{:.1f}'.format(x) for x in df_plot['true_pi12'].values*scale_coef]);
    plt.locator_params(axis='y', nbins=4)

def plot_simu_bgmg_rg_or_rho12(df_plot, do_rg, do_title=True):  # do_rg: plot rg or rho_beta?
    h2_index=df_plot.true_h2.iloc[0]
    spow_index=df_plot.spow.iloc[0]
    pi1u=df_plot.true_pi1u.iloc[0]
    pi2u=df_plot.true_pi2u.iloc[0]
    scale_coef, scale_text = (1e4, '10^{-4}') if pi2u==3e-4 else (1e3, '10^{-3}')

    bars1 = df_plot['rg_mean'].values if do_rg else df_plot['rho12_mean'].values
    bars2 = float(true_rg)*df_plot['true_pi12'].values / np.sqrt(df_plot['true_pi1u']*df_plot['true_pi2u']) if do_rg else float(true_rg)*np.ones(df_plot['true_pi12'].shape)
    yer1 = df_plot['rg_se_mean'].values if do_rg else df_plot['rho12_se_mean'].values
    yer2 = df_plot['rg_std'].values if do_rg else df_plot['rho12_std'].values

    barWidth = 0.3
    r1 = np.arange(len(bars1)); r2 = [x + barWidth for x in r1] # # The x position of bars

    if do_rg==False: yer2[0]=np.nan
    pltbar = plt.bar(r1, bars1, width = barWidth, color = 'blue', edgecolor = 'black', yerr=yer2, capsize=7, label='estimated')
    if do_rg==False: pltbar[0].set_visible(False)
    pltbar = plt.bar(r2, bars2, width = barWidth, color = 'cyan', edgecolor = 'black',            capsize=7, label='true')
    if do_rg==False: pltbar[0].set_visible(False)

    plt.ylabel('estimated $rg$' if do_rg else 'estimated $\\rho_{12}$')
    plt.xlabel(     'true $\pi_{12}\;(\\times ' + scale_text + '$)' )
    plt.gca().set_xticks([r + 0.5*barWidth for r in range(len(bars1))]);
    plt.gca().set_xticklabels(['{:.1f}'.format(x) for x in df_plot['true_pi12'].values*scale_coef]);
    #plt.legend()
    #if do_title: plt.title('$h2=${}, $\pi^u_1=${:.0e}, $\pi^u_2=${:.0e}'.format(h2_index, pi1u, pi2u).replace('e-0', 'e-'))
    if do_title: plt.title('${}=${}, $\pi^u_1=${:.0e}, $\pi^u_2=${:.0e}'.format('S' if (spow_index != 0) else h2, spow_index if (spow_index != 0) else h2_index, pi1u, pi2u).replace('e-0', 'e-'))
    if df_plot.true_rg.iloc[0]=='0.0':
        a=max(np.abs(plt.gca().get_ylim())); plt.gca().set_ylim([-a, a]);

for censored in [False]:
    if DO_QQ_MODEL_DATA_vs_NULL:
        fig = plt.figure(figsize=(15, 9+9+4.5), dpi=80)
        for index, trait in enumerate(traits14_ordered):
            plt.subplot(5,3,1+index)
            fname = format_string_UGMG_test_censored.format(trait) if censored else format_string_UGMG_test.format(trait)
            data = json_loads(fname)
            qq = data['result']['univariate'][0]['qq_plot_data']
            make_qq_plot(qq)
            plt.ylim(0, 7.3); plt.xlim(0, 7.3)
            plt.title(traits14[trait])
            #plt.legend(['data', 'model', 'null'], loc='lower right')
            #plt.xlabel('Expected(-$log_{10}$ p-value)')
            #plt.ylabel('Observed(-$log_{10}$ p-value)')
        fig.subplots_adjust(hspace=0.35)
        fig.text(0.5, 0.075, 'Expected(-$log_{10}$ p-value)', ha='center', fontsize=24)
        fig.text(0.075, 0.5, 'Observed(-$log_{10}$ p-value)', va='center', fontsize=24, rotation='vertical')
        savefig(figures_folder, 'QQ_MODEL_DATA_vs_NULL'+(' (censored)' if censored else ''))

    if DO_QQ_MODEL_vs_DATA:
        plt.figure(figsize=(12, 9), dpi=80)
        leg_labels = []
        for index, trait in enumerate(traits_ordered):
            fname = format_string_UGMG_test_censored.format(trait) if censored else format_string_UGMG_test.format(trait)
            data = json.loads(open(fname).read())
            qq = data['result']['univariate'][0]['qq_plot_data']
            #hModelvsData = plt.plot(qq['model_logpvec'], qq['data_logpvec'])
            hModelvsData = plt.plot(interp1d(qq['data_logpvec'],qq['hv_logp'],bounds_error=False)(qq['hv_logp']),interp1d(qq['model_logpvec'],qq['hv_logp'],bounds_error=False)(qq['hv_logp']))
            plt.ylim(0, 7.3)
            plt.xlim(0, 7.3)
            plt.title('')
            plt.legend(['model-vs-data', 'null'], loc='lower right')
            leg_labels.append(traits[trait])
        hNull = plt.plot(qq['hv_logp'], qq['hv_logp'], 'k--')
        leg_labels.append('null')
        plt.legend(leg_labels, loc='lower right')
        plt.title('QQ plot, model vs data')
        plt.xlabel('Observed(-$log_{10}$ p-value, model)')
        plt.ylabel('Observed(-$log_{10}$ p-value, data)')
        
        savefig(figures_folder, 'QQ_MODEL_vs_DATA'+(' (censored)' if censored else ''))

if DO_QQ_BINS_MODEL_DATA_vs_NULL:
    for index, trait in enumerate(traits14_ordered):
        fig = plt.figure(figsize=(16, 12), dpi=80)
        fig.suptitle(traits14[trait]) #, fontsize=16)
        fname = format_string_UGMG_test.format(trait)
        data = json_loads(fname)
        for subplot_index in range(0,9):
            axis = plt.subplot(3,3,1+subplot_index)
            qq=data['result']['univariate'][0]['qq_plot_bins_data'][subplot_index]
            make_qq_plot(qq)

            #plt.title(str(subplot_index) + qq["options"]["title"].replace('$$', '$')) #.split('$$  $$')
            qq["options"]["title"] = qq["options"]["title"].replace('L \\in', 'LDscore \\in$\n$').replace('maf', 'MAF').replace(',', ', ').replace('-Infinity', '0').replace('Infinity', 'Inf')

            if True:
                axis.xaxis.set_label_position('top')
                axis.yaxis.set_label_position('right') 
                label_indices = {'y':[2, 5, 8], 'x':[0, 1, 2]}
            else:
                label_indices = {'y':[0, 3, 6], 'x':[6, 7, 8]}
            
            if subplot_index in label_indices['y']: plt.ylabel(qq["options"]["title"].split('$  $')[1].replace('$$', '$'))
            if subplot_index in label_indices['x']: plt.xlabel(qq["options"]["title"].split('$  $')[0].replace('$$', '$'))
            plt.ylim(0, 7.3); plt.xlim(0, 7.3)
            #plt.tight_layout()
        fig.text(0.5, 0.075, 'Expected(-$log_{10}$ p-value)', ha='center')
        fig.text(0.075, 0.5, 'Observed(-$log_{10}$ p-value)', va='center', rotation='vertical')
        savefig(figures_folder, 'QQ_BINS_MODEL_DATA_vs_NULL_{}'.format(trait))

if DO_POWER_PLOT:
    for repeat in [1, 0]:
        traits_this_run_ordered = traits_ordered if (repeat==0) else traitsAM_ordered + traitsIM_ordered
        fig = plt.figure(figsize=(12, 9), dpi=80)
        leg_labels = []
        current_n = []; current_s = []
        for index, trait in enumerate(traits_this_run_ordered):
            fname = format_string_UGMG_test.format(trait)
            data = json_loads(fname)
            pdata = data['result']['univariate'][0]['power_plot_data']
            plt.plot(np.log10(pdata['power_nvec']), pdata['power_svec'])

            current_n.append(data['result']['options']['trait1_nval'])
            current_s.append(interp1d(np.log10(pdata['power_nvec']),pdata['power_svec'])(np.log10(current_n[-1])))

            leg_labels.append('{0} ({1:.1f}%)'.format(traits14_short[trait],100* current_s[-1]))
        plt.plot([-1], [-1], '*k', markersize=8)
        leg_labels.append('Current N')
        for index,trait in enumerate(traits_this_run_ordered):
            plt.plot([np.log10(current_n[index])], [current_s[index]], '*k', markersize=8)

        plt.legend(leg_labels, loc='upper left',frameon=False, numpoints=1)
        plt.title('Power plot')
        plt.xlabel('$log_{10}$ N')
        plt.ylabel('Fraction of h2 explained by\ngenome-wide significant loci')
        plt.xlim([3, 8])
        plt.ylim([0, 1])

        savefig(figures_folder, 'POWER_PLOT_{}'.format('psych' if (repeat==0) else 'antro_immuno'))

if DO_ADJ_TRAJECTORY:
    fig = plt.figure(figsize=(18+9, 12+6), dpi=80)
    for index, trait in enumerate(traits_ordered):
        plt.subplot(2+1,2+1,1+index)
        fname = format_string_UGMG_fit.format(trait)
        data = json.loads(open(fname).read())
        lat = data['result']['univariate'][0]['loglike_adj_trajectory']
        hData = plt.plot(lat['pivec'], lat['cost'], '.-')
        plt.title(traits[trait], loc='right')
        #plt.legend(['data', 'model', 'null'], loc='lower right')
        plt.xlabel('$\pi_1$')
        plt.ylabel('log likelihood')
        plt.locator_params(axis='x', nbins=5)
        
        params_pi = data['result']['univariate'][0]['params']['pi_vec']
        params_cost = interp1d(lat['pivec'],lat['cost'])(params_pi)
        plt.plot([params_pi],[params_cost], '*k', markersize=8)
        fmtr = mpl.ticker.ScalarFormatter(useOffset=math.floor(min(lat['cost'])))
        plt.gca().get_yaxis().set_major_formatter(fmtr)
    
    fig.subplots_adjust(hspace=0.25)
    savefig(figures_folder, 'ADJ_TRAJECTORY')

table_index += 2
if DO_UGMG_TABLE:
    df_data = {}
    def insert_key_to_dictionary_as_list(key, value):
        if key not in df_data:
            df_data[key] = []
        df_data[key].append(value)

    for censored in [False, True]:
        for index, trait in enumerate(traits14_ordered):
            fname = format_string_UGMG_fit_censored.format(trait) if censored else format_string_UGMG_fit.format(trait)
            data = json_loads(fname)
            formats = {'h2':'{:.4f}', 'sig2_zero':'{:.4f}', 'pi_vec':'{:.2e}', 'sig2_beta':'{:.2e}'}
            ci = data['result']['univariate'][0]['ci']
            insert_key_to_dictionary_as_list('trait', traits14[trait] + (' (censored)' if censored else ''))
            for key in ci:
                if key == 'sig2_zero_minus1': continue
                column_name = key

                insert_key_to_dictionary_as_list(column_name, formats[key].format(ci[key]['point_estimate']))
                insert_key_to_dictionary_as_list(column_name + '_se',  formats[key].format(ci[key]['se']))
            ldsr_data = {}
            for line in open(format_string_LDSR_h2.format(trait)).read().split('\n'):
                if 'Total Observed scale h2: ' in line: ldsr_data['LDSR(2)_h2'], ldsr_data['LDSR(2)_h2_se'] = tuple(line.split(': ')[1].replace('(', '').replace(')', '').split(' '))
                if 'Intercept: ' in line: ldsr_data['LDSR(2)_intercept'], ldsr_data['LDSR(2)_intercept_se'] = tuple(line.split(': ')[1].replace('(', '').replace(')', '').split(' '))
            for line in open(format_string_LDSR_phase3_h2.format(trait)).read().split('\n'):
                if 'Total Observed scale h2: ' in line: ldsr_data['LDSR_h2'], ldsr_data['LDSR_h2_se'] = tuple(line.split(': ')[1].replace('(', '').replace(')', '').split(' '))
                if 'Intercept: ' in line: ldsr_data['LDSR_intercept'], ldsr_data['LDSR_intercept_se'] = tuple(line.split(': ')[1].replace('(', '').replace(')', '').split(' '))
            if len(ldsr_data) != 8: raise('LDSR data error')
            for key in ldsr_data: insert_key_to_dictionary_as_list(key, float(ldsr_data[key]))

        df = pd.DataFrame(df_data)
        for c in [c for c in df.columns if c.startswith('pi_')]: df[c.replace('pi_', 'nc@p9_').replace('_vec', '')] = ['{:.1f}'.format(x*NCKoef/1000) for x in df[c].astype(float).values*nsnps_LDSR]
        for c in [c for c in df.columns if c.startswith('pi_')]: df[c.replace('pi_', 'nc_'    ).replace('_vec', '')] = ['{:.1f}'.format(x/1000) for x in df[c].astype(float).values*nsnps_LDSR]
        df = df[['trait'] +  [x for x in df.columns if x != 'trait']]
        df = df[[x for x in df.columns if not x.startswith('LDSR')] +  [x for x in df.columns if x.startswith('LDSR')]]
        df.to_csv(os.path.join(tables_folder, 'UGMG_TABLE.csv'), sep='\t',index=False)

        if tables_writer is None: tables_writer = pd.ExcelWriter(os.path.join(tables_folder, 'all_tables.xlsx'))
        concat_se(df).to_excel(tables_writer,'table{}'.format(table_index-1), index=False)
        
        df_main = df[~df['trait'].str.contains('censored') & df['trait'].isin([traits4[k] for k in traits4])][['trait', 'nc@p9', 'nc@p9_se', 'h2', 'h2_se', 'LDSR_h2', 'LDSR_h2_se']].copy()
        concat_se(df_main).to_excel(tables_writer,'table{}'.format(table_index), index=False)

table_index += 2
if DO_BGMG_TABLE:
    bgmg_table_df_censored = None
    for censored in [False, True]:
        df_data = {}
        def insert_key_to_dictionary_as_list(key, value):
            if key not in df_data:
                df_data[key] = []
            df_data[key].append(value)

        for trait1, trait2 in list(itertools.combinations(traits_ordered, 2)) + list(itertools.combinations(traitsAM_ordered, 2)) + list(itertools.combinations(traitsIM_ordered, 2)):
            format_string = format_string_BGMG_fit_censored if censored else format_string_BGMG_fit 
            fname1 = format_string.format(trait1, trait2)
            fname2 = format_string.format(trait2, trait1)
            keys_rg = {'rg':'rg', 'rho_beta':'rho12'}
            keys_pifrac = {'pi12_over_min_piXu': 'pi12frac'}

            if os.path.exists(fname2):
                data = json_loads(fname2)
                keys_pivec = {'pi_vec_C1':'nc2','pi_vec_C2':'nc1','pi_vec_C3':'nc12'}
                #trait1, trait2 = trait2, trait1
                flip_data = True
            elif os.path.exists(fname1):
                data = json_loads(fname1)
                keys_pivec = {'pi_vec_C1':'nc1','pi_vec_C2':'nc2','pi_vec_C3':'nc12'}
                flip_data = False
            else:
                print('missing: {} vs {}'.format(trait1, trait2))
                continue

            insert_key_to_dictionary_as_list('trait1', traits14[trait1] +(' (censored)' if censored else ''))
            insert_key_to_dictionary_as_list('trait2', traits14[trait2] +(' (censored)' if censored else ''))
            insert_key_to_dictionary_as_list('trait1_ID', trait1)
            insert_key_to_dictionary_as_list('trait2_ID', trait2)

            ci = data['result']['bivariate']['ci']
            formats = {'pi_vec_C1':'{:.1f}', 'pi_vec_C2':'{:.1f}', 'pi_vec_C3':'{:.1f}','pi12_over_min_piXu':'{:.4f}','rg':'{:.4f}','rho_beta':'{:.4f}'}

            for key in keys_rg:
                insert_key_to_dictionary_as_list(keys_rg[key], formats[key].format(ci[key]['point_estimate']))
                insert_key_to_dictionary_as_list(keys_rg[key] + '_se', formats[key].format(ci[key]['se']))
            for key in keys_pivec:
                insert_key_to_dictionary_as_list(keys_pivec[key], formats[key].format(ci[key]['point_estimate']*nsnps_LDSR))
                insert_key_to_dictionary_as_list(keys_pivec[key] + '_se', formats[key].format(ci[key]['se']*nsnps_LDSR))
                insert_key_to_dictionary_as_list(keys_pivec[key]+'@p9', '{:.2f}'.format(ci[key]['point_estimate']*nsnps_LDSR*NCKoef))
                insert_key_to_dictionary_as_list(keys_pivec[key]+'@p9_se', '{:.2f}'.format(ci[key]['se']*nsnps_LDSR*NCKoef))
            #for key in keys_pifrac:
            #    insert_key_to_dictionary_as_list(keys_pifrac[key], formats[key].format(ci[key]['point_estimate']))
            #    insert_key_to_dictionary_as_list(keys_pifrac[key] + '_se', formats[key].format(ci[key]['se']))

            ldsr_header = 'p1 p2 rg se z p h2_obs h2_obs_se  h2_int h2_int_se   gcov_int gcov_int_se'
            has_ldsr_data = 0
            for ldsr_formatter, LDSRtag in [(format_string_LDSR_rg, 'LDSR(2)_'), (format_string_LDSR_phase3_rg, 'LDSR_')]:
                for line in open(ldsr_formatter).read().split('\n'):
                    if '--rg' in line: continue
                    if 'Reading summary statistics from' in line: continue
                    if (trait1 not in line) or (trait2 not in line): continue
                    for key, value in zip(ldsr_header.split(), line.split()):
                        #if key not in 'rg se p p1 p2'.split(): continue   # good itea to take a look at p1 p2 to make sure we are picking correct traits
                        if key not in 'rg se'.split(): continue
                        insert_key_to_dictionary_as_list(LDSRtag + {'rg':'rg', 'se':'rg_se', 'p':'p'}[key], {'p1':'{}', 'p2':'{}', 'p':'{:.3e}', 'rg':'{:.4f}', 'se':'{:.4f}'}[key].format(value if key in ['p1', 'p2'] else float(value)))
                    has_ldsr_data += 1
                    break
            if (has_ldsr_data != 2): raise('error, LDSR rg data for {} vs {} is not available'.format(trait1, trait2))

        bgmg_table_df = pd.DataFrame(df_data)
        #bgmg_table_df.loc[(bgmg_table_df.trait2_ID == 'PGC_MDD_2018_no23andMe') & (bgmg_table_df.trait1_ID=='PGC_SCZ_2014_EUR'),'rho12_se']=0
        #bgmg_table_df.loc[(bgmg_table_df.trait2_ID == 'PGC_MDD_2018_no23andMe') & (bgmg_table_df.trait1_ID=='PGC_SCZ_2014_EUR'),'rg_se']=0
        bgmg_table_df = bgmg_table_df[['trait1', 'trait2'] +  [x for x in bgmg_table_df.columns if x not in ['trait1', 'trait2']]]
        bgmg_table_df = bgmg_table_df[[x for x in bgmg_table_df.columns if not x.startswith('LDSR')] +  [x for x in bgmg_table_df.columns if x.startswith('LDSR')]]
        for col in bgmg_table_df.columns: bgmg_table_df[col] = pd.to_numeric(bgmg_table_df[col],errors='ignore')

        # save bgmg_table_df to dictionary, and clear bgmg_table_df variable 
        # this is to make sure we don't confuse censored and non-censored results in venn diagram
        if ('bgmg_table_df_censored' not in locals()) or (bgmg_table_df_censored is None):
            bgmg_table_df_censored = {}
        bgmg_table_df_censored[censored] = bgmg_table_df
    
    bgmg_table_df = pd.concat([bgmg_table_df_censored[k] for k in bgmg_table_df_censored])
    df = bgmg_table_df[[x for x in bgmg_table_df.columns if x not in ['trait1_ID', 'trait2_ID']]].copy()
    for c in df.columns:
        if c.startswith('nc'): df[c] =  ['{:.2f}'.format(x/1000) for x in df[c].values]
    df.to_csv(os.path.join(tables_folder, 'BGMG_TABLE'+(' (censored)' if censored else '')+'.csv'), sep='\t',index=False)

    if tables_writer is None: tables_writer = pd.ExcelWriter(os.path.join(tables_folder, 'all_tables.xlsx'))
    concat_se(df)[['trait1','trait2', 'nc12@p9 (se)', 'nc1@p9 (se)','nc2@p9 (se)','nc12 (se)','nc1 (se)','nc2 (se)','rho12 (se)','rg (se)','LDSR_rg (se)','LDSR(2)_rg (se)']].to_excel(tables_writer,'table{}'.format(table_index), index=False)
    df_main = df[~df['trait1'].str.contains('censored') & ~df['trait2'].str.contains('censored') &
                 df['trait1'].isin([traits4[k] for k in traits4]) & 
                 df['trait2'].isin([traits4[k] for k in traits4])].copy()
    df_main = df_main[['trait1', 'trait2', 'nc12@p9', 'nc12@p9_se', 'nc1@p9', 'nc1@p9_se', 'nc2@p9', 'nc2@p9_se',
                       'rho12', 'rho12_se', 'rg', 'rg_se', 'LDSR_rg', 'LDSR_rg_se']].copy()    
    for c in [c for c in df_main.columns if (('rg' in c) or ('rho' in c))]: df_main[c] = ['{:.3f}'.format(x) for x in df_main[c].astype(float).values]
    concat_se(df_main, zeros_as_nan=False).to_excel(tables_writer,'table{}'.format(table_index-1), index=False)
     
    bgmg_table_df = None
        

if DO_VENN_DIAGRAMS:
    for censored in [False, True]:
        bgmg_table_df = bgmg_table_df_censored[censored]
        fig=plt.figure(figsize=(12, 9), dpi=80)
        cm = plt.cm.get_cmap('tab10')
        for index, row in bgmg_table_df[bgmg_table_df.trait1_ID.isin(traits4_ordered) & bgmg_table_df.trait2_ID.isin(traits4_ordered)].reset_index(drop=True).iterrows():
            plt.subplot(2, 3, 1+index)
            n1 = row['nc1@p9']/1000; n1_se = row['nc1@p9_se']/1000
            n2 = row['nc2@p9']/1000; n2_se = row['nc2@p9_se']/1000;
            n12 = row['nc12@p9']/1000; n12_se = row['nc12@p9_se']/1000
            max_size = 12
            v = venn2(subsets = (n1, n2, n12), normalize_to=(n1+n2+n12)/max_size, set_labels = ('',''))
            v.get_patch_by_id('100').set_color(cm.colors[trait4_index_map[row.trait1_ID]])
            v.get_patch_by_id('010').set_color(cm.colors[trait4_index_map[row.trait2_ID]])    
            v.get_patch_by_id('110').set_color(cm.colors[7])    
            v.get_label_by_id('100').set_text('{:.1f}\n({:.1f})'.format(n1, n1_se))
            v.get_label_by_id('010').set_text('{:.1f}\n({:.1f})'.format(n2, n2_se))
            v.get_label_by_id('110').set_text('{:.1f}\n({:.1f})'.format(n12, n12_se))
            #print(plt.xlim(), plt.ylim())bot
            plt.xlim([-0.75, 0.75]), plt.ylim([-0.6, 0.6])
            plt.title(traits4_short[row.trait1_ID] +' & ' + traits4_short[row.trait2_ID])
        fig.subplots_adjust(hspace=-0.6, wspace=-0.3)
        plt.tight_layout()
        savefig(figures_folder, 'VENN_DIAGRAMS'+(' (censored)' if censored else ''))

        for trait_index_map, trait_group, venn_figsize in [(traitPS_index_map, 'psych', (30, 28)), (traitIM_index_map, 'immuno', (15, 14)), (traitAM_index_map, 'antro', (15, 14))]:
            plt.figure(figsize=venn_figsize, dpi=80)
            cm = plt.cm.get_cmap('tab10')
            for index, row in bgmg_table_df.iterrows():
                if (row.trait1_ID not in trait_index_map) or (row.trait2_ID not in trait_index_map): continue
                i1 = trait_index_map[row.trait1_ID]
                i2 = trait_index_map[row.trait2_ID]
                n1 = row['nc1@p9']/1000; n1_se = row['nc1@p9_se']/1000
                n2 = row['nc2@p9']/1000; n2_se = row['nc2@p9_se']/1000;
                n12 = row['nc12@p9']/1000; n12_se = row['nc12@p9_se']/1000
                t1 = row.trait1_ID
                t2 = row.trait2_ID
                #if i1<i2: i1, i2 = i2, i1
                f = lambda x: x if x < 7 else x+1
                for i in range(1,3):
                    if i1>i2:
                        plt.subplot(len(trait_index_map)-1, len(trait_index_map)-1, 1+(i1-1)+(len(trait_index_map)-1)*i2)
                        v = venn2(subsets = (n1, n2, n12), set_labels = (traits14_short[t1], traits14_short[t2]))
                        v.get_patch_by_id('100').set_color(cm.colors[f(trait_index_map[t1])])
                        v.get_patch_by_id('010').set_color(cm.colors[f(trait_index_map[t2])])    
                        v.get_patch_by_id('110').set_color(cm.colors[7])   
                        format_numers = '{:.1f}\n({:.1f})' if trait_group != 'immuno' else '{:.2f}\n({:.2f})'
                        v.get_label_by_id('100').set_text(format_numers.format(n1, n1_se))
                        v.get_label_by_id('010').set_text(format_numers.format(n2, n2_se))
                        v.get_label_by_id('110').set_text(format_numers.format(n12, n12_se))
                    i1, i2 = i2, i1
                    t1, t2 = t2, t1
                    n1, n2 = n2, n1; n1_se, n2_se = n2_se, n1_se

            fig.subplots_adjust(hspace=0.10, wspace=0.10)
            plt.tight_layout()
            savefig(figures_folder, 'VENN_DIAGRAMS_{}_sym'.format(trait_group) +(' (censored)' if censored else ''))

if DO_VENN_DIAGRAMS_SUPPL:
    cm = plt.cm.get_cmap('tab10')
    for censored in [True, False]:
        bgmg_table_df = bgmg_table_df_censored[censored]

        venn_diag_suppl_tasks = []
        for trait1 in traits_ordered:
            venn_diag_suppl_tasks.append((traits_short[trait1], 17, traitPS_index_map, ([(trait1, trait2) for trait2 in traits_ordered if trait2 != trait1])))
        venn_diag_suppl_tasks.append(('immuno', 0.5, traitIM_index_map, list(itertools.combinations(traitsIM_ordered, 2))))
        venn_diag_suppl_tasks.append(('antro', 6, traitAM_index_map, list(itertools.combinations(traitsAM_ordered, 2))))

        for task_name, max_size, traitXX_index_map, venn_diag_suppl_sub_tasks in venn_diag_suppl_tasks:
            fig=plt.figure(figsize=(16, 28), dpi=80)
            for index, (trait1_orig, trait2_orig) in enumerate(venn_diag_suppl_sub_tasks):
                trait1 = trait1_orig; trait2=trait2_orig
                row = bgmg_table_df[bgmg_table_df.trait1_ID.isin([trait1, trait2]) & bgmg_table_df.trait2_ID.isin([trait1, trait2])].reset_index(drop=True).iloc[0]
                flip_data = (row.trait1_ID != trait1)

                n1 = row['nc1@p9']/1000; n1_se = row['nc1@p9_se']/1000
                n2 = row['nc2@p9']/1000; n2_se = row['nc2@p9_se']/1000;
                n12 = row['nc12@p9']/1000; n12_se = row['nc12@p9_se']/1000
                t1 = row.trait1_ID; t2 = row.trait2_ID
                if flip_data: t1, t2 = t2, t1;n1, n2 = n2, n1; n1_se, n2_se = n2_se, n1_se
                f = lambda x: x if x < 7 else x+1

                plt.subplot(len(traitPS_index_map)-1, 3, 3*index + 1)
                v = venn2(subsets = (n1, n2, n12), normalize_to=(n1+n2+n12)/max_size, set_labels = ("", ""))
                v.get_patch_by_id('100').set_color(cm.colors[f(traitXX_index_map[t1])])
                v.get_patch_by_id('010').set_color(cm.colors[f(traitXX_index_map[t2])])    
                v.get_patch_by_id('110').set_color(cm.colors[7])   
                formatter = '{:.1f}\n({:.1f})' if (task_name != 'immuno') else '{:.2f}\n({:.2f})'
                v.get_label_by_id('100').set_text(formatter.format(n1, n1_se))
                v.get_label_by_id('010').set_text(formatter.format(n2, n2_se))
                v.get_label_by_id('110').set_text(formatter.format(n12, n12_se))

                plt.xlim([-0.75, 0.75]), plt.ylim([-0.6, 0.6])
                if flip_data: plt.title(traits14_short[row.trait2_ID] +' & ' + traits14_short[row.trait1_ID])
                else: plt.title(traits14_short[row.trait1_ID] +' & ' + traits14_short[row.trait2_ID])

                fname1 = (format_string_BGMG_test_censored if censored else format_string_BGMG_test).format(trait1, trait2)
                fname2 = (format_string_BGMG_test_censored if censored else format_string_BGMG_test).format(trait2, trait1)
                if os.path.exists(fname2):
                    data = json_loads(fname2)
                    flip_data = True
                elif os.path.exists(fname1):
                    data = json_loads(fname1)
                    flip_data = False
                else:
                    print('missing: {} vs {}'.format(trait1, trait2))
                    continue

                for repeat in range(0,2):
                    plt.subplot(len(traitPS_index_map)-1, 3, 3*index + 2 + ((1-repeat) if flip_data else repeat))
                    sqq = data['result']['bivariate']['stratified_qq_plot_fit_data']
                    for index2, qq in enumerate(sqq['trait2' if repeat==0 else 'trait1']):
                        hData = plt.plot(qq['data_logpvec'], qq['hv_logp'], color=cm.colors[index2], linestyle='solid')
                    hNull = plt.plot(qq['hv_logp'], qq['hv_logp'], 'k--')
                    plt.legend(['All SNPs'] + '$P{0}\leq0.1$ $P{0}\leq0.01$ $P{0}\leq0.001$'.format('_{' + traits14_short[trait2]+'}').split(), loc='lower right', fontsize=14, borderpad=0.2, frameon=False, borderaxespad=0.2, labelspacing=0.1)
                    for index2, qq in enumerate(sqq['trait2' if repeat==0 else 'trait1']):
                        hModel = plt.plot(qq['model_logpvec'], qq['hv_logp'], color=cm.colors[index2], linestyle='dashed')
                    plt.ylim(0, 7.3); plt.xlim(0, 7.3)
                    if flip_data: plt.title('{} | {}'.format(traits14_short[trait2], traits14_short[trait1]))
                    else: plt.title('{} | {}'.format(traits14_short[trait1], traits14_short[trait2]))
                    trait1, trait2 = trait2, trait1  # swap traits
            fig.subplots_adjust(hspace=0.25, wspace=0.35)

            savefig(figures_folder, 'SF_VENN_and_stratQQ_{}'.format(task_name)  +(' (censored)' if censored else '') )


if DO_ADJ_TRAJECTORY_BGMG:
    fig = plt.figure(figsize=(2*16, 3*16), dpi=80)
    for index, (trait1, trait2) in enumerate(list(itertools.combinations(traits_ordered, 2))):
        fname1 = format_string_BGMG_fit.format(trait1, trait2)
        fname2 = format_string_BGMG_fit.format(trait2, trait1)
        if os.path.exists(fname2):
            data = json.loads(open(fname2).read())
            trait1, trait2 = trait2, trait1
        elif os.path.exists(fname1):
            data = json.loads(open(fname1).read())        
        else:
            print('missing: {} vs {}'.format(trait1, trait2))
            continue
        plt.subplot(6, 4, 1+index)

        lat = data['result']['bivariate']['loglike_adj_trajectory']
        fmtr = mpl.ticker.ScalarFormatter(useOffset=math.floor(min(lat['cost'])))
        plt.gca().get_yaxis().set_major_formatter(fmtr)

        pi12vector = [x[2] for x in lat['pivec']]
        plt.plot(pi12vector, lat['cost'], '.-')
        plt.xlabel('$\pi_{12}$')
        plt.ylabel('log likelihood')
        plt.locator_params(axis='x', nbins=6)

        params_pi12 = data['result']['bivariate']['params']['pi_vec'][2]
        params_cost = interp1d(pi12vector,lat['cost'])(max(pi12vector[0], min(pi12vector[-1], params_pi12)))
        plt.plot([params_pi12],[params_cost], '*k', markersize=8)

        plt.title('{} vs {}'.format(traits_short[trait1], traits_short[trait2]), loc='right')
    fig.subplots_adjust(hspace=0.35, wspace=0.25)
    savefig(figures_folder, 'ADJ_TRAJECTORY_BGMG')

if DO_STRATIFIED_QQ:
    fig = plt.figure(figsize=(18, 18), dpi=80)
    cm = plt.cm.get_cmap('tab10')
    censored = False
    
    for index, trait in enumerate(traits4_ordered):
        plt.subplot(4,4,1+index+4*index)
        fname = format_string_UGMG_test.format(trait)
        data = json_loads(fname)
        qq = data['result']['univariate'][0]['qq_plot_data']
        hData = plt.plot(qq['data_logpvec'], qq['hv_logp'], color=cm.colors[0], linestyle='solid')
        hModel = plt.plot(qq['model_logpvec'], qq['hv_logp'], color=cm.colors[0], linestyle='dashed')
        hNull = plt.plot(qq['hv_logp'], qq['hv_logp'], 'k--')
        plt.ylim(0, 7.3); plt.xlim(0, 7.3)
        plt.title(traits4_short[trait])
        plt.legend(['data', 'model', 'null'], loc='lower right', fontsize=14, borderpad=0.2, frameon=False, borderaxespad=0.2, labelspacing=0.1)

    for index, (trait1, trait2) in enumerate(list(itertools.combinations(traits4_ordered, 2))):
        fname1 = (format_string_BGMG_test_censored if censored else format_string_BGMG_test).format(trait1, trait2)
        fname2 = (format_string_BGMG_test_censored if censored else format_string_BGMG_test).format(trait2, trait1)
        if os.path.exists(fname2):
            data = json_loads(fname2)
            trait1, trait2 = trait2, trait1
        elif os.path.exists(fname1):
            data = json_loads(fname1)        
        else:
            print('missing: {} vs {}'.format(trait1, trait2))
            continue
        
        for repeat in range(0,2):
            ti1 = trait4_index_map[trait1]
            ti2 = trait4_index_map[trait2]
            plt.subplot(4,4,1 + ti1 + 4*ti2)

            sqq = data['result']['bivariate']['stratified_qq_plot_fit_data']
            for index, qq in enumerate(sqq['trait2' if repeat==0 else 'trait1']):
                hData = plt.plot(qq['data_logpvec'], qq['hv_logp'], color=cm.colors[index], linestyle='solid')
            hNull = plt.plot(qq['hv_logp'], qq['hv_logp'], 'k--')
            plt.legend(['All SNPs'] + '$P{0}\leq0.1$ $P{0}\leq0.01$ $P{0}\leq0.001$'.format('_{' + traits4_short[trait2]+'}').split(), loc='lower right', fontsize=14, borderpad=0.2, frameon=False, borderaxespad=0.2, labelspacing=0.1)
            for index, qq in enumerate(sqq['trait2' if repeat==0 else 'trait1']):
                hModel = plt.plot(qq['model_logpvec'], qq['hv_logp'], color=cm.colors[index], linestyle='dashed')
            plt.ylim(0, 7.3); plt.xlim(0, 7.3)
            plt.title('{} | {}'.format(traits4_short[trait1], traits4_short[trait2]))
            trait1, trait2 = trait2, trait1  # swap traits

    fig.subplots_adjust(hspace=0.35)
    fig.text(0.5, 0.075, 'Expected(-$log_{10}$ p-value)', ha='center')
    fig.text(0.075, 0.5, 'Observed(-$log_{10}$ p-value)', va='center', rotation='vertical')
   
    savefig(figures_folder, 'STRATIFIED_QQ')

if DO_SIMU_QQ:
    fig = plt.figure(figsize=(24, 32), dpi=80)
    cm = plt.cm.get_cmap('tab10')

    for h2_index, h2 in enumerate(['0.1', '0.4', '0.7']):
        for pi_index, pi in enumerate(['3.0000e-03', '3.0000e-04']):
            ylim_value = 0
            plt.subplot(3,2,1+2*h2_index+ pi_index)
            for repi in range(1, 5):
                for po_index, polyOverlap in enumerate([False, True]):
                    pi12=(pi if polyOverlap else '0.0000e+00')
                    fname2=r'C:\work\SIMU_BGMG_11pifrac_simple2\simu_h2={h2}_rg=0.0_pi1u={pi}_pi2u={pi}_pi12={pi12}_rep={rep}_tag1=customPolygenicOverlapAt{po}_tag2=evenPolygenicity.bgmg.test.json'.format(h2=h2,pi=pi,pi12=pi12, rep=repi, po=('1p0' if polyOverlap else '0p0'))
                    if not os.path.isfile(fname2): continue
                    data = json.loads(open(fname2).read())

                    for repeat in range(0,2): # T1|T2, T2|T1
                        sqq = data['result']['bivariate']['stratified_qq_plot_fit_data']
                        for index, qq in enumerate(sqq['trait2' if repeat==0 else 'trait1']):
                            hData = plt.plot(qq['data_logpvec'], qq['hv_logp'], color=cm.colors[index], linestyle='solid')
                            ylim_value = max(ylim_value, np.max([y for x, y in zip(qq['data_logpvec'], qq['hv_logp']) if np.isfinite(x)]))
                            break
                        hNull = plt.plot(qq['hv_logp'], qq['hv_logp'], 'k--')
                        for index, qq in enumerate(sqq['trait2' if repeat==0 else 'trait1']):
                            hModel = plt.plot(qq['model_logpvec'], qq['hv_logp'], color=cm.colors[index+1], linestyle='solid')
                            break
            plt.xlim(0, 7.3);#plt.ylim(0, 20); 
            plt.ylim(0, ylim_value*1.05 + 1); 
            if h2_index==2: plt.xlabel('Expected(-$log_{10}$ p-value)')
            if pi_index==0: plt.ylabel('Observed(-$log_{10}$ p-value)')
            plt.title('h2={}, pi1u={}'.format(h2, pi.replace('.0000e-0', 'e-')))
    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    savefig(figures_folder, 'SIMU_QQ')

if DO_SIMU_QQ_BINS:
    for pi_index, pi in enumerate(['3.0000e-03', '3.0000e-04']):
        fig = plt.figure(figsize=(18, 18), dpi=80)
        cm = plt.cm.get_cmap('tab10')

        for repi in range(1, 11):
            fname2=r'C:\work\SIMU_BGMG_11pifrac_simple2\simu_h2=0.4_rg=0.0_pi1u={pi}_pi2u={pi}_pi12=0.0000e+00_rep={rep}_tag1=customPolygenicOverlapAt0p0_tag2=evenPolygenicity.trait1.test.json'.format(pi=pi, rep=repi)
            if not os.path.isfile(fname2): continue
            data = json.loads(open(fname2).read())
            for subplot_index in range(0,9):
                axis = plt.subplot(3,3,1+subplot_index)
                qq = data['result']['univariate'][0]['qq_plot_bins_data'][subplot_index]
                hData = plt.plot(qq['data_logpvec'], qq['hv_logp'], color=cm.colors[0])
                hModel = plt.plot(qq['model_logpvec'], qq['hv_logp'], color=cm.colors[1])
                hNull = plt.plot(qq['hv_logp'], qq['hv_logp'], 'k--')

                #plt.title(str(subplot_index) + qq["options"]["title"].replace('$$', '$')) #.split('$$  $$')
                qq["options"]["title"] = qq["options"]["title"].replace('L \\in', 'LDscore \\in$\n$').replace('maf', 'MAF').replace(',', ', ').replace('Infinity', 'Inf')

                if True:
                    label_indices = {'y':[2, 5, 8], 'x':[0, 1, 2]}
                    axis_indices = {'y':[0, 3, 6], 'x':[6, 7, 8]}
                else:
                    label_indices = {'y':[0, 3, 6], 'x':[6, 7, 8]}
                    axis_indices = {'y':[2, 5, 8], 'x':[0, 1, 2]}

                if subplot_index in label_indices['y']: axis.yaxis.set_label_position('right'), plt.ylabel(qq["options"]["title"].split('$  $')[1].replace('$$', '$'))
                if subplot_index in axis_indices['y']: plt.ylabel('Observed(-$log_{10}$ p-value)')
                if subplot_index in label_indices['x']:  plt.gca().xaxis.set_label_position('top'), plt.xlabel(qq["options"]["title"].split('$  $')[0].replace('$$', '$'))
                if subplot_index in axis_indices['x']: plt.xlabel('Expected(-$log_{10}$ p-value)')
                plt.ylim(0, 80 if (pi_index==1) else 30); 
                plt.xlim(0, 7.3)

        fig.subplots_adjust(hspace=0.25, wspace=0.25)
        savefig(figures_folder, 'SIMU_QQ_BINS_{}'.format(pi))


if DO_SIMU_STRATIFIED_QQ:
    fig = plt.figure(figsize=(24, 32), dpi=80)
    cm = plt.cm.get_cmap('tab10')

    for repi in range(6, 7):
        for h2_index, h2 in enumerate(['0.1', '0.4', '0.7']):
            for pi_index, pi in enumerate(['3.0000e-03', '3.0000e-04']):
                for po_index, polyOverlap in enumerate([True, False]):
                    pi12=(pi if polyOverlap else '0.0000e+00')
                    fname2=r'C:\work\SIMU_BGMG_11pifrac_simple2\simu_h2={h2}_rg=0.0_pi1u={pi}_pi2u={pi}_pi12={pi12}_rep={rep}_tag1=customPolygenicOverlapAt{po}_tag2=evenPolygenicity.bgmg.test.json'.format(h2=h2,pi=pi,pi12=pi12, rep=repi, po=('1p0' if polyOverlap else '0p0'))
                    if not os.path.isfile(fname2): continue
                    data = json.loads(open(fname2).read())

                    for repeat in range(0,1): # T1|T2, T2|T1
                        plt.subplot(4,3,1+h2_index+ 3*(pi_index + 2*(1-polyOverlap)))
                        sqq = data['result']['bivariate']['stratified_qq_plot_fit_data']
                        ylim_value = 0
                        for index, qq in enumerate(sqq['trait2' if repeat==0 else 'trait1']):
                            hData = plt.plot(qq['data_logpvec'], qq['hv_logp'], color=cm.colors[index], linestyle='solid')
                            ylim_value = max(ylim_value, np.max([y for x, y in zip(qq['data_logpvec'], qq['hv_logp']) if np.isfinite(x)]))
                        hNull = plt.plot(qq['hv_logp'], qq['hv_logp'], 'k--')
                        for index, qq in enumerate(sqq['trait2' if repeat==0 else 'trait1']):
                            hModel = plt.plot(qq['model_logpvec'], qq['hv_logp'], color=cm.colors[index], linestyle='dashed')
                        plt.xlim(0, 7.3);plt.ylim(0, 20); 
                        if (pi_index==1) and (po_index==1): plt.xlabel('Expected(-$log_{10}$ p-value)')
                        if h2_index==0: plt.ylabel('Observed(-$log_{10}$ p-value)')
                        #plt.ylim(0, ylim_value + 1); 
                        plt.title('h2={}, pi1u={}, pi12={}'.format(h2, pi.replace('.0000e-0', 'e-'), pi12.replace('0.0000e+00', '0').replace('.0000e-0', 'e-')))
    fig.subplots_adjust(hspace=0.25, wspace=0.25)
    savefig(figures_folder, 'SIMU_STRATIFIED_QQ')

if DO_SIMU_UGMG:
    for DO_SPOW, figtag in [(False, ''), (True, '.spow')]:

        df_data = {}
        colnames = {'h2':'true_h2', 'rg':'true_rg', 'pi1u':'true_pi1u', 'pi2u':'true_pi2u', 'pi12':'true_pi12'}

        def insert_key_to_dictionary_as_list(key, value):
            if key not in df_data:
                df_data[key] = []
            df_data[key].append(value)

        if DO_SPOW:
            files = glob.glob(globpat_SIMU_BGMG_spow2) 
            figure_group_field = 'spow'
            figure_group_values = [-0.25, -0.5, -0.75]
        else:
            files = glob.glob(globpat_SIMU_BGMG_11pifrac)
            figure_group_field = 'true_h2'
            figure_group_values = [0.1, 0.4, 0.7]

        #if DO_GCTA: files = glob.glob(folder_SIMU_BGMG_GCTA+r'/*.trait*.fit.short.json')
        for fname in files:
            if 'h2=0.8' in fname: continue
            data = json_loads(fname)
            for trait in ['1', '2']:
                # general info about the run ['h2', 'rg', 'pi1u', 'pi2u', 'pi12', 'rep', 'tag1', 'tag2', 'outtag']
                rep = {'.bgmg':'', '.short':'', '.json':'', 'outtag=run1_outtag':'outtag'}
                for k in rep: fname = fname.replace(k, rep[k])

                for key, value in [tuple(x.split('=')) for x in os.path.basename(fname).split('_') if ('=' in x)]:
                    insert_key_to_dictionary_as_list(colnames[key] if key in key in colnames else key, value)

                if 'spow' not in fname:
                    insert_key_to_dictionary_as_list('spow', '0.0')

                # polygenicity and heritability estimates
                ci = data['result']['bivariate']['ci']
                insert_key_to_dictionary_as_list('h2', ci['h2_T{}'.format(trait)]['point_estimate'])
                insert_key_to_dictionary_as_list('h2_se', ci['h2_T{}'.format(trait)]['se'])
                insert_key_to_dictionary_as_list('pi_vec', ci['pi{}u'.format(trait)]['point_estimate'])
                insert_key_to_dictionary_as_list('pi_vec_se', ci['pi{}u'.format(trait)]['se'])

        df = pd.DataFrame(df_data)
        colorder = ['tag1', 'tag2', 'outtag', 'true_h2', 'true_rg', 'true_pi1u', 'true_pi2u', 'true_pi12', 'rep']
        df = df[[c for c in colorder if c in df] + [x for x in df if x not in colorder]]
        df = df[df['pi_vec']<0.95].copy() # drop 1 outlier

        df_agg = df[['true_pi1u', 'true_h2', 'spow', 'pi_vec', 'pi_vec_se', 'h2', 'h2_se']].groupby(['true_pi1u', 'true_h2', 'spow']).agg({'pi_vec':['mean', 'std'], 'pi_vec_se':['mean'], 'h2':['mean', 'std'], 'h2_se':['mean']})
        df_agg.columns = ['_'.join(col).strip() for col in df_agg.columns.values]
        df_agg = df_agg.reset_index()

        # https://python-graph-gallery.com/8-add-confidence-interval-on-barplot/
        # width of the bars
        barWidth = 0.3

        fig = plt.figure(figsize=(12, 16), dpi=80)

        if DO_SPOW:
            figure_group_vec = [('-0.25', '3.0000e-03'), ('-0.5', '3.0000e-03'), ('-0.75', '3.0000e-03'), ('-0.25', '3.0000e-04'), ('-0.5', '3.0000e-04'), ('-0.75', '3.0000e-04')]
        else:
            figure_group_vec = [('0.1', '3.0000e-03'), ('0.4', '3.0000e-03'), ('0.7', '3.0000e-03'), ('0.1', '3.0000e-04'), ('0.4', '3.0000e-04'), ('0.7', '3.0000e-04')]
        
        for index in range(0, 2):
            plt.subplot(2,1,1+index)
            bars1 = df_agg['h2_mean' if index==0 else 'pi_vec_mean'].values 
            bars2 = df_agg['true_h2' if index==0 else 'true_pi1u'].astype(float).values
            yer1 = df_agg['h2_se_mean' if index==0 else 'pi_vec_se_mean'].values
            yer2 = df_agg['h2_std' if index==0 else 'pi_vec_std']

            r1 = np.arange(len(bars1)); r2 = [x + barWidth for x in r1] # # The x position of bars

            plt.bar(r1, bars1, width = barWidth, color = 'blue', edgecolor = 'black', yerr=yer2, capsize=7, label='estimated')
            plt.bar(r2, bars2, width = barWidth, color = 'cyan', edgecolor = 'black',            capsize=7, label='true')
            if index==1: plt.gca().set_yscale('log')

            plt.xticks([r + barWidth for r in range(len(bars1))], ['{}={}\n$\pi_1$={:.0e}'.format('S' if DO_SPOW else 'h2', h2_or_spow, float(pi)).replace('e-0', 'e-') for h2_or_spow, pi in figure_group_vec])
            plt.ylabel('h2' if index==0 else '$\pi_1$' )
            plt.legend()
            if index==1: plt.ylim([1e-4, 1e-2])

            fig.subplots_adjust(hspace=0.35)
            savefig(figures_folder, 'SIMU_UGMG_pi_and_h2' + ('.spow' if DO_SPOW else ''))
                   
table_index += 1
if DO_SIMU_UGMG and DO_SIMU_UGMG_TABLE:
    df_agg['true_h2'] = df_agg['true_h2'].astype(float)
    df_agg['true_pi1u'] = df_agg['true_pi1u'].astype(float)
    df_agg.to_csv(os.path.join(tables_folder, 'SIMU_UGMG_TABLE.csv'), sep='\t',index=False)
    
    if tables_writer is None: tables_writer = pd.ExcelWriter(os.path.join(tables_folder, 'all_tables.xlsx'))
    df=df_agg.copy()
    for c in [a for a in df.columns if 'pi_vec' in a]: df[c] = ['{:.2e}'.format(x) for x in df[c]]
    for c in [a for a in df.columns if 'h2' in a]: df[c] = ['{:.3f}'.format(float(x)) for x in df[c]]
    concat_se_mean_std(df).to_excel(tables_writer,'table{}'.format(table_index), index=False)

plt.rcParams.update({'mathtext.default': 'regular', 'font.size': 16 })

table_index += 3
if DO_SIMU_BGMG or DO_SIMU_BGMG_TABLE:
    for DO_SPOW, figtag in [(False, 'h2'), (True, 'spow')]:
        df_data = {}
        colnames = {'h2':'true_h2', 'rg':'true_rg', 'pi1u':'true_pi1u', 'pi2u':'true_pi2u', 'pi12':'true_pi12'}
        
        def insert_key_to_dictionary_as_list(key, value):
            if key not in df_data:
                df_data[key] = []
            df_data[key].append(value)

        if DO_SPOW:
            files = glob.glob(globpat_SIMU_BGMG_spow2) 
            figure_group_field = 'spow'
            figure_group_values = [-0.25, -0.5, -0.75]
        else:
            files = glob.glob(globpat_SIMU_BGMG_11pifrac) + glob.glob(globpat_SIMU_BGMG_11pifrac_wave2)
            figure_group_field = 'true_h2'
            figure_group_values = [0.1, 0.4, 0.7]
        #if DO_GCTA: files = glob.glob(folder_SIMU_BGMG_GCTA + r'\*.bgmg.fit.short.json')

        for fname in files:
            if 'h2=0.8' in fname: continue
            data = json_loads(fname)

            # general info about the run ['h2', 'rg', 'pi1u', 'pi2u', 'pi12', 'rep', 'tag1', 'tag2', 'outtag']
            rep = {'.bgmg':'', '.short':'', '.json':'', 'outtag=run1_outtag':'outtag'}
            for k in rep: fname = fname.replace(k, rep[k])
            for key, value in [tuple(x.split('=')) for x in os.path.basename(fname).split('_') if ('=' in x)]:
                insert_key_to_dictionary_as_list(colnames[key] if key in key in colnames else key, value)
            if 'spow' not in fname:
                insert_key_to_dictionary_as_list('spow', '0.0')
                
            # polygenicity and heritability estimates
            ci = data['result']['bivariate']['ci']
            extract = [('pi12', 'pi_vec_C3'), ('pi1u', 'pi1u'), ('pi2u', 'pi2u'), ('rho12', 'rho_beta'), ('rg', 'rg')]
            for a,b in extract:
                insert_key_to_dictionary_as_list(a+'', ci[b]['point_estimate'])
                insert_key_to_dictionary_as_list(a+'_se', ci[b]['se'])

        df = pd.DataFrame(df_data)
        colorder = ['tag1', 'tag2', 'outtag', 'true_h2', 'spow', 'true_rg', 'true_pi1u', 'true_pi2u', 'true_pi12', 'rep']
        df = df[[c for c in colorder if c in df] + [x for x in df if x not in colorder]]
        df = df[(df['pi1u']<0.95) & (df['pi2u']<0.95)].copy() # drop 1 outlier

        df_agg = df.groupby(['true_pi1u', 'true_pi2u', 'true_pi12' , 'true_h2', 'true_rg', 'spow']).agg({'pi12':['mean', 'std'], 'pi12_se':['mean'], 'rg':['mean', 'std'], 'rg_se':['mean'], 'rho12':['mean', 'std'], 'rho12_se':['mean']})
        df_agg.columns = ['_'.join(col).strip() for col in df_agg.columns.values]
        df_agg = df_agg.reset_index()
        for col in ['true_pi1u','true_pi2u','true_pi12','true_h2', 'spow']: df_agg[col] = df_agg[col].astype(float)

        if DO_SIMU_BGMG and not DO_SPOW:
            fig = plt.figure(figsize=(24, 24), dpi=80)
            for true_rg_index, true_rg in enumerate(['0.0', '0.5']):
                for index, (pi1u, pi2u) in enumerate([(3e-03, 3e-03), (3e-4, 3e-04), (3e-03, 3e-04)]):
                    scale_coef, scale_text = (1e4, '10^{-4}') if pi2u==3e-4 else (1e3, '10^{-3}')
                    plt.subplot(3, 2, 2*index+1+true_rg_index)
                    for h2_index in [0.1, 0.4, 0.7]:
                        df_plot = df_agg[(df_agg.true_rg==true_rg) & (df_agg.true_pi1u == pi1u) & (df_agg.true_pi2u == pi2u) & (df_agg.true_h2 == h2_index)].sort_values('true_pi12')
                        plt.plot(df_plot.true_pi12.values*scale_coef, df_plot.pi12_mean.values*scale_coef, '.-')
                    plt.plot(df_plot.true_pi12.values*scale_coef, df_plot.true_pi12.values*scale_coef, 'k--')
                    plt.ylabel('estimated $\pi_{12}\;(\\times ' + scale_text + '$)' )  # @10k
                    plt.xlabel(     'true $\pi_{12}\;(\\times ' + scale_text + '$)' )
                    plt.legend(['h2=0.1', 'h2=0.4', 'h2=0.7', 'true'])
                    plt.title('$\pi_1$ = {:.1e}, $\pi_2$ = {:.1e}'.format(pi1u, pi2u) + ', $\\rho_{12}$=' + true_rg)
                    plt.locator_params(axis='x', nbins=7)
                    plt.ticklabel_format(style='sci', axis='both')
            fig.subplots_adjust(hspace=0.30, wspace=0.50)
            savefig(figures_folder, 'SIMU_BGMG_{}_rg={}_old_style_pi12'.format(figtag, true_rg))

        for true_rg in ['0.0', '0.5']:
            if not DO_SIMU_BGMG: break
            fig = plt.figure(figsize=(18, 16), dpi=80)
            for index, (pi1u, pi2u) in enumerate([(3e-03, 3e-03), (3e-4, 3e-04), (3e-03, 3e-04)]):
                if DO_SPOW and (pi1u!=pi2u): continue
                for jndex, h2_index in enumerate(figure_group_values):
                    df_plot = df_agg[(df_agg.true_rg==true_rg) & (df_agg.true_pi1u == pi1u) & (df_agg.true_pi2u == pi2u) & (df_agg[figure_group_field] == h2_index)].sort_values('true_pi12')
                    if df_plot.empty: continue
                    ax = plt.subplot(3,3,index+3*jndex+1)
                    plot_simu_bgmg_pi12(df_plot)

            fig.subplots_adjust(hspace=0.50, wspace=0.25)
            savefig(figures_folder, '{}_{}_pi12_rg={}'.format('SIMU_BGMG', figtag, true_rg))

        for do_rg in [True, False]:  # rg or rho_beta?
            if not DO_SIMU_BGMG: break
            for true_rg in ['0.0', '0.5']:
                fig = plt.figure(figsize=(18, 16), dpi=80)
                for index, (pi1u, pi2u) in enumerate([(3e-03, 3e-03), (3e-4, 3e-04), (3e-03, 3e-04)]):
                    if DO_SPOW and (pi1u!=pi2u): continue
                    for jndex, h2_index in enumerate(figure_group_values):
                        df_plot = df_agg[(df_agg.true_rg==true_rg) & (df_agg.true_pi1u == pi1u) & (df_agg.true_pi2u == pi2u) & (df_agg[figure_group_field] == h2_index)].sort_values('true_pi12')
                        if df_plot.empty: continue
                        plt.subplot(3,3,index+3*jndex+1)
                        plot_simu_bgmg_rg_or_rho12(df_plot, do_rg)

                fig.subplots_adjust(hspace=0.30, wspace=0.50)                    
                savefig(figures_folder, '{}_{}_{}_rg={}.'.format('SIMU_BGMG', figtag, 'rg' if do_rg else 'rho12', true_rg))

        if not DO_SPOW:
            # aggregate together selected simulations for the main text figure
            true_rg = '0.5'; pi1u=3e-04; pi2u=3e-04; h2=0.4;
            df_plot = df_agg[(df_agg.true_rg==true_rg) & (df_agg.true_pi1u == pi1u) & (df_agg.true_pi2u == pi2u) & (df_agg.true_h2 == h2_index)].sort_values('true_pi12')

            fig = plt.figure(figsize=(18, 4.5), dpi=80)
            ax = plt.subplot(1,3,1)
            plot_simu_bgmg_pi12(df_plot, do_title=False)
            plt.title('A', loc='left')

            ax = plt.subplot(1,3,2)   
            plot_simu_bgmg_rg_or_rho12(df_plot, do_rg=False, do_title=False)
            plt.title('B', loc='left')

            ax = plt.subplot(1,3,3)   
            plot_simu_bgmg_rg_or_rho12(df_plot, do_rg=True, do_title=False)
            plt.title('C', loc='left')

            fig.subplots_adjust(hspace=0.30, wspace=0.50)                    
            savefig(figures_folder, 'SIMU_BGMG_main_figure2')

        if DO_SIMU_BGMG_TABLE:
            df_agg.to_csv(os.path.join(tables_folder, 'SIMU_BGMG_TABLE_{}.csv'.format(figtag)), sep='\t',index=False)
            if tables_writer is None: tables_writer = pd.ExcelWriter(os.path.join(tables_folder, 'all_tables.xlsx'))
            df=df_agg.copy()
            for c in [a for a in df.columns if 'pi' in a]: df[c] = ['{:.2e}'.format(x) for x in df[c]]
            for c in [a for a in df.columns if 'rg' in a]: df[c] = ['{:.3f}'.format(float(x)) for x in df[c]]
            for c in [a for a in df.columns if 'rho12' in a]: df[c] = ['{:.3f}'.format(float(x)) for x in df[c]]
            df = concat_se_mean_std(df)
            df[(df['true_pi1u'] == '3.00e-03') & (df['true_pi2u'] == '3.00e-03')].to_excel(tables_writer,'table{}'.format(table_index-2), index=False)
            df[(df['true_pi1u'] == '3.00e-04') & (df['true_pi2u'] == '3.00e-04')].to_excel(tables_writer,'table{}'.format(table_index-1), index=False)
            df[(df['true_pi1u'] == '3.00e-03') & (df['true_pi2u'] == '3.00e-04')].to_excel(tables_writer,'table{}'.format(table_index-0), index=False)

    table_index += 1
    if DO_GWAS_DATA_TABLE:
        if tables_writer is None: tables_writer = pd.ExcelWriter(os.path.join(tables_folder, 'all_tables.xlsx'))
        df_gwas_data = pd.read_table('gwas_data.csv',sep='\t')
        del df_gwas_data['ID']
        df_gwas_data.to_excel(tables_writer,'table{}'.format(table_index), index=False)

if DO_SIMU_UGMG_ANNOTENRICH:
    for mafSuffix in ['bgmgMAF']: # ['']
        df_data = {}
        colnames = {'h2':'true_h2', 'rg':'true_rg', 'pi1u':'true_pi1u', 'pi2u':'true_pi2u', 'pi12':'true_pi12'}

        def insert_key_to_dictionary_as_list(key, value):
            if key not in df_data:
                df_data[key] = []
            df_data[key].append(value)

        files = glob.glob(folder_SIMU_BGMG_annotenrich + '\\*' + mafSuffix + '*bgmg.short.json')
        for fname in files:
            if (mafSuffix=='bgmgMAF') and ('300K' in fname): continue
            data = json.loads(open(fname).read())

            for trait_index in range(0, 2):
                # general info about the run ['h2', 'rg', 'pi1u', 'pi2u', 'pi12', 'rep', 'tag1', 'tag2', 'outtag']
                rep = {'.bgmg':'', '.short':'', '.json':''}
                for k in rep: fname = fname.replace(k, rep[k])
                for key, value in [tuple(x.split('=')) for x in os.path.basename(fname).split('_') if ('=' in x)]:
                    insert_key_to_dictionary_as_list(colnames[key] if key in key in colnames else key, value)

                # polygenicity and heritability estimates
                ci = data['result']['univariate'][trait_index]['ci']
                insert_key_to_dictionary_as_list('h2', ci['h2']['point_estimate'])
                insert_key_to_dictionary_as_list('h2_se', ci['h2']['se'])
                insert_key_to_dictionary_as_list('pi_vec', ci['pi_vec']['point_estimate'])
                insert_key_to_dictionary_as_list('pi_vec_se', ci['pi_vec']['se'])

        df = pd.DataFrame(df_data)
        colorder = ['tag1', 'tag2', 'outtag', 'true_h2', 'true_rg', 'true_pi1u', 'true_pi2u', 'true_pi12', 'rep']
        df = df[[c for c in colorder if c in df] + [x for x in df if x not in colorder]]
        df = df[df['pi_vec']<0.95].copy() # drop 1 outlier
        

        df_agg = df[['piTag', 'true_h2', 'pi_vec', 'pi_vec_se', 'h2', 'h2_se']].groupby(['piTag', 'true_h2']).agg({'pi_vec':['mean', 'std'], 'pi_vec_se':['mean'], 'h2':['mean', 'std'], 'h2_se':['mean']})
        df_agg.columns = ['_'.join(col).strip() for col in df_agg.columns.values]
        df_agg = df_agg.reset_index()
        df_agg['true_pi1u'] = df_agg['piTag'].str.replace('K', '').astype(float)*1e3/11015833
        # https://python-graph-gallery.com/8-add-confidence-interval-on-barplot/
        # width of the bars
        barWidth = 0.3

        fig = plt.figure(figsize=(18, 16), dpi=80)

        h2_pi_vec = [('0.1', '300K'), ('0.4', '300K'), ('0.7', '300K'), ('0.1', '30K'), ('0.4', '30K'), ('0.7', '30K'), ('0.1', '3K'), ('0.4', '3K'), ('0.7', '3K')]
        if mafSuffix=='bgmgMAF': h2_pi_vec = [('0.1', '30K'), ('0.4', '30K'), ('0.7', '30K'), ('0.1', '3K'), ('0.4', '3K'), ('0.7', '3K')]
        for index in range(0, 2):
            plt.subplot(2,1,1+index)
            bars1 = df_agg['h2_mean' if index==0 else 'pi_vec_mean'].values 
            bars2 = df_agg['true_h2' if index==0 else 'true_pi1u'].astype(float).values
            yer1 = df_agg['h2_se_mean' if index==0 else 'pi_vec_se_mean'].values
            yer2 = df_agg['h2_std' if index==0 else 'pi_vec_std']

            r1 = np.arange(len(bars1)); r2 = [x + barWidth for x in r1] # # The x position of bars

            plt.bar(r1, bars1, width = barWidth, color = 'blue', edgecolor = 'black', yerr=yer2, capsize=7, label='estimated')
            plt.bar(r2, bars2, width = barWidth, color = 'cyan', edgecolor = 'black',            capsize=7, label='true')
            if index==1: plt.gca().set_yscale('log')

            plt.xticks([r + barWidth for r in range(len(bars1))], ['h2={}\n$nc_1$={}'.format(h2, piTag).replace('e-0', 'e-') for h2, piTag in h2_pi_vec])
            plt.ylabel('h2' if index==0 else '$\pi_1$' )
            if (index==1) and (mafSuffix=='bgmgMAF'): plt.ylim([1e-4, 1e-2])

            plt.legend()

            fig.subplots_adjust(hspace=0.35)
            savefig(figures_folder, 'SIMU_UGMG_ANNOTENRICH{}_pi_and_h2'.format(mafSuffix))

for mafSuffix in ['_bgmgMAF', '']:
    table_index += 1
    if DO_SIMU_BGMG_ANNOTENRICH_TABLE:
        df_data = {}
        colnames = {'h2':'true_h2', 'rg':'true_rg', 'pi1u':'true_pi1u', 'pi2u':'true_pi2u', 'pi12':'true_pi12'}

        df={'rep':[], 'piTag':[], 'true_n12':[]}
        for rep in range(1, 21):
            for poly in ['3K', '30K', '300K']:
                if (mafSuffix == '_bgmgMAF') and (poly=='300K'): continue
                data = [open('C:\work\SIMU_BGMG_annotenrich\simu{}_{}_enriched.rep={}.{}.snps'.format(mafSuffix, poly, rep, trait), 'r').read().split() for trait in ['trait1', 'trait2']]
                df['true_n12'].append(len(set(data[0]).intersection(set(data[1]))))
                df['piTag'].append(poly)
                df['rep'].append(rep)
        df_annotenrich_true_n12 = pd.DataFrame(df)

        def insert_key_to_dictionary_as_list(key, value):
            if key not in df_data:
                df_data[key] = []
            df_data[key].append(value)

        files = glob.glob(folder_SIMU_BGMG_annotenrich + '\\*' +mafSuffix+'*bgmg.short.json')
        for fname in files:
            if (mafSuffix=='_bgmgMAF') and ('300K' in fname): continue
            data = json.loads(open(fname).read())

            # general info about the run ['h2', 'rg', 'pi1u', 'pi2u', 'pi12', 'rep', 'tag1', 'tag2', 'outtag']
            rep = {'.bgmg':'', '.short':'', '.json':''}
            for k in rep: fname = fname.replace(k, rep[k])
            repi = None; piTag = None;
            for key, value in [tuple(x.split('=')) for x in os.path.basename(fname).split('_') if ('=' in x)]:
                insert_key_to_dictionary_as_list(colnames[key] if key in key in colnames else key, value)

            # polygenicity and heritability estimates
            ci = data['result']['bivariate']['ci']
            insert_key_to_dictionary_as_list('pi12', ci['pi_vec_C3']['point_estimate'])
            insert_key_to_dictionary_as_list('pi12_se', ci['pi_vec_C3']['se'])
            insert_key_to_dictionary_as_list('pi1u', ci['pi1u']['point_estimate'])
            insert_key_to_dictionary_as_list('pi1u_se', ci['pi1u']['se'])
            insert_key_to_dictionary_as_list('pi2u', ci['pi2u']['point_estimate'])
            insert_key_to_dictionary_as_list('pi2u_se', ci['pi2u']['se'])
            insert_key_to_dictionary_as_list('pifrac', ci['pi12_over_pi1u']['point_estimate'])
            insert_key_to_dictionary_as_list('pifrac_se', ci['pi12_over_pi1u']['se'])
            insert_key_to_dictionary_as_list('h2', ci['h2_T1']['point_estimate'])
            insert_key_to_dictionary_as_list('h2_se', ci['h2_T1']['se'])

        df = pd.DataFrame(df_data)
        colorder = ['tag1', 'tag2', 'outtag', 'true_h2', 'true_rg', 'true_pi1u', 'true_pi2u', 'true_pi12', 'rep']
        df = df[[c for c in colorder if c in df] + [x for x in df if x not in colorder]]
        df = df[(df['pi1u']<0.95) & (df['pi2u']<0.95)].copy() # drop 1 outlier

        df_annotenrich_true_n12['rep']=df_annotenrich_true_n12['rep'].astype(str)
        df = pd.merge(df, df_annotenrich_true_n12, how='left',on=['piTag', 'rep'])

        df_agg = df[['piTag', 'true_h2', 'pi12', 'pi12_se', 'pi1u', 'pi1u_se', 'pifrac', 'pifrac_se', 'h2', 'h2_se', 'true_n12']].groupby(['piTag', 'true_h2']).agg({'pi12':['mean', 'std'], 'pi12_se':['mean'], 'pi1u':['mean', 'std'], 'pi1u_se':['mean'], 'pifrac':['mean', 'std'], 'pifrac_se':['mean'], 'h2':['mean', 'std'], 'h2_se':['mean'], 'true_n12':['mean','std']})
        df_agg.columns = ['_'.join(col).strip() for col in df_agg.columns.values]
        df_agg = df_agg.reset_index()
        df_agg['true_pi1u'] = df_agg['piTag'].str.replace('K', '').astype(float)*1e3/11015833

        #	b Standard	errors	(the	first	number	is	the	mean	of	theoretical	standard	errors	derived	from	variance	formulas	and	the	second	one	is	the	empirical	standard	errors	across	100	replications).
        koef=1e-03;  # if change koef one must also change the logic that calculates expected true_n1U. Currently it is taken from piTag.
        df_agg['N12 (x1000)'] = ['{:.2f} ({:.2f}/{:.2f})'.format(v*koef*nsnps_HAPGEN,se*koef*nsnps_HAPGEN,sd*koef*nsnps_HAPGEN) for v,se,sd in zip(df_agg['pi12_mean'],df_agg['pi12_se_mean'],df_agg['pi12_std'])]
        df_agg['N1u (x1000)'] = ['{:.2f} ({:.2f}/{:.2f})'.format(v*koef*nsnps_HAPGEN,se*koef*nsnps_HAPGEN,sd*koef*nsnps_HAPGEN) for v,se,sd in zip(df_agg['pi1u_mean'],df_agg['pi1u_se_mean'],df_agg['pi1u_std'])]
        df_agg['h2']   = ['{:.2f} ({:.2f}/{:.2f})'.format(v,se,sd) for v,se,sd in zip(df_agg['h2_mean'],df_agg['h2_se_mean'],df_agg['h2_std'])]
        df_agg['N12/N1u']   = ['{:.2f} ({:.2f}/{:.2f})'.format(v,se,sd) for v,se,sd in zip(df_agg['pifrac_mean'],df_agg['pifrac_se_mean'],df_agg['pifrac_std'])]
        df_agg['true_n12']   = ['{:.2f} ({:.2f})'.format(v*koef,se*koef) for v,se in zip(df_agg['true_n12_mean'],df_agg['true_n12_std'])]
        df_agg['true_N1u (x1000)'] = [x.replace('K','') for x in df_agg['piTag']]
        df_agg = df_agg[['true_N1u (x1000)', 'true_h2', 'true_n12', 'N1u (x1000)', 'N12 (x1000)', 'N12/N1u','h2']].copy()

        df_agg.to_csv(os.path.join(tables_folder, 'SIMU_BGMG_ANNOTENRICH{}_TABLE.csv'.format(mafSuffix)), sep='\t',index=False)

        if tables_writer is None: tables_writer = pd.ExcelWriter(os.path.join(tables_folder, 'all_tables.xlsx'))
        df_agg.to_excel(tables_writer,'table{}'.format(table_index), index=False)

plt.rcParams.update({'mathtext.default': 'regular', 'font.size': 20 })

table_index += 1
if DO_SIMU_UGMG_SUBREF:
    df_data = {}
    colnames = {'h2':'true_h2', 'rg':'true_rg', 'pi1u':'true_pi1u', 'pi2u':'true_pi2u', 'pi12':'true_pi12'}

    def insert_key_to_dictionary_as_list(key, value):
        if key not in df_data:
            df_data[key] = []
        df_data[key].append(value)

    files = glob.glob(folder_SIMU_BGMG_subref + r'\*run2.fit.short.json')
    for fname in files:
        data = json.loads(open(fname).read())
        # general info about the run ['h2', 'rg', 'pi1u', 'pi2u', 'pi12', 'rep', 'tag1', 'tag2', 'outtag']
        rep = {'.bgmg':'', '.short':'', '.json':''}
        for k in rep: fname = fname.replace(k, rep[k])
        for key, value in [tuple(x.split('=')) for x in os.path.basename(fname).split('_') if ('=' in x)]:
            insert_key_to_dictionary_as_list(colnames[key] if key in key in colnames else key, value)

        # polygenicity and heritability estimates
        ci = data['result']['univariate'][0]['ci']
        insert_key_to_dictionary_as_list('h2', ci['h2']['point_estimate'])
        insert_key_to_dictionary_as_list('h2_se', ci['h2']['se'])
        insert_key_to_dictionary_as_list('pi_vec', ci['pi_vec']['point_estimate'])
        insert_key_to_dictionary_as_list('pi_vec_se', ci['pi_vec']['se'])
        insert_key_to_dictionary_as_list('sig2_beta', ci['sig2_beta']['point_estimate'])
        insert_key_to_dictionary_as_list('sig2_beta_se', ci['sig2_beta']['se'])

    df = pd.DataFrame(df_data)
    for col in df.columns: df[col] = pd.to_numeric(df[col],errors='ignore')
    df['hat_ncausal'] = nsnps_HAPGEN*df['frac'].values*df['pi_vec'].values
    df['hat_ncausal_se'] = nsnps_HAPGEN*df['frac'].values*df['pi_vec_se'].values
    df['true_ncausal'] = nsnps_HAPGEN*df['true_pi1u']
    for col in ['h2', 'h2_se']: df[col] = ['{:.3f}'.format(x) for x in df[col].values]
    for col in ['pi_vec', 'pi_vec_se']: df[col] = ['{:.3e}'.format(x) for x in df[col].values]
    for col in ['sig2_beta', 'sig2_beta_se']: df[col] = ['{:.3e}'.format(x) for x in df[col].values]
    for col in ['hat_ncausal', 'hat_ncausal_se', 'true_ncausal']: df[col] = ['{:.1f}'.format(x/1000) for x in df[col].values]
    df_agg = concat_se(df)[['true_h2', 'true_ncausal','frac', 'h2 (se)', 'hat_ncausal (se)', 'sig2_beta (se)', 'pi_vec (se)']].copy()

    df_agg.to_csv(os.path.join(tables_folder, 'SIMU_UGMG_SUBREF_TABLE.csv'), sep='\t',index=False)
    if tables_writer is None: tables_writer = pd.ExcelWriter(os.path.vjoin(tables_folder, 'all_tables.xlsx'))
    df_agg.to_excel(tables_writer,'table{}'.format(table_index), index=False)

if tables_writer: tables_writer.save()

    