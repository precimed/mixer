import os
import os.path
import sys
import glob
import itertools

dry_run=False
num_submit = int(sys.argv[1])

#out_prefix = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/plsa_mixer/real_run4_ukb/'
#out_prefix = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/plsa_mixer/mixer_storm/'
out_prefix = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/plsa_mixer/mixer_without_hm3/'
#out_prefix = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/alexeas_power_test/'

ugmg_pattern = out_prefix+'{}.fit.rep{}.json'
bgmg_pattern = out_prefix+'{}_vs_{}.fit.rep{}.json'

g0 = ['SSGAC_EDU_2018_no23andMe', 'PGC_SCZ_2014_EUR', 'PGC_BIP_2016', 'GIANT_HEIGHT_2018_UKB'] # 'PGC_MDD_2018_with23andMe', 'PGC_ADHD_2017_EUR', 'PGC_ASD_2017_iPSYCH',  
g1 = ['CTG_COG_2018', 'SSGAC_EDU_2018_no23andMe']
g2 = ['PGC_SCZ_2014_EUR', 'PGC_SCZ_0518_EUR', 'PGC_BIP_2016', 'PGC_MDD_2018_Howard_no23andMe', 'PGC_ADHD_2017_EUR', 'PGC_ASD_2017_iPSYCH']
g3 = ['GSCAN_DRINK_2019_DrinksPerWeek', 'CTG_NEUR_2018_no23andMe', 'ICC_CANNABIS_2018_UKB', 'UKB_LONELY_2018_Loneliness', 'SSGAC_RISK_2019',  'UKB_CHRONOTYPE_2018', 'UKB_SLEEP_2018', 'CTG_INSOMNIA_2018']
g3r =['GSCAN_DRINK_2019_DrinksPerWeek',                            'ICC_CANNABIS_2018_UKB', 'UKB_LONELY_2018_Loneliness', 'SSGAC_RISK_2019']
g5 = ['PGC_SCZ_2014_EUR', 'PGC_SCZ_0518_EUR', 'PGC_BIP_2016', 'PGC_MDD_2018_with23andMe', 'PGC_MDD_2018_Howard_with23andMe']
g5r= ['PGC_SCZ_2014_EUR',                     'PGC_BIP_2016',                             'PGC_MDD_2018_Howard_with23andMe']
g6 = ['UKB_CHRONOTYPE_2018', 'UKB_SLEEP_2018', 'CTG_INSOMNIA_2018']
g7 = ['UKB_HEIGHT_2018_irnt']
g8 = ['IHGC_MIG_2016_with23andMe', 'PGC_SCZ_0518_EUR', 'CLOZUK_SCZ_2018_withPGC', 'PGC_BIP_2016','PGC_MDD_2018_Howard_with23andMe', 'PGC_MDD_2018_with23andMe', 'PGC_ADHD_2017_EUR']
g9 = ['PGC_SCZ_2014_EUR', 'PGC_BIP_2016', 'PGC_MDD_2018_Howard_no23andMe', 'PGC_ASD_2017_iPSYCH', 'PGC_ADHD_2017_EUR', 'CTG_COG_2018', 'SSGAC_EDU_2018_no23andMe', 
      'CTG_NEUR_2018_no23andMe', 'UKB_LONELY_2018_Loneliness',
      'ICC_CANNABIS_2018_UKB', 'GSCAN_SMOKE_2019_CigarettesPerDay', 'GSCAN_DRINK_2019_DrinksPerWeek', 'EGG_BIRTHWEIGHT_2016', 'GIANT_HEIGHT_2018_UKB',
      'CARDIOGRAM_CAD_2015', 'GIANT_BMI_2018_UKB_v2', 'DIAGRAM_T2D_2017'] 
g9b=[ 'PGC_SCZ_2014_EUR', 'PGC_BIP_2016', 'PGC_MDD_2018_Howard_no23andMe', 'PGC_ASD_2017_iPSYCH', 'PGC_ADHD_2017_EUR', 'CTG_COG_2018', 'SSGAC_EDU_2018_no23andMe', ]
g9r=[ 'CTG_NEUR_2018_no23andMe', 'UKB_LONELY_2018_Loneliness',
      'ICC_CANNABIS_2018_UKB', 'GSCAN_SMOKE_2019_CigarettesPerDay', 'GSCAN_DRINK_2019_DrinksPerWeek', 'EGG_BIRTHWEIGHT_2016', 'GIANT_HEIGHT_2018_UKB',
      'CARDIOGRAM_CAD_2015', 'GIANT_BMI_2018_UKB_v2', 'DIAGRAM_T2D_2017'] 
g10=['PGC_SCZ_0518_EUR', 'PGC_SCZ_0418b', 'PGC_SCZ_2014_EUR']
g11=['CLOZUK_SCZ_2018_withPGC']
g12=['PGC_MDD_2018_no23andMe', 'SSGAC_EDU_2018_no23andMe', 'PGC_BIP_2019_wave3', 'PGC_SCZ_0518_EUR', 'PGC_SCZ_2014_EUR', 'PGC_BIP_2016', 'PGC_MDD_2018_Howard_no23andMe', 'PGC_ASD_2017_iPSYCH', 'PGC_ADHD_2017_EUR']
g13=['23andMe_SWB_2018', '23andMe_EXTRA_2016', '23andMe_AGREE_2016', '23andMe_CONSC_2016', '23andMe_NEUR_2016', '23andMe_OPEN_2016']
g14=['SSGAC_DEPRESSIVE_2016', 'SSGAC_SWB_2016', 'SSGAC_NEUROTICISM_2016', 'CTG_NEUR_2018_no23andMe']


g15 = ['PGC_SCZ_2014_EUR', 'PGC_SCZ_0518_EUR', 'PGC_BIP_2016', 'PGC_BIP_2019_wave3', 'PGC_MDD_2018_with23andMe', 'PGC_MDD_2018_Howard_with23andMe', 'GIANT_HEIGHT_2018_UKB', 'CTG_COG_2018', 'SSGAC_EDU_2018_no23andMe', 'PGC_ASD_2017_iPSYCH', 'PGC_ADHD_2017_EUR'] 

g20 = ['assoc.qt.h2_0.41.pi_0.0032.run_{}.320K.all.gwas'.format(run) for run in range(10, 20)]
g21 = ['assoc.qt.h2_0.41.pi_0.00032.run_{}.320K.all.gwas'.format(run) for run in range(1, 11)]

g22 = 'PGC_ADHD_2017_EUR PGC_ASD_2017_iPSYCH PGC_BIP_2016 PGC_BIP_2019_wave3 CLOZUK_SCZ_2018_withPGC PGC_SCZ_0518_EUR PGC_MDD_2018_no23andMe PGC_MDD_2018_with23andMe'.split()
g23 = 'CTG_COG_2018 SSGAC_EDU_2018_no23andMe'.split()
g24 = 'CTG_NEUR_2018_no23andMe 23andMe_SWB_2018 SSGAC_RISK_2019 UKB_MOOD_2019 CTG_NEUR_2018_onlyUKB_with23andMe'.split()
g25 = ['GIANT_HEIGHT_2018_UKB']

g30 = 'PGC_SCZ_2014_EUR PGC_MDD_2018_Howard_no23andMe CLOZUK_SCZ_2018_withPGC PGC_BIP_2019_wave3 PGC_BIP_2016 PGC_ADHD_2017_EUR PGC_ASD_2017_iPSYCH SSGAC_EDU_2018_no23andMe CTG_COG_2018 GIANT_HEIGHT_2018_UKB UKB_HEIGHT_2018_irnt'.split()
g31 =  'UKB_MOOD_2019 SSGAC_RISK_2019 UKB_LONELY_2018_Loneliness IHGC_MIG_2016_with23andMe GIANT_BMI_2018_UKB_v2'.split()
#g30 = 'PGC_BIP_2016 PGC_ADHD_2017_EUR PGC_ASD_2017_iPSYCH PGC_MDD_2018_Howard_no23andMe UKB_HEIGHT_2018_irnt GIANT_HEIGHT_2018_UKB'.split()
#g30 = 'PGC_BIP_2016 PGC_ADHD_2017_EUR PGC_ASD_2017_iPSYCH UKB_HEIGHT_2018_irnt GIANT_HEIGHT_2018_UKB'.split()

ugmg_traits = []
bgmg_traits = []

ugmg_traits = g30 + g31 + g6

#ugmg_traits=g21
#ugmg_traits = g22 + g23 + g24 + g25
#ugmg_traits = g0 + g1 + g2
#ugmg_traits = g0 + g1 + g2 + g5 + g7 + g12 + g13 + g14

#bgmg_traits.extend(              list(itertools.combinations(set(g15), 2)) )
#bgmg_traits.extend(              list( itertools.combinations(set(g13+g14), 2)))
#bgmg_traits.extend(              list( itertools.product(set(g15), set(g13+g14))))
#bgmg_traits.extend(              list( itertools.product(set(['PGC_SCZ_2014_EUR', 'PGC_BIP_2019_wave3', 'CTG_COG_2018', 'SSGAC_EDU_2018_no23andMe']), set(g13+g14))))

#bgmg_traits.extend ( list ( itertools.product ( set(g22 + g23 + g24 + g25), set(g22 + g23 + g24 + g25))))


#bgmg_traits = list(itertools.combinations(set(g0+g1+g2), 2))
#bgmg_traits = list(itertools.combinations(set(g0 + g1 + g2), 2))
#bgmg_traits = list(itertools.combinations(g9b + ['PGC_SCZ_0518_EUR'], 2))
#bgmg_traits =  list(itertools.product(g11, g12+g9r+g9b))
#bgmg_traits =  list(itertools.product(['CLOZUK_SCZ_2018_withPGC'], g30))
#bgmg_traits = list(itertools.combinations(set(g0), 2))


######################################################################
template_head='''
#$ -cwd
#$ -l h_vmem=120G
#$ -l h_rt=36:00:00
#$ -pe smp 3

set MIXER_ROOT=/home/oleksandr/github/mixer 
set SUMSTAT=/space/syn03/1/data/GWAS/SUMSTAT
set SUMSTATnomhc=/space/syn03/1/data/GWAS/SUMSTAT/TMP/nomhc
#set SUMSTATldsr=/space/syn03/1/data/GWAS/SUMSTAT/TMP/ldsr
set SUMSTATldsr=/space/syn03/1/data/GWAS/SUMSTAT/TMP/nomhc
#set SUMSTATnomhc=/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/alexeas_power_test
#set SUMSTATldsr=/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/alexeas_power_test
set PYTHON=/home/oleksandr/miniconda3/bin/python3
set LDFILE=/space/syn03/1/data/GWAS/SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld
set BIMFILE=/space/syn03/1/data/GWAS/SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim
#      --bim-file /space/gwas-syn1/1/data/GWAS/UKBioBank/projects/plsa_mixer/ukb_genetics_qc/ukb_bed/ukb_imp_chr@_v3_qc.bim \\
#      --ld-file /space/gwas-syn1/1/data/GWAS/UKBioBank/projects/plsa_mixer/ukb_genetics_qc/ukb_bed/ukb_imp_chr@_v3_qc.run1.ld \\
'''

template_ugmg='''
set OUT={out}
set TRAIT={trait}
set EXTRACT=/space/syn03/1/data/GWAS/SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep{rep}.snps

$PYTHON $MIXER_ROOT/precimed/mixer.py fit1 \\
      --trait1-file $SUMSTATldsr/$TRAIT.sumstats.gz \\
      --out ${{OUT}}$TRAIT.fit.rep{rep} \\
      --bim-file $BIMFILE --ld-file $LDFILE \\
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so --threads 12 \\
      --extract $EXTRACT \\

#      --chr2use 21 --fit-sequence diffevo-fast neldermead-fast --diffevo-fast-repeats 2
#      --extract $SUMSTAT/LDSR/w_hm3.justrs \\

$PYTHON $MIXER_ROOT/precimed/mixer.py test1 \\
      --trait1-file $SUMSTATnomhc/$TRAIT.sumstats.gz \\
      --load-params-file $OUT/$TRAIT.fit.rep{rep}.json \\
      --out ${{OUT}}$TRAIT.test.rep{rep} \\
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so --threads 12 \\
      --bim-file $BIMFILE --ld-file $LDFILE \\

#      --chr2use 21
'''

template_bgmg='''
set OUT={out}
set TRAIT1={trait1}
set TRAIT2={trait2}
set EXTRACT=/space/syn03/1/data/GWAS/SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep{rep}.snps

$PYTHON $MIXER_ROOT/precimed/mixer.py fit2 \\
      --trait1-file $SUMSTATldsr/$TRAIT1.sumstats.gz \\
      --trait2-file $SUMSTATldsr/$TRAIT2.sumstats.gz \\
      --trait1-params-file ${{OUT}}$TRAIT1.fit.rep{rep}.json \\
      --trait2-params-file ${{OUT}}$TRAIT2.fit.rep{rep}.json \\
      --out ${{OUT}}${{TRAIT1}}_vs_${{TRAIT2}}.fit.rep{rep} \\
      --extract $EXTRACT \\
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so --threads 12 \\
      --bim-file $BIMFILE --ld-file $LDFILE \\

#      --chr2use 21 --fit-sequence diffevo-fast neldermead-fast --diffevo-fast-repeats 2
#      --extract $SUMSTAT/LDSR/w_hm3.justrs \\

$PYTHON $MIXER_ROOT/precimed/mixer.py test2 \\
      --trait1-file $SUMSTAT/TMP/nomhc/$TRAIT1.sumstats.gz \\
      --trait2-file $SUMSTAT/TMP/nomhc/$TRAIT2.sumstats.gz \\
      --load-params-file ${{OUT}}${{TRAIT1}}_vs_${{TRAIT2}}.fit.rep{rep}.json \\
      --out ${{OUT}}${{TRAIT1}}_vs_${{TRAIT2}}.test.rep{rep} \\
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so --threads 12 \\
      --bim-file $BIMFILE --ld-file $LDFILE \\

#      --chr2use 21
'''
##########################################################

numrep=21

ugmg_traits =list(sorted(set(ugmg_traits)))
#bgmg_traits = []
#bgmg_traits = list(sorted(set([tuple(sorted(x)) for x in bgmg_traits])))
#bgmg_traits = list(sorted(set([x for x in bgmg_traits])))

skipped = 0; missing=0; 
processed = set()  # exclude duplicates

cmd_list = []
touch_list = []

for trait, rep in itertools.product(ugmg_traits, list(range(1, numrep))):
    if not os.path.exists('/space/syn03/1/data/GWAS/SUMSTAT/TMP/nomhc/{}.sumstats.gz'.format(trait)): raise ValueError('file note exits: {}'.format(trait))
    out_name = ugmg_pattern.format(trait, rep)
    if os.path.exists(out_name): skipped+=1; print('exists {}'.format(out_name)); continue;
    processed.add(out_name)

    cmd = template_head + template_ugmg.format(out=out_prefix, trait=trait, rep=rep)
    cmd_list.append(cmd)
    touch_list.append(out_name)

for (trait1, trait2), rep in itertools.product(bgmg_traits, list(range(1, numrep))):
    if trait1==trait2: continue
    out_name = bgmg_pattern.format(trait1, trait2, rep)
    if os.path.exists(out_name): skipped+=1; print('exists {}'.format(out_name)); continue;
    if os.stat(ugmg_pattern.format(trait1, rep)).st_size == 0: missing+=1; continue; # print('skip {}'.format(out_name)); continue
    if os.stat(ugmg_pattern.format(trait2, rep)).st_size == 0: missing+=1; continue; # print('skip {}'.format(out_name)); continue
    #print('tbd {}'.format(out_name))
    processed.add(out_name)

    cmd = template_head + template_bgmg.format(out=out_prefix, trait1=trait1, trait2=trait2, rep=rep)
    cmd_list.append(cmd)
    touch_list.append(out_name)

print('submission:')
for cmd, touch in zip(cmd_list, touch_list):
    with open('run_script.sh', 'w') as f: f.write(cmd)
    if num_submit > 0:
        if not dry_run: os.system('qsub run_script.sh')
        if not dry_run: os.system('touch {}'.format(touch))
        print('qsub {}'.format(touch))
    num_submit -= 1
    if num_submit <= 0: break

print('processed: {}, skipped: {}, missing input, total: {}:  {}'.format(len(processed), skipped,  missing, len(cmd_list)))

