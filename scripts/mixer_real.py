import os.path
import sys
import glob
import itertools

dry_run=True
num_submit = int(sys.argv[1])

out_prefix = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/plsa_mixer/real_run1/mixer'
ugmg_pattern = out_prefix+'_{}.fit.json'
bgmg_pattern = out_prefix+'_{}_vs_{}.fit.json'

g0 = ['SSGAC_EDU_2018_no23andMe', 'PGC_SCZ_2014_EUR', 'PGC_BIP_2016', 'GIANT_HEIGHT_2018_UKB'] # 'PGC_MDD_2018_with23andMe', 'PGC_ADHD_2017_EUR', 'PGC_ASD_2017_iPSYCH',  
g1 = ['CTG_COG_2018', 'SSGAC_EDU_2018_no23andMe']
g2 = ['PGC_SCZ_2014_EUR', 'PGC_SCZ_0518_EUR', 'PGC_BIP_2016', 'PGC_MDD_2018_Howard_no23andMe', 'PGC_ADHD_2017_EUR', 'PGC_ASD_2017_iPSYCH']
g3 = ['GSCAN_DRINK_2019_DrinksPerWeek', 'CTG_NEUR_2018_no23andMe', 'ICC_CANNABIS_2018_UKB', 'UKB_LONELY_2018_Loneliness', 'SSGAC_RISK_2019',  'UKB_CHRONOTYPE_2018', 'UKB_SLEEP_2018', 'CTG_INSOMNIA_2018']
g3r =['GSCAN_DRINK_2019_DrinksPerWeek',                            'ICC_CANNABIS_2018_UKB', 'UKB_LONELY_2018_Loneliness', 'SSGAC_RISK_2019']
g4 = ['UKB_HS_2019', 'GIANT_HEIGHT_2018_UKB', 'UKB_HEIGHT_2018_irnt', 'GIANT_BMI_2018_UKB_v2', 'EGG_BIRTHWEIGHT_2016', 'GIANT_BMI_2015_EUR', 'GIANT_HEIGHT_2014', 'GIANT_WHR_2015_EUR']
g4r= [                                                                'GIANT_BMI_2018_UKB_v2',                         'GIANT_BMI_2015_EUR',                      'GIANT_WHR_2015_EUR']
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
g12=['PGC_MDD_2018_no23andMe', 'SSGAC_EDU_2018_no23andMe', 'PGC_SCZ_2014_EUR', 'PGC_BIP_2016', 'PGC_MDD_2018_Howard_no23andMe', 'PGC_ASD_2017_iPSYCH', 'PGC_ADHD_2017_EUR']

ugmg_traits = g10
bgmg_traits = list(itertools.combinations(g9b + ['PGC_SCZ_0518_EUR'], 2))
#bgmg_traits =  list(itertools.product(g11, g12+g9r+g9b))


######################################################################
template_head='''
#!/bin/bash
#$ -cwd
#$ -l h_vmem=120G
#$ -l h_rt=36:00:00

set MIXER_ROOT=/home/oleksandr/github/mixer 
set SUMSTAT=/space/syn03/1/data/GWAS/SUMSTAT
set PYTHON=/home/oleksandr/miniconda3/bin/python3
'''

template_ugmg='''
set OUT_DIR={out}
set TRAIT={trait}

$PYTHON $MIXER_ROOT/precimed/mixer.py fit \
      --trait1-file $SUMSTAT/TMP/ldsr/$TRAIT.sumstats.gz \
      --out $OUT_DIR/$TRAIT.fit \
      --extract $SUMSTAT/LDSR/w_hm3.justrs --ci-alpha 0.05 \
      --bim-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so \

$PYTHON $MIXER_ROOT/precimed/mixer.py fit \
      --trait1-file $SUMSTAT/TMP/nomhc/$TRAIT.sumstats.gz \
      --load-params-file $OUT_DIR/$TRAIT.fit.json \
      --out $OUT_DIR/$TRAIT.test \
      --fit-sequence load inflation --power-curve --qq-plots --kmax 100 \
      --bim-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so \
'''

template_bgmg='''
set OUT_DIR={out}
set TRAIT1={trait1}
set TRAIT2={trait2}

$PYTHON $MIXER_ROOT/precimed/mixer.py fit \
      --trait1-file $SUMSTAT/TMP/ldsr/$TRAIT1.sumstats.gz \
      --trait2-file $SUMSTAT/TMP/ldsr/$TRAIT2.sumstats.gz \
      --trait1-params-file $OUT_DIR/$TRAIT1.fit.json \
      --trait2-params-file $OUT_DIR/$TRAIT2.fit.json \
      --out $OUT_DIR/${{TRAIT1}}_vs_${{TRAIT2}}.fit \
      --extract $SUMSTAT/LDSR/w_hm3.justrs --ci-alpha 0.05 \
      --bim-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so \

$PYTHON $MIXER_ROOT/precimed/mixer.py fit \
      --trait1-file $SUMSTAT/TMP/nomhc/$TRAIT1.sumstats.gz \
      --trait2-file $SUMSTAT/TMP/nomhc/$TRAIT2.sumstats.gz \
      --load-params-file $OUT_DIR/${{TRAIT1}}_vs_${{TRAIT2}}.fit.json \
      --out $OUT_DIR/${{TRAIT1}}_vs_${{TRAIT2}}.test \
      --fit-sequence load inflation --qq-plots --kmax 100 \
      --bim-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file $SUMSTAT/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld \
      --lib $MIXER_ROOT/src/build/lib/libbgmg.so \
'''
##########################################################

skipped = 0; duplicated = 0; processed = set()  # exclude duplicates

cmd_list = []
touch_list = []

for trait in ugmg_traits:
    out_name = ugmg_pattern.format(trait)
    if os.path.exists(out_name): skipped+=1; continue
    if out_name in processed: duplicated+=1; continue
    processed.add(out_name)

    cmd = template_head + template_ugmg.format(out=out_prefix, trait=trait)
    cmd_list.append(cmd)
    touch_list.append(out_name)

for trait1, trait2 in bgmg_traits:
    out_name = bgmg_pattern.format(trait1, trait2)
    if os.path.exists(out_name): skipped+=1; continue
    if out_name in processed: duplicated+=1; continue
    processed.add(out_name)

    cmd = template_head + template_bgmg.format(out=out_prefix, trait1=trait1, trait2=trait2)
    cmd_list.append(cmd)
    touch_list.append(out_name)

for cmd, touch in zip(cmd_list, touch_list):
    with open('run_script.sh', 'w') as f: f.write(cmd)
    if num_submit > 0:
        if not dry_run: os.system('qsub run_script.sh')
        if not dry_run: os.system('touch {}'.format(touch))
        print('qsub {}'.format(cmd))
    num_submit -= 1
    if num_submit <= 0: break

print('skipped: {}, duplicated: {}'.format(skipped, duplicated))
