DATA=/space/syn03/1/data/GWAS/SUMSTAT/LDSR/MATLAB_Data

TRAITS=\
	PGC_SCZ_2014_EUR_qc\
	PGC_SCZ_2014\
	PGC_SCZ_0917\
	PGC_SCZ_0917_trios_asia\
	PGC_BIP_2012_lift\
	PGC_BIP_2016_qc\
	GIANT_HEIGHT_2014_lift\
	GIANT_BMI_2015_EUR_lift\
	GIANT_WHR_2015_EUR\
	IIBDGC_CD_2015_EUR\
	LIPIDS_TG_2013\
	IIBDGC_UC_2015_EUR\
	ReproGen_MENARCHE_2017_qc_lift\

#	SSGAC_DEPRESSIVE_2016\
#	SSGAC_NEUROTICISM_2016\
#	SSGAC_SWB_2016\


TRAITS2=$(TRAITS)

#	PGC_ADHD_2017_EUR\
#	IGAP_AD_2013_noAPOE\


TRAITS3=\
	23andMe_AGREE_2016\
	23andMe_CONSC_2016\
	23andMe_EXTRA_2016\
	23andMe_NEUR_2016\
	23andMe_OPEN_2016\
	BROADABC_ASB_2017_qc_lift\
	CTG_INSOMNIA_2017\
	CTG_INTELLIGENCE_2017\
	EAGLE_ADHD_2016_noGC_lift\
	ENIGMA_HV_2016_noPGC_noGC\
	ENIGMA_ICV_2016_noPGC_noGC_lift\
	ENIGMA_PALL_2016_noPGC_noGC\
	ENIGMA_PUT_2016_noPGC_noGC\
	GAMEON_BREAST_2013_BCAC\
	GAMEON_COLON_2015_CORECT\
	GAMEON_LUNG_2014_TRICL_6study\
	GAMEON_OVARIAN_2013_FOCI\
	GIANT_BMI_2015_EUR_lift\
	GIANT_HEIGHT_2014_lift\
	GIANT_WHR_2015_EUR\
	IGAP_AD_2013_noAPOE\
	IIBDGC_CD_2015_EUR\
	IIBDGC_IBD_2015_EUR\
	IIBDGC_UC_2015_EUR\
	LIPIDS_HDL_2013\
	LIPIDS_LDL_2013\
	LIPIDS_TG_2013\
	OKADA_RA_2014_EUR\
	PGC_ADHD_2017_EUR\
	PGC_BIP_2012_lift\
	PGC_BIP_2016_qc\
	PGC_SCZ_0917\
	PGC_SCZ_0917_trios_asia\
	PGC_SCZ_2014\
	PGC_SCZ_2014_EUR_qc\
	ReproGen_MENARCHE_2014_lift\
	ReproGen_MENARCHE_2017_qc_lift\
	ReproGen_MENOPAUSE_2015_lift\
	SSGAC_DEPRESSIVE_2016\
	SSGAC_EDU_2016\
	SSGAC_NEUROTICISM_2016\
	SSGAC_SWB_2016\
	TAG_CIGPERDAY_2010_lift\
	TAG_EVERSMOKE_2010_lift\
	TAG_FORMERSMOKE_2010_lift\
	TAG_SMOKEONSET_2010_lift\
	UKB_CHRONOTYPE_2016\
	UKB_COLLEGE_2016\
	UKB_MEMORY_2016\
	UKB_RT_2016\
	UKB_SLEEP_2016\
	UKB_VNR_2016\
	XXX_CRP_2009\

TRAITS_noMHC=$(addsuffix _noMHC, $(TRAITS))
TRAITS2_noMHC=$(addsuffix _noMHC, $(TRAITS2))

# https://stackoverflow.com/questions/15102267/nested-loops-in-makefile-compatible-with-j-n
jobs := $(foreach i,$(TRAITS_noMHC),$(foreach j,$(TRAITS2_noMHC),BGMG_run_$i-$j.txt))
#jobs := $(foreach i,$(TRAITS_noMHC),$(foreach j,$(TRAITS_noMHC),$(if $(findstring $i, $(word 1, $(sort $i $j))), $(if $(findstring $i, $j), ,BGMG_run_$i-$j.txt))))

.PHONY: all
all: $(jobs) ; echo $@ Success

i = $(firstword $(subst -, ,$*))
j = $(lastword $(subst -, ,$*))

UGMG: $(addsuffix .txt, $(addprefix UGMG_run_, $(TRAITS_noMHC)))

UGMG_run_%.txt: $(DATA)/%.mat
	matlab -nodisplay -nosplash -nodesktop -r "trait1='$<'; data_path=''; BGMG_run; exit;"

$(jobs): BGMG_run_%.txt:
	echo $@ 
	matlab -singleCompThread -nodisplay -nosplash -nodesktop -r "trait1='$(DATA)/$i.mat'; trait2='$(DATA)/$j.mat'; data_path=''; BGMG_run; exit;"
