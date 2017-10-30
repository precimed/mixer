# Makefile to run BGMG analysis for each pairs of phenotypes (e.g. each from TRAITS vs each from TRAITS2, aka rectangular block)

DATA=/space/syn03/1/data/GWAS/SUMSTAT/LDSR/MATLAB_Data

TRAITS=\
	PGC_SCZ_2014\
	PGC_BIP_2012_lift\

TRAITS2=\
	PGC_SCZ_2014\
	PGC_BIP_2012_lift\
	LIPIDS_TG_2013\
	IIBDGC_CD_2015_EUR\

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
	matlab -singleCompThread -nodisplay -nosplash -nodesktop -r "options.use_legacy_impl=1;trait1='$(DATA)/$i.mat'; trait2='$(DATA)/$j.mat'; data_path=''; BGMG_run; exit;"
