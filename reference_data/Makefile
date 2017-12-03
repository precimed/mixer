LDSR=/mnt/h/NORSTORE/MMIL/SUMSTAT/LDSR
BINARY_ANNOT=$(LDSR)/LDSR_Annot/1000G_EUR_Phase3_baseline_split

WORK_FOLDER_NAME=bgmg-annot-1000G-Phase3-EUR-dec2017
WORK=/mnt/h/work/$(WORK_FOLDER_NAME)
PYTHON=python
GITHUB=/mnt/c/Users/Oleksandr/Documents/GitHub
STD=/mnt/h/NORSTORE/MMIL/SUMSTAT/STD

GENOTYPE=$(LDSR)/LDSR_Annot/1000G_EUR_Phase3_plink_with_cm
REF_10M=$(GENOTYPE)/snplist.gz
REF_1M=$(LDSR)/1m.ref

# NB: ldsc.py must support --l4, --bias-r2, --r2-min, --r2-max flags (https://github.com/ofrei/ldsc)
LDSC_PY=$(GITHUB)/ldsc/ldsc.py
SUMSTATS_PY=$(GITHUB)/python_convert/sumstats.py

ANNOTS=\
	base

# NB! Rules for good make file
# 1. Target must be a file - not a fake label
# seq -s ' ' 0 0.05 0.95
R2_BINS_MIN = 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95
R2_BINS_MAX = 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00

R2_BINS = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
#R2_BINS = 1 10

CHRS = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
#CHRS = 21 22

jobs := $(foreach i,$(ANNOTS),$(foreach j,$(CHRS),$(foreach k,$(R2_BINS),$i-chr$j-r2bin$k)))

jobs2 := $(foreach i,$(ANNOTS),$(foreach j,$(CHRS),$i-chr$j))

jobsIK := $(foreach i,$(ANNOTS),$(foreach k,$(R2_BINS),$i-r2bin$k))

i = $(word 1,$(subst -, ,$*))
j = $(subst chr,,$(word 2,$(subst -, ,$*)))
k = $(subst r2bin,,$(word 3,$(subst -, ,$*)))
K = $(subst r2bin,,$(word 2,$(subst -, ,$*)))

# mat depends on all but this is not stated in dependency graph - you have to call it manually
all: LDr2hist biased_LDr2hist biased_LDr4hist
mat: LDr2hist_mat biased_LDr2hist_mat biased_LDr4hist_mat $(WORK)/w_ld.mat $(WORK)/ref_10m.mat $(WORK)/ref_1m.mat

LDr2hist: $(foreach job,$(jobs),$(WORK)/LDr2hist/$(job).l2.ldscore.gz)
LDr2hist_mat: $(foreach job,$(jobsIK),$(WORK)/LDr2hist/$(job).l2.mat)
biased_LDr2hist: $(foreach job,$(jobs),$(WORK)/biased_LDr2hist/$(job).l2.ldscore.gz)
biased_LDr2hist_mat: $(foreach job,$(jobsIK),$(WORK)/biased_LDr2hist/$(job).l2.mat)
biased_LDr4hist: $(foreach job,$(jobs),$(WORK)/biased_LDr4hist/$(job).l4.ldscore.gz)
biased_LDr4hist_mat: $(foreach job,$(jobsIK),$(WORK)/biased_LDr4hist/$(job).l4.mat)

input: | $(REF_1M) $(REF_10M) $(LDSR)/w_hm3.snplist $(LDSR)/MATLAB_Annot/infomat.mat $(GENOTYPE) $(LDSR)/LDSR_Annot/1000G_Phase3_frq $(WORK)

# ================  Folder structure ==================
$(WORK):
	mkdir $@

$(WORK)/w_ld: | $(WORK)
	mkdir $@

$(WORK)/LDr2hist: | $(WORK)
	mkdir $@
	
$(WORK)/biased_LDr2hist: | $(WORK)
	mkdir $@
	
$(WORK)/biased_LDr4hist: | $(WORK)
	mkdir $@

# ================ ref_10m and ref_1m ==================
$(WORK)/ref_10m.mat: $(REF_10M)
	$(PYTHON) $(SUMSTATS_PY) ref-to-mat --ref $(REF_10M) --out $@ --force --numeric-only

$(WORK)/ref_1m.mat: $(REF_1M)
	$(PYTHON) $(SUMSTATS_PY) ref-to-mat --ref $(REF_1M) --out $@ --force --numeric-only

# ======================== w_ld ========================

$(WORK)/w_ld/chr%.l2.ldscore.gz: | input $(WORK)/w_ld
	$(PYTHON) $(LDSC_PY) --l2 --ld-wind-cm 1 --bfile $(GENOTYPE)/1000G.EUR.QC.$* --out $(WORK)/w_ld/chr$* --extract $(LDSR)/w_hm3.snplist

$(WORK)/w_ld.mat: $(addprefix $(WORK)/w_ld/chr, $(addsuffix .l2.ldscore.gz, $(CHRS))) 
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(REF_1M) --out $@ --force \
		--ldscore $(WORK)/w_ld/chr@.l2.ldscore.gz \
		--M $(WORK)/w_ld/chr@.l2.M \
		--M-5-50 $(WORK)/w_ld/chr@.l2.M_5_50 \

# ======================== mafvec ========================

$(WORK)/mafvec_1m.mat: | input
	$(PYTHON) $(SUMSTATS_PY) frq-to-mat --ref $(REF_1M) --out $@ --force \
		--frq $(LDSR)/LDSR_Annot/1000G_Phase3_frq/1000G.EUR.QC.@.frq

$(WORK)/mafvec_10m.mat: | input
	$(PYTHON) $(SUMSTATS_PY) frq-to-mat --ref $(REF_10M) --out $@ --force \
		--frq $(LDSR)/LDSR_Annot/1000G_Phase3_frq/1000G.EUR.QC.@.frq

# ========================  LDr2hist =====================

$(WORK)/LDr2hist/%.l2.ldscore.gz: | $(WORK)/LDr2hist
	$(PYTHON) $(LDSC_PY) --l2 --ld-wind-cm 1\
		--bfile $(GENOTYPE)/1000G.EUR.QC.$j\
		--annot $(BINARY_ANNOT)/$i/$j.annot.gz\
		--r2-min $(word $k,$(R2_BINS_MIN)) \
		--r2-max $(word $k,$(R2_BINS_MAX)) \
		--out $(basename $(basename $(basename $@)))

$(WORK)/LDr2hist/%.l2.mat:
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(REF_10M) --force \
		--ldscore $(WORK)/LDr2hist/$i-chr@-r2bin$K.l2.ldscore.gz \
		--M $(WORK)/LDr2hist/$i-chr@-r2bin$K.l2.M \
		--M-5-50 $(WORK)/LDr2hist/$i-chr@-r2bin$K.l2.M_5_50 \
		--out $@

$(WORK)/biased_LDr4hist/%.l4.ldscore.gz: | $(WORK)/biased_LDr4hist
	$(PYTHON) $(LDSC_PY) --l4 --ld-wind-cm 1\
		--bfile $(GENOTYPE)/1000G.EUR.QC.$j\
		--annot $(BINARY_ANNOT)/$i/$j.annot.gz\
		--r2-min $(word $k,$(R2_BINS_MIN)) \
		--r2-max $(word $k,$(R2_BINS_MAX)) \
		--out $(basename $(basename $(basename $@)))

$(WORK)/biased_LDr4hist/%.l4.mat:
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(REF_10M) --force \
		--ldscore $(WORK)/biased_LDr4hist/$i-chr@-r2bin$K.l4.ldscore.gz \
		--M $(WORK)/biased_LDr4hist/$i-chr@-r2bin$K.l4.M \
		--M-5-50 $(WORK)/biased_LDr4hist/$i-chr@-r2bin$K.l4.M_5_50 \
		--out $@

$(WORK)/biased_LDr2hist/%.l2.ldscore.gz: | $(WORK)/biased_LDr2hist
	$(PYTHON) $(LDSC_PY) --l2 --ld-wind-cm 1\
		--bias-r2\
		--bfile $(GENOTYPE)/1000G.EUR.QC.$j\
		--annot $(BINARY_ANNOT)/$i/$j.annot.gz\
		--r2-min $(word $k,$(R2_BINS_MIN)) \
		--r2-max $(word $k,$(R2_BINS_MAX)) \
		--out $(basename $(basename $(basename $@)))

$(WORK)/biased_LDr2hist/%.l2.mat:
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(REF_10M) --force \
		--ldscore $(WORK)/biased_LDr2hist/$i-chr@-r2bin$K.l2.ldscore.gz \
		--M $(WORK)/biased_LDr2hist/$i-chr@-r2bin$K.l2.M \
		--M-5-50 $(WORK)/biased_LDr2hist/$i-chr@-r2bin$K.l2.M_5_50 \
		--out $@

# ======================== upload =====================

upload:	
	rsync -avzP $(WORK) oleksanf@login.nird.sigma2.no:/nird/home/oleksanf/NS9114K/users/oleksanf/$(WORK_FOLDER_NAME)
