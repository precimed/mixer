LDSR=/mnt/h/NORSTORE/MMIL/SUMSTAT/LDSR
WORK=/mnt/h/work/bgmg-annot-1000G-Phase3-EUR
PYTHON=python
GITHUB=/mnt/c/Users/Oleksandr/Documents/GitHub

# NB: ldsc.py must support --l4 flag, implemented by ofrei
LDSC_PY=$(GITHUB)/ldsc/ldsc.py
SUMSTATS_PY=$(GITHUB)/python_convert/sumstats.py

# NB! Rules for good make file
# 1. Target must be a file - not a fake label

CHRS = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
#CHRS = 21 22

TARGETS = ref_l2.mat ref_l4.mat biased_ref_l2.mat biased_ref_l4.mat w_ld.mat mafvec.mat mafvec_all_snps_1000G_Phase3_frq.mat

all: $(addprefix $(WORK)/, $(TARGETS))

biased_ref_l2: $(addprefix $(WORK)/biased_ref_l2/chr, $(addsuffix .l2.ldscore.gz, $(CHRS))) | input $(WORK)
biased_ref_l4: $(addprefix $(WORK)/biased_ref_l4/chr, $(addsuffix .l4.ldscore.gz, $(CHRS))) | input $(WORK)

input: | $(LDSR)/1m.ref $(LDSR)/w_hm3.snplist $(LDSR)/MATLAB_Annot/infomat.mat $(LDSR)/LDSR_Annot/1000G_EUR_Phase3_plink_with_cm $(LDSR)/LDSR_Annot/1000G_Phase3_frq $(WORK)

# Folder structure
$(WORK):
	mkdir $@

$(WORK)/ref_l2: | $(WORK)
	mkdir $@

$(WORK)/ref_l4: | $(WORK)
	mkdir $@

$(WORK)/biased_ref_l2: | $(WORK)
	mkdir $@

$(WORK)/biased_ref_l4: | $(WORK)
	mkdir $@

$(WORK)/w_ld: | $(WORK)
	mkdir $@

$(WORK)/ref_l2/chr%.l2.ldscore.gz: | input $(WORK)/ref_l2 
	$(PYTHON) $(LDSC_PY) --l2 --ld-wind-cm 1 --bfile $(LDSR)/LDSR_Annot/1000G_EUR_Phase3_plink_with_cm/1000G.EUR.QC.$* --out $(WORK)/ref_l2/chr$*

$(WORK)/ref_l4/chr%.l4.ldscore.gz: | input $(WORK)/ref_l4
	$(PYTHON) $(LDSC_PY) --l4 --ld-wind-cm 1 --bfile $(LDSR)/LDSR_Annot/1000G_EUR_Phase3_plink_with_cm/1000G.EUR.QC.$* --out $(WORK)/ref_l4/chr$*

$(WORK)/biased_ref_l2/chr%.l2.ldscore.gz: | input $(WORK)/biased_ref_l2
	$(PYTHON) $(LDSC_PY) --l2 --ld-wind-cm 1 --bfile $(LDSR)/LDSR_Annot/1000G_EUR_Phase3_plink_with_cm/1000G.EUR.QC.$* --out $(WORK)/biased_ref_l2/chr$* --bias-r2 --r2-min 0.05

$(WORK)/biased_ref_l4/chr%.l4.ldscore.gz: | input $(WORK)/biased_ref_l4
	$(PYTHON) $(LDSC_PY) --l4 --ld-wind-cm 1 --bfile $(LDSR)/LDSR_Annot/1000G_EUR_Phase3_plink_with_cm/1000G.EUR.QC.$* --out $(WORK)/biased_ref_l4/chr$* --bias-r2 --r2-min 0.05

$(WORK)/w_ld/chr%.l2.ldscore.gz: | input $(WORK)/w_ld
	$(PYTHON) $(LDSC_PY) --l2 --ld-wind-cm 1 --bfile $(LDSR)/LDSR_Annot/1000G_EUR_Phase3_plink_with_cm/1000G.EUR.QC.$* --out $(WORK)/w_ld/chr$* --extract $(LDSR)/w_hm3.snplist

$(WORK)/ref_l2.mat: $(addprefix $(WORK)/ref_l2/chr, $(addsuffix .l2.ldscore.gz, $(CHRS))) 
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(LDSR)/1m.ref --out $@ --force \
		--ldscore $(WORK)/ref_l2/chr@.l2.ldscore.gz \
		--M $(WORK)/ref_l2/chr@.l2.M \
		--M-5-50 $(WORK)/ref_l2/chr@.l2.M_5_50 \

$(WORK)/ref_l4.mat: $(addprefix $(WORK)/ref_l4/chr, $(addsuffix .l4.ldscore.gz, $(CHRS))) 
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(LDSR)/1m.ref --out $@ --force \
		--ldscore $(WORK)/ref_l4/chr@.l4.ldscore.gz \

$(WORK)/biased_ref_l2.mat: $(addprefix $(WORK)/biased_ref_l2/chr, $(addsuffix .l2.ldscore.gz, $(CHRS))) 
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(LDSR)/1m.ref --out $@ --force \
		--ldscore $(WORK)/biased_ref_l2/chr@.l2.ldscore.gz \
		--M $(WORK)/biased_ref_l2/chr@.l2.M \
		--M-5-50 $(WORK)/biased_ref_l2/chr@.l2.M_5_50 \

$(WORK)/biased_ref_l4.mat: $(addprefix $(WORK)/biased_ref_l4/chr, $(addsuffix .l4.ldscore.gz, $(CHRS))) 
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(LDSR)/1m.ref --out $@ --force \
		--ldscore $(WORK)/biased_ref_l4/chr@.l4.ldscore.gz \


$(WORK)/w_ld.mat: $(addprefix $(WORK)/w_ld/chr, $(addsuffix .l2.ldscore.gz, $(CHRS))) 
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(LDSR)/1m.ref --out $@ --force \
		--ldscore $(WORK)/w_ld/chr@.l2.ldscore.gz \
		--M $(WORK)/w_ld/chr@.l2.M \
		--M-5-50 $(WORK)/w_ld/chr@.l2.M_5_50 \

$(WORK)/mafvec.mat: | input
	$(PYTHON) $(SUMSTATS_PY) frq-to-mat --ref $(LDSR)/1m.ref --out $@ --force \
		--frq $(LDSR)/LDSR_Annot/1000G_Phase3_frq/1000G.EUR.QC.@.frq

$(WORK)/mafvec_all_snps_1000G_Phase3_frq.mat: | input
	$(PYTHON) $(SUMSTATS_PY) frq-to-mat --out $@ --force \
		--frq $(LDSR)/LDSR_Annot/1000G_Phase3_frq/1000G.EUR.QC.@.frq
	
upload:	
	rsync -avzP $(WORK) oleksanf@login.nird.sigma2.no:/nird/home/oleksanf/NS9114K/users/oleksanf/bgmg-annot
