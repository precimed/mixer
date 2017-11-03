LDSR=/mnt/h/NORSTORE/MMIL/SUMSTAT/LDSR
WORK=/mnt/h/work/bgmg-annot-1000G-Phase3-EUR
PYTHON=python
GITHUB=/mnt/c/Users/Oleksandr/Documents/GitHub
ALL_REFFILE=$(GITHUB)/BGMG/reference_data/10m.ref.gz
STD=/mnt/h/NORSTORE/MMIL/SUMSTAT/STD

TRAITS=PGC_SCZ_2014 PGC_BIP_2012_lift PGC_BIP_2016_qc PGC_SCZ_0917 PGC_SCZ_0917_trios_asia

# NB: ldsc.py must support --l4 flag, implemented by ofrei
LDSC_PY=$(GITHUB)/ldsc/ldsc.py
SUMSTATS_PY=$(GITHUB)/python_convert/sumstats.py
MAKE_LD_MATRIX_PY=$(GITHUB)/python_convert/make_ld_matrix/make_ld_matrix.py

# NB! Rules for good make file
# 1. Target must be a file - not a fake label

CHRS = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
#CHRS = 21 22

#10m.ref.gz produced as follows:
#	import pandas as pd
#	df = pd.concat([pd.read_table(r'H:\NORSTORE\MMIL\SUMSTAT\LDSR\LDSR_Annot\1000G_EUR_Phase3_plink_with_cm\1000G.EUR.QC.{}.bim'.format(chri), sep='\t', header=None) for chri in range(1, 23)])
#	df.columns=['CHR', 'SNP', 'GP', 'BP', 'A1', 'A2']
#	df.to_csv(r'C:\Users\Oleksandr\Documents\GitHub\BGMG\reference_data\10m.ref.gz', compression='gzip', sep='\t', index=False)


# Manually merge matlab LD matrices
#outname = 'reference_data_10m/binary_ld_p8.mat'; for chri=1:22, df = load(sprintf('H:\\work\\bgmg-annot-1000G-Phase3-EUR\\binary_ld_p8\\%i.mat', chri), 'id1', 'id2', 'nsnp'); if chri==1, id1=df.id1; id2=df.id2; nsnp=df.nsnp; else id1 = [id1; df.id1]; id2 = [id2; df.id2]; end; end; LDmat = sparse(double(id1),double(id2),true,double(nsnp),double(nsnp));LDmat = LDmat | speye(double(nsnp));LDmat = LDmat | (LDmat - LDmat');save(outname, 'LDmat', '-v7.3');clear;
#outname = 'reference_data_10m/binary_ld_p2.mat'; for chri=1:22, df = load(sprintf('H:\\work\\bgmg-annot-1000G-Phase3-EUR\\binary_ld_p2\\%i.mat', chri), 'id1', 'id2', 'nsnp'); if chri==1, id1=df.id1; id2=df.id2; nsnp=df.nsnp; else id1 = [id1; df.id1]; id2 = [id2; df.id2]; end; end; LDmat = sparse(double(id1),double(id2),true,double(nsnp),double(nsnp));LDmat = LDmat | speye(double(nsnp));LDmat = LDmat | (LDmat - LDmat');save(outname, 'LDmat', '-v7.3');clear;
#outname = 'reference_data_10m/binary_ld_p1.mat'; for chri=1:22, df = load(sprintf('H:\\work\\bgmg-annot-1000G-Phase3-EUR\\binary_ld_p1\\%i.mat', chri), 'id1', 'id2', 'nsnp'); if chri==1, id1=df.id1; id2=df.id2; nsnp=df.nsnp; else id1 = [id1; df.id1]; id2 = [id2; df.id2]; end; end; LDmat = sparse(double(id1),double(id2),true,double(nsnp),double(nsnp));LDmat = LDmat | speye(double(nsnp));LDmat = LDmat | (LDmat - LDmat');save(outname, 'LDmat', '-v7.3');clear;

TARGETS = ref_l2.mat ref_l4.mat biased_ref_l2.mat biased_ref_l4.mat w_ld.mat ref_mafvec.mat all_l2.mat all_l4.mat biased_all_l2.mat biased_all_l4.mat all_mafvec.mat all.mat

all: $(addprefix $(WORK)/, $(TARGETS))

biased_ref_l2: $(addprefix $(WORK)/biased_ref_l2/chr, $(addsuffix .l2.ldscore.gz, $(CHRS))) | input $(WORK)
biased_ref_l4: $(addprefix $(WORK)/biased_ref_l4/chr, $(addsuffix .l4.ldscore.gz, $(CHRS))) | input $(WORK)

input: | $(LDSR)/1m.ref $(LDSR)/w_hm3.snplist $(LDSR)/MATLAB_Annot/infomat.mat $(LDSR)/LDSR_Annot/1000G_EUR_Phase3_plink_with_cm $(LDSR)/LDSR_Annot/1000G_Phase3_frq $(WORK)

binary_ld_p1: $(addprefix $(WORK)/binary_ld_p1/, $(addsuffix .ld.gz, $(CHRS))) | input $(WORK)
binary_ld_p2: $(addprefix $(WORK)/binary_ld_p2/, $(addsuffix .ld.gz, $(CHRS))) | input $(WORK)
binary_ld_p8: $(addprefix $(WORK)/binary_ld_p8/, $(addsuffix .ld.gz, $(CHRS))) | input $(WORK)
binary_ld_mat_p1: $(addprefix $(WORK)/binary_ld_p1/, $(addsuffix .mat, $(CHRS))) | input $(WORK)
binary_ld_mat_p2: $(addprefix $(WORK)/binary_ld_p2/, $(addsuffix .mat, $(CHRS))) | input $(WORK)
binary_ld_mat_p8: $(addprefix $(WORK)/binary_ld_p8/, $(addsuffix .mat, $(CHRS))) | input $(WORK)

binary_ld: binary_ld_p1 binary_ld_p2 binary_ld_p8
binary_ld_mat: binary_ld_mat_p1 binary_ld_mat_p2 binary_ld_mat_p8

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

$(WORK)/binary_ld_p1: | $(WORK)
	mkdir $@

$(WORK)/binary_ld_p2: | $(WORK)
	mkdir $@

$(WORK)/binary_ld_p8: | $(WORK)
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

$(WORK)/binary_ld_p1/%.ld.gz: | $(WORK)/binary_ld_p1
	plink --memory 4096 --threads 2 --bfile $(LDSR)/LDSR_Annot/1000G_EUR_Phase3_plink_with_cm/1000G.EUR.QC.$* --r2 gz --ld-window-kb 10000 --ld-window 999999 --ld-window-r2 0.1 --out $(basename $(basename $@))

$(WORK)/binary_ld_p2/%.ld.gz: | $(WORK)/binary_ld_p2
	plink --memory 4096 --threads 2 --bfile $(LDSR)/LDSR_Annot/1000G_EUR_Phase3_plink_with_cm/1000G.EUR.QC.$* --r2 gz --ld-window-kb 10000 --ld-window 999999 --ld-window-r2 0.2 --out $(basename $(basename $@))

$(WORK)/binary_ld_p8/%.ld.gz: | $(WORK)/binary_ld_p8
	plink --memory 4096 --threads 2 --bfile $(LDSR)/LDSR_Annot/1000G_EUR_Phase3_plink_with_cm/1000G.EUR.QC.$* --r2 gz --ld-window-kb 10000 --ld-window 999999 --ld-window-r2 0.8 --out $(basename $(basename $@))

$(WORK)/binary_ld_p8/%.mat: $(WORK)/binary_ld_p8/%.ld.gz | $(WORK)/binary_ld_p8
	$(PYTHON) $(MAKE_LD_MATRIX_PY) --ldfile $< --ref $(ALL_REFFILE) --savemat $@

$(WORK)/binary_ld_p2/%.mat: $(WORK)/binary_ld_p2/%.ld.gz | $(WORK)/binary_ld_p2
	$(PYTHON) $(MAKE_LD_MATRIX_PY) --ldfile $< --ref $(ALL_REFFILE) --savemat $@

$(WORK)/binary_ld_p1/%.mat: $(WORK)/binary_ld_p1/%.ld.gz | $(WORK)/binary_ld_p1
	$(PYTHON) $(MAKE_LD_MATRIX_PY) --ldfile $< --ref $(ALL_REFFILE) --savemat $@

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

$(WORK)/ref_mafvec.mat: | input
	$(PYTHON) $(SUMSTATS_PY) frq-to-mat --ref $(LDSR)/1m.ref --out $@ --force \
		--frq $(LDSR)/LDSR_Annot/1000G_Phase3_frq/1000G.EUR.QC.@.frq


$(WORK)/all_mafvec.mat: | input
	$(PYTHON) $(SUMSTATS_PY) frq-to-mat --ref $(ALL_REFFILE) --out $@ --force \
		--frq $(LDSR)/LDSR_Annot/1000G_Phase3_frq/1000G.EUR.QC.@.frq

$(WORK)/all_l2.mat: $(addprefix $(WORK)/ref_l2/chr, $(addsuffix .l2.ldscore.gz, $(CHRS))) 
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(ALL_REFFILE) --out $@ --force \
		--ldscore $(WORK)/ref_l2/chr@.l2.ldscore.gz \
		--M $(WORK)/ref_l2/chr@.l2.M \
		--M-5-50 $(WORK)/ref_l2/chr@.l2.M_5_50 \

$(WORK)/all_l4.mat: $(addprefix $(WORK)/ref_l4/chr, $(addsuffix .l4.ldscore.gz, $(CHRS))) 
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(ALL_REFFILE) --out $@ --force \
		--ldscore $(WORK)/ref_l4/chr@.l4.ldscore.gz \

$(WORK)/biased_all_l2.mat: $(addprefix $(WORK)/biased_ref_l2/chr, $(addsuffix .l2.ldscore.gz, $(CHRS))) 
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(ALL_REFFILE) --out $@ --force \
		--ldscore $(WORK)/biased_ref_l2/chr@.l2.ldscore.gz \
		--M $(WORK)/biased_ref_l2/chr@.l2.M \
		--M-5-50 $(WORK)/biased_ref_l2/chr@.l2.M_5_50 \

$(WORK)/biased_all_l4.mat: $(addprefix $(WORK)/biased_ref_l4/chr, $(addsuffix .l4.ldscore.gz, $(CHRS))) 
	$(PYTHON) $(SUMSTATS_PY) ldsc-to-mat --ref $(ALL_REFFILE) --out $@ --force \
		--ldscore $(WORK)/biased_ref_l4/chr@.l4.ldscore.gz \

$(WORK)/all.mat: $(ALL_REFFILE)
	$(PYTHON) $(SUMSTATS_PY) ref-to-mat --ref $(ALL_REFFILE) --out $@ --force

mat: $(addsuffix .mat, $(addprefix $(WORK)/,$(TRAITS)))

$(WORK)/%.mat : $(STD)/%.csv.gz
	$(PYTHON) $(SUMSTATS_PY) mat --force --ref $(ALL_REFFILE) --sumstats $< --out $@
	
upload:	
	rsync -avzP $(WORK) oleksanf@login.nird.sigma2.no:/nird/home/oleksanf/NS9114K/users/oleksanf/bgmg-annot
