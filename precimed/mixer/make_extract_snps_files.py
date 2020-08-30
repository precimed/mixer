import sys
sys.path.append('/home/oleksanf/github/mixer')

import precimed
import precimed.mixer
import precimed.mixer.libbgmg
import precimed.mixer.utils

from precimed.mixer.utils import AnnotUnivariateParams
from precimed.mixer.utils import AnnotUnivariateParametrization

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import random

libbgmg_so = '/home/oleksanf/github/mixer/src/build/lib/libbgmg.so'
logfile = '/home/oleksanf/github/mixer/mixer.log'
bim_file = '/home/oleksanf/vmshare/data/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim'
ld_file = '/home/oleksanf/vmshare/data/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld'
#bim_file = '/home/oleksanf/vmshare/data/bfile_merged/chr@.bim'
#ld_file = '/home/oleksanf/vmshare/data/bfile_merged/chr@.run4.ld'

trait1_file = ''
trait2_file = ''
extract = ''
exclude = ''
#annot_file = '/home/oleksanf/vmshare/data/bfile_merged/baseline.chr@.annot.gz'

libbgmg = precimed.mixer.libbgmg.LibBgmg(libbgmg_so, dispose=True)
libbgmg.init_log(logfile)

chr_labels= list(range(1, 23))
libbgmg.init(bim_file, "", chr_labels, trait1_file, trait2_file, exclude, extract)

for chr_label in chr_labels: 
    libbgmg.set_ld_r2_coo_from_file(int(chr_label), ld_file.replace('@', str(chr_label)))
    libbgmg.set_ld_r2_csr(int(chr_label))

mafvec = np.minimum(libbgmg.mafvec, 1-libbgmg.mafvec)
ref=pd.concat([pd.read_csv(bim_file.replace('@', str(chr_label)), sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split()) for chr_label in chr_labels])

snps_goodMAF = np.sum(mafvec>=maf_thresh)
maf_thresh = 0.05
r2_thresh = 0.8
subset = 2000000 # int(snps_goodMAF/5)
seed = 123
out_file = '/home/oleksanf/vmshare/data/LDSR/1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep@.snps'
repeats = 20

if seed is not None: np.random.seed(seed)
    
sets = []
for repeat in range(repeats):
    print(repeat)
    # step0 - generate random values for clumping (this is a way to implement random pruning)
    buf = np.random.rand(libbgmg.num_tag, 1)
    
    # step1 - filter SNPs below MAF threshold
    buf[mafvec<maf_thresh] = np.nan
    
    # step2 - select a random subset of SNPs
    indices = list(np.where(np.isfinite(buf))[0])
    sample = random.sample(indices, min(subset, len(indices)))
    mask = np.ones(buf.shape); mask[sample] = 0
    buf[mask == 1] = np.nan

    # step3 - further prune SNPs at certain r2 threshold
    buf_pruned = libbgmg.perform_ld_clump(r2_thresh, buf.flatten()) 
    ref.SNP[np.isfinite(buf_pruned)].to_csv(out_file.replace('@', str(repeat+1)), index=False, header=False)
    sets.append(set(ref.SNP[np.isfinite(buf_pruned)].values))

if False:
    overlaps = np.zeros((repeats, repeats))
    for i in range(repeats):
        for j in range(repeats):
            overlaps[i][j] = len(sets[i] & sets[j])

    overlaps.astype(np.int32)
