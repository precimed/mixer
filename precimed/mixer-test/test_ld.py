import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
from mixer import libbgmg

ld=np.loadtxt('mixer-test/tiny/chr21qc.ld.gz'); r2=np.multiply(ld, ld); r4=np.multiply(r2, r2)
frq=pd.read_csv('mixer-test/tiny/chr21qc.frq', delim_whitespace=True)  # nsubj=100, nchr=200 => 3 digits precision is accurate
bim=pd.read_csv('mixer-test/tiny/chr21qc.bim', delim_whitespace=True, header=None, names='CHR SNP GP BP A1 A2'.split())

lib=libbgmg.LibBgmg('../src/build/lib/libbgmg.so', init_log='mixer-test/test.log', context_id=0, dispose=True)
lib.defvec = bim['CHR'].notnull()
lib.chrnumvec = bim['CHR'].astype(int).values

r2min = 0.0499 # 0.05
ldscore_r2min = 0.001 # 0.001

lib.calc_ld_matrix(bfile='mixer-test/tiny/chr21qc', outfile='mixer-test/tiny/chr21qc.mixer.ld', r2min=r2min, ldscore_r2min=ldscore_r2min)
lib.set_ld_r2_coo_from_file(21, 'mixer-test/tiny/chr21qc.mixer.ld')
lib.set_ld_r2_csr()
lib_hetvec = 2*np.multiply(lib.mafvec, 1-lib.mafvec)
hetvec = 2*np.multiply(frq.MAF, 1-frq.MAF)
r2_het = np.multiply(r2, np.tile(hetvec, (len(hetvec), 1)))
r4_het = np.multiply(r2_het, r2_het)

# check that MAF is OK, except that libbgmg MAF is opposite (e.g. about A2, not about A1)
print('max MAF difference:', np.max(np.abs(1-lib.mafvec) - frq.MAF))
print('max HET difference:', np.max(np.abs(lib_hetvec - hetvec)))

lib.set_option('retrieve_ld_tag_type', 0); ld_tag_r2_sum_above = lib.ld_tag_r2_sum
lib.set_option('retrieve_ld_tag_type', 1); ld_tag_r2_sum_below = lib.ld_tag_r2_sum
lib.set_option('retrieve_ld_tag_type', 2); ld_tag_r2_sum_adjust_for_hvec_above = lib.ld_tag_r2_sum
lib.set_option('retrieve_ld_tag_type', 3); ld_tag_r2_sum_adjust_for_hvec_below = lib.ld_tag_r2_sum
lib.set_option('retrieve_ld_tag_type', 0); ld_tag_r4_sum_above = lib.ld_tag_r4_sum
lib.set_option('retrieve_ld_tag_type', 2); ld_tag_r4_sum_adjust_for_hvec_above = lib.ld_tag_r4_sum

print('max diff in r2_sum_above: ', np.max(np.abs(np.sum(np.multiply(r2, r2 >= r2min), 1) - ld_tag_r2_sum_above)))
print('max diff in r2_sum_below: ', np.max(np.abs(np.sum(np.multiply(r2, (r2 < r2min) & (r2 >= ldscore_r2min) ), 1) - ld_tag_r2_sum_below)))
print('max diff in r2_sum_total: ', np.max(np.abs(np.sum(np.multiply(r2, (r2 >= ldscore_r2min) ), 1) - (ld_tag_r2_sum_above + ld_tag_r2_sum_below) )))
print('max diff in r4_sum_above: ', np.max(np.abs(np.sum(np.multiply(r4, r2 >= r2min), 1) - ld_tag_r4_sum_above)))

print('max diff in r2_sum_het_above: ', np.max(np.abs(np.sum(np.multiply(r2_het, r2 >= r2min), 1) - ld_tag_r2_sum_adjust_for_hvec_above)))
print('max diff in r2_sum_het_below: ', np.max(np.abs(np.sum(np.multiply(r2_het, (r2 < r2min) & (r2 >= ldscore_r2min) ), 1) - ld_tag_r2_sum_adjust_for_hvec_below)))
print('max diff in r2_sum_het_total: ', np.max(np.abs(np.sum(np.multiply(r2_het, (r2 >= ldscore_r2min) ), 1) - (ld_tag_r2_sum_adjust_for_hvec_above + ld_tag_r2_sum_adjust_for_hvec_below) )))
print('max diff in r4_sum_het_above: ', np.max(np.abs(np.sum(np.multiply(r4_het, r2 >= r2min), 1) - ld_tag_r4_sum_adjust_for_hvec_above)))

[snp, tag, lib_r] = lib.get_ld_r2_chr(21)
lib_r_mat = coo_matrix((lib_r, (snp, tag)), shape=(lib.num_snp, lib.num_tag)).toarray()
lib_r2 = np.multiply(lib_r_mat, lib_r_mat)

r2_above = np.multiply(ld, ld); r2_above = np.multiply(r2_above, r2_above >= r2min);
print('max diff in r2 sparse matrix: ', np.max(np.abs(lib_r2 - r2_above)))

"""
max MAF difference: 2.8610229518832853e-08
max HET difference: 5.4162740702190515e-08
max diff in r2_sum_above:  0.0009558380511123232
max diff in r2_sum_below:  0.0009984926233030933
max diff in r2_sum_total:  0.0012759102995758553
max diff in r4_sum_above:  0.0014510578888824455
max diff in r2_sum_het_above:  0.00021826482092990318
max diff in r2_sum_het_below:  0.00047987763431756036
max diff in r2_sum_het_total:  0.0005588952533823743
max diff in r4_sum_het_above:  0.0001568508609040009
max diff in r2 sparse matrix:  3.027609337069581e-05
"""
