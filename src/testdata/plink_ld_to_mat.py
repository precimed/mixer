import pandas as pd
import scipy.io as sio
import sys

if __name__ == '__main__':
	in_name = sys.argv[1]  # 'bfile_merged_10K_ldmat_p01_SNPwind50k_chr22.ld.gz'
	out_name = sys.argv[2]  # 'bfile_merged_10K_ldmat_p01_SNPwind50k_chr22.ld.mat'
	ref_name = sys.argv[3] # r"H:\Dropbox\shared\BGMG\all_chromosomes.ref.gz"

	print('args: {}'.format(sys.argv[1:]))
	
	print('reading {}...'.format(in_name))
	df=pd.read_table(in_name,delim_whitespace=True)

	print('reading {}...'.format(ref_name))
	ref = pd.read_table(ref_name,delim_whitespace=True)
	ref['index_A']=ref.index; ref['index_B']=ref.index; ref['SNP_A']=ref['SNP']; ref['SNP_B']=ref['SNP']

	print('merging 1...')
	df_with_pos = df;
	df_with_pos = pd.merge(df_with_pos, ref[['SNP_A', 'index_A']], how='left', on='SNP_A')

	print('merging 2...')
	df_with_pos = pd.merge(df_with_pos, ref[['SNP_B', 'index_B']], how='left', on='SNP_B')

	print('saving {}'.format(out_name))
	sio.savemat(out_name,
				{'index_A': df_with_pos['index_A'].values,
				 'index_B': df_with_pos['index_B'].values,
				 'r2' : df_with_pos['R2'].values}, format='5', do_compression=False, oned_as='column', appendmat=False)

	print('done')
