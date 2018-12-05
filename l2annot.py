import datetime
import pandas as pd
import numpy as np
from scipy import io as sio


def get_snps(bim_f, verbose=False):
    """
    Get SNPs with their chromosomes and positions from bim_f.
    """
    bim_df = pd.read_table(bim_f, header=None, names=["CHR","SNP","BP"],
           usecols=[0,1,3], index_col="SNP", na_filter=False)
    # drop duplicated SNPs from annot file
    i = bim_df.index.duplicated(keep='first')
    if i.any():
        bim_df = bim_df.loc[~i,:]
    if verbose: print(f"{len(bim_df)} SNPs in {bim_f}")
    return bim_df


def get_annot(annot_f, bim_f=None, bim_df=None, verbose=False):
    """
    Get annotations for SNPs in bim_f/bim_df.
    Return:
        pd.DataFrame(index=snp_id, columns=annot_id, data=binary_annot,
            dtype='bool')
    """
    if bim_f is None and bim_df is None:
        raise ValueError("in get_annot: either bim_f or bim_df must be not None")
    elif bim_df is None:
        bim_df = get_snps(bim_f, verbose)

    annot_df = pd.read_table(annot_f, index_col=0, na_filter=False)
    if verbose: print(f"{len(annot_df.columns)} annot in {annot_f}")

    # drop duplicated SNPs from annot file
    i = annot_df.index.duplicated(keep='first')
    if i.any():
        annot_df = annot_df.loc[~i,:]
    # take only SNPs from bim file in the corresponding order
    annot_df = annot_df.loc[bim_df.index,:]
    return annot_df.astype("bool")



def get_l2annot(ld_f, bim_f=None, annot_f=None, bim_df=None, annot_df=None,
        verbose=False):
    """
    Get l2 values for SNPs in bim_f/bim_df for each annotation category in
    annot_f/annot_df.
    l2(snp_i, category_k) = sum_{for all snp_j from category_k}(r^2(snp_i, snp_j))
    Return:
        pd.DataFrame(index=snp_id, columns=annot_id, data=l2_matrix, dtype='float')
    """
    if bim_f is None and bim_df is None:
        raise ValueError("in get_l2annot: either bim_f or bim_df must be not None")
    elif bim_df is None:
        bim_df = get_snps(bim_f, verbose)
    if annot_f is None and annot_df is None:
        raise ValueError("in get_l2annot: either annot_f or annot_df must be not None")
    elif annot_df is None:
        annot_df = get_annot(annot_f, bim_f, bim_df, verbose)

    if not (bim_df.index == annot_df.index).all():
        raise ValueError("in get_l2annot: bim_df and annot_df must have the same SNPs in the same order")

    snp_ind = {s:i for i,s in enumerate(bim_df.index)}
    snp2ind = lambda s: snp_ind[s]

    annot_ind = {a:i for i,a in enumerate(annot_df.columns)}
    annot2ind = lambda a: annot_ind[a]

    # take only SNPs from bim file in the corresponding order
    annot_mat = annot_df.values

    l2annot_mat = np.zeros((len(bim_df), len(annot_df.columns)))
    l2annot_mat[annot_mat] += 1 # each SNP is in LD with itself
   
    if verbose: print(f"processing {ld_f}")
    chunk_size = 1000000
    ld_df_reader = pd.read_table(ld_f, na_filter=False, delim_whitespace=True,
            usecols=["SNP_A","SNP_B","R2"], dtype={"R2":float},
            converters={"SNP_A":snp2ind,"SNP_B":snp2ind}, chunksize=chunk_size)
    for chunk_i, ld_df_chunk in enumerate(ld_df_reader):
        # add r2 for SNP_A
        snp_b_i, snp_b_annot_i = np.nonzero(annot_mat[ld_df_chunk.SNP_B])
        snp_a_i = ld_df_chunk.SNP_A.values[snp_b_i]
        r2 = ld_df_chunk.R2.values[snp_b_i]
        # use np.add.at as some elements can be incremented multiple times
        np.add.at(l2annot_mat, (snp_a_i, snp_b_annot_i), r2)

        # add r2 for SNP_B
        snp_a_i, snp_a_annot_i = np.nonzero(annot_mat[ld_df_chunk.SNP_A])
        snp_b_i = ld_df_chunk.SNP_B.values[snp_a_i]
        r2 = ld_df_chunk.R2.values[snp_a_i]
        np.add.at(l2annot_mat, (snp_b_i, snp_a_annot_i), r2)

        if verbose: print(f"{(chunk_i+1)*chunk_size} lines processed")
    return pd.DataFrame(index=bim_df.index, columns=annot_df.columns,
            data=l2annot_mat)


def get_maf(frq_f, bim_f=None, bim_df=None, verbose=False):
    """
    Get maf values for SNPs in bim_f/bim_df.
    Return:
        pd.Series(index=snp_id, data=maf)
    """
    if bim_f is None and bim_df is None:
        raise ValueError("in get_l2annot: either bim_f or bim_df must be not None")
    elif bim_df is None:
        bim_df = get_snps(bim_f, verbose)
    maf = pd.read_table(frq_f, usecols=["SNP","MAF"], index_col="SNP",
            squeeze=True, delim_whitespace=True)
    if verbose: print(f"{len(maf)} SNPs in {frq_f}")
    maf = maf[bim_df.index]
    return maf 


if __name__ == "__main__":
    print(f"l2annot.py started at {datetime.datetime.now()}")

    annot_f = "/mnt/seagate10/genotypes/ldsc489eur10m/annot/template.ldsc489eur10m.sorted.complete_annot_hg19.annomat.uniq.txt.gz"
    chromosomes = [str(i) for i in range(1,23)]
    bim_files = [f"/mnt/seagate10/genotypes/ldsc489eur10m/1000G.EUR.QC.{c}.bim" for c in chromosomes]
    ld_files = [f"/mnt/seagate10/genotypes/ldsc489eur10m/1000G.EUR.QC.{c}.r2.ld.gz" for c in chromosomes]
    frq_files = [f"/mnt/seagate10/genotypes/ldsc489eur10m/1000G.EUR.QC.{c}.maf.frq" for c in chromosomes]
    verbose = True
    
    save_mat = False
    save_hdf = False
    save_csv = False
    # extension of the output file (mat/h5/csv) will be added automatically
    out_files = [f"/mnt/seagate10/genotypes/ldsc489eur10m/annot/l2annot/schork.chr{c}" for c in chromosomes]


    make_ld_informed_annot = True
    add_intergenic = True
    out_files_ld_informed = [f"/mnt/seagate10/genotypes/ldsc489eur10m/annot/schork/ld_informed.chr{c}.csv" for c in chromosomes]
    ld_files_schork = [f"/mnt/seagate10/genotypes/ldsc489eur10m/schork/1000G.EUR.QC.{c}.schork.r2.ld.gz" for c in chromosomes]
    annot2use = ["5UTR", "3UTR", "Exon", "Intron", "1kUp", "1kDown", "10kUp", "10kDown"]
    aux_annot = ["NoncodingTranscript", "100kUp", "100kDown", "mirna", "tfbs"]

    ld_informed_i = 0
    for bim_f,ld_f,frq_f,out_f in zip(bim_files,ld_files,frq_files,out_files):
        bim_df = get_snps(bim_f, verbose=verbose)
        annot_df = get_annot(annot_f, bim_df=bim_df, verbose=verbose)
        l2annot_df = get_l2annot(ld_f, bim_df=bim_df, annot_df=annot_df,
                verbose=verbose)
        maf = get_maf(frq_f, bim_df=bim_df, verbose=verbose)

        if save_mat:
            out_f_mat = f"{out_f}.mat"
            sio.savemat(out_f_mat,{"l2_annot":l2annot_df.values,
                "bin_annot":annot_df.values, "maf":maf.values,
                "annot":annot_df.columns.values, "chr":bim_df.CHR.values,
                "bp":bim_df.BP.values})
            print(f"saved to {out_f_mat}")
        if save_hdf or save_csv:
            l2annot_df_columns_backup = l2annot_df.columns.copy()
            l2annot_df.columns = [f"{c}_l2" for c in l2annot_df.columns]
            out_df = pd.concat([bim_df, maf, annot_df, l2annot_df], axis=1)
            out_df.reset_index(inplace=True)
            if save_hdf:
                out_f_hdf = f"{out_f}.h5"
                out_df.to_hdf(out_f_hdf, key="df" , mode='w', format="fixed")
                print(f"saved to {out_f_hdf}")
            if save_csv:
                out_f_csv = f"{out_f}.csv"
                out_df.to_csv(out_f_csv, index=False)
                print(f"saved to {out_f_csv}")
            l2annot_df.columns = l2annot_df_columns_backup
            
        if make_ld_informed_annot:
            ld_f_schork = ld_files_schork[ld_informed_i]
            l2annot_schork_df = get_l2annot(ld_f_schork, bim_df=bim_df, annot_df=annot_df,
                    verbose=verbose)
            ld_informed_annot_df = (l2annot_schork_df.loc[:,annot2use] >= 1)
            if add_intergenic:
                df_aux_annot = annot_df.loc[:,aux_annot]
                intergenic = ~(ld_informed_annot_df.any(1) | df_aux_annot.any(1))
                ld_informed_annot_df.loc[:,"Intergenic"] = intergenic
            out_f_ld_informed = out_files_ld_informed[ld_informed_i]
            ld_informed_annot_df.to_csv(out_f_ld_informed)
            ld_informed_i += 1
            print(f"saved to {out_f_ld_informed}")


            

    print(f"ld2annot.py finished at {datetime.datetime.now()}")
