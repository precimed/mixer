import configparser
import datetime
import sys
import os
import pandas as pd
from multiprocessing import Pool


def read_bim(bim_f, verbose=False):
    """
    Get SNPs with their chromosomes and positions from bim_f.
    Return:
        pd.DataFrame(index=snp_id, data=[chr, bp, a1, a2])
    """
    bim_df = pd.read_table(bim_f, header=None, usecols=[0,1,3,4,5],
           names=["CHR","SNP","BP","A1","A2"], index_col="SNP", na_filter=False)
    if verbose: print(f"{bim_df.shape[0]} SNPs in {bim_f}")
    return bim_df


def read_snp(bim_f, verbose=False):
    """
    Get SNPs with their chromosomes and positions from bim_f.
    Return:
        pd.Series(index=1..N, data=snp_id)
    """
    snp = pd.read_table(bim_f, header=None, names=["SNP"], usecols=[1],
            na_filter=False, squeeze=True)
    if verbose: print(f"{snp.size} SNPs in {bim_f}")
    return snp


def read_frq(frq_f, verbose=False):
    """
    Get frq values for SNPs.
    Return:
        pd.Series(index=snp_id, data=frq)
    """
    frq = pd.read_table(frq_f, usecols=["SNP","MAF"], index_col="SNP",
            na_filter=False, squeeze=True, delim_whitespace=True)
    if verbose: print(f"{frq.size} SNPs in {frq_f}")
    return frq


def read_ld(ld_f, snp, verbose=False):
    """
    Get r2 values from ld_f for SNPs in snp.
    Return:
        pd.DataFrame(index=1..N, data=[snp_a_ind, snp_b_ind, r2])
        where snp_a_ind/snp_b_ind are indices in bim_f of snp_a and snp_b from
        ld_f correspondingly
    """
    snp_ind = {s:i for i,s in enumerate(snp)}
    snp2ind = lambda s: snp_ind[s]
    # use 32-bit float type ('f4') for r2 valeus
    ld_df = pd.read_table(ld_f, na_filter=False, delim_whitespace=True,
            usecols=["SNP_A","SNP_B","R2"], dtype={"R2":"f4"},
            converters={"SNP_A":snp2ind,"SNP_B":snp2ind})
    if verbose: print(f"{ld_df.shape[0]} r2 values in {ld_f}")
    return ld_df


def make_ld_pkl(bim_f, ld_f, ld_pkl_f, verbose=False):
    snp = read_snp(bim_f, verbose=verbose)
    ld_df = read_ld(ld_f, snp, verbose=verbose)
    ld_df.to_pickle(ld_pkl_f)
    if verbose: print(f"LD information pickled to {ld_pkl_f}")



if __name__ == "__main__":
    print(f"makeinput.py started at {datetime.datetime.now()}")
    print(f"Reading config from {sys.argv[1]}")
    cfg = configparser.ConfigParser()
    cfg.read(sys.argv[1])

    n_proc = cfg["general"].getint("n_proc")
    verbose = cfg["general"].getboolean("verbose")

    if verbose: print(f"Using {n_proc} processors")

    template_dir = cfg["template"].get("template_dir")
    chr2use = cfg["template"].get("chr2use").split()
    bim_prefix = cfg["template"].get("bim_prefix")
    bim_suffix = cfg["template"].get("bim_suffix")
    ld_prefix = cfg["template"].get("ld_prefix")
    ld_suffix = cfg["template"].get("ld_suffix")
    out_dir = cfg["output"].get("out_dir")
    if verbose: print(f"Using {len(chr2use)} chromosomes: {', '.join(chr2use)}")

    bim_fiels = [os.path.join(template_dir,"".join((bim_prefix,c,bim_suffix)))
            for c in chr2use]
    ld_files = [os.path.join(template_dir,"".join((ld_prefix,c,ld_suffix)))
            for c in chr2use]
    ld_pkl_files = [os.path.join(out_dir, f"chr{c}.ld.pkl") for c in chr2use]

    for f in (bim_fiels + ld_files):
        if not os.path.isfile(f):
            raise ValueError(f"{f} file does not exist")
    if not os.path.isdir(out_dir):
        raise ValueError(f"{out_dir} dir does not exist")

    v = [verbose]*len(bim_fiels)
    with Pool(processes=n_proc) as pool:
        pool.starmap(make_ld_pkl, zip(bim_fiels,ld_files,ld_pkl_files,v))

    print(f"makeinput.py finished at {datetime.datetime.now()}")
