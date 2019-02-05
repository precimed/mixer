import configparser
import argparse
import logging
import sys
import numpy as np
import pandas as pd
from scipy import sparse
from multiprocessing import Pool


def _create_logger():
    logger = logging.getLogger("makeinput")
    logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # create formatter
    formatter = logging.Formatter('[%(asctime)s | %(name)s | %(levelname)s] : %(message)s')
    # add formatter to ch
    ch.setFormatter(formatter)
    # add ch to logger
    logger.addHandler(ch)
    return logger

# create and configure logger
LOGGER = _create_logger()


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Prepare input data for optimization.")

    parser.add_argument("config_file", help="configuration file")
    parser.add_argument("--ld2npz", action="store_true",
        help="make LD npz files")
    parser.add_argument("--rprune", action="store_true",
        help="make random pruned list of SNPs")
    parser.add_argument("--rpruneweights", action="store_true",
        help="estimate SNP weights based on random pruning iterations")
    
    return parser.parse_args(args)



def ld2npz(ld_f, bim_f, out_npz_f):
    """
    Produces npz file with ld block size, i2 and r2 data taking SNP indices
    according to the bim_f. All SNP ids from plink ld file must be in bim_f.
    Args:
        ld_f: plink ld file
        bim_f: plink bim file
        out_npz_f: output npz file name
    Output:
        out_npz_f npz file with 3 fields:
            bs: size of ld block, for SNPs from bim_f
            i2: index of the second SNP (assuming i1 are ordered)
            r2: r2 correlation values between i1 and i2
    Example:
    i1: 1 1 1 1 2 2 2 3 3 3 3 3 4 4 5 5 5   # len(i1) = 2*len(ld_f) + len(bim_f)
    bs: 4       3     5         2   3       # len(bs) = len(bim_f)
    i2: 2 3 5 1 1 3 2 1 2 4 5 3 3 4 1 3 5   # len(i2) = len(i1)
    r2: x y z 1 a s 1 q w e r 1 t 1 d f 1   # x, y ... from [0,1], len(r2) = len(i1)
    """
    LOGGER.info(f"Processing {ld_f}")
    snps = np.loadtxt(bim_f, usecols=1, dtype=bytes)

    ind = np.arange(len(snps), dtype='u4')
    snp_ind_dict = dict(zip(snps, ind))
    if len(snps) != len(snp_ind_dict):
        raise ValueError(f"Duplicated variant ids in {bim_f}. Must be unique.")

    snp2ind = lambda snp: snp_ind_dict[snp]

    #rm ld_file = os.path.join(template_dir, f"chr{c}.r2.ld.gz")
    i1, i2, r2 = np.loadtxt(ld_f, skiprows=1, usecols=[2,5,6], 
        dtype={'names':("i1","i2","r2"),'formats':('u4','u4','f4')},
        converters={2:snp2ind,5:snp2ind}, unpack=True)

    i12 = np.concatenate([i1, i2, ind])
    i12_sort_idx = np.argsort(i12)
    i12 = i12[i12_sort_idx]
    # bs[i] = number of SNPs in ld with i-th variant from bim_f
    bs = np.searchsorted(i12, ind, side='right')
    bs[1:] -= bs[:-1]

    i21 = np.concatenate([i2, i1, ind])
    i21 = i21[i12_sort_idx]

    r2 = np.concatenate([r2, r2, np.ones(len(ind), dtype='f4')])
    r2 = r2[i12_sort_idx]

    np.savez(out_npz_f,bs=bs, i2=i21, r2=r2)
    LOGGER.info(f"{out_npz_f} created")



def rprune(ld_npz_f, bim_f, ld_thresh, out_f, snps2overlap=None):
    """
    Generate a out_f with SNP ids which survives a sungle iteration of random
    pruning with ld_thresh.
    Overlap SNP ids from bim_f with SNP ids from snps2overlap_f and
    Wrapper for _rprune which works with ld_npz_f produced by ld2npz.
    Args:
        ld_npz_f: a file generated with ld2npz
        bim_f: plink bim file
        ld_thresh: threshold for r2 values (r2 >= ld_thresh are pruned)
        out_f: output file
        snps2overlap: SNP ids to overlap with SNPs from bim_f
    Output:
        out_f: a file with SNP ids which survived random pruning
    """
    LOGGER.info(f"pruning {ld_npz_f}")
    ld_csr = _ld_npz2ld_csr(ld_npz_f)
    snps = pd.read_table(bim_f, sep='\t', header=None, usecols=[1],
            dtype='str', squeeze=True, engine="c", na_filter=False)

    idx2use = None if snps2overlap is None else snps.isin(snps2overlap)

    survived_idx = _rprune(ld_csr, ld_thresh, idx2use)
    snps[survived_idx].to_csv(out_f, index=False, header=False)
    LOGGER.info(f"{out_f} created")



def rpruneweights(ld_npz_f, bim_f, ld_thresh, n_iter, out_f, snps2overlap=None):
    """
    Wrapper for _rpruneweights which works with ld_npz_f produced by ld2npz.
    """
    LOGGER.info(f"estimating pruning weights for {ld_npz_f} from {n_iter} iterations")
    ld_csr = _ld_npz2ld_csr(ld_npz_f)

    snps = pd.read_table(bim_f, sep='\t', header=None, usecols=[1],
            dtype='str', squeeze=True, engine="c", na_filter=False)

    idx2use = None if snps2overlap is None else snps.isin(snps2overlap)

    snp_weights = _rpruneweights(ld_csr, ld_thresh, n_iter, idx2use)
    pd.Series(snp_weights).to_csv(out_f, index=False, header=False)
    LOGGER.info(f"{out_f} created")



def _ld_npz2ld_csr(ld_npz_f):
    """
    Make scipy.sparse.csr matrix from ld_npz_f produced by ld2npz.
    """
    ld_npz = np.load(ld_npz_f)
    bs = ld_npz.get("bs")
    i2 = ld_npz.get("i2")
    r2 = ld_npz.get("r2")
    i1 = np.repeat(np.arange(bs.size), bs)
    return sparse.coo_matrix((r2, (i1, i2)), shape=(bs.size, bs.size), dtype='f4').tocsr()



def _rprune(ld_csr, ld_thresh, idx2use):
    """
    Get random pruned indices corresponding to ld_thresh considering only SNPs
    which are True in idx2use. SNPs which have False in idx2use will not
    survive random pruning.
    Args:
        ld_csr: npz file created with ld2npz
        ld_thresh: threshold for LD r2 values, r2 >= ld_thresh are pruned
        idx2use: indices of SNPs to use
    Output:
        np.array(len=#(SNPs in ld_npz_f), dtype=bool)
        Bool array with True values corresponding to variants that survived
        random pruning, and False values corresponding to non-survivals.
    """
    n_snps = ld_csr.shape[0]
    idx = np.arange(n_snps)
    np.random.shuffle(idx)
    if idx2use is None:
        survived = np.ones(n_snps, dtype='bool')
    else:
        survived = np.zeros(n_snps, dtype='bool')
        survived[idx2use] = True
    for i in idx:
        if survived[i]:
            _, j, r2 = sparse.find(ld_csr[i,:])
            survived[j[r2>=ld_thresh]] = False
            survived[i] = True
    return survived



def _rpruneweights(ld_csr, ld_thresh, n_iter, idx2use=None):
    """
    Make weights induced by random pruning with ld_thresh threshold taking into
    consideration only SNPs which are True in idx2use. SNPs which are False in
    idx2use will have the weight = 0.
    """
    n_snps = ld_csr.shape[0]
    rprune_w = np.zeros(n_snps)
    for i in range(n_iter):
        survived = _rprune(ld_csr, ld_thresh, idx2use)
        rprune_w += survived
    rprune_w /= n_iter
    return rprune_w



if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    cfg = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    cfg.read(args.config_file)
    # add file handler (writing to log file) to the LOGGER
    log_f = cfg["misc"].get("log_file")
    fh = logging.FileHandler(log_f)
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('[%(asctime)s | %(name)s | %(levelname)s] : %(message)s')
    fh.setFormatter(formatter)
    LOGGER.addHandler(fh)
    LOGGER.info(f"Reading config from {args.config_file}")

    n_proc = cfg["misc"].getint("n_proc", fallback=1)

    if args.ld2npz:
        LOGGER.info("Running ld2npz")
        ld_files = cfg["ld2npz"].get("ld_files").split()
        bim_files = cfg["ld2npz"].get("bim_files").split()
        out_npz_files = cfg["ld2npz"].get("out_npz_files").split()

        if not (len(ld_files) == len(bim_files) == len(out_npz_files)):
            raise ValueError(
                f"Number of bim, ld and out_npz files in {args.config_file} must be equal")
        n_proc = min(n_proc, len(ld_files))
        LOGGER.info(f"Using {n_proc} cores")
        with Pool(processes=n_proc) as pool:
            pool.starmap(ld2npz, zip(ld_files, bim_files, out_npz_files))

    if args.rprune:
        LOGGER.info("Running rprune")
        ldnpz_files = cfg["rprune"].get("ldnpz_files").split()
        bim_files = cfg["rprune"].get("bim_files").split()
        ld_threshold = cfg["rprune"].getfloat("ld_threshold")
        out_files = cfg["rprune"].get("out_files").split()
        snps2overlap_f = cfg["rprune"].get("snps2overlap_f")

        if not (len(ldnpz_files) == len(bim_files) == len(out_files)):
            raise ValueError(
                f"Number of bim, ldnpz and out files in {args.config_file} must be equal")
        
        if snps2overlap_f is None:
            snps2overlap = None
        else:
            LOGGER.info(f"reading snps2overlap from {snps2overlap_f}")
            snps2overlap = pd.read_table(snps2overlap_f, header=None, usecols=[0],
                    dtype='str', squeeze=True, engine="c", na_filter=False)

        for ldnpz_f, bim_f, out_f in zip(ldnpz_files, bim_files, out_files):
            rprune(ldnpz_f, bim_f, ld_threshold, out_f, snps2overlap)

    if args.rpruneweights:
        LOGGER.info("Running rpruneweights")
        ldnpz_files = cfg["rpruneweights"].get("ldnpz_files").split()
        bim_files = cfg["rpruneweights"].get("bim_files").split()
        ld_threshold = cfg["rpruneweights"].getfloat("ld_threshold")
        n_iter = cfg["rpruneweights"].getint("n_iter")
        out_files = cfg["rpruneweights"].get("out_files").split()
        snps2overlap_f = cfg["rpruneweights"].get("snps2overlap_f")

        if not (len(ldnpz_files) == len(bim_files) == len(out_files)):
            raise ValueError(
                f"Number of bim, ldnpz and out files in {args.config_file} must be equal")
        
        if snps2overlap_f is None:
            snps2overlap = None
        else:
            LOGGER.info(f"reading snps2overlap from {snps2overlap_f}")
            snps2overlap = pd.read_table(snps2overlap_f, header=None, usecols=[0],
                    dtype='str', squeeze=True, engine="c", na_filter=False)

        for ldnpz_f, bim_f, out_f in zip(ldnpz_files, bim_files, out_files):
            rpruneweights(ldnpz_f, bim_f, ld_threshold, n_iter, out_f, snps2overlap)

    LOGGER.info("Finished")
