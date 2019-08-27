import configparser
import datetime
import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import norm
from multiprocessing import Pool


def process_template_chr(c, template_dir, bfile_prefix, out_dir):
    print(f"processing template chr {c}")
    bim_file = os.path.join(template_dir, f"{bfile_prefix}{c}.bim")
    snp = np.loadtxt(bim_file, usecols=1, dtype=bytes)
    ind = np.arange(len(snp), dtype='u4')
    snp_ind_dict = dict(zip(snp, ind))
    assert len(snp) == len(snp_ind_dict), f"Duplicated variant ids"
    # print(f"   {len(snp)} variants")

    snp2ind = lambda snp: snp_ind_dict[snp]

    ld_file = os.path.join(template_dir, f"{bfile_prefix}{c}.r2.ld.gz")
    r2_1i, r2_2i, r2 = np.loadtxt(ld_file, skiprows=1, usecols=[2,5,6], 
        dtype={'names':("r2_1i","r2_2i","r2"),'formats':('u4','u4','f4')},
        converters={2:snp2ind,5:snp2ind}, unpack=True)
    # print(f"   {len(r2)} r2 coefficients")

    r2_12i = np.concatenate([r2_1i, r2_2i, ind])
    r2_21i = np.concatenate([r2_2i, r2_1i, ind])
    r2 = np.concatenate([r2, r2, np.ones(len(ind), dtype='f4')])
    i = np.argsort(r2_12i)
    r2_12i = r2_12i[i]
    r2_21i = r2_21i[i]
    r2 = r2[i]

    # ld_bs = LD block size, ld_bs[i] = number of variants in ld with i-th variant from snp (template) array
    ld_bs = np.searchsorted(r2_12i, ind, side='right')
    ld_bs[1:] -= ld_bs[:-1]

    frq_file = os.path.join(template_dir, f"{bfile_prefix}{c}.maf.frq")
    si, maf = np.loadtxt(frq_file, skiprows=1, usecols=[1,4], 
        dtype={'names':("si","maf"),'formats':('u4','f4')},
        converters={1:snp2ind}, unpack=True)
    assert len(maf) == len(snp), "Length of maf and snp vectors should be equal"
    maf = maf[np.argsort(si)] # just to be sure that we have correct order
    het = 2*maf*(1 - maf)

    snp = snp.astype(str) # convert to unicode in python 3
    out_file = os.path.join(out_dir, f"template.chr{c}.npz")
    # May be use numpy.savez_compressed instead of uncompressed savez
    # r2_21i: [u4] indices of variants in ld blocks corresponding to r2 array, r2 = corr^2(r2_12i, r2_21i)
    # r2_21i is required to get het and annot indices for s2 and annot_s2 arrays correspondingly
    np.savez(out_file, snp=snp, ind=ind, het=het, ld_bs=ld_bs, r2_21i=r2_21i, r2=r2)
    print(f"{out_file} created")
    return len(snp)


def process_sumstats(sumstats_file, sumstats_out_id, out_dir, snp_col, p_col,
    effect_col, effect_baseline, n_col=None, ncase_col=None, ncontrol_col=None):
    print(f"Processing sumstats {sumstats_file}")
    usecols = [snp_col, p_col, effect_col]
    if (n_col is None) and (not ncase_col is None) and (not ncontrol_col is None):
        usecols += [ncase_col, ncontrol_col]
    elif (not n_col is None):
        usecols.append(n_col)
    else:
        raise ValueError(
            "Either 'n_col' or both 'ncase_col' and 'ncontrol_col' must be not None")
    df = pd.read_table(sumstats_file, delim_whitespace=True, usecols=usecols,
        index_col=snp_col)
    print(f"{len(df)} variants")
    df = df.loc[np.isfinite(df[p_col]),:]
    print(f"{len(df)} variants with defined p-value")
    df = df.loc[df[p_col]>0,:]
    print(f"{len(df)} variants with p-value > 0")
    i = df.index.duplicated(keep='first')
    if i.any():
        print(f"{i.sum()} duplicated variant ids")
        print("Only the first row with duplicated id will be retained")
        df = df.loc[~i,:]
    positive_effect = df[effect_col] > effect_baseline
    df["Z"] = norm.ppf(0.5*df[p_col]) # here all Z are < 0
    df.loc[positive_effect,"Z"] *= -1

    if n_col is None:
        ssize = 4*df[ncase_col]*df[ncontrol_col]/(df[ncase_col]+df[ncontrol_col])
    else:
        ssize = df[n_col]
    ssize = ssize.values
    
    out_file = os.path.join(out_dir, f"sumstats.{sumstats_out_id}.npz")
    np.savez(out_file, snp=df.index.values, z=df.Z.values, ssize=ssize)
    print(f"{out_file} created")


def process_annot(annot_file, annot_out_id, out_dir, snp_col, annot_col):
    print(f"Processing annot {annot_file}")
    df = pd.read_table(annot_file, usecols=[snp_col, annot_col],
        index_col=snp_col, dtype={annot_col:str})
    print(f"{len(df)} variants")
    i = df.index.duplicated(keep='first')
    if i.any():
        print(f"{i.sum()} duplicated variant ids")
        print("Only the first row with duplicated id will be retained")
        df = df.loc[~i,:]
    unique_annot = np.sort(df[annot_col].unique())
    print(f"{len(unique_annot)} unique annotation categories")

    assert len(unique_annot) < 256, "Maximum 255 annotation categories are supported"

    unique_annot_ind = np.arange(len(unique_annot), dtype='u1')
    annot_ind_dict = dict(zip(unique_annot, unique_annot_ind))
    annot = np.array([annot_ind_dict[s] for s in df[annot_col]], dtype='u1')

    out_file = os.path.join(out_dir, f"annot.{annot_out_id}.npz")
    np.savez(out_file, snp=df.index.values, annot=annot, categories=unique_annot,
        category_ind=unique_annot_ind)
    print(f"{out_file} created")


def process_qq_annot(qq_annot_file, qq_annot_out_id, out_dir):
    """
    Produce npz file with annotations for qq plot. Each SNP in qq
    plot annotaions can belong to several catigories (in contrast
    to template annotations used for optimization, which must be
    unique).
    qq_annot matrix in the output npz file should have dtype=np.bool
    Args:
        qq_annot_file: a file with annotation for qq plot.
            Must have a header line. 
            1st column must contain SNP ids.
            2nd and further columns form a binary matrix indicating
            whether corresponding SNP belongs to the annotation
            category (1) or not (0).
        qq_annot_out_id: string id which will be used to name output
            npz file (qq_annot.{qq_annot_out_id},npz).
        out_dir: directory where output file will be placed.
    """
    print(f"Processing qq annot file {qq_annot_file}")
    df = pd.read_table(qq_annot_file, index_col=0)
    print(f"{len(df)} variants")
    print(f"{len(df.columns)} annotation categories")
    i = df.index.duplicated(keep='first')
    if i.any():
        print(f"{i.sum()} duplicated variant ids")
        print("Only the first row with duplicated id will be retained")
        df = df.loc[~i,:]
    out_file = os.path.join(out_dir, f"qq_annot.{qq_annot_out_id}.npz")
    np.savez(out_file, snp=df.index.values,
        qq_annot=df.values.astype('u1'), categories=df.columns)
    print(f"{out_file} created")




if __name__ == "__main__":
    print(f"makeinput.py started at {datetime.datetime.now()}")
    print("Preparing input data for cmm analysis")
    print(f"Reading config from {sys.argv[1]}")
    cfg = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    cfg.read(sys.argv[1])

    make_template = cfg["general"].getboolean("make_template")
    make_annot = cfg["general"].getboolean("make_annot")
    make_qq_annot = cfg["general"].getboolean("make_qq_annot")
    make_sumstats = cfg["general"].getboolean("make_sumstats")
    out_dir = cfg["general"].get("out_dir")


    if make_template:
        # process template
        template_dir = cfg["template"].get("template_dir")
        bfile_prefix = cfg["template"].get("bfile_prefix")
        out_dir = out_dir
        chr2use = [c for c in cfg["template"].get("chr2use").split()]
        chr2use.sort(key=lambda x: int(x))
        n_proc = cfg["template"].getint("n_proc")

        with Pool(processes=n_proc) as pool:
            n = sum(pool.starmap(process_template_chr,
                ((c, template_dir, bfile_prefix, out_dir) for c in chr2use)))
        print(f"{n} variants in template")


    if make_annot:
        # process annot
        annot_file = cfg["annot"].get("annot_file")
        annot_out_id = cfg["annot"].get("annot_out_id")
        out_dir = out_dir
        snp_col = cfg["annot"].get("snp_col")
        annot_col = cfg["annot"].get("annot_col")
        process_annot(annot_file, annot_out_id, out_dir, snp_col, annot_col)

    if make_qq_annot:
        # process qq annot
        qq_annot_file = cfg["qq_annot"].get("qq_annot_file")
        qq_annot_out_id = cfg["qq_annot"].get("qq_annot_out_id")
        out_dir = out_dir
        process_qq_annot(qq_annot_file, qq_annot_out_id, out_dir)


    if make_sumstats:
        #process sumstats
        sumstats_file = cfg["sumstats"].get("sumstats_file")
        sumstats_out_id = cfg["sumstats"].get("sumstats_out_id")
        out_dir = out_dir
        snp_col = cfg["sumstats"].get("snp_col")
        p_col = cfg["sumstats"].get("p_col")
        effect_col = cfg["sumstats"].get("effect_col")
        effect_baseline = cfg["sumstats"].getfloat("effect_baseline") # e.g. 0 for BETA, 1 for OR
        ncase_col = cfg["sumstats"].get("ncase_col", fallback=None)
        ncontrol_col = cfg["sumstats"].get("ncontrol_col", fallback=None)
        n_col = cfg["sumstats"].get("n_col", fallback=None)
        process_sumstats(sumstats_file, sumstats_out_id, out_dir, snp_col, p_col,
            effect_col, effect_baseline, n_col=n_col, ncase_col=ncase_col,
            ncontrol_col=ncontrol_col)

    print(f"makeinput.py finished at {datetime.datetime.now()}")
