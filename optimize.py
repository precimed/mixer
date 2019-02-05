import numpy as np
import pandas as pd
import configparser
import sys
import time
import datetime
import os
import re
import cmmcost_omp.cmmcost_omp as cmmcost_omp
import test.cmmutils as cmmutils
from scipy.optimize import minimize
import scipy.stats as sstats


def _drop_duplicates(df):
    i = df.index.duplicated(keep='first')
    if i.any(): df = df.loc[~i,:]
    return df


def load_data(bim_files, frq_files, ldnpz_files, sumstats_f, annot_f,
    annot2use=None, snp2use_f=None, subset_size=None, subset_seed=None):
    """
    Args:
        bim_files   [K] (str): list of template bim files
        frq_files   [K] (str): list of template frq files
        ldnpz_files [K] (str): list of files produced by makeinput.ld2npz
        sumstats_f      (str): sumstats file, must contain columns: SNP, Z, N.
                               Only SNPs with z-score in sumstats are used for
                               fitting
        annot_f         (str): annotation file, must contain header with SNP
                               column and at least one 1/0 annotation column.
                               Only SNPs with annotation are used fo fitting
        annot2use   [M] (str): list of annotation categories to use, if None all
                               except SNP columns form annot_f are interpreted
                               as annotations
        snp2use_f       (str): no-headr single-column file with SNP ids to use.
                               If None all SNPs survived so far are used for
                               fitting
        subset_size     (int): number of randomly selected SNPs. This applies
                               after overlapping template with sumstats, annot
                               and snp2use files. If None all SNPs survived so
                               far are used for fitting
        subset_seed     (int): seed for random subsampling
    Return:
        snp         [N]   (str): SNP ids to use for fitting
        snp_z       [N]   (f8) : z-scores of SNPs in snp array
        s2          [L]   (f4) : vector of het*r2*ssize, L = is2[-1]
        is2         [N+1] (u8) : vector of LD-block indices in s2 array,
                                 is2[i] = start of LD-block for snp[i]
        annot_s2    [L]   (u1) : annotation codes for SNPs in s2 array
                                 0 <= i <= M corresponds to annot[i]
        annot_names [M]   (str): annotation categories
    """
    print("Loading data")

    # process sumstats
    sumstats_df = pd.read_csv(sumstats_f, sep="\t", usecols=["SNP", "Z", "N"],
        index_col="SNP", dtype={"SNP":str, "Z":'f8', "N":'f8'}, na_filter=False)
    print(f"{sumstats_df.shape[0]} SNPs in {sumstats_f}")
    sumstats_df = _drop_duplicates(sumstats_df)
    print(f"{sumstats_df.shape[0]} non-duplicated SNPs in {sumstats_f}")

    # process annot
    annot_usecols = None if (annot2use is None) else (["SNP"] + annot2use)
    annot_df = pd.read_csv(annot_f, sep="\t", usecols=annot_usecols,
        index_col="SNP", dtype={"SNP":str}, na_filter=False)
    annot_df = annot_df.astype(bool)
    annot_names = annot_df.columns.to_numpy()
    print(f"{annot_df.shape[0]} SNPs in {annot_f}")
    print(f"{annot_names.size} annotation categories in {annot_f}")
    print(f"annotation categories: {', '.join(annot_names)}")
    annot_df = _drop_duplicates(annot_df)
    print(f"{annot_df.shape[0]} non-duplicated SNPs in {annot_f}")

    commonsnp = sumstats_df.index.intersection(annot_df.index)
    print(f"{commonsnp.size} common SNPs in sumstats and annot")

    # process snp2use
    if not (snp2use_f is None):
        snp2use = pd.read_csv(snp2use_f, header=None, names=["SNP"],
            squeeze=True, na_filter=False)
        print(f"{snp2use.size} SNPs in {snp2use_f}")
        snp2use.drop_duplicates(keep='first', inplace=True)
        print(f"{snp2use.size} non-duplicated SNPs in {snp2use_f}")
        snp2use = pd.Index(snp2use)
        commonsnp = commonsnp.intersection(snp2use)
        print(f"{commonsnp.size} common SNPs in sumstats, annot and snp2use")

    # read SNPs from template and overlap with commonsnp
    template_snp = [pd.read_csv(f, sep="\t", usecols=[1], header=None,
                                squeeze=True, na_filter=False, dtype=str)
                    for f in bim_files]
    template_snp = pd.concat(template_snp)
    print(f"{template_snp.size} SNPs in template")
    commonsnp = commonsnp.intersection(template_snp)
    print(f"{commonsnp.size} common SNPs in sumstats, annot, snp2use and template")

    # get annotation codes for template SNPs
    annot_df = annot_df.loc[template_snp,:]
    #TODO: improve the following check
    # check whether each template SNP belongs to one and only one annotation category
    if not (annot_df.sum(axis=1) == 1).all():
        wrong_annot = (annot_df.sum(axis=1) != 1)
        wrong_annot_5snp = annot_df.index[wrong_annot][:10] # take max first 5
        raise ValueError("All template SNPs must belong to one and only one " +
            f"annotation category! Check these SNPs: {', '.join(wrong_annot_5snp)}")
    #get annotation categories for template SNPs coded as index of 1/True in annot_df
    annot = pd.Series(annot_df.to_numpy().nonzero()[1], index=annot_df.index, dtype='u1')

    # take a random subset from commonsnp
    if not (subset_size is None):
        if not (subset_seed is None):
            np.random.seed(subset_seed)
        commonsnp = commonsnp[np.random.choice(commonsnp.size, subset_size, replace=False)]
        print(f"{commonsnp.size} random SNPs are taken")
    print(f"{commonsnp.size} SNPs will be used for fitting")

    # take further only required SNPs from sumstats
    sumstats_df = sumstats_df.loc[commonsnp,:]

    # estimate len of s2 and get indices of common SNPs
    common_idx = [] # common_idx = [idx_chr_1, idx_chr_2, ...],
                    # idx_chr_i =  indices of commonsnp in template
    bs = [] # bs = [bs_chr_1, bs_chr_2, ...],
            # bs_chr_i = LD block sizes of SNPs in bim_f[i]
    for bim_f, ldnpz_f in zip(bim_files, ldnpz_files):
        # assuming no duplicates in bim file
        bim_snp = pd.read_csv(bim_f, sep="\t", usecols=[1], header=None,
                squeeze=True, na_filter=False)
        i = bim_snp.isin(commonsnp).to_numpy().nonzero()[0].astype('u4')
        common_idx.append(i)

        ldnpz = np.load(ldnpz_f)
        bs.append(ldnpz.get("bs"))
    len_s2 = sum(b[common_idx[i]].sum() for i,b in enumerate(bs))
    print(f"{len_s2} elements in s2 array")
    
    bs_all = np.concatenate(bs)
    print(f"{bs_all.mean():.2f} SNPs in LD block on average")
    print(f"{bs_all.max()} SNPs in the largest LD block")
    is2 = np.zeros(template_snp.size+1, dtype='u8')
    is2[1:] = bs_all.cumsum(dtype='u8') 

    snp = []
    snp_z = []
    s2 = np.zeros(len_s2, dtype='f4')
    annot_s2 = np.zeros(len_s2, dtype='u1')

    s2_2use_i_start = 0

    for i, (frq_f, ldnpz_f) in enumerate(zip(frq_files, ldnpz_files)):
        print(f"processing {ldnpz_f}")
        # assuming the same order of SNPs in bim and frq files
        ldnpz = np.load(ldnpz_f)
        maf = pd.read_csv(frq_f, delim_whitespace=True, usecols=["SNP", "MAF"],
            index_col="SNP", squeeze=True, na_filter=False,
            dtype={"SNP":str, "MAF":'f8'})
        commonsnp_chri = maf.index[common_idx[i]] # SNP ids for i-th chr in the corresponding order
        print(f"{commonsnp_chri.size} common SNPs")
        snp.append(commonsnp_chri.to_numpy())
        snp_z.append(sumstats_df.loc[commonsnp_chri,"Z"].to_numpy())

        het_i = (2*maf*(1-maf)).to_numpy()
        annot_i = annot.loc[maf.index].to_numpy()
        is2_i = np.zeros(het_i.size+1, dtype='u8')
        is2_i[1:] = bs[i].cumsum(dtype='u8')

        r2_i = ldnpz.get("r2")
        i2_i = ldnpz.get("i2")

        ssize_i = sumstats_df.loc[commonsnp_chri, "N"].to_numpy()
        ind_2use_i = common_idx[i]
        is2_2use_i = np.zeros(ind_2use_i.size+1, dtype='u8')
        is2_2use_i[1:] = bs[i][ind_2use_i].cumsum(dtype='u8')

        s2_2use_i_end = int(s2_2use_i_start + is2_i[-1]) # w/o int() result of + has float type ??
        s2_2use_i = s2[s2_2use_i_start:s2_2use_i_end]
        annot_s2_2use_i = annot_s2[s2_2use_i_start:s2_2use_i_end]

        cmmutils.fill_s2(het_i, annot_i, is2_i, r2_i, i2_i, ind_2use_i,
                ssize_i, is2_2use_i, s2_2use_i, annot_s2_2use_i)

        s2_2use_i_start = s2_2use_i_end

    snp = np.concatenate(snp)
    snp_z = np.concatenate(snp_z)

    s2_min, s2_mean, s2_max = cmmcost_omp.get_nonzero_min_mean_max(s2)
    print(f"{s2_min:.2f} = min(s2)")
    print(f"{s2_mean:.2f} = mean(s2)")
    print(f"{s2_max:.2f} = max(s2)")

    return snp, snp_z, s2, is2, annot_s2, annot_names


def process_input(sumstats_f, annot_f, template_dir):
    """
    Args:
        sumstats_f: sumstats npz file produced by makeinput.py
        annot_f: annotation npz file produced by makeinput.py
        template_dir: dir with template.chr*.npz files produced by makeinput.py

    Get input in the format produced by makeinput.py and produce input in the
    format required for optimization, i.e.:
    template_snp: array of SNP ids
    template_z: [f8], vector of z-scores for the template (0 if not in template)
    template_snp_in_sumstats: [u1], vector indicating whether snp from template
        is also in sumstats (i.e. has z-score). 1 if in sumstats 0 if not.
    template_annot: [u1], vector of annotations  for the template
    s2: [f4], vector of het*r2*ssize
    is2: [u8], vector of ld-block indices in s2 values,
        s2[is2[i]:is2[i+1]] = s2 values for i-th variant from template
    annot_s2: [u1], vector of annotations for s2 values
    n_categories: int, number of annotation categories
         = len(unique(template_annot))
    """
    print("Getting input data")
    # define template
    chr2use = []
    p = re.compile('template.chr([0-9]{1,2}).npz')
    for n in os.scandir(template_dir):
        m = re.match(p, n.name)
        if m:
            chr2use.append(m.group(1))
    chr2use.sort(key=lambda x: int(x))
    print(f"    {len(chr2use)} chromosomes in {template_dir}")

    # read sumstats
    sumstats = np.load(sumstats_f)
    snp_z_dict = dict(zip(sumstats["snp"], sumstats["z"]))
    snp_ssize_dict = dict(zip(sumstats["snp"], sumstats["ssize"]))

    # read annotations
    annotations = np.load(annot_f)
    annot_categories = annotations.get("categories")
    n_categories = len(annot_categories)
    snp_annot_dict = dict(zip(annotations["snp"], annotations["annot"]))

    # get the number of elements in s2/annot_s2 arrays
    len_s2 = 0
    for c in chr2use:
        chr_template_file = os.path.join(template_dir, f"template.chr{c}.npz")
        template = np.load(chr_template_file)
        len_s2 += len(template["r2"])

    # create template_snp, template_snp_in_sumstats, z, is2, s2, annot_s2 and template_annot arrays
    template_snp = []
    template_snp_in_sumstats = []
    template_z = []
    template_annot = []
    is2 = [[0]]
    s2 = np.zeros(len_s2, dtype='f4')
    annot_s2 = np.zeros(len_s2, dtype='u1')
    start_i = 0
    for c in chr2use:
        chr_template_file = os.path.join(template_dir, f"template.chr{c}.npz")
        template = np.load(chr_template_file)
        snp = template.get("snp")
        ld_bs = template.get("ld_bs")
        r2 = template.get("r2")
        het = template.get("het")
        r2_21i = template.get("r2_21i")
        template_snp.append(snp)
        # can store r2_12i directly in template.npz file
        r2_12i = np.repeat(np.arange(len(snp), dtype='u4'), ld_bs)
        ssize = np.array([snp_ssize_dict.get(s, 0) for s in snp])
        is2.append(ld_bs)
        template_snp_in_sumstats += [s in snp_z_dict for s in snp]
        template_z += [snp_z_dict.get(s, 0) for s in snp]
        try:
            template_annot_chr = np.array([snp_annot_dict[s] for s in snp], dtype='u1')
        except KeyError as ke:
            print(f"Some SNPs from template are not annotated: {ke}")
            raise(ke)
        template_annot.append(template_annot_chr)
        end_i = start_i + len(r2)
        s2[start_i:end_i] = (ssize[r2_12i]*het[r2_21i])*r2
        annot_s2[start_i:end_i] = template_annot_chr[r2_21i]
        start_i = end_i

    template_snp = np.concatenate(template_snp)
    template_annot = np.concatenate(template_annot)
    template_snp_in_sumstats = np.array(template_snp_in_sumstats, dtype='u1')
    template_z = np.array(template_z, dtype='f8')
    is2 = np.concatenate(is2)
    mean_bs = is2[1:].mean()
    max_bs = is2[1:].max()
    is2 = is2.cumsum(dtype='u8') # if not specified np.cumsum automatically switches dtype to 'i8'

    s2_min, s2_mean, s2_max = cmmcost_omp.get_nonzero_min_mean_max(s2)

    print(f"    {len(template_snp)} variants in template")
    print(f"    {len_s2} variants in all LD neighbourhoods = len(s2)")
    print(f"    {mean_bs:.2f} variants in LD neighbourhood on average")
    print(f"    {max_bs} variants in the largest LD neighbourhood")
    print(f"    {s2_max} max element in s2 = ssize*r2*het")
    print(f"    {s2_mean:.9f} mean element in s2")
    print(f"    {s2_min} min element in s2")
    print(f"    {n_categories} annotation categories:")
    print(f"        {', '.join(annot_categories)}")
    print(f"    {len(snp_z_dict)} variants in sumstats")
    print(f"    {template_snp_in_sumstats.sum()} variants in the overlap of the template and sumstats")

    # ensure that types are correct
    assert template_z.dtype == 'f8'
    assert template_snp_in_sumstats.dtype == 'u1'
    assert s2.dtype == 'f4'
    assert is2.dtype == 'u8'
    assert template_annot.dtype == 'u1'
    assert annot_s2.dtype == 'u1'

    return template_snp, template_z, template_snp_in_sumstats, template_annot, s2, is2, annot_s2, annot_categories 


def get_qq_annot4snp(snp, qq_annot_f):
    """
    Get qq nnotations for SNPs. Raise exception if some SNPs are
    not annotated.
    Args:
        snp: list of SNP ids
        qq_annot_f: qq annotation npz file produced by makeinput.py
    Return:
        [np.bool] matrix NxM, where N = number of SNPs, M = number of categories
        names of annotation categories
    """
    print(f"Getting qq annotations from {qq_annot_f}")
    qq_annot = np.load(qq_annot_f)
    qq_df = pd.DataFrame(index=qq_annot["snp"], columns=qq_annot["categories"],
            data=qq_annot["qq_annot"])
    qq_df = qq_df.loc[snp,:]
    return qq_df.values, qq_df.columns


def get_z_cdf_2tails(z_grid, template_snp_in_sumstats, n_samples, p, sb2, s02, s2, is2, annot_s2,
        qq_template_annot):
    # Function to produce model data for QQ plot
    print("Getting z-score CDF")
    print(f"{template_snp_in_sumstats.sum()} SNPs will be used")

    z_cdf_total = np.zeros(len(z_grid))
    n_categories = qq_template_annot.shape[1]
    z_cdf_annot = np.zeros((n_categories, len(z_grid)))

    # uniform weights
    total_c = template_snp_in_sumstats.sum()
    annot_c = qq_template_annot.sum(axis=0).reshape((n_categories,1)) # sum of SNPs in each annot

    max_block_size = max([(is2[i+1]-is2[i]) for i in range(len(template_snp_in_sumstats))])
    print(f"Largest LD block size: {max_block_size}")
    rg = np.random.rand(n_samples, max_block_size)
    for i_templ, (istart, iend) in enumerate(zip(is2[:-1], is2[1:])):
        if i_templ%10000 == 0: print(f"{i_templ} variants processed")
        if template_snp_in_sumstats[i_templ]:
            p_in_ld = p[annot_s2[istart:iend]]
            sb2_in_ld = sb2[annot_s2[istart:iend]]*s2[istart:iend]
            curr_ld_block_size = iend-istart
            sigma2 = (sb2_in_ld*(rg[:,:curr_ld_block_size] < p_in_ld)).sum(axis=1) + s02
            if len(sigma2) > 0:
                sigma2 = sigma2.reshape((n_samples, 1))
                z_cdf = sstats.norm.cdf(z_grid, 0, np.sqrt(sigma2)).sum(axis=0)/len(sigma2) # estimate z_cdf as a cdf of equally weighted mixture of n_samples elements with sigma2 variance
                # z_cdf_total += z_cdf*weight_total[i_templ]
                # z_cdf_annot[qq_template_annot[i_templ]] += z_cdf*weight_annot[i_templ]
                z_cdf_total += z_cdf
                z_cdf_annot[qq_template_annot[i_templ]] += z_cdf
    # if weighted, z_cdf_total and z_cdf_annot should be normalized by the sum of all weights in the corresponding annotation, i.e. z_cdf_total /= sum(weights_all), z_cdf_annot[i] /= sum(weights_annot[i])
    # see cmma.jl code for more details
    # multiply by 2 since we want to have 2 tails
    return 2*z_cdf_total/total_c, 2*z_cdf_annot/annot_c


def get_xy_from_p(p, p_weights=None, nbins=200):
    """
    Thins function is taken from qq.py (excluding top_as_dot argument)
    """
    if p_weights is None:
        p_weights = np.ones(len(p))
    p_weights /= sum(p_weights) # normalize weights

    i = np.argsort(p)
    p = p[i]
    p_weights = p_weights[i]
    cum_p_weights = np.cumsum(p_weights)

    y = np.logspace(np.log10(p[-1]), np.log10(p[0]), nbins)
    y_i = np.searchsorted(p, y, side='left')
    y_i[0] = len(p) - 1  # last index in cum_p_weights
    y_i[-1] = 0
    p_cdf = cum_p_weights[y_i]
    x = -np.log10(p_cdf)
    y = -np.log10(y)
    return x, y


def get_params(opt_result_file):
    print(f"Getting parameters from {opt_result_file}")
    res = np.load(opt_result_file)
    p = res.get("p_opt")
    sb2 = res.get("sb2_opt")
    s02 = res.get("s02_opt")
    print(f"p: {p}")
    print(f"sb2: {sb2}")
    print(f"s02: {s02}")
    return p, sb2, s02


def overlap(z2use, template_snp, snp2use_f):
    """
    The function doesn't modify z2use array but creates a copy and modifies it
    Args:
        z2use: ['u1'] array indicating which snp to use for optimization
        template_snp: array of snp ids from template
        snp2use_f: a file with SNP ids (a single column with one id per line).
            An overlap between SNPs from template and SNPs from this file
            will be marked as "1" in z2use array
    Return:
        z2use: new z2use array, where all snps from template which are not in
            snp2use_f are marked with 0
    """
    print(f"Overlapping template snps with variants from {snp2use_f}")
    z2use_overlap = z2use.copy()
    with open(snp2use_f) as f:
        snpids = set(map(str.rstrip, f))
        print(f"    {len(snpids)} SNPs in {snp2use_f}")
    snp_overlap = np.array([s in snpids for s in template_snp])
    print(f"    {snp_overlap.sum()} SNPs in the overlap of template and {snp2use_f}")
    z2use_overlap[~snp_overlap] = 0
    print(f"    {z2use_overlap.sum()} SNPs in z2use after overlap with {snp2use_f}")
    return z2use_overlap


def randsubset(z2use, subset_size=100000, seed=None):
    """
    The function doesn't modify z2use but creates new array and modifies it.
    Return z2use array with subset_size True values. If z2use[i] = False it
    stays False in the returned array.
    """
    # get indices of True in z2use
    print("Getting random subset of z2use")
    assert z2use.sum() >= subset_size
    print(f"    subset size: {subset_size}")
    if not seed is None:
        print(f"    setting np.random.seed to {seed}")
        np.random.seed(seed)
    nz_i = np.nonzero(z2use)[0]
    # get subset of current True indices
    nz_i = np.random.choice(nz_i, subset_size, replace=False)
    z2use_subset = np.zeros(len(z2use), dtype=z2use.dtype)
    z2use_subset[nz_i] = 1
    return z2use_subset


def rand_initial_guess(n_categories, same_pi, same_sb2, seed, min_p=1e-4, max_p=1e-2,
    min_sb2=1e-5, max_sb2=1e-3, min_s02=1, max_s02=1.4):
    print("Generating random initial guess")
    if not seed is None:
        print(f"    setting np.random.seed to {seed}")
        np.random.seed(seed)
    p0 = np.random.uniform(min_p, max_p, n_categories)
    if same_pi:
        p0 = p0[:1]
    sb20 = np.random.uniform(min_sb2, max_sb2, n_categories)
    if same_sb2:
        sb20 = sb20[:1]
    s020 = np.random.uniform(min_s02, max_s02)
    print(f"    p_0: {p0}")
    print(f"    sb2_0: {sb20}")
    print(f"    s02_0: {s020}")
    return p0, sb20, s020



def process_idump(dump_input_file):
    print(f"Loading dumped imnput from {dump_input_file}")
    idump = np.load(dump_input_file)
    template_snp = idump.get("template_snp") 
    template_z = idump.get("template_z")
    template_snp_in_sumstats = idump.get("template_snp_in_sumstats")
    template_annot = idump.get("template_annot")
    s2 = idump.get("s2")
    is2 = idump.get("is2")
    annot_s2 = idump.get("annot_s2")
    annot_categories = idump.get("annot_categories")
    return template_snp, template_z, template_snp_in_sumstats, template_annot, s2, is2, annot_s2, annot_categories


def logistic(x):
    """ Logistic function. Maps [-∞; ∞] -> [0; 1].
    """
    return 1/(1 + np.exp(-x))


def logit(x):
    """ Inverse logistic function (logit). Maps [0; 1] -> [-∞; ∞].
    """
    return np.log(x/(1 - x))


def objective_func(x, z, z2use, s2, is2, annot_s2, n_categories):
    assert len(x) == 2*n_categories + 1, "Wrong x input in objective function"
    p = logistic(x[:n_categories])
    sb2 = np.exp(x[n_categories:2*n_categories])
    s02 = np.exp(x[-1])

    print(f"pi: {p}, sb2: {sb2}, s02: {s02}")

    cost = cmmcost_omp.get_cost(z, z2use, s2, is2, p, sb2, s02, annot_s2)
    print(f"cost: {cost}")
    return cost


def objective_func_same_sb2(x, z, z2use, s2, is2, annot_s2, n_categories):
    assert len(x) == n_categories + 2, "Wrong x input in objective function"
    p = logistic(x[:n_categories])
    sb2 = np.exp(np.repeat(x[-2],n_categories))
    s02 = np.exp(x[-1])

    print(f"pi: {p}, sb2: {sb2}, s02: {s02}")

    cost = cmmcost_omp.get_cost(z, z2use, s2, is2, p, sb2, s02, annot_s2)
    print(f"cost: {cost}")
    return cost


def objective_func_same_pi(x, z, z2use, s2, is2, annot_s2, n_categories):
    assert len(x) == n_categories + 2, "Wrong x input in objective function"
    p = logistic(np.repeat(x[0], n_categories))
    sb2 = np.exp(x[1:-1])
    s02 = np.exp(x[-1])

    print(f"pi: {p}, sb2: {sb2}, s02: {s02}")

    cost = cmmcost_omp.get_cost(z, z2use, s2, is2, p, sb2, s02, annot_s2)
    print(f"cost: {cost}")
    return cost


def objective_func_same_pi_same_sb2(x, z, z2use, s2, is2, annot_s2, n_categories):
    assert len(x) == 3, "Wrong x input in objective function"
    p = logistic(np.repeat(x[0], n_categories))
    sb2 = np.exp(np.repeat(x[1],n_categories))
    s02 = np.exp(x[-1])

    print(f"pi: {p}, sb2: {sb2}, s02: {s02}")

    cost = cmmcost_omp.get_cost(z, z2use, s2, is2, p, sb2, s02, annot_s2)
    print(f"cost: {cost}")
    return cost


def run_optimization(p0, sb20, s020, z, z2use, s2, is2, annot_s2, annot_categories,
    adaptive, maxiter, same_pi, same_sb2, opt_result_file):
    """
    Start optimization with initial guess given by p0, sb20 and s020 arguments.

    More info:
    https://docs.scipy.org/doc/scipy/reference/optimize.minimize-neldermead.html
    https://github.com/scipy/scipy/pull/5205
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
    """
    print("Starting optimization")
    x0 = np.concatenate([logit(p0), np.log(sb20), [np.log(s020)]]) # initial guess in opt space
    n_categories = len(annot_categories)
    args = (z, z2use, s2, is2, annot_s2, n_categories)
    print(f"Use adaptive strategy: {adaptive}")
    print(f"Maximum number of iterations: {maxiter}")

    print(f"Use same pi for all categories: {same_pi}")
    print(f"Use same sb2 for all categories: {same_sb2}")

    # choose proper objective function
    if same_pi and same_sb2:
        ofunc = objective_func_same_pi_same_sb2
    elif same_pi:
        ofunc = objective_func_same_pi
    elif same_sb2:
        ofunc = objective_func_same_sb2
    else:
        ofunc = objective_func

    res = minimize(ofunc, x0, args, method='Nelder-Mead',
        options={'fatol':1e-7, 'xatol':1e-5, 'adaptive':adaptive, 'maxiter':maxiter})

    if same_pi and same_sb2:
        p_opt = logistic(np.repeat(res.x[0],n_categories))
        sb2_opt = np.exp(np.repeat(res.x[1],n_categories))
    elif same_pi:
        p_opt = logistic(np.repeat(res.x[0],n_categories))
        sb2_opt = np.exp(res.x[1:-1])
    elif same_sb2:
        p_opt = logistic(res.x[:n_categories])
        sb2_opt = np.exp(np.repeat(res.x[-2],n_categories))
    else:
        p_opt = logistic(res.x[:n_categories])
        sb2_opt = np.exp(res.x[n_categories:2*n_categories])
    s02_opt = np.exp(res.x[-1])
    print(f"Opt success: {res.success}")
    print(f"Solution: {res.x}")
    print(f"p_opt: {p_opt}")
    print(f"sb2_opt: {sb2_opt}")
    print(f"s02_opt: {s02_opt}")
    print(f"Final obj func value: {res.fun}")
    print(f"Termination message: {res.message}")
    print(f"Number of iterations: {res.nit}")
    print(f"Number of function evaluations: {res.nfev}")
    print(f"Final simplex: {res.final_simplex}")

    np.savez(opt_result_file, success=res.success, x=res.x, x0=x0, p_opt=p_opt,
        sb2_opt=sb2_opt, s02_opt=s02_opt, fun=res.fun, nit=res.nit, nfev=res.nfev,
        final_simplex_x=res.final_simplex[0],final_simplex_fun=res.final_simplex[1],
        annot_categories=annot_categories)
    print(f"Optimization result saved to {opt_result_file}")


if __name__ == "__main__":
    template_dir = "/mnt/seagate10/genotypes/1000genomes503eur9m"
    chr2use = [str(i) for i in range(1,23)]

    bim_files = [os.path.join(template_dir, f"chr{c}.bim") for c in chr2use]
    frq_files = [os.path.join(template_dir, f"chr{c}.maf.frq") for c in chr2use]
    ldnpz_files = [os.path.join(template_dir, "ld2npz", f"chr{c}.ld.npz") for c in chr2use]
    # /mnt/seagate10/sumstats/DIAGRAM_T2D_2018_adjBMI/DIAGRAM_T2D_2018_adjBMI.sumstats.gz
    # /mnt/seagate10/sumstats/GIANT_HEIGHT_2018_UKB/GIANT_HEIGHT_2018_UKB.sumstats.gz
    sumstats_f = "/mnt/seagate10/sumstats/GIANT_HEIGHT_2018_UKB/GIANT_HEIGHT_2018_UKB.sumstats.gz"
    annot_f = "/mnt/seagate10/genotypes/1000genomes503eur9m/annot/1000genomes503eur9m.utr5utr3exon_other.binary.txt"
    annot2use = None
    snp2use_f = os.path.join(template_dir, "snps.no_mhc.txt")
    subset_size = 500000
    subset_seed = 1
    load_data(bim_files, frq_files, ldnpz_files, sumstats_f, annot_f,
        annot2use, snp2use_f, subset_size, subset_seed)


# if __name__ == "__main__":
#     print(f"optimize.py started at {datetime.datetime.now()}")
#     print(f"Reading config from {sys.argv[1]}")
#     cfg = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
#     cfg.read(sys.argv[1])

#     os.environ["OMP_NUM_THREADS"] = cfg["omp"].get("OMP_NUM_THREADS")
#     print(f"OMP_NUM_THREADS is set to {os.environ['OMP_NUM_THREADS']}")

#     # input parameters
#     template_dir = cfg["general"].get("template_dir")
#     sumstats_f = cfg["general"].get("sumstats_f")
#     annot_f = cfg["general"].get("annot_f")
#     opt_result_file = cfg["general"].get("opt_result_file")

#     load_idump = cfg["dump"].getboolean("load_idump")    
#     dump_input = cfg["dump"].getboolean("dump_input")
#     dump_input_file = cfg["dump"].get("dump_input_file")
#     if load_idump:
#         template_snp, template_z, template_snp_in_sumstats, template_annot, s2, is2, annot_s2, annot_categories = process_idump(dump_input_file)
#     else:
#         template_snp, template_z, template_snp_in_sumstats, template_annot, s2, is2, annot_s2, annot_categories = process_input(sumstats_f, annot_f, template_dir)
#     if dump_input:
#         print(f"Dumping input to {dump_input_file}")
#         np.savez(dump_input_file, template_snp=template_snp, template_z=template_z, template_snp_in_sumstats=template_snp_in_sumstats, template_annot=template_annot,
#             s2=s2, is2=is2, annot_s2=annot_s2, annot_categories=annot_categories)

#     run_opt = cfg["optimization"].getboolean("run_opt")
#     if run_opt:
#         print("Preparing optimization")
#         adaptive = cfg["optimization"].getboolean("adaptive")
#         maxiter = cfg["optimization"].getint("maxiter")
#         same_pi = cfg["optimization"].getboolean("same_pi")
#         same_sb2 = cfg["optimization"].getboolean("same_sb2")
#         snp2use_f = cfg["optimization"].get("snp2use_f", fallback=None)
#         subset_size = cfg["optimization"].getint("subset_size", fallback=None)
#         subset_seed = cfg["optimization"].getint("subset_seed", fallback=None)
#         rand_init_seed = cfg["optimization"].getint("rand_init_seed", fallback=None)

#         n_categories = len(annot_categories)
#         p, sb2, s02 = rand_initial_guess(n_categories, same_pi, same_sb2, rand_init_seed)
#         z2use = template_snp_in_sumstats.copy()
#         if snp2use_f:
#             z2use = overlap(z2use, template_snp, snp2use_f)
#         if subset_size:
#             z2use = randsubset(z2use, subset_size, subset_seed)
#         run_optimization(p, sb2, s02, template_z, z2use, s2, is2, annot_s2, annot_categories,
#                          adaptive, maxiter, same_pi, same_sb2, opt_result_file)


#     make_qq = cfg["qq"].getboolean("make_qq")
#     if make_qq:
#         print("Making QQ plot")
#         modelqq_out_file = cfg["qq"].get("modelqq_out_file")
#         template_qq_annot_file = cfg["qq"].get("template_annot_file")
#         n_samples = cfg["qq"].getint("n_samples")
#         opt_result_file = cfg["qq"].get("opt_result_file")
#         # get template annotations
#         qq_template_annot, annot_names = get_qq_annot4snp(template_snp, template_qq_annot_file)
#         # parameters derived from data
#         p, sb2, s02 = get_params(opt_result_file)
#         p_experimental = 2*sstats.norm.cdf(-np.abs(template_z))
#         y_max = min(50, max(-np.log10(p_experimental)))

#         # get model data
#         n_grid = 150
#         p_grid = np.logspace(-y_max,0,n_grid)
#         # multiply by 0.5 since we want to have two tailed quantiles    
#         z_grid = sstats.norm.ppf(0.5*p_grid)

#         # qq plot should be only for the SNPs which present in sumstats data, i.e. template_snp_in_sumstats
#         # TODO: allow taking only a subset of SNPs for the qq plot. This subset must be used both for modeled and experemental plot.
#         # z_cdf_total, z_cdf_annot = get_z_cdf_2tails(z_grid, template_snp_in_sumstats, n_samples,
#         #     p, sb2, s02, s2, is2, annot_s2, qq_template_annot)

#         z_cdf_total, z_cdf_annot = cmmcost_omp.get_cdfsampling(z_grid, template_snp_in_sumstats, s2, is2, p,
#                 sb2, s02, annot_s2, qq_template_annot.astype('u1'), n_samples)

#         model_total_x = -np.log10(z_cdf_total)
#         model_annot_x = -np.log10(z_cdf_annot)
#         model_y =  -np.log10(p_grid)

#         # get experimental data
#         data_total_x, data_total_y = get_xy_from_p(p_experimental[template_snp_in_sumstats.astype('bool')])
#         # take only SNPs which are in sumstats
#         data_total_x = data_total_x
#         data_total_y = data_total_y
#         data_annot_x = []
#         data_annot_y = []
#         for i in range(len(annot_names)):
#             annot_i = qq_template_annot[:,i].astype('bool')
#             # take only SNPs which are in sumstats
#             annot_i &= template_snp_in_sumstats.astype('bool')
#             p_experimental_annot = p_experimental[annot_i]
#             annot_x, annot_y = get_xy_from_p(p_experimental_annot)
#             data_annot_x.append(annot_x)
#             data_annot_y.append(annot_y)
#         data_annot_x = np.array(data_annot_x)
#         data_annot_y = np.array(data_annot_y)        

#         # save results
#         np.savez(modelqq_out_file, annot_names=annot_names, data_total_x=data_total_x,
#             data_total_y=data_total_y, data_annot_x=data_annot_x, data_annot_y=data_annot_y,
#             model_total_x=model_total_x, model_annot_x=model_annot_x, model_y=model_y)

#         print(f"QQ plot data saved to {modelqq_out_file}")


#     run_single = cfg["test"].getboolean("run_single")
#     if run_single:
#         print("Estimating test cost")
#         # estimate for a given set of parameters
#         n_categories = len(annot_categories)
#         p = np.array([0.005]*n_categories, dtype='f8') # [0.02, 0.005, 0.05, 0.001] [0.005]*n_categories
#         sb2 = np.array([0.0002]*n_categories, dtype='f8') # [5.786066041591348e-05, 5.837030778549846e-05, 5.703057621611015e-05, 5.6426373881400386e-05] [0.0002]*n_categories
#         s02 = 1.0 # 1.0
#         start = time.time()
#         cost = cmmcost_omp.get_cost(template_z, z2use, s2, is2, p, sb2, s02, annot_s2)
#         end = time.time()
#         print(f"{end-start} seconds per single cost function evaluation")
#         print(f"Cost is: {cost}")

#     print(f"optimize.py finished at {datetime.datetime.now()}")
