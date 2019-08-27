import numpy as np
import pandas as pd
import configparser
import sys
import time
import datetime
import os
import re
import c.cmmcost_omp as cmmcost_omp
from scipy.optimize import minimize, basinhopping, differential_evolution
import scipy.stats as sstats
import random


def process_input(template_dir, annot_f, sumstats_kwargs):
    """
    Args:
        template_dir: dir with template.chr*.npz files produced by makeinput.py
        annot_f: annotation npz file produced by makeinput.py
        sumstats_kwargs: "sumstats" section from the config file converted to dict

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
    sumstats_kwargs["effect_baseline"] = float(sumstats_kwargs["effect_baseline"])
    ss_snp, ss_z, ss_ssize = process_sumstats(**sumstats_kwargs)
    snp_z_dict = dict(zip(ss_snp, ss_z))
    snp_ssize_dict = dict(zip(ss_snp, ss_ssize))

    # read annotations
    annotations = np.load(annot_f, allow_pickle=True)
    annot_categories = annotations.get("categories")
    n_categories = len(annot_categories)
    snp_annot_dict = dict(zip(annotations["snp"], annotations["annot"]))

    # get the number of elements in s2/annot_s2 arrays
    len_s2 = 0
    for c in chr2use:
        chr_template_file = os.path.join(template_dir, f"template.chr{c}.npz")
        template = np.load(chr_template_file, allow_pickle=True)
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
        template = np.load(chr_template_file, allow_pickle=True)
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


def process_sumstats(sumstats_file, snp_col, p_col, effect_col, effect_baseline,
    n_col=None, ncase_col=None, ncontrol_col=None):
    print(f"Processing sumstats {sumstats_file}")
    usecols = [snp_col, p_col, effect_col]
    if (n_col is None) and (not ncase_col is None) and (not ncontrol_col is None):
        usecols += [ncase_col, ncontrol_col]
    elif (not n_col is None):
        usecols.append(n_col)
    else:
        raise ValueError(
            "Either 'n_col' or both 'ncase_col' and 'ncontrol_col' must be not None")
    df = pd.read_csv(sumstats_file, delim_whitespace=True, usecols=usecols, index_col=snp_col)
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
    df["Z"] = sstats.norm.ppf(0.5*df[p_col]) # here all Z are < 0
    df.loc[positive_effect,"Z"] *= -1

    if n_col is None:
        ssize = 4*df[ncase_col]*df[ncontrol_col]/(df[ncase_col]+df[ncontrol_col])
    else:
        ssize = df[n_col]
    ssize = ssize.values

    ss_snp = df.index.values
    ss_z = df.Z.values
    ss_ssize = ssize
    return ss_snp, ss_z, ss_ssize


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
    qq_annot = np.load(qq_annot_f, allow_pickle=True)
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


def get_xy_from_p(p, p_weights=None, nbins=500):
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
    res = np.load(opt_result_file, allow_pickle=True)
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


def rand_initial_guess(n_categories, same_pi, same_sb2, seed, min_p=1e-4, max_p=1e-1,
    min_sb2=1e-5, max_sb2=1e-2, min_s02=1, max_s02=1.5):
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
    print(f"    rand p_0: {p0}")
    print(f"    rand sb2_0: {sb20}")
    print(f"    rand s02_0: {s020}")
    return p0, sb20, s020


def estimate_qq_slope(z):
    """
    Initialize s02 as a b^2 where b is a slope of the linear regression for the
    QQ plot of theoretical standard normal quantiles vs. quantiles of z scores
    from the summary statistics. Middle 50% of the distribution is used, i.e.
    quantiles from 25% to 75% (with the step = 1%).
    This procedure returns s02_init which is always greater than the actual s02.
    Args:
        z: z-scores from summary statistics
    Return:
        s02_init: initial guess for s02 parameter
    """
    print("Estimating initila guess for s02 as a QQ plot slope")
    print(f"{len(z)} z-scores are used for the s02 initialization")
    quantiles = np.linspace(0.25,0.75,26)
    theor_q = sstats.norm.ppf(quantiles)
    sample_q = np.percentile(z, 100*quantiles)
    slope, intercept, r_value, p_value, std_err = sstats.linregress(theor_q, sample_q)
    print(f"slope: {slope:.4f}    intercept: {intercept:.4f}    p-value: {p_value:.2E}")
    if p_value > 1.E-20:
        print(f"WARNING. Fairly high p-value ({p_value}) may indicate that linear approximation doesn't work good.")
    s02_init = slope*slope
    print(f"s02_init = {s02_init:.4f}")
    return s02_init



def process_idump(dump_input_file):
    print(f"Loading dumped imnput from {dump_input_file}")
    idump = np.load(dump_input_file, allow_pickle=True)
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

    cost = cmmcost_omp.get_cost(z, z2use, s2, is2, p, sb2, s02, annot_s2, n_categories)
    print(f"cost: {cost}")
    return cost


def objective_func_grad(x, z, z2use, s2, is2, annot_s2, n_categories):
    assert len(x) == 2*n_categories + 1, "Wrong x input in objective function"
    p = logistic(x[:n_categories])
    sb2 = np.exp(x[n_categories:2*n_categories])
    s02 = np.exp(x[-1])

    print(f"pi: {p}, sb2: {sb2}, s02: {s02}")

    grad = cmmcost_omp.get_cost_der(z, z2use, s2, is2, p, sb2, s02, annot_s2, n_categories)

    # derivative of coordinate transformations
    # [1/(1 + exp(-x))]' = exp(-x)/(1 + exp(-x))**2
    p_transform_der = []
    for y in x[:n_categories]:
        if y < -100: p_transform_der.append(1)
        elif y > 100: p_transform_der.append(0)
        else: p_transform_der.append(np.exp(-y)/(1 + np.exp(-y))**2)

    # p_transform_der = [np.exp(-y)/(1 + np.exp(-y))**2 for y in x[:n_categories]]
    s2_transform_der = [np.exp(y) for y in x[n_categories:-1]]
    s02_transform_der = [np.exp(x[-1])]
    transform_der = p_transform_der + s2_transform_der + s02_transform_der
    grad = np.array([g*d for g,d in zip(grad, transform_der)])

    print(f"gradient: [{', '.join('%.8f' % e for e in grad)}]")
    return grad


def objective_func_same_sb2(x, z, z2use, s2, is2, annot_s2, n_categories):
    assert len(x) == n_categories + 2, "Wrong x input in objective function"
    p = logistic(x[:n_categories])
    sb2 = np.exp(np.repeat(x[-2],n_categories))
    s02 = np.exp(x[-1])

    print(f"pi: {p}, sb2: {sb2}, s02: {s02}")

    cost = cmmcost_omp.get_cost(z, z2use, s2, is2, p, sb2, s02, annot_s2, n_categories)
    print(f"cost: {cost}")
    return cost


def objective_func_same_pi(x, z, z2use, s2, is2, annot_s2, n_categories):
    assert len(x) == n_categories + 2, "Wrong x input in objective function"
    p = logistic(np.repeat(x[0], n_categories))
    sb2 = np.exp(x[1:-1])
    s02 = np.exp(x[-1])

    print(f"pi: {p}, sb2: {sb2}, s02: {s02}")

    cost = cmmcost_omp.get_cost(z, z2use, s2, is2, p, sb2, s02, annot_s2, n_categories)
    print(f"cost: {cost}")
    return cost


def objective_func_same_pi_same_sb2(x, z, z2use, s2, is2, annot_s2, n_categories):
    assert len(x) == 3, "Wrong x input in objective function"
    p = logistic(np.repeat(x[0], n_categories))
    sb2 = np.exp(np.repeat(x[1],n_categories))
    s02 = np.exp(x[-1])

    print(f"pi: {p}, sb2: {sb2}, s02: {s02}")

    cost = cmmcost_omp.get_cost(z, z2use, s2, is2, p, sb2, s02, annot_s2, n_categories)
    print(f"cost: {cost}")
    return cost


def run_optimization(p0, sb20, s020, z, z2use, s2, is2, annot_s2, annot_categories,
    local_opt_method, local_opt_options, global_opt, global_opt_options, same_pi,
    same_sb2, opt_result_file):
    """
    Start optimization with initial guess given by p0, sb20 and s020 arguments.

    More info:
    https://docs.scipy.org/doc/scipy/reference/optimize.minimize-neldermead.html
    https://github.com/scipy/scipy/pull/5205
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
    """
    # np.seterr(all="raise")
    print("Starting optimization")
    print(f"{z2use.sum()} variants are used for optimization")
    n_categories = len(annot_categories)
    args = (z, z2use, s2, is2, annot_s2, n_categories)
    print(f"Perform initial global optimization: {global_opt}")
    if global_opt:
        print("Global optimization options:")
        for key, val in global_opt_options.items():
            print(f"   {key}: {val}")

    print(f"Using {local_opt_method} local optimization method with options:")
    for key, val in local_opt_options.items(): print(f"   {key}: {val}")

    print(f"Use same pi for all categories: {same_pi}")
    print(f"Use same sb2 for all categories: {same_sb2}")

    # choose proper objective function
    if same_pi and same_sb2:
        if local_opt_method == "BFGS":
            raise NotImplementedError(f"BFGS is only implemented for same_pi = False and same_sb2 = False")
        ofunc = objective_func_same_pi_same_sb2
    elif same_pi:
        if local_opt_method == "BFGS":
            raise NotImplementedError(f"BFGS is only implemented for same_pi = False and same_sb2 = False")
        ofunc = objective_func_same_pi
    elif same_sb2:
        if local_opt_method == "BFGS":
            raise NotImplementedError(f"BFGS is only implemented for same_pi = False and same_sb2 = False")
        ofunc = objective_func_same_sb2
    else:
        ofunc = objective_func

    if global_opt:
        print("Starting global optimization")
        # bounds in opt space for preliminary scan
        p_bounds = (logit(5.E-5), logit(5.E-1))
        s2_bounds = (np.log(5.E-6), np.log(5.E-2))
        s02_bounds = (np.log(0.9), np.log(2.5))
        # p_bounds = (logit(1.E-3), logit(1.E-1))
        # s2_bounds = (np.log(1.E-4), np.log(1.E-2))
        # s02_bounds = (np.log(1.), np.log(2.5))
        bounds4opt = []
        if same_pi: bounds4opt.append(p_bounds)
        else: bounds4opt += [p_bounds]*n_categories
        if same_sb2: bounds4opt.append(s2_bounds)
        else: bounds4opt += [s2_bounds]*n_categories
        bounds4opt.append(s02_bounds)

        res_global = differential_evolution(ofunc, bounds4opt, args=args, tol=0.01,
            mutation=(0.5, 1), recombination=0.7, atol=0, updating='immediate',
            polish=False, workers=1, **global_opt_options)
        x0 = res_global.x
    else:
        x0 = np.concatenate([logit(p0), np.log(sb20), [np.log(s020)]]) # initial guess in opt space

    print("Starting local optimization")
    if local_opt_method == 'BGFS':
        res = minimize(ofunc, x0, args, method='BFGS', jac=objective_func_grad,
            options=local_opt_options)
    else:
        # initial_simplex = np.array([[1E-3,1E-3,1E-5,1E-5,1],
        #                             [1E-5,1E-3,1E-5,1E-5,1],
        #                             [1E-3,1E-5,1E-5,1E-5,1],
        #                             [1E-3,1E-3,1E-1,1E-5,1],
        #                             [1E-3,1E-3,1E-5,1E-1,1],
        #                             [1E-3,1E-3,1E-5,1E-5,2.5]])
        # for j in range(initial_simplex.shape[0]):
        #     s = initial_simplex[j]
        #     for i in (0,1):
        #         s[i] = logit(s[i])
        #     for i in (2,3):
        #         s[i] = np.log(s[i])
        #     s[-1] = np.log(s[-1])
        #     initial_simplex[j] = s
        # local_opt_options["initial_simplex"] = initial_simplex
        # print(f"initial simplex = {initial_simplex}")

        res = minimize(ofunc, x0, args, method=local_opt_method,
            options=local_opt_options)

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

    data2save = {"success":res.success, "x":res.x, "x0":x0, "p_opt":p_opt,
            "sb2_opt":sb2_opt, "s02_opt":s02_opt, "fun":res.fun, "nit":res.nit,
            "nfev":res.nfev, "annot_categories":annot_categories}

    if local_opt_method == 'Nelder-Mead':
        data2save.update({"final_simplex_x":res.final_simplex[0],
            "final_simplex_fun":res.final_simplex[1]})
    elif local_opt_method == 'BGFS':
        data2save.update({"hess_inv":res.hess_inv, "njev":res.njev})

    if global_opt:
        data2save.update({"global_fun":res_global.fun, "global_nfev":res_global.nfev,
            "global_nit":res_global.nit, "global_success":res_global.success,
            "global_x":res_global.x})

    np.savez(opt_result_file, **data2save)

    for k, v in data2save.items():
        print(f"{k}: {v}")

    print(f"Optimization result saved to {opt_result_file}")



if __name__ == "__main__":
    print(f"optimize.py started at {datetime.datetime.now()}")
    print(f"Reading config from {sys.argv[1]}")
    cfg = configparser.ConfigParser(interpolation=configparser.ExtendedInterpolation())
    cfg.read(sys.argv[1])

    os.environ["OMP_NUM_THREADS"] = cfg["omp"].get("OMP_NUM_THREADS")
    print(f"OMP_NUM_THREADS is set to {os.environ['OMP_NUM_THREADS']}")

    # input parameters
    template_dir = cfg["general"].get("template_dir")
    sumstats_kwargs = dict(cfg["sumstats"])
    annot_f = cfg["general"].get("annot_f")
    opt_result_file = cfg["general"].get("opt_result_file")

    template_snp, template_z, template_snp_in_sumstats, template_annot, s2, is2, annot_s2, annot_categories = process_input(template_dir, annot_f, sumstats_kwargs)

    run_opt = cfg["optimization"].getboolean("run_opt")
    if run_opt:
        print("Preparing optimization")
        local_opt_method = cfg["optimization"].get("local_opt_method")
        assert local_opt_method in ('Nelder-Mead', 'Powell', 'BGFS'), f"{local_opt_method} method is not supported"
        local_opt_options = dict(cfg["local_opt"])
        global_opt = cfg["optimization"].getboolean("global_opt")
        global_opt_options = dict(cfg["global_opt"])
        abs_z_threshold = cfg["optimization"].getfloat("abs_z_threshold", fallback=None)
        adaptive = cfg["optimization"].getboolean("adaptive")
        maxiter = cfg["optimization"].getint("maxiter")
        same_pi = cfg["optimization"].getboolean("same_pi")
        same_sb2 = cfg["optimization"].getboolean("same_sb2")
        init_s02_with_qq_slope = cfg["optimization"].getboolean("init_s02_with_qq_slope")
        snp2use_f = cfg["optimization"].get("snp2use_f", fallback=None)
        subset_size = cfg["optimization"].getint("subset_size", fallback=None)
        subset_seed = cfg["optimization"].getint("subset_seed", fallback=None)
        p_0 = cfg["optimization"].get("p_0", fallback=None)
        if not (p_0 is None):
            p_0 = np.array(list(map(float, p_0.split())))
        sb2_0 = cfg["optimization"].get("sb2_0", fallback=None)
        if not (sb2_0 is None):
            sb2_0 = np.array(list(map(float, sb2_0.split())))
        s02_0 = cfg["optimization"].getfloat("s02_0", fallback=None)
        rand_init_seed = cfg["optimization"].getint("rand_init_seed", fallback=None)

        # convert options from strings to corresponding values
        for key, val in global_opt_options.items():
            try: global_opt_options[key] = eval(val) 
            except NameError: pass
        for key, val in local_opt_options.items():
            try: local_opt_options[key] = eval(val) 
            except NameError: pass

        z2use = template_snp_in_sumstats.copy()
        if not abs_z_threshold is None:
            large_z_i = np.abs(template_z) > abs_z_threshold
            print(f"{large_z_i.sum()} variants have |z| > {abs_z_threshold:.4f}")
            z2use[large_z_i] = 0
        if snp2use_f:
            z2use = overlap(z2use, template_snp, snp2use_f)
        if subset_size:
            z2use = randsubset(z2use, subset_size, subset_seed)

        n_categories = len(annot_categories)
        p, sb2, s02 = rand_initial_guess(n_categories, same_pi, same_sb2, rand_init_seed)
        if not (p_0 is None):
            p = p_0
        if not (sb2_0 is None):
            sb2 = sb2_0

        if not (s02_0 is None):
            s02 = s02_0
        elif init_s02_with_qq_slope:
            s02 = estimate_qq_slope(template_z[z2use.astype("bool")])
            print(f"s02 initialized as {s02}")

        print(f"p_0: {p}")
        print(f"sb2_0: {sb2}")
        print(f"s02_0: {s02}")

        run_optimization(p, sb2, s02, template_z, z2use, s2, is2, annot_s2, annot_categories,
                         local_opt_method, local_opt_options, global_opt, global_opt_options,
                         same_pi, same_sb2, opt_result_file)

    make_qq = cfg["qq"].getboolean("make_qq")
    if make_qq:
        print("Making QQ plot")

	# make qq for pruned set of SNPs
        if cfg["qq"].getboolean("use_pruned_snps"):
            snp2use_f = cfg["optimization"].get("snp2use_f", fallback=None)
            if snp2use_f:
                template_snp_in_sumstats = overlap(template_snp_in_sumstats, template_snp, snp2use_f)


        modelqq_out_file = cfg["qq"].get("modelqq_out_file")
        template_qq_annot_file = cfg["qq"].get("template_annot_file")
        n_samples = cfg["qq"].getint("n_samples")
        opt_result_file = cfg["qq"].get("opt_result_file")
        # get template annotations
        qq_template_annot, annot_names = get_qq_annot4snp(template_snp, template_qq_annot_file)
        # parameters derived from data
        p, sb2, s02 = get_params(opt_result_file)
        p_experimental = 2*sstats.norm.cdf(-np.abs(template_z))
        y_max = min(50, max(-np.log10(p_experimental)))

        # get model data
        n_grid = 150
        p_grid = np.logspace(-y_max,0,n_grid)
        # multiply by 0.5 since we want to have two tailed quantiles    
        z_grid = sstats.norm.ppf(0.5*p_grid)

        # qq plot should be only for the SNPs which present in sumstats data, i.e. template_snp_in_sumstats
        # TODO: allow taking only a subset of SNPs for the qq plot. This subset must be used both for modeled and experemental plot.
        # z_cdf_total, z_cdf_annot = get_z_cdf_2tails(z_grid, template_snp_in_sumstats, n_samples,
        #     p, sb2, s02, s2, is2, annot_s2, qq_template_annot)

        z_cdf_total, z_cdf_annot = cmmcost_omp.get_cdfsampling(z_grid, template_snp_in_sumstats, s2, is2, p,
                sb2, s02, annot_s2, qq_template_annot.astype('u1'), n_samples)

        model_total_x = -np.log10(z_cdf_total)
        model_annot_x = -np.log10(z_cdf_annot)
        model_y =  -np.log10(p_grid)

        # get experimental data
        data_total_x, data_total_y = get_xy_from_p(p_experimental[template_snp_in_sumstats.astype('bool')])
        # take only SNPs which are in sumstats
        data_total_x = data_total_x
        data_total_y = data_total_y
        data_annot_x = []
        data_annot_y = []
        for i in range(len(annot_names)):
            annot_i = qq_template_annot[:,i].astype('bool')
            # take only SNPs which are in sumstats
            annot_i &= template_snp_in_sumstats.astype('bool')
            p_experimental_annot = p_experimental[annot_i]
            annot_x, annot_y = get_xy_from_p(p_experimental_annot)
            data_annot_x.append(annot_x)
            data_annot_y.append(annot_y)
        data_annot_x = np.array(data_annot_x)
        data_annot_y = np.array(data_annot_y)        

        # save results
        np.savez(modelqq_out_file, annot_names=annot_names, data_total_x=data_total_x,
            data_total_y=data_total_y, data_annot_x=data_annot_x, data_annot_y=data_annot_y,
            model_total_x=model_total_x, model_annot_x=model_annot_x, model_y=model_y)

        print(f"QQ plot data saved to {modelqq_out_file}")


    run_single = cfg["test"].getboolean("run_single")
    if run_single:
        print("Running 'run_single'")
        # estimate for a given set of parameters
        n_categories = len(annot_categories)
        n_samples = 10000
        z2use = np.ones(len(template_z), dtype='u1')
        if sum(z2use) > 200000:
            z2use = randsubset(z2use, 300000, 1)
        p = np.array([0.005]*n_categories, dtype='f8') # [0.02, 0.005, 0.05, 0.001] [0.005]*n_categories
        sb2 = np.array([0.0002]*n_categories, dtype='f8') # [5.786066041591348e-05, 5.837030778549846e-05, 5.703057621611015e-05, 5.6426373881400386e-05] [0.0002]*n_categories
        s02 = 1.0 # 1.0

        print("Estimating direct cost")
        start = time.time()
        cost = cmmcost_omp.get_cost(template_z, z2use, s2, is2, p, sb2, s02, annot_s2, n_categories)
        end = time.time()
        print(f"   {cost:.8f} = cost direct is")
        print(f"   {end-start:.8f} seconds per single cost function evaluation")

        print("Estimating sampling cost")
        start = time.time()
        cost = cmmcost_omp.get_costsampling(template_z, z2use, s2, is2, p, sb2, s02, annot_s2, n_samples)
        end = time.time()
        print(f"   {cost:.8f} = cost sampling is")
        print(f"   {end-start:.8f} seconds per single cost function evaluation")
        

        print("Estimating gradient of direct cost")
        start = time.time()
        cost_grad = cmmcost_omp.get_cost_der(template_z, z2use, s2, is2, p, sb2, s02, annot_s2, n_categories)
        end = time.time()
        print(f"   [{', '.join('%.4f' % e for e in cost_grad)}] = gradient")
        print(f"   {end-start:.8f} seconds per single gradient evaluation")
        

    print(f"optimize.py finished at {datetime.datetime.now()}")
