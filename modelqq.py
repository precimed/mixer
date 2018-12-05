import numpy as np
import datetime
import scipy.stats as sstats


def get_z_cdf_2tails(z_grid, snps2use, n_samples, p, sb2, s02, s2, is2, annot_s2, template_annot,
    n_categories):
    print("Getting z-score CDF")

    z_cdf_total = np.zeros(len(z_grid))
    z_cdf_annot = np.zeros((n_categories, len(z_grid)))

    # uniform weights
    annot_i, annot_c = np.unique(template_annot, return_counts=True)
    total_c = annot_c.sum()
    annot_c = annot_c.reshape((len(annot_c),1))

    max_block_size = max([(is2[i+1]-is2[i]) for i in range(len(template_annot))])
    print(f"Largest LD block size: {max_block_size}")
    rg = np.random.rand(n_samples, max_block_size)
    for i_templ, (istart, iend) in enumerate(zip(is2[:-1], is2[1:])):
        if i_templ%10000 == 0: print(f"{i_templ} variants processed")
        if snps2use[i_templ]:
            p_in_ld = p[annot_s2[istart:iend]]
            sb2_in_ld = sb2[annot_s2[istart:iend]]*s2[istart:iend]
            curr_ld_block_size = iend-istart
            sigma2 = (sb2_in_ld*(rg[:,:curr_ld_block_size] < p_in_ld)).sum(axis=1) + s02
            if len(sigma2) > 0:
                sigma2 = sigma2.reshape((n_samples, 1))
                z_cdf = sstats.norm.cdf(z_grid, 0, np.sqrt(sigma2)).sum(axis=0)/len(sigma2) # estimate z_cdf as a cdf of equally weighted mixture of n_samples elements with sigma2 variance
                # z_cdf_total += z_cdf*weight_total[i_templ]
                # z_cdf_annot[template_annot[i_templ]] += z_cdf*weight_annot[i_templ]
                z_cdf_total += z_cdf
                z_cdf_annot[template_annot[i_templ]] += z_cdf
    # if weighted, z_cdf_total and z_cdf_annot should be normalized by the sum of all weights in the corresponding annotation, i.e. z_cdf_total /= sum(weights_all), z_cdf_annot[i] /= sum(weights_annot[i])
    # see cmma.jl code for more details
    # multiply by 2 since we want to have 2 tails
    return 2*z_cdf_total/total_c, 2*z_cdf_annot/annot_c
        

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


def load_idump(dump_f):
    print(f"Loading dumped imnput from {dump_f}")
    idump = np.load(dump_f)
    z = idump.get("z")
    z2use = idump.get("z2use")
    template_annot = idump.get("template_annot")
    s2 = idump.get("s2")
    is2 = idump.get("is2")
    annot_s2 = idump.get("annot_s2")
    n_categories = len(np.unique(template_annot))
    return z, z2use, template_annot, s2, is2, annot_s2, n_categories


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



if __name__ == '__main__':
    print(f"modelqq.py started at {datetime.datetime.now()}")

    dump_input_file = "/mnt/seagate10/projects/cmm/experiments/hapgen10k11m_chrs_21_22.chrs_21_22.fixed_or_free_pi_and_sigma/idump.hapgen10k11m_chrs_21_22.chrs_21_22.sumstats.hapgen10k11m_chrs_21_22.chrs_21_22.chr_21_p001.chr_22_p0001.npz"
    opt_result_file = "/mnt/seagate10/projects/cmm/experiments/hapgen10k11m_chrs_21_22.chrs_21_22.fixed_or_free_pi_and_sigma/results/optimize.hapgen10k11m_chrs_21_22.chrs_21_22.sumstats.hapgen10k11m_chrs_21_22.chrs_21_22.chr_21_p001.chr_22_p0001.test.run_1.npz"
    modelqq_out_file = "/mnt/seagate10/projects/cmm/experiments/hapgen10k11m_chrs_21_22.chrs_21_22.fixed_or_free_pi_and_sigma/results/qqtest.hapgen10k11m_chrs_21_22.chrs_21_22.sumstats.hapgen10k11m_chrs_21_22.chrs_21_22.chr_21_p001.chr_22_p0001.same_pi.run_1.test.npz"
    # parameters for model qq plot
    p, sb2, s02 = get_params(opt_result_file)

    z, z2use, template_annot, s2, is2, annot_s2, n_categories = load_idump(dump_input_file)
    n_samples = 100

    # parameters derived from data
    p_experimental = 2*sstats.norm.cdf(-np.abs(z))
    y_max = min(50, max(-np.log10(p_experimental)))

    # get model data
    n_grid = 150
    p_grid = np.logspace(-y_max,0,n_grid)
    # multiply by 0.5 since we want to have two tailed quantiles    
    z_grid = sstats.norm.ppf(0.5*p_grid)
    z_cdf_total, z_cdf_annot = get_z_cdf_2tails(z_grid, z2use, n_samples, p, sb2, s02, s2, is2, annot_s2, template_annot, n_categories)
    model_total_x = -np.log10(z_cdf_total)
    model_annot_x = -np.log10(z_cdf_annot)
    model_y =  -np.log10(p_grid)

    # get experimental data
    data_total_x, data_total_y = get_xy_from_p(p_experimental)
    annot_categories = np.unique(template_annot)
    data_annot_x = []
    data_annot_y = []
    for a in annot_categories:
        i = (template_annot == a)
        p_experimental_annot = p_experimental[i]
        annot_x, annot_y = get_xy_from_p(p_experimental_annot)
        data_annot_x.append(annot_x)
        data_annot_y.append(annot_y)
    data_annot_x = np.array(data_annot_x)
    data_annot_y = np.array(data_annot_y)


    # save results
    np.savez(modelqq_out_file, data_total_x=data_total_x, data_total_y=data_total_y,
        data_annot_x=data_annot_x, data_annot_y=data_annot_y,
        model_total_x=model_total_x, model_annot_x=model_annot_x, model_y=model_y)


    print(f"Output saved to {modelqq_out_file}")
