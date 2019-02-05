import scipy.stats as ss
import numpy as np
import cmmcost_omp as cmmcost
import argparse
import sys
import time


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--no-precision", action="store_true",
            help="Don't measure precision")
    parser.add_argument("--no-speed", action="store_true",
            help="Don't measure speed")
    parser.add_argument("--no-cdf", action="store_true",
            help="Don't measure cdf precision")
    args = parser.parse_args(sys.argv[1:])


    z2use = np.array([True, True, False], dtype='u1')
    s2 = np.array([1.593, 0.934, 2.463, 0.719, 1.847, 3.012, 1.927, 0.896, 1.401], dtype='f4')
    z = np.array([1.275, -0.496, 2.983])
    is2 = np.array([0,3,5,9], dtype='u8')
    p = np.array([1.,1.,1.])
    sb2 = np.array([1.17, 2.03, 0.954])
    s02 = 1.03
    annot=np.array([0,1,0,2,1,1,0,0,2], dtype='u1')

    n_samples = 100000

    # additional parameters for cdfsampling function
    n_zgrid = 4
    n_qq_annot = 3
    zgrid = np.array([0.0, -0.1, -0.01, -0.001]);
    # qq_template_annot : [nz x n_qq_annot] 2D array
    qq_template_annot = np.array([[True, True, True],
                                  [False, True, True],
                                  [False, False, True]], dtype='u1')


    # SNP   S2_eff
    # 0     1.593*1.17 + 0.934*2.03 + 2.463*1.17 + 1.03 = 7.671539999999999
    # 1     0.719*0.954 + 1.847*2.03 + 1.03 = 5.465336
    # 2     3.012*2.03 + 1.927*1.17 + 0.896*1.17 + 1.401*0.954 + 1.03 = 11.783824

    # z2use = {true, true, false}
    # qq_template_annot = { {true,  true,  true}, {false, true,  true}, {false,  false, true} }
    # zgrid = {0.0, -0.1, -0.01, -0.001}
    # zcdf = {1.0, 0.9685399893037775, 0.9968531741508169, 0.9996853165901003}
    # zcdf_qq_annot = { 1.0                 1.0                 1.0
    #                   0.971199207323047   0.9685399893037775  0.9685399893037775
    #                   0.9971193012706704  0.9968531741508169  0.9968531741508169
    #                   0.9997119295074843  0.9996853165901003  0.9996853165901003 }


    if not args.no_precision:
        print("Running precision test")
        # get exact cost
        sigma1 = s2[0]*sb2[annot[0]] + s2[1]*sb2[annot[1]] + s2[2]*sb2[annot[2]] + s02
        sigma2 = s2[3]*sb2[annot[3]] + s2[4]*sb2[annot[4]] + s02
        exact_cost = 0
        exact_cost += -np.log(ss.norm.pdf(z[0], 0, np.sqrt(sigma1)))
        exact_cost += -np.log(ss.norm.pdf(z[1], 0, np.sqrt(sigma2)))
        exact_cost *= 0.5

        thresh = 1.E-12

        print(f"   exact cost:     {exact_cost}")
        print("")

        # get estimated direct cost
        estimated_cost = cmmcost.get_cost(z, z2use, s2, is2, p, sb2, s02, annot)
        deviation = abs(estimated_cost - exact_cost)
        print(f"   estimated direct cost: {estimated_cost}")
        print(f"   deviation:             {deviation}")

        if deviation < thresh:
            print("   precision test PASSED")
        else:
            print(f"   deviation is larger than {thresh}")
            print("   precision test FAILED")
        print("")

        # get estimated sampling cost
        estimated_cost = cmmcost.get_costsampling(z, z2use, s2, is2, p, sb2, s02, annot, n_samples)
        deviation = abs(estimated_cost - exact_cost)
        print(f"   estimated sampling cost: {estimated_cost}")
        print(f"   deviation:               {deviation}")

        if deviation < thresh:
            print("   precision test PASSED")
        else:
            print(f"   deviation is larger than {thresh}")
            print("   precision test FAILED")
        print("")

    if not args.no_speed:
        print("Running speed test")
        niter = 100000
        start = time.time()
        [cmmcost.get_cost(z, z2use, s2, is2, p, sb2, s02, annot) for _ in range(niter)]
        end = time.time()
        print("   direct cost:")
        print(f"      time per {niter} iterations (seconds): {end-start}")

        start = time.time()
        [cmmcost.get_costsampling(z, z2use, s2, is2, p, sb2, s02, annot, n_samples)]
        end = time.time()
        print("   sampling cost:")
        print(f"      time per {niter} iterations (seconds): {end-start}")
        print("")

    if not args.no_cdf:
        print("Running CDF test")
        # get exact cost
        sigma1 = s2[0]*sb2[annot[0]] + s2[1]*sb2[annot[1]] + s2[2]*sb2[annot[2]] + s02
        sigma2 = s2[3]*sb2[annot[3]] + s2[4]*sb2[annot[4]] + s02
        sigmas = np.sqrt([sigma1, sigma2])
        zcdf_exact = (2/len(sigmas)) * np.array([ss.norm.cdf(z, 0, sigmas).sum() for z in zgrid])

        annot_sigmas = [np.sqrt([sigma1]), sigmas, sigmas]
        zcdf_qq_annot_exact = [(2/len(annot_s))*ss.norm.cdf(z, 0, annot_s).sum() for annot_s in annot_sigmas for z in zgrid]
        zcdf_qq_annot_exact = np.array(zcdf_qq_annot_exact).reshape((4,3), order='F')

        thresh = 1.E-10

        zcdf, zcdf_qq_annot = cmmcost.get_cdfsampling(zgrid, z2use, s2, is2, p,
                sb2, s02, annot, qq_template_annot, n_samples)
        zcdf_deviation = np.max(np.abs(zcdf - zcdf_exact))
        zcdf_qq_annot_deviation = np.max(np.abs(zcdf_qq_annot - zcdf_qq_annot_exact))

        if zcdf_deviation < thresh:
            print("   zcdf precision test PASSED")
        else:
            print(f"   deviation is larger than {thresh}")
            print("   zcdf precision test FAILED")

        if zcdf_qq_annot_deviation < thresh:
            print("   zcdf_qq_annot precision test PASSED")
        else:
            print(f"   deviation is larger than {thresh}")
            print("   zcdf_qq_annot precision test FAILED")



