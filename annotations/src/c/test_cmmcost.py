import scipy.stats as ss
import numpy as np
from scipy.stats import norm
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
    parser.add_argument("--no-derivative", action="store_true",
            help="Don't measure derivative precision")
    args = parser.parse_args(sys.argv[1:])


    z2use = np.array([True, True, False], dtype='u1')
    s2 = np.array([1.593, 0.934, 2.463, 0.719, 1.847, 3.012, 1.927, 0.896, 1.401], dtype='f4')
    z = np.array([1.275, -0.496, 2.983])
    is2 = np.array([0,3,5,9], dtype='u8')
    p = np.array([1.,1.,1.])
    sb2 = np.array([1.17, 2.03, 0.954])
    s02 = 1.03
    annot=np.array([0,1,0,2,1,1,0,0,2], dtype='u1')
    n_annot = 3

    n_samples = 100000

    # additional parameters for cdfsampling function
    n_zgrid = 4
    n_qq_annot = 3
    zgrid = np.array([0.0, -0.1, -0.01, -0.001]);
    # qq_template_annot : [nz x n_qq_annot] 2D array
    qq_template_annot = np.array([[True, True, True],
                                  [False, True, True],
                                  [False, False, True]], dtype='u1')

    print(f"z2use = [{', '.join('true' if z else 'false' for z in z2use)}]")
    print("")

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
        estimated_cost = cmmcost.get_cost(z, z2use, s2, is2, p, sb2, s02, annot, n_annot)
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
        [cmmcost.get_cost(z, z2use, s2, is2, p, sb2, s02, annot, n_annot) for _ in range(niter)]
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
        print("")

    if not args.no_derivative:
        print("Running cost function derivative test")
        thresh = 1.E-8

        def normal_pdf_sb2_derivative(x, s2):
            return (np.exp(-x*x/(2*s2))/np.sqrt(2*np.pi*s2))*(x*x/(2*s2*s2) - 1/(2*s2))

        z_vec = [1.275, -0.496, 2.983]
        nnz = sum(z2use)
        # SNP   S2_eff (see _cmmcost_omp.c, main())
        # 0     1.593*1.17 + 0.934*2.03 + 2.463*1.17 + 1.03 = 7.671539999999999
        # 1     0.719*0.954 + 1.847*2.03 + 1.03 = 5.465336
        # 2     3.012*2.03 + 1.927*1.17 + 0.896*1.17 + 1.401*0.954 + 1.03 = 11.783824
        #
        # sb2 = {1.17, 2.03, 0.954}
        s2_eff_vec = [7.671539999999999, 5.465336, 11.783824]

        # estimate true (analytical values)
        likelihood_der_sb2_0 = 0
        if z2use[0]: likelihood_der_sb2_0 += -(1.593 + 2.463)*normal_pdf_sb2_derivative(z_vec[0], s2_eff_vec[0])/norm.pdf(z_vec[0],0,np.sqrt(s2_eff_vec[0]))
        if z2use[2]: likelihood_der_sb2_0 += -(1.927 + 0.896)*normal_pdf_sb2_derivative(z_vec[2], s2_eff_vec[2])/norm.pdf(z_vec[2],0,np.sqrt(s2_eff_vec[2]))
        likelihood_der_sb2_0 /= nnz

        likelihood_der_sb2_1 = 0
        if z2use[0]: likelihood_der_sb2_1 += -(0.934)*normal_pdf_sb2_derivative(z_vec[0], s2_eff_vec[0])/norm.pdf(z_vec[0],0,np.sqrt(s2_eff_vec[0]))
        if z2use[1]: likelihood_der_sb2_1 += -(1.847)*normal_pdf_sb2_derivative(z_vec[1], s2_eff_vec[1])/norm.pdf(z_vec[1],0,np.sqrt(s2_eff_vec[1]))
        if z2use[2]: likelihood_der_sb2_1 += -(3.012)*normal_pdf_sb2_derivative(z_vec[2], s2_eff_vec[2])/norm.pdf(z_vec[2],0,np.sqrt(s2_eff_vec[2]))
        likelihood_der_sb2_1 /= nnz

        likelihood_der_sb2_2 = 0
        if z2use[1]: likelihood_der_sb2_2 += -(0.719)*normal_pdf_sb2_derivative(z_vec[1], s2_eff_vec[1])/norm.pdf(z_vec[1],0,np.sqrt(s2_eff_vec[1]))
        if z2use[2]: likelihood_der_sb2_2 += -(1.401)*normal_pdf_sb2_derivative(z_vec[2], s2_eff_vec[2])/norm.pdf(z_vec[2],0,np.sqrt(s2_eff_vec[2]))
        likelihood_der_sb2_2 /= nnz

        likelihood_der_s02 = 0
        if z2use[0]: likelihood_der_s02 += -normal_pdf_sb2_derivative(z_vec[0], s2_eff_vec[0])/norm.pdf(z_vec[0],0,np.sqrt(s2_eff_vec[0]))
        if z2use[1]: likelihood_der_s02 += -normal_pdf_sb2_derivative(z_vec[1], s2_eff_vec[1])/norm.pdf(z_vec[1],0,np.sqrt(s2_eff_vec[1]))
        if z2use[2]: likelihood_der_s02 += -normal_pdf_sb2_derivative(z_vec[2], s2_eff_vec[2])/norm.pdf(z_vec[2],0,np.sqrt(s2_eff_vec[2]))
        likelihood_der_s02 /= nnz

        exact_sb2_s02_gradient = np.array([0, 0, 0, likelihood_der_sb2_0, likelihood_der_sb2_1, likelihood_der_sb2_2, likelihood_der_s02]) # set dummy values vor p derivatives
        gradient = cmmcost.get_cost_der(z, z2use, s2, is2, p, sb2, s02, annot, n_annot)
        print(f"   gradient = [{', '.join('%.6f' % g for g in gradient)}]")
        gradient_sb2_s02_deviation = np.max(np.abs(exact_sb2_s02_gradient[n_annot:] - gradient[n_annot:]))

        if gradient_sb2_s02_deviation < thresh:
            print("   gradient precision test (for sb2 and s02) PASSED")
        else:
            print(f"   deviation is larger than {thresh}")
            print("   gradient precision test FAILED")




