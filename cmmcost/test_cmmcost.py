import scipy.stats as ss
import numpy as np
import cmmcost
import argparse
import sys
import time


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--no-precision", action="store_true",
            help="Don't measure precision")
    parser.add_argument("--no-speed", action="store_true",
            help="Don't measure speed")
    args = parser.parse_args(sys.argv[1:])


    z2use = np.array([True, True, False], dtype='u1')
    s2 = np.array([1.593, 0.934, 2.463, 0.719, 1.847, 3.012, 1.927, 0.896, 1.401], dtype='f4')
    z = np.array([1.275, -0.496, 2.983])
    is2 = np.array([0,3,5,9], dtype='u8')
    p = np.array([1.,1.,1.])
    sb2 = np.array([1.17, 2.03, 0.954])
    s02 = 1.03
    annot=np.array([0,1,0,2,1,1,0,0,2], dtype='u1')

    if not args.no_precision:
        print("Running precision test")
        # get estimated cost
        estimated_cost = cmmcost.get_cost(z, z2use, s2, is2, p, sb2, s02, annot)

        # get exact cost
        sigma1 = s2[0]*sb2[annot[0]] + s2[1]*sb2[annot[1]] + s2[2]*sb2[annot[2]] + s02
        sigma2 = s2[3]*sb2[annot[3]] + s2[4]*sb2[annot[4]] + s02
        exact_cost = 0
        exact_cost += -np.log(ss.norm.pdf(z[0], 0, np.sqrt(sigma1)))
        exact_cost += -np.log(ss.norm.pdf(z[1], 0, np.sqrt(sigma2)))
        exact_cost *= 0.5

        thresh = 1.E-12
        deviation = abs(estimated_cost - exact_cost)

        print(f"   estimated cost: {estimated_cost}")
        print(f"   exact cost:     {exact_cost}")
        print(f"   deviation:      {deviation}")

        if deviation < thresh:
            print("   precision test PASSED")
        else:
            print(f"   deviation is larger than {thresh}")
            print("   precision test FAILED")

    if not args.no_speed:
        print("Running speed test")
        niter = 100000
        start = time.process_time()
        [cmmcost.get_cost(z, z2use, s2, is2, p, sb2, s02, annot) for _ in range(niter)]
        end = time.process_time()
        print(f"   time per {niter} iterations (seconds): {end-start}")

