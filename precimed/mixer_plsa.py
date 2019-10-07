#!/usr/bin/env python
'''
(c) 2018-2019 Oleksandr Frei, Alexey A. Shadrin
MiXeR software: Univariate and Bivariate Causal Mixture for GWAS
NB! PLSA does not refer to Probabilistic latent semantic analysis, rather
    P <- causal mixture model (pi) 
    L <- LD-dependent architectures
    S <- MAF-dependent architectures
    A <- annotation-informed
'''

import logging
import sys

from mixer.plsa import parse_args
from mixer.plsa import log_header
from mixer.libbgmg import LibBgmg

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    args = parse_args(sys.argv[1:])

    libbgmg = LibBgmg(args.lib)
    libbgmg.init_log(args.log if args.log else args.out + '.log')
    log_header(args, sys.argv[1], libbgmg)
    libbgmg.dispose()

    args.func(args)
