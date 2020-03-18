#!/usr/bin/env python
'''
(c) 2016-2020 Oleksandr Frei, Alexey A. Shadrin, Dominic Holland
MiXeR software: Univariate and Bivariate Causal Mixture for GWAS
'''

import logging
import sys

from mixer.cli import parse_args
from mixer.cli import log_header
from mixer.libbgmg import LibBgmg

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    args = parse_args(sys.argv[1:])

    libbgmg = LibBgmg(args.lib)
    libbgmg.init_log(args.log if args.log else args.out + '.log')
    log_header(args, sys.argv[1], libbgmg)
    libbgmg.dispose()

    args.func(args)
