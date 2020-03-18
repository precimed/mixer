#!/usr/bin/env python
'''
(c) 2016-2020 Oleksandr Frei, Alexey A. Shadrin, Dominic Holland
MiXeR software: Univariate and Bivariate Causal Mixture for GWAS
'''

import logging
import sys

from mixer.figures import parse_args

if __name__ == "__main__":
    logging.getLogger().setLevel(logging.INFO)
    args = parse_args(sys.argv[1:])
    args.func(args)
