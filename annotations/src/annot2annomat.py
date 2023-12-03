import sys
import pandas as pd
from collections import defaultdict
import numpy as np
import argparse


def parseArgs(args):
    parser = argparse.ArgumentParser(
        description="Creates binary annotation matrix from the output of knownGene2annot.py")
    parser.add_argument("annot_file", help="A file with snp id in column 4 and annotation category name in column 8")
    parser.add_argument("template_file", help="Template bed file")
    parser.add_argument("out_file", help="Output file name")
    return parser.parse_args(args)


if __name__ == "__main__":
    args = parseArgs(sys.argv[1:])
    print("Processing %s" % args.annot_file)
    df = pd.read_csv(args.annot_file, header=None, sep='\t')
    dd = defaultdict(set)

    for t in df.itertuples():
        dd[t[4]].add(t[8]) # t[0] = index

    template_snps = pd.read_csv(args.template_file,sep='\t',usecols=[3],
                                header=None,names=["SNP"],squeeze=True)
    # hardcoded list of categories supported by knownGene2annot.py
    l = ["100kDown", "10kDown", "1kDown", "100kUp", "10kUp", "1kUp", "3UTR",
        "5UTR", "Exon", "Intron", "ProteinCoding", "NoncodingTranscript",
        "mirna", "tfbs"]

    data = np.zeros((len(template_snps), len(l)), dtype=int)
    for i,s in enumerate(template_snps):
        for j, k in enumerate(l):
            if k in dd[s]:
                data[i,j] = 1

    dfa = pd.DataFrame(index=template_snps, columns=l, data=data)
    print("Writing results to %s" % args.out_file)
    # out_file = "../data/9m_template_knownGene_ucsc_hg19_annomat.txt.gz"
    dfa.to_csv(args.out_file, sep='\t',
        compression='gzip' if args.out_file.endswith('.gz') else None)
    print("Done")

