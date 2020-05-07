import pandas as pd
import numpy as np
import sys
import argparse

"""
Contains a logic for making unique binary annotation categories from overlapping
binary annotations produced by annot2annomat.py.

All annotation categories:
100kDown  10kDown  1kDown  100kUp  10kUp  1kUp  3UTR  5UTR  Exon
Intron  ProteinCoding  NoncodingTranscript  mirna  tfbs

Categories used for annotations:
3UTR  5UTR  Exon  Intron  10kDown  1kDown  10kUp  1kUp

Only variants in protein coding transcripts are considered for annotation:
NoncodingTranscript == 0
Priority order is given by the annot2use list below.

Returns a file where all categories form annot2use are mutually exclusive,
categories from auxiliary_annot are directly copied from the input table.
"""



annot2use = ["5UTR", "3UTR", "Exon", "Intron", "1kUp", "1kDown", "10kUp", "10kDown"]
auxiliary_annot = ["NoncodingTranscript", "100kUp", "100kDown", "mirna", "tfbs"]


def parseArgs(args):
    parser = argparse.ArgumentParser(
        description=("Creates nonoverlapping annotations for 'main' categories "
            "(annot2use). Auxiliary categories (auxiliary_annot) may overlap."))
    parser.add_argument("annomat_file", help="An output of annot2annomat.py")
    parser.add_argument("out_file",
        help="Output file name (will be gzipped if ends with '.gz')")
    return parser.parse_args(args)


if __name__ == "__main__":
    args = parseArgs(sys.argv[1:])
    print("Processing %s" % args.annomat_file)
    # annomat_file = "../data/test_1000genomes_chr22/1kgenomes_chr22_complete_annot_hg19_annomat.txt.gz"
    # out_file = "../data/test_1000genomes_chr22/1kgenomes_chr22_complete_annot_hg19_nonoverlapping_and_aux_annot.txt.gz"
    df = pd.read_table(args.annomat_file, index_col="SNP")

    assert len(set(annot2use) - set(df.columns)) == 0, "Some annot2use categories are absent in the input table"

    df_tmp = df[df["NoncodingTranscript"] == 0].copy()

    ddf = pd.DataFrame(index=df.index, columns=annot2use+auxiliary_annot,
        data=np.zeros((len(df), len(annot2use+auxiliary_annot)), dtype=int))
    for c in annot2use:
        snps_in_c = list(df_tmp[df_tmp[c] == 1].index)
        print("%d variants in %s category" % (len(snps_in_c), c))
        ddf.loc[snps_in_c,c] = 1
        df_tmp.drop(snps_in_c, inplace=True)

    for c in auxiliary_annot:
        ddf[c] = df[c]
    print("Writing output to %s" % args.out_file)
    compression = 'gzip' if args.out_file.endswith(".gz") else  None
    ddf.to_csv(args.out_file, sep='\t', compression=compression)

    print("Done")