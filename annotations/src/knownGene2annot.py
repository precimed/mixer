import gzip
import sys
import argparse

SHOW_WARNINGS = False

def parseArgs(args):
    parser = argparse.ArgumentParser(
        description="Create annotation file from UCSC knownGene file.")
    parser.add_argument("known_gene_file", help="UCSC knownGene file")
    parser.add_argument("out_file", help="Output file name")
    parser.add_argument("--show-warns", action="store_true",
        help="Suppress warnings")
    return parser.parse_args(args)

def myOpen(f_name, mode='r'):
    if f_name.endswith(".gz"):
        return gzip.open(f_name, mode)
    else:
        return open(f_name, mode)


def write2file(columnList, f):
    line = '\t'.join(columnList)
    if columnList[-2] == columnList[-1]:
        msg = "Warning generating annotation:\n%s\n" % line
        msg += "Empty segment generated. It will be ignored."
        if SHOW_WARNINGS: print(msg)
    else:
        f.write("%s\n" % line)


def parseLine(line, outFile):
    """
    Parse a line from UCSC knownGene file. Write annotated segments to outFile.
    UCSC knownGene file for hg19 can be downloaded from here:
        http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/
    knownGene format description is here:
        http://genome.ucsc.edu/cgi-bin/hgTables
    Supported categories:
        5UTR, 3UTR, Exon, ProteinConding, Intron, 1kUp, 10kUp, 100kUp, 1kDown,
        10kDown, 100kDown
        for coding transcripts and a single category NoncodingTranscript for
        noncoding transcripts (to be consistent wit classification used in 
        "All SNPs are not vreated equal" paper by A. Schork).
    Output file contains 7 tab-delimited columns:
        gene name, chromosome, strand, protein ID, category, start, end
    Algorithm (assuming the gene is on "+" strand):
                                 Gene
                start                             end
     -------------[--Exon--|----Intron----|--Exon--]-------------
           [-1kUp-]                                [-1kDown-]
      [---10kUp---]                                [---10kUp---]
    All Exons are taken as they are reported in the knownGene file, all segments
    between consequent exons are taken as Introns.
    If cdsStart == cdsEnd: no coding sequence, no further processing
    else:
        each exon is compared to cds start/end as follows (6 different cases):
                s                  e
     -----------|-------Exon-------|-----------
        cs   ce
     1: [-cds-]                                  (s,e)   - 3UTR
     2: [----cds----]                            (s,ce)  - ProtCod
                                                 (ce,e)  - 3UTR
     3: [----------------cds----------------]    (s,e)   - ProtCod
     4:             [---cds---]                  (s,cs)  - 5UTR
                                                 (cs,ce) - ProtCod
                                                 (ce,e)  - 3UTR
     5:             [----------cds----------]    (s,cs)  - 5UTR
                                                 (cs,e)  - ProtCod
     6:                               [-cds-]    (s,e)   - 5UTR

     additional special sub-cases are when cds limits are equal to exon limits.
    """
    spltLine = line.split("\t")
    (name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount,
        exonStarts, exonEnds, proteinID, alignID) = spltLine
    txStart, txEnd, cdsStart, cdsEnd, exonCount = map(int, spltLine[3:8])
    exonStarts = exonStarts.split(",")[:-1]
    exonEnds = exonEnds.split(",")[:-1]
    columnList = [name, chrom, strand, proteinID]
    if cdsStart == cdsEnd:
        # noncoding transcript
        s,e = str(txStart), str(txEnd)
        write2file(columnList + ["NoncodingTranscript", s, e], outFile)
    else:
        for s,e in zip(exonEnds[:-1], exonStarts[1:]):
            write2file(columnList + ["Intron", s, e], outFile)
        for s,e in zip(exonStarts, exonEnds):
            write2file(columnList + ["Exon", s, e], outFile)
        exonStarts = map(int,exonStarts)
        exonEnds = map(int, exonEnds)
        if strand == "+":
            s, e = max(0, txStart-1000), txStart
            write2file(columnList + ["1kUp", "%d" % s, "%d" % e], outFile)
            s, e = max(0, txStart-10000), txStart
            write2file(columnList + ["10kUp", "%d" % s, "%d" % e], outFile)
            s, e = max(0, txStart-100000), txStart
            write2file(columnList + ["100kUp", "%d" % s, "%d" % e], outFile)
            s, e = txEnd, txEnd+1000
            write2file(columnList + ["1kDown", "%d" % s, "%d" % e], outFile)
            s, e = txEnd, txEnd+10000
            write2file(columnList + ["10kDown", "%d" % s, "%d" % e], outFile)
            s, e = txEnd, txEnd+100000
            write2file(columnList + ["100kDown", "%d" % s, "%d" % e], outFile)
            if cdsStart < cdsEnd: # otherwise this gene is non-coding
                # consider 6 different cases described above
                for s,e in zip(exonStarts, exonEnds):
                    if cdsEnd <= s: # case 1
                        write2file(columnList + ["3UTR", "%d" % s, "%d" % e], outFile)
                    elif e <= cdsStart: # case 6
                        write2file(columnList + ["5UTR", "%d" % s, "%d" % e], outFile)
                    elif cdsStart < s and cdsEnd < e: # case 2
                        write2file(columnList + ["ProteinCoding", "%d" % s, "%d" % cdsEnd], outFile)
                        write2file(columnList + ["3UTR", "%d" % cdsEnd, "%d" % e], outFile)
                    elif cdsStart <= s and cdsEnd >= e: # case 3
                        write2file(columnList + ["ProteinCoding", "%d" % s, "%d" % e], outFile)
                    elif s <= cdsStart and cdsEnd <= e: # case 4
                        write2file(columnList + ["ProteinCoding", "%d" % cdsStart, "%d" % cdsEnd], outFile)
                        if s < cdsStart:
                            write2file(columnList + ["5UTR", "%d" % s, "%d" % cdsStart], outFile)
                        if cdsEnd < e:
                            write2file(columnList + ["3UTR", "%d" % cdsEnd, "%d" % e], outFile)
                    elif s < cdsStart and cdsEnd > e: # case 5
                        write2file(columnList + ["ProteinCoding", "%d" % cdsStart, "%d" % e], outFile)
                        write2file(columnList + ["5UTR", "%d" % s, "%d" % cdsStart], outFile)
                    else:
                        msg = "Error processing line:\n%s\n" % line
                        msg += "unclassified exon: (start, end) = (%d,%d)" % (s, e)
                        raise ValueError(msg)
        elif strand == "-":
            s, e = max(0, txStart-1000), txStart
            write2file(columnList + ["1kDown", "%d" % s, "%d" % e], outFile)
            s, e = max(0, txStart-10000), txStart
            write2file(columnList + ["10kDown", "%d" % s, "%d" % e], outFile)
            s, e = max(0, txStart-100000), txStart
            write2file(columnList + ["100kDown", "%d" % s, "%d" % e], outFile)
            s, e = txEnd, txEnd+1000
            write2file(columnList + ["1kUp", "%d" % s, "%d" % e], outFile)
            s, e = txEnd, txEnd+10000
            write2file(columnList + ["10kUp", "%d" % s, "%d" % e], outFile)
            s, e = txEnd, txEnd+100000
            write2file(columnList + ["100kUp", "%d" % s, "%d" % e], outFile)
            if cdsStart < cdsEnd: # otherwise this gene is non-coding
                # consider 6 different cases described above, but for "-" strand
                for s,e in zip(exonStarts, exonEnds):
                    if cdsEnd <= s: # case 1
                        write2file(columnList + ["5UTR", "%d" % s, "%d" % e], outFile)
                    elif e <= cdsStart: # case 6
                        write2file(columnList + ["3UTR", "%d" % s, "%d" % e], outFile)
                    elif cdsStart < s and cdsEnd < e: # case 2
                        write2file(columnList + ["ProteinCoding", "%d" % s, "%d" % cdsEnd], outFile)
                        write2file(columnList + ["5UTR", "%d" % cdsEnd, "%d" % e], outFile)
                    elif cdsStart <= s and cdsEnd >= e: # case 3
                        write2file(columnList + ["ProteinCoding", "%d" % s, "%d" % e], outFile)
                    elif s <= cdsStart and cdsEnd <= e: # case 4
                        write2file(columnList + ["ProteinCoding", "%d" % cdsStart, "%d" % cdsEnd], outFile)
                        if s < cdsStart:
                            write2file(columnList + ["3UTR", "%d" % s, "%d" % cdsStart], outFile)
                        if cdsEnd < e:
                            write2file(columnList + ["5UTR", "%d" % cdsEnd, "%d" % e], outFile)
                    elif s < cdsStart and cdsEnd > e: # case 5
                        write2file(columnList + ["ProteinCoding", "%d" % cdsStart, "%d" % e], outFile)
                        write2file(columnList + ["3UTR", "%d" % s, "%d" % cdsStart], outFile)
                    else:
                        msg = "Error processing line:\n%s\n" % line
                        msg += "unclassified exon: (start, end) = (%d,%d)" % (s, e)
                        raise ValueError(msg)
        else:
            msg = "Error processing line:\n%s\n" % line
            msg += "unknown strand: '%s'" % strand
            raise ValueError(msg)



if __name__ == "__main__":
    args = parseArgs(sys.argv[1:])
    SHOW_WARNINGS = args.show_warns
    with myOpen(args.known_gene_file, 'rt') as f, myOpen(args.out_file, 'wt') as of:
        print("Processing UCSC knownGene file %s" % args.known_gene_file)
        print("Writing to %s" % args.out_file)
        for i, l in enumerate(f):
            parseLine(l, of)
            if (i+1)%10000 == 0:
                print("%d lines processed" % (i+1))
        print("%d lines processed in total" % (i+1))
    print("Completed")
