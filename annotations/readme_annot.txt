1.1. Download KnownGene table
wget -O /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.txt.gz ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz



1.2. Download MiRNA target regions and regulatory features from Ensembl BioMart.
MiRNA:
wget -O /mnt/seagate10/projects/annotations/data/test/mirna_targets.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_mirna_target_feature" interface = "default" ><Attribute name = "chromosome_name" /><Attribute name = "chromosome_start" /><Attribute name = "chromosome_end" /><Attribute name = "accession" /></Dataset></Query>'

Regulatory features:
wget -O /mnt/seagate10/projects/annotations/data/test/regulatory_features.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_regulatory_feature" interface = "default" ><Attribute name = "chromosome_name" /><Attribute name = "chromosome_start" /><Attribute name = "chromosome_end" /><Attribute name = "regulatory_stable_id" /><Attribute name = "feature_type_name" /></Dataset></Query>'

http://grch37.ensembl.org/biomart/martview/891d411f1bfeb1cd8d2a0f46fe62cff2?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_mirna_target_feature.default.mirna_target_feature.chromosome_name|hsapiens_mirna_target_feature.default.mirna_target_feature.chromosome_start|hsapiens_mirna_target_feature.default.mirna_target_feature.chromosome_end|hsapiens_mirna_target_feature.default.mirna_target_feature.accession&FILTERS=&VISIBLEPANEL=attributepanel



2. Process UCSC annotations with knownGene2annot.py.
python knownGene2annot.py /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.txt.gz /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.annot.txt.gz



3.1. Convert gene, mirna and regulatury feature files into bed format (mark features accordingly in the 4-th column of bed file).
zcat /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.annot.txt.gz | awk -F$'\t' 'BEGIN{OFS="\t"} {if($2 ~ "chr[0-9]+") print($2,$6,$7,$5)}' | gzip -c > /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.annot.bed.gz

awk -F$'\t' 'BEGIN{OFS="\t"} {print("chr"$1,$2-1,$3,"mirna"}' /mnt/seagate10/projects/annotations/data/test/mirna_targets.txt | gzip -c > /mnt/seagate10/projects/annotations/data/test/mirna_targets.bed.gz
awk -F$'\t' 'BEGIN{OFS="\t"} {print("chr"$1,$2-1,$3,"tfbs"}' /mnt/seagate10/projects/annotations/data/test/regulatory_features.txt | gzip -c > /mnt/seagate10/projects/annotations/data/test/regulatory_features.bed.gz



3.2. Merge and sort all bed files for known genes, mirna and tfbs.
zcat /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.annot.bed.gz /mnt/seagate10/projects/annotations/data/test/mirna_targets.bed.gz /mnt/seagate10/projects/annotations/data/test/regulatory_features.bed.gz | sort -k1,1 -k2,2n | gzip -c > /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.annot.complete.sorted.bed.gz



4. Create sorted template bed file.
cat /mnt/seagate10/genotypes/1000genomes503eur9m/chr[0-9]*.bim | awk 'BEGIN{OFS="\t";} {print("chr"$1, $4-1, $4-1+length($5), $2)}' | sort -k1,1 -k2,2n | gzip -c > /mnt/seagate10/projects/annotations/data/test/template.1000genomes503eur9m.sorted.bed.gz



5. Intersect annotations bed file with template bed file using bedtools.
/mnt/seagate10/projects/github/bedtools2/bedtools intersect -a /mnt/seagate10/projects/annotations/data/test/template.1000genomes503eur9m.sorted.bed.gz -b /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.annot.complete.sorted.bed.gz -wa -wb -sorted | gzip -c > /mnt/seagate10/projects/annotations/data/test/template.1000genomes503eur9m.complete_annot_hg19.intersect.txt.gz



6. Create binary annotations (this step is slow).
python annot2annomat.py /mnt/seagate10/projects/annotations/data/test/template.1000genomes503eur9m.complete_annot_hg19.intersect.txt.gz /mnt/seagate10/projects/annotations/data/test/template.1000genomes503eur9m.sorted.bed.gz /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.annomat.txt.gz



7. Create unique annotations (many variants belong to multiple annotation categories, this step assigns a single category to each variant, prioritizing categoryes according to "All SNPs ...").
python uniq_annot.py /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.annomat.txt.gz /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.annomat.uniq.txt.gz



8. Calculate LD r2 coefficients with parameters from "All SNPs are not created equal" supplement.
for ((i=1;i<23;i++)); do plink --bfile /mnt/seagate10/genotypes/1000genomes503eur9m/chr${i} --r2 inter-chr gz yes-really --ld-window-r2 0.2 --out /mnt/seagate10/genotypes/1000genomes503eur9m/schork/chr${i}.schork.r2; done



9. Create ld-induced categories for the experiment (this step is very slow).
python ld_informed_annot.py
python ld_informed_annot.py /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.annomat.uniq.txt.gz /mnt/seagate10/genotypes/1000genomes503eur9m/schork/ /mnt/seagate10/projects/annotations/data/test/knownGene_ucsc_hg19.annot.ld_informed.txt.gz

