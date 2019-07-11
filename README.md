
# Introduction

Suppose you would like to use different feature boundaries than what 10x
Genomics reports. This repository contains an example on how to count
alternative user-defined features for 10x Genomics scATAC-seq data. The code uses `samtools`, `pysam`, `Subread`'s `featureCounts` software, with a few linux utilities.

# Minimal working example

Suppose you are in the outs/ directory produced by `cellranger-atac count`. For the sake of the example, down sample to 1% of the aligned reads.

```bash
# subsample to 1% of aligned reads
samtools view -bs 42.01 possorted_bam.bam > subsampled.bam
# index the bam file
samtools index subsampled.bam
```


Run the python script `cb2rg.py`, which replaces the RG tag with the cell
barcode (CB) tag) and then creates a header file 'header.sam'.

```bash
python cb2rg.py -i subsampled.bam -o subsampled.clean.bam
# insert header into modified bam.
samtools reheader -P -i header.sam subsampled.clean.bam > subsampled.clean.reheader.bam
```

Run featureCounts with the user-defined features in AnnotationFile.txt.
Importantly, quantify with the `--byReadGroup` option.

```bash
# fix the path to featureCounts
FEATCOUNTS=/path/to/Subread/1.6.3/bin/featureCounts
$FEATCOUNTS -T 10 \
  --byReadGroup \
  --readExtension5 1000 --readExtension3 1000
  -a AnnotationFile.txt -o subsampled.featurecounts.txt \
  -F SAF \
  subsampled.clean.reheader.bam
# Drop the file prefix from Subread's columns
sed 's/subsampled.clean.reheader.bam://g' \
  subsampled.featurecounts.txt > subsampled.featurecounts.shortname.txt
# drop intermediate file
rm subsampled.featurecounts.txt
```

Subset only the GEMs labeled as filtered cells

```bash
python subset_cells.py -i subsampled.featurecounts.shortname.txt \
  -b filtered_peak_bc_matrix/barcodes.tsv \
  -o filtered_counts.txt
# drop intermediate file it is huge.
rm subsampled.featurecounts.shortname.txt
```
