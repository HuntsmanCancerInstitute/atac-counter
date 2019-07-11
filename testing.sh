
#!/bin/bash

# Assume you are in the outs/ directory of cellranger-atac count.

# subsample to 1% of aligned reads
samtools view -bs 42.01 possorted_bam.bam > subsampled.bam
# index the bam file
samtools index subsampled.bam
# add the RG tags and create a header file 'header.sam'
python cb2rg.py -i subsampled.bam -o subsampled.clean.bam
# insert header into file.
samtools reheader -P -i header.sam subsampled.clean.bam > subsampled.clean.reheader.bam
# run feature counts <fix the path>
FEATCOUNTS=/path/to/Subread/1.6.3/bin/featureCounts
$FEATCOUNTS -T 10 --readExtension5 1000 --readExtension3 1000 --byReadGroup \
  -a AnnotationFile.txt -o subsampled.featurecounts.txt \
  -F SAF \
  subsampled.clean.reheader.bam
# Drop the file prefix from Subread's columns
sed 's/subsampled.clean.reheader.bam://g' subsampled.featurecounts.txt > subsampled.featurecounts.shortname.txt
# subset only the GEMS labeled as filtered cells
python subset_cells.py -i subsampled.featurecounts.shortname.txt -b filtered_peak_bc_matrix/barcodes.tsv
