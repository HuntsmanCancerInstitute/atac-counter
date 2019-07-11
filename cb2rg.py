#!/usr/bin/env python

import os
import argparse
import pysam

# parse input
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_bam', type=str,
        help='BAM file from cellranger-atac count',
        required = True)
parser.add_argument('-o','--output_bam', type=str,
        help='BAM file with RG tags added ',
        required = True)
parser.add_argument('-d','--header_file', type=str,
            default = 'header.sam',
            help='Name of modified sam header including RG tags.',
        required = False)
args = parser.parse_args()

# sanitize inputs
if not os.path.isfile(args.input_bam):
    raise Exception('Input file does not exist. Check the argument --input_bam.')
if os.path.isfile(args.output_bam):
    raise Exception('Output file already exists. Specify different name for --output_bam.')

# open connection to input and output files.
input_bam = pysam.AlignmentFile(args.input_bam, "rb")
output_bam = pysam.AlignmentFile(args.output_bam, "wb", template=input_bam)

#get SAM header
header = input_bam.header.to_dict()
# initialize read groups set
read_groups = set()
# keep track of alignments with a cell barcode
nmiss = 0
nkeep = 0

print("Parsing BAM file.")
for read in input_bam.fetch():
    if read.has_tag('CB'):
        cell_barcode = read.get_tag('CB')
        read_groups.add(cell_barcode)
        # replace RG tag with the cell barcode
        read.set_tag('RG', cell_barcode, 'Z', True)
        nkeep += 1
        output_bam.write(read)
    else:
        nmiss += 1
input_bam.close()
output_bam.close()

print('Reads missing CB tag (not reported): ' + str(nmiss))
print('Reads with CB tag (reported): ' + str(nkeep))

print("Writing header.sam with updated RG tags.")
# write out RG tags to header
rg_header = [{'ID' : id } for id in read_groups]
header['RG'] = rg_header
with pysam.AlignmentFile(args.header_file, "wh", header=header) as outf:
    a = pysam.AlignedSegment()
    outf.write(a)

# re-header the output_bam file in place
pysam.reheader("-P", "-i", args.header_file, args.output_bam)
