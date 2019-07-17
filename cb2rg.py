#!/usr/bin/env python

import os
import argparse
import pysam

# parse input
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_bam', type=str,
                    help='BAM file from cellranger-atac count',
                    required=True)
parser.add_argument('-o', '--output_bam', type=str,
                    help='BAM file with RG tags added ',
                    required=True)
parser.add_argument('-d', '--header_file', type=str,
                    default='header.sam',
                    help='Name of modified sam header including RG tags.',
                    required=False)
parser.add_argument('-b', '--barcodes', type=str,
                    help='Path to Cell Ranger file filtered_peak_bc_matrix/barcodes.tsv',
                    required=True)
args = parser.parse_args()

# sanitize inputs
if not os.path.isfile(args.input_bam):
    raise Exception(
        'Input file does not exist. Check the argument --input_bam.')
if os.path.isfile(args.output_bam):
    raise Exception(
        'Output file already exists. Specify different name for --output_bam.')
if not os.path.isfile(args.barcodes):
    raise Exception(
        'Input file does not exist. Check the argument --barcodes.')


# open connection to input and output files.
input_bam = pysam.AlignmentFile(args.input_bam, "rb")
output_bam = pysam.AlignmentFile(args.output_bam, "wb", template=input_bam)
fname_barcodes = args.barcodes
file_barcodes = open(fname_barcodes, "r")

print('Reading in filtered cell barcodes...')
filtered_barcodes = [line.rstrip('\n') for line in file_barcodes]
print('Found ' + str(len(filtered_barcodes)) + ' filtered barcodes.')

# get SAM header
header = input_bam.header.to_dict()
# initialize read groups set
read_groups = set()
# keep track of alignments with a cell barcode
nmiss = 0
nkeep = 0
nfiltered = 0

print("Parsing BAM file.")
for read in input_bam.fetch():
    if read.has_tag('CB'):
        cell_barcode = read.get_tag('CB')
        if cell_barcode in filtered_barcodes:
            read_groups.add(cell_barcode)
            # replace RG tag with the cell barcode
            read.set_tag('RG', cell_barcode, 'Z', True)
            nkeep += 1
            output_bam.write(read)
        else:
            nfiltered += 1
    else:
        nmiss += 1
input_bam.close()
output_bam.close()

print('\n')
print('No. of alignments missing CB tag: ' + str(nmiss))
print('No. of alignments contains filtered out CB tag: ' + str(nmiss))
print('No. of Reads with passing-filter CB tag: ' + str(nfiltered))
print('No. of read groups: ' + str(len(read_groups)))
print('\n')

print("Writing header.sam with updated RG tags.")
# write out RG tags to header
rg_header = [{'ID': id} for id in read_groups]
header['RG'] = rg_header
with pysam.AlignmentFile(args.header_file, "wh", header=header) as outf:
    a = pysam.AlignedSegment()
    outf.write(a)

# re-header the output_bam file in place
pysam.reheader("-P", "-i", args.header_file, args.output_bam)
