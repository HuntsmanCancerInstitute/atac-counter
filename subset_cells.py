#!/usr/bin/env python

import re
import os
import argparse

# parse input
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_featurecounts', type=str,
        help='Path to file output from featureCounts program.',
        required = True)
parser.add_argument('-b','--barcodes', type=str,
        help='Path to Cell Ranger file filtered_peak_bc_matrix/barcodes.tsv',
        required = True)
parser.add_argument('-o','--outfile', type=str,
        default = 'filtered_counts.txt',
        help='Outfile name',
        required = True)
args = parser.parse_args()

# fname_featcounts = 'subsampled.featurecounts.shortname.txt'
fname_featcounts = args.input_featurecounts
fname_barcodes = args.barcodes
file_featcounts = open(fname_featcounts, "r")
file_barcodes = open(fname_barcodes, "r")

# get 'featureCounts' command
command_featureCounts = file_featcounts.readline()
# get column names
print('Reading in columnn names...')
raw_col_names = file_featcounts.readline()
col_names = re.split(r'\t+', raw_col_names.rstrip('\n'))

print('Reading in filtered cell barcodes...')
filtered_barcodes = [line.rstrip('\n') for line in file_barcodes]
print('Found ' + str(len(filtered_barcodes)) + ' filtered barcodes.')
# get indices (use base 1 b/c of cut)
print("Finding column indices of filtered cell barcodes...")
idx_filtered_barcodes = [str(idx) for idx,col in enumerate(col_names,1) if col in filtered_barcodes]
print('Found ' + str(len(idx_filtered_barcodes)) + ' barcodes for subsetting.')

if len(filtered_barcodes) > len(idx_filtered_barcodes):
    raise Exception("""Did not find all the barcodes in featureCounts file.
    Exiting program.""")

barcode_idxs = ','.join(idx_filtered_barcodes)
anno_idxs = ','.join([str(i) for i in range(1,7)])
idxs = ','.join([anno_idxs,barcode_idxs])
# idxs = barcode_idxs
print('Applying linux cut program to subset indices...')
fname_out = os.path.basename(fname_featcounts)
cmd_cut = 'cut -f ' + idxs + r' < ' + fname_featcounts + r' > ' + args.outfile
exitval = os.system(cmd_cut)
if exitval != 0:
    raise Exception('The cut command failed with error code: ' + exitval)
else:
    print('Success! The subsetted featureCounts file is named ' + args.outfile)
