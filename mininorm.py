#!/bin/env python3

from collections import Counter
from statistics import median
from Fastq import FastqReader

import sys
import argparse

# Set up appropriate argument parsing

parser = argparse.ArgumentParser(description="Digitally normalise long-read DNA sequence read files using k-mer minimisers")
parser.add_argument("input-file", nargs=1, type=str, help="FASTQ file of long-read DNA sequences")
parser.add_argument("-o", "--outfile", type=str, help="name of FASTQ file to store the downsampled reads (use '-' for stdout)", default="-")
parser.add_argument("-r", "--reject", type=str, help="name of FASTQ file to store the reads rejected as liekly duplicate, if not specifed, rejected reads are not collected", default="", metavar="rejects-file")
parser.add_argument("-w", "--window-size", type=int, help="window size (default=20)", default=20, metavar="w")
parser.add_argument("-k", "--kmer-size", type=int, help="k-mer size (default=20)", default=20, metavar="k")
parser.add_argument("-c", "--coverage-threshold", type=int, help="coverage threshold.  Median minimiser count above which a read will be discarded as a likely duplicate (default=20)", default=20, metavar="coverage")
parser.add_argument("-s", "--stats", type=str, help="filename in which to place details of each sequences as it is analysed, if not specified, stats are not collected", default="", metavar="stats-file")
parser.add_argument("-n", "--counts", type=str, help="filename to store counts of all the minimisers, which will be very large", default="", metavar="counts-file")

args = parser.parse_args()

# Process the arguments

if args.outfile == "-":
    outfile = sys.stdout
else:
    outfile = open(args.outfile, 'w')

# Optional output files for rejected reads, read statistics and minimiser counts

rejects_file = False
if args.reject != "":
    rejects_file = open(args.reject, 'w')

stats_file = False
if args.stats != "":
    stats_file = open(args.stats, 'w')
    print("Sequence length\tNum minimisers\tCumulative minimisers\tNew minimisers\tMedian minimiser count", file=stats_file)

counts_file = False
if args.counts != "":
    counts_file = open(args.counts, 'w')
    print("Count")

# Core parameters for window size, kmer size and threshold for coverage

window_size = args.window_size
kmer_size = args.kmer_size
coverage_threshold = args.coverage_threshold

num_seqs = 0

minimiser_counts = Counter()
last_minimisers = 0

for read in FastqReader(args.input_file):
    num_seqs += 1
    seq = read.seq
    revcomp = read.reverse_complement()
    minimisers = set()

# Extract the (overlapping) kmers from the sequence and its reverse complement, and reverse the order
# of the RC kmers, so that each kmer and its reverse complement are in the corresponding elements of the
# two lists

    kmers = [seq[i:i+kmer_size] for i in range(len(seq)-kmer_size)]
    rkmers = [revcomp[i:i+kmer_size] for i in range(len(revcomp)-kmer_size)]
    rkmers.reverse()

# For each window along the sequence, add the minimal kmer to the list of minimisers, regardless of
# whether it came from the forward or RC sequence. Note: profiling shows that almost half of the time
# is spent in these calls to the "min" function, probably becase we do it so many times per sequence

    for i in range(len(kmers) - window_size):
        minimisers.add(min(kmers[i:i+window_size] + rkmers[i:i+window_size]))

# We could add end minimisers here, but these are really only required for strict end-to-end overlap detection, which is not 
# what we are doing here...

# Add to the counts for the minimisers we found in this sequence

    for minimiser in minimisers:
        minimiser_counts[minimiser] += 1

# Calculate the median count for the minimisers from this sequence

    median_minimiser_count = median([minimiser_counts[minimiser] for minimiser in minimisers])

# Decide whether to accept or reject the read

    if media_minimiser_count <= coverage_threshold:
        print(read.to_fastq(), file=outfile)
    else:
        if rejects_file:
            print(read.to_fastq(), file=rejects_file)

# Write the stats for this sequence out to the stats file

    all_minimisers = len(minimiser_counts)
    if stats_file:
        print("%d\t%d\t%d\t%d\t%d" % 
            (len(seq), len(minimisers), all_minimisers, all_minimisers-last_minimisers, median_minimiser_count), 
            file=stats_file)
    last_minimisers = all_minimisers
    if stats_file and num_seqs % 1000 == 0:
        stats_file.flush()

outfile.close()
stats_file.close()

# Output the counts, if required

if counts_file:
    print("Count", file=counts_file)
    for minimiser, count in minimiser_counts.items():
        print(count, file=counts_file)
    counts_file.close()
