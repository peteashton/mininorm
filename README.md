# Mininorm - minimiser-based digital normalisation of large long-read DNA sequencing datasets

k-mer based normalisation has been used to perform digital normalisation of large DNA sequence datasets.  
khmer utilises a k-mer counting method to perform such normalisation in a single pass through a file of 
short-read sequences.  Recently minimiser and minhash approaches have been used to perform fast analyses 
of long-read data, and this project attempts to use minimisers to normalise long-read sequencing data in 
a similar manner to khmer's approach.  We anticipate this will be particularly useful in e.g., metagenomics
applications, where several subsets of the DNA in a sample are present at different abundances.  In normalising 
the data, we expect that the representation of the less abundant DNAs in the samples will be improved in 
subsequent assemblies of the data.

```
usage: mininorm.py [-h] [-o OUTFILE] [-r rejects-file] [-w w] [-k k]
                   [-c coverage] [-s stats-file] [-n counts-file]
                   inputfile

Digitally normalise long-read DNA sequence read files using k-mer minimisers

positional arguments:
  inputfile            FASTQ file of long-read DNA sequences

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        name of FASTQ file to store the downsampled reads (use
                        '-' for stdout)
  -r rejects-file, --reject rejects-file
                        name of FASTQ file to store the reads rejected as
                        liekly duplicate, if not specifed, rejected reads are
                        not collected
  -w w, --window-size w
                        window size (default=20)
  -k k, --kmer-size k   k-mer size (default=20)
  -c coverage, --coverage-threshold coverage
                        coverage threshold. Median minimiser count above which
                        a read will be discarded as a likely duplicate
                        (default=20)
  -s stats-file, --stats stats-file
                        filename in which to place details of each sequences
                        as it is analysed, if not specified, stats are not
                        collected
  -n counts-file, --counts counts-file
                        filename to store counts of all the minimisers, which
                        will be very large
```
