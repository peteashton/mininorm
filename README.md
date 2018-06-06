# Mininorm - minimiser-based digital normalisation of large long-read DNA sequencing datasets

k-mer based normalisation has been used to perform digital normalisation of large DNA sequence datasets.  
khmer utilises a k-mer counting method to perform such normalisation in a single pass through a file of 
short-read sequences.  Recently minimiser and minhash approaches have been used to perform fast analyses 
of long-read data, and this project attempts to use minimisers to normalise long-read sequencing data in 
a similar manner to khmer's approach.  We anticipate this will be particularly useful in e.g., metagenomics
applications, where several subsets of the DNA in a sample are present at different abundances.  In normalising 
the data, we expect that the representation of the less abundant DNAs in the samples will be improved in 
subsequent assemblies of the data.

