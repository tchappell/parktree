# parktree
parallel k-tree clustering for biological sequences using topological signatures 

## Usage
`./ParKTree (options) [fasta input]`

## Options
* -sw [signature width (default = 256)]
* -k [kmer length (default = 5)]
* -d [signature density (default = 0.0476..)]
* -o [tree order (default = 10)]
* -c [capacity (default = 1000000)]
* --fasta-output

## Requirements

A version of gcc (or compatible compiler) with support for OpenMP and `__builtin_popcountll()`. 

## Installation

Type `make` to install the software. If make or gcc is not available, this software only consists of a single source file which can be compiled manually. If OpenMP is not available the `#pragma` directives can be ignored to build a single-threaded version of the software. If `__builtin_popcountll()` is not available, it can be replaced with whatever builtin is needed to emit a 64-bit `POPCNT` instruction with your compiler and architecture.

## Operation

Run ParKTree, passing it the fasta file containing the sequences to be clustered and any options needed. It will then produce a series of clusters as output to standard output, which can then be redirected as necessary. Note that both the sequences and the signatures generated for each sequence are stored in memory during clustering.

ParKTree is multithreaded with OpenMP. The number of threads used can hence be controlled with OpenMP environment variables such as `OMP_NUM_THREADS`.

## Detailed description of options

### -sw [signature width]

This is the number of bits to use to store the signatures that are generated to represent each sequence. This value must be a multiple of 64.

The signature width affects the representational capacity of the signatures, and when sequences are long and/or long k-mer lengths are used to represent the signatures, a larger signature width can be beneficial. This comes with a cost to signature generation and clustering speed.

### -k [kmer length]

This is the number of characters making up each k-mer. k-mers of the given length are created for every position within each sequence to represent that sequence for the purposes of signature construction. For example, with a k-mer length of 5, the sequence `ACAAGATGCCATTG` results in the following k-mers: `ACAAG CAAGA AAGAT AGATG GATGC ATGCC TGCCA GCCAT CCATT CATTG`.

The k-mer length plays an important role in both the representational capacity of the generated signatures and the fragility of those sequences, and should only be tweaked with testing to ensure that the new value improves performance.

### -d [signature density]

When signatures are created from sequences, each k-mer is hashed into a k-mer signature of the density provided to this parameter. The k-mer signatures are then combined to create a signature for the sequence. The signature density determines the portion of bits set in each k-mer signature. This value may need to be tweaked based on the length of the signatures and/or the k-mer length.

### -o [tree order]

The number of child nodes each node of the tree can reach before splitting.

### -c [tree capacity]

The number of nodes that can be allocated to the tree during construction.

### --fasta-output

By default ParKTree will produce a two-column CSV consisting of the sequence ID and cluster ID of each sequence. An alternative output is available by passing in this parameter; instead, ParKTree will produce a fasta-format file containing the same sequences passed in, but with the name of each sequence replaced with the cluster number that sequence is a part of.
