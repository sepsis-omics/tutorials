#  Variants: Background

[from LSCC docs]

## Identifying SNPs

There are many methods of identifying SNPs, but all rely on the per-base evidence provided by all the reads that have mapped to particular position in the sequence. So, it's useful to aggregate the evidence from all reads that relate to a particular base in the sequence. One method is to generate a pileup: a summary of sequence information from the entire set of reads across each relevant base in the reference sequence along with quality metrics relating to base quality and the mapping quality for each/all reads.

[see Pileup file information](../../data/formats/pileup.md) FIXME: can't link
