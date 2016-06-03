# File formats

FIXME: link to images?


File format | File extensions | Description
----------  | --------------  | -----------
FASTA       |                 | <tx>>SequenceID <br> sequence </tx> <br> [More information.](https://en.wikipedia.org/wiki/FASTA_format)
            | .fasta  .fas  .fa | generic FASTA file
            | .fna            | nucleotides
            | .ffn            | nucleotides for coding regions only?
            | .faa            | amino acids
            | .frn            | nucleotides? for non-coding RNA regions
FASTQ       | .fastq          | <tx>@SequenceID <br> sequence <br> + <br> nucleotide quality scores</tx>
SAM         | .sam            | tab-delimited text file of reads aligned to a reference. e.g. mapped position, sequence, quality scores. [More](../sam.md) [FIXME: can't link]
BAM         | .bam            | Compressed version of SAM file. [More](../bam.md) [FIXME: can't link]
GBK         | .gbk            | Genbank format. Sequence information, features, protein translations, DNA sequence.
GFF         | .gff            | General Feature Format. Tab-delimited. Each line is a feature. also known as GTF? mention gff3? e.g. it needs to be version 3? [More information.](https://en.wikipedia.org/wiki/General_feature_format)
VCF         | .vcf            | Variant Call Format. Tab-delimited. Header, then one line per identified variant.
pileup      |                 | Each line is a nucleotide in the sequence, with information on how all the reads are mapped to that position. [More](../pileup.md) [FIXME: can't link]

there are some nice images etc in here https://docs.google.com/document/pub?id=1fouC29Lq0CXxQQCpuojrR5RXbdzMdxRf8ZID01XYNqI#h.18e90b8fc68f
