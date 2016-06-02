# File formats


File format | File extensions | Description
----------  | --------------  | -----------
FASTA       |                 | <tx>>SequenceID <br> sequence </tx> <br> [More information.](https://en.wikipedia.org/wiki/FASTA_format)
            | .fasta  .fas  .fa | generic FASTA file
            | .fna            | nucleotides
            | .ffn            | nucleotides for coding regions only?
            | .faa            | amino acids
            | .frn            | nucleotides? for non-coding RNA regions
FASTQ       | .fastq          | <tx>@SequenceID <br> sequence <br> + <br> nucleotide quality scores</tx>
SAM         | .sam            | tab-delimited text file of reads aligned to a reference. e.g. mapped position, sequence, quality scores. [More information.](https://samtools.github.io/hts-specs/SAMv1.pdf)
BAM         | .bam            | Compressed version of SAM file.
GBK         | .gbk            | Genbank format. Sequence information, features, protein translations, DNA sequence.
GFF         | .gff            | General Feature Format. Tab-delimited. Each line is a feature. also known as GTF? mention gff3? e.g. it needs to be version 3? [More information.](https://en.wikipedia.org/wiki/General_feature_format)
VCF         | .vcf            | Variant Call Format. Tab-delimited. Header, then one line per identified variant. 
