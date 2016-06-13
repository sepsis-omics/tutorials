# SAM File

[from LSCC docs]

Sequence Alignment/Map (SAM) format records all information relevant to how a set of reads aligns to a reference genome. A SAM file has an optional set of header lines describing the context of the alignment, then one line per read, with the following format:

11 mandatory fields (+ variable number of optional fields)
1    QNAME: Query name of the read
2    FLAG
3    RNAME: Reference sequence name
4    POS: Position of alignment in reference sequence
5    MAPQ: Mapping quality (Phred-scaled)
6    CIGAR: String that describes the specifics of the alignment against the reference
7    MRNM
8    MPOS
9    ISIZE
10  SEQQuery: Sequence on the same strand as the reference
11  QUAL: Query quality (ASCII-33=Phred base quality)


## SAM example
SRR017937.312 16   chr20 43108717 37   76M   *    0    0
TGAGCCTCCGGGCTATGTGTGCTCACTGACAGAAGACCTGGTCACCAAAGCCCGGGAAGAGCTGCAGGAAAAGCCG
?,@A=A<5=,@==A:BB@=B9(.;A@B;>@ABBB@@9BB@:@5<BBBB9)>BBB2<BBB@BBB?;;BABBBBBBB@

For this example:

QNAME = SRR017937.312 - this is the name of this read
FLAG = 16 - see the format description below
RNAME = chr20 - this read aligns to chromosome 20
POS = 43108717 - this read aligns the sequence on chr20 at position 43108717
MAPQ = 37 - this is quite a high quality score for the alignment (b/w 0 and 90)
CIGAR = 76M - this read aligns to the reference segment across all bases (76 Matches means no deletions or insertions. Note that 'aligns' can mean 'aligns with mismatches' - mismatches that don't affect the alignment are not recorded in this field)
MRNM = * - see the format description below
MPOS = 0 as there is no mate for this read - the sequenced DNA library was single ended, not mate paired*.
ISIZE = 0 as there is no mate for this read
SEQQuery = the 76bp sequence of the reference segment
QUAL = per-base quality scores for each position on the alignment. This is just a copy of what is in the FASTQ file

NOTE: reads are shown mapped to the "sense" strand of the reference, and bases are listed in 5' -> 3' order. This is important because an actual read might be from the other strand of DNA. The alignment tool will try to map the read as it is, and also the reverse compliment. If it was on the other strand then the reverse compliment is shown in the SAM file, rather than the original read itself

[More information.](https://samtools.github.io/hts-specs/SAMv1.pdf)

![SAM file](sam.png) [FIXME link]
