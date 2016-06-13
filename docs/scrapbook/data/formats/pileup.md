# Pileup format

[from LSCC docs]


## Pileup Format
A pileup file has as many lines as there are bases in the reference sequence that are aligned with reads in the SAM/BAM file. Each line contains information about every base found in the sequence reads that corresponds to the reference base on that line. The format of a pileup record is:

- ReferenceSeq [string] - name of the reference sequence
- Coordinate [integer] - position in the reference sequence
- ReferenceBase [A/C/G/T/N] - reference base at that position
- Num. Reads [integer] - number of reads aligning to that base
- ReadBases [variable length string, see below]
- BaseQualities [variable length string, Phred encoded]

ReadBases:

- Each separate read that covers the base is represented here. The more reads that cover this base, the longer this string
- . = match on forward strand for that base
- , = match on reverse strand
- ACGTN = mismatch on forward
- acgtn = mismatch on reverse
- +[0-9]+[ACGTNacgtn]+' = insertion between this reference position and the next
- [0-9]+[ACGTNacgtn]+' = deletion between this reference position and the next
- ^ = start of read
- $ = end of read

BaseQualities = one character per base in ReadBases, ASCII encoded Phred scores

Example:

- chr1 272 T 24 ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
- chr1 273 T 23 ,.....,,.,.,...,,,.,..A <<<;<<<<<<<<<3<=<<<;<<+
- chr1 274 T 23 ,.$....,,.,.,...,,,.,... 7<7;<;<<<<<<<<<=<;<;<<6
- chr1 275 A 23 ,$....,,.,.,...,,,.,...^l. <+;9*<<<<<<<<<=<<:;<<<<
- chr1 276 G 22 ...T,,.,.,...,,,.,.... 33;+<<7=7<<7<&<<1;<<6<

In this example there are 5 chromosomal positions represented, with between 22 and 24 reads aligning to each of the positions. There are two mismatches: an 'A' in position 273 and a 'T' in position 276. That is, only a single read contained that mismatch in either case.
More information on pileup format here: <http://samtools.sourceforge.net/pileup.shtml>

## Pileup file in Galaxy

![Pileup image](./images/pileup.png)
