# Spades - command line

FIXME: add general info about spades (copy from galaxy-spades)

## Start

- open your mGVL (e.g. via ssh in terminal)
- check spades is installed: `which spades`
- the path to spades should appear
- FIXME: it's already loaded? seems to work
- make a new directory for these analyses: `mkdir spades`

- copy the read files into your mGVL
    - e.g. R1.fastq and R2.fastq (e.g. download from the galaxy workshops)
    - in terminal, navigate to the folder containing these files
    - `scp <file names> <user@mGVL:/path to put data>`

## Run spades

```bash
spades.py -1 R1.fastq -2 R2.fastq -k 33,55,77 -o output
```

## Output

- Go to the folder called output
    - contigs.fasta
    - scaffolds.fasta
    - assembly_graph.fastg
    - subfolder: corrected reads [FIXME explain]

## What next?

- view assembly_graph.fasta in Bandage

- generate summary statistics with Quast
