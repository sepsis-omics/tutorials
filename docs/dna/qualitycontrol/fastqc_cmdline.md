# FastQC - commandline

## Start
- on your local machine, make sure XQuartz is installed (but it doesn't have to be open - it will open automatically later).
- in terminal, ssh to your virtual machine with -X and -Y, e.g. `ssh -X -Y ubuntu@111.111.111.111` (the -X and -Y means it will use your local XQuartz to display some files).
- `module load fastqc_dist_0_10_1` [FIXME: this may have a different name/ or be already loaded)
- navigate to where you want to make a FastQC analysis folder.
- Make a folder: `mkdir fastqc_analyses`
- Move to that folder: `cd fastqc_analyses`

## Input
- [FIXME: get the data into this folder]

## Run
- `fastqc R1reads.fastq` [runs fastqc]
- type `fastqc --help` to see settings that you can change, and defaults
- FIXME: any to change
- FIXME: repeat for R2reads?

## Output
- R1reads_fastqc: folder containing the output, e.g. <fn>fastqc_report.html</fn>
- to view this, type: `firefox fastqc_report.html` - firefox should open and display the report
(you may get an error message in terminal but ignore this)
