# Prokka on commandline

- ssh to your mGVL
- wget data - scaffolds of .fna for several species
- (alternatively - have already put data on your mGVL)
- type in:
```
prokka --outdir [name of output folder for 1 sample] --locustag [tag eg sample number] [fna filename]
```
- make a new directory for the gff files e.g. "gff_files"
- mv -v */*.gff gff_files/  (means move any .gff files into that folder)

- FIXME: prokka options see manual
