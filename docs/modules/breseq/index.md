# Breseq

## Overview

Breseq is a tool to find variants, by comparing sequencing reads against a reference.

http://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing
https://www.ncbi.nlm.nih.gov/pubmed?Db=pubmed&Cmd=ShowDetailView&TermToSearch=24838886

## Get data

ref sequence in .gbk
R1.fq
R2.fq


## Run

```
breseq -j 16 -r ref.gbk R1.fq R2.fq -o output
```


-j number of cores
-o output dir


## Results

Navigate to output directory.

use readlink command

go into own mac terminal

rsync -av (output)

then open in browser address bar
