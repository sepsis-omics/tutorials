# Variant calling with Snippy

## Background

Variant calling is the process of identifying differences between to genome samples.
Usually differences are limited to single nucleotide polymorphisms (SNPs) and small insertions and deletions (indels). Larger structural variation such as inversions, duplications and large deletions are not typically covered by "variant calling".

## Learning Objectives

At the end of this tutorial you should be able to:

1. Find variants between a reference genome and a set of reads
2. Determine the effect of those variants on genomic features
3. Understand if the SNP is potentially affecting the phenotype

<!---
## Experiment

FIXME: You are working on a bug and you notice one of them is forming smaller colonies than normal. You want to find out why this small colony variant (SCV) is doing at the DNA level.
--->

## Prepare reference

We will use the same data that we used in the [Assembly with Spades tutorial.](../spades/index.md) This should still be in your current galaxy history. If not, re-import the data into a new history using the instructions in that tutorial.

For variant calling, we need a reference genome that is of the same strain as the input sequence reads. This can be the plain nucelotide sequence in FASTA format, or annotated in GBK or GFF3 format.

- upload <fn>wildtpe.gbk</fn> file

Our reference is the <fn>wildtype.gbk</fn> file and our reads are <fn>mutant_R1.fastq</fn> and <fn>mutant_R2.fastq</fn>.

- screenshot of files

<!--
!!! note
    Please make sure your reference genome includes all chromosomes and plasmids
-->

## Call variants with Snippy

- Go to <ss>Tools &rarr; NGS Analysis &rarr; NGS: Variant Analysis &rarr; snippy</ss>
- For <ss>Reference type</ss> select *Genbank*.
- Then for <ss>Reference Genbank</ss> choose the <fn>wildtype.gbk</fn> file.
- For <ss>Single or Paired-end reads</ss> choose *Paired*.
- Then choose the first set of reads, <fn>mutant_R1.fastq</fn> and second set of reads, <fn>mutant_R2.fastq</fn>.
- For <ss>Cleanup the non-snp output files</ss> select *No*.
- Click <ss>Execute</ss>.

to do:

- look in the output tab file to see list of SNPs
- what types, what do they mean, snpeff predictions
- any other snippy output to look at?


- go to Graph/Display data
- JBrowse
- Fasta seq: choose wildtype.fna
- Produce standalone- yes
- genetic code - 11
- insert track group - name it
- insert annotation track
   - for snps gff, vcf, bam  (autogenerate snp track - yes)
- this makes a file - click on eye to view
- select tracks from the left
- go to position x to see a snp
- go back to snippy tab output to see the snpeff result for this snp





<!---
FIXME: talk about multimapping reads?

## Filter variants
FIXME:  vcffilter? something else?  mindepth, homozygous?

## Annotate consequencs

-- how to get ref genome in?

FIXME: snpEff - but it is hard to add a genome

!!! hint
    Just use Snippy and all this will happen magically?
--->


## What next?

<!--
* SNPs can be used to build [phylogentic trees](/trees/index.md).
-->
