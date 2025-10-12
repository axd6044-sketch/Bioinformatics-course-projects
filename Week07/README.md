## Week 7 Assignment - Reusable Makefile, using bigwig and bam file visualization

A makefile is given in this repository. The aims of the makefile are:
1. Downloads the genome 
2. Downloads fastq reads by taking SRR number in parameters
3. Indexes the reference genome
4. Aligns the single and paired end reads to the reference genome
5. Converts the bam files generated to bigwig for IGV visualization

## Makefile options: 

### get_genome 
Use the get_genome target to download both the FASTA (genome sequence) and GFF (annotation) files for a specific NCBI accession.
```
make get_genome genome=GCF_000882815.3
```
These files will be unzipped and renamed for convenient use in later steps.

### 2. get_fastq and get_pairedfastq
Download single end fastq for Illumina NestSeq dataset and Paired-end fastq data for Illumina MiSeq dataset.
```
make get_fastq fastq=SRR3194430
make get_pairedfastq fastq=SRR3191544
```
### 3. index
```
make index genome_fa=ref/genome.fa
```
### 4. align
```
make align fastq=SRR3194430
make align fastq=SRR3191544
```
### 5. bigwig
```
make bigwig fastq=SRR3194430 genome_fa=ref/genome.fa
make bigwigfastq=SRR3191544 genome_fa=ref/genome.fa
```
## IGV visualization of bam files 

### First : Illumina NextSeq for single end dataset 

### Second : Illumina MiSeq for paired end dataset

## Briefly describe the differences between the alignment in both files.

## Briefly compare the statistics for the two BAM files.

## How many primary alignments dows each of your BAM files contain?

## What coordinate has the largest observed coverage?

## Select a gene of interest. How many alignments on a forward strand cover the gene?
