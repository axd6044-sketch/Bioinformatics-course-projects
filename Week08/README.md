## Assignment for Week08 
the assignment requirements:
Makefile — produces multiple BAMs named by sample, plus FastQC and BigWig
design.csv — connects SRR → sample name → layout (SINGLE/PAIRED)
README.md — explains how to run the pipeline including GNU Parallel from design.csv

## Requirements
```
mamba activate bioinfo
```

### First finding out what is the PRJNA number 
PRJNA313294

### The Genome file we need is of Zika Virus 
This will be downloaded once in the beginning 
GCF = GCF_000882815.3 (Zika MR766)

### Creating the design file with all the samples metadata 
The design.csv file contains all the SRR numbers and their meta data
```
$ bio search -H --csv PRJNA313294 > design.csv
```
Now some of the samples will have single end data and some paired -end and we need to download fastq separately for them therefore 
### we will first see if we have single end and paired end files: 
```
$ head -n 1 design.csv | tr ',' '\n' | nl
$ awk -F',' 'NR>1 {print $1, $14}' design.csv 
```
Output: 
     1	run_accession
     2	sample_accession
     3	sample_alias
     4	sample_description
     5	first_public
     6	country
     7	scientific_name
     8	fastq_bytes
     9	base_count
    10	read_count
    11	library_name
    12	library_strategy
    13	library_source
    14	library_layout
    15	instrument_platform
    16	instrument_model
    17	study_title
    18	fastq_url
    19	info

SRR3191542 PAIRED
SRR3191545 PAIRED
SRR3194429 SINGLE
SRR3191543 PAIRED
SRR3191544 PAIRED
SRR3194430 SINGLE
SRR3194431 SINGLE
SRR3194428 SINGLE


Now, this SRR numbers are required for running the makefile 
the makefile which was used in week07 will be used here 
the makefile is present in this repository too 
### the things this makefile can do : 
1. Downloads the genome 
2. Downloads fastq reads by taking SRR number in parameters
3. Indexes the reference genome
4. Aligns the single and paired end reads to the reference genome
5. Converts the bam files generated to bigwig for IGV visualization

the make all runs :  get_genome	get_fastq get_pairedfastq index alignsingle alignpaired bigwigse bigwigpe
For more details on how to run this makefile, Week07 README file can be referred. 

We will now run the single end files and paired end files using data from design file and input those during running makefile - how to do that is shown below

### Getting the genome which is common and has to be done once for all the samples 
```
make get_genome genome=GCF_000882815.3
```

### First dry run 
```
awk -F',' 'NR>1 && $14=="SINGLE" {print $1}' design.csv | \
parallel --dry-run "
  make get_fastq fastq={};
  make alignsingle fastq={};
  make bigwigse fastq={} genome_fa=ref/genome.fa
"
awk -F',' 'NR>1 && $14=="PAIRED" {print $1}' design.csv | \
parallel --dry-run "
  make get_pairedfastq fastq={};
  make alignpaired fastq={};
  make bigwigpe fastq={} genome_fa=ref/genome.fa
"


```
### Now real run 
```
awk -F',' 'NR>1 && $14=="SINGLE" {print $1}' design.csv | \
parallel --jobs 2 --bar "
  make get_fastq fastq={};
  make alignsingle fastq={};
  make bigwigse fastq={} genome_fa=ref/genome.fa
"

awk -F',' 'NR>1 && $14=="PAIRED" {print $1}' design.csv | \
parallel --jobs 2 --bar "
  make get_pairedfastq fastq={};
  make alignpaired fastq={};
  make bigwigpe fastq={} genome_fa=ref/genome.fa
"

```
### Outputs per sample:
FASTQ: reads/<sample>.fastq or reads/<sample>_R1.fastq, _R2.fastq
FastQC: reads/fastqc_reports/*_fastqc.html
BAM + BAI: alignments/<sample>.sorted.bam(.bai)
BigWig: coverage/<sample>.bw
Stats: stats/<sample>_flagstat.txt