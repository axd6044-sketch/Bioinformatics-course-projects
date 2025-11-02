

## Prerequisites

```bash
mamba activate bioinfo
```

## Pipeline Overview

This pipeline automates the processing of RNA-seq data using a Makefile and design.csv file. The pipeline handles both single-end and paired-end reads, producing aligned BAM files and visualization files.

## Running the Pipeline

### 1. Prepare design file 
```bash
make design PRJNA=PRJNA313294
```
### 2. Download the reference genome 
Download and index the reference genome (only needed once):

```bash
# download human chromosome 22 
make get_genome species=human genome=chr22 REF=ref
# download whole human genome
make get_genome species=human genome=hg38_ucsc REF=ref
# download zika genome
make get_genome species=zika genome=GCF_000882815.3 REF=ref

#index only human chromosome 22
make genome_index species=human REF=ref/human_chr22.fa
# for whole human genome
make genome_index species=human REF=ref/human_hg38.fa
# or for zika
make genome_index species=zika REF=ref/ZikaVirus_genome.fa

```
## For processing single sample (single or paired-end), run:

```bash
make get_fastq srr=<SRR_ID> reads=reads fastqcreports=reads/fastqc_reports
make alignreads genome_fa=ref/hmanchr22 reads=reads srr=SRR3191544 bam=bam
make bigwig srr=SRR3191544 bam=bam genome_fa=ref/hmanchr22
make call_variants species=human REF=ref/human_hg38.fa bam=bam srr=SRR3194430 vcf=vcf

```
## Parallel processing of multiple samples 
To process all samples in parallel from design.csv, make sure GNU Parallel can create temporary files. On macOS the default temp dir sometimes isn't writable from conda environments — create a per-user tmpdir and pass it with `--tmpdir`.

Generate the design file
```
make design
```
 
```
# Create a safe tmpdir (do this once per session)
mkdir -p ~/parallel_tmp
chmod 700 ~/parallel_tmp

# Download FASTQ files for all human samples
awk -F',' 'NR>1 {print $1}' design.csv \
  | parallel --tmpdir ~/parallel_tmp --jobs 2 --bar \
      'make get_fastq srr={} READS=reads fastqcreports=reads/fastqc_reports'

# Align all samples to human chr22 (hg38 build)
awk -F',' 'NR>1 {print $1}' design.csv \
  | parallel --tmpdir ~/parallel_tmp --jobs 2 --bar \
      'make alignreads species=human REF=ref/human_chr22 READS=reads srr={} bam=bam'

# Generate BigWig files for visualization
awk -F',' 'NR>1 {print $1}' design.csv \
  | parallel --tmpdir ~/parallel_tmp --jobs 2 --bar \
      'make bigwig species=human srr={} bam=bam REF=ref/human_chr22.fa'

# Call variants for each human sample (per-sample)
awk -F',' 'NR>1 {print $$1}' design.csv | while read srr; do
    make call_variants species=human REF=ref/human_chr22.fa bam=bam srr=$$srr vcf=vcf
done

#Merge all individual VCFs into one multi-sample file
# creating multisample VCF
bcftools merge vcf/*.vcf.gz -O z -o vcf/multisample_merged.vcf.gz
bcftools index -t vcf/multisample_merged.vcf.gz

```

Visualize in IGV or JBrowse
# ========================================================
# Load these in IGV:
#   - Reference genome: ref/human_chr22.fa
#   - Annotation file:  ref/human_hg38.gtf.gz
#   - Aligned BAMs:     bam/SRR*_human.bam
#   - Variants:         vcf/multisample_merged.vcf.gz
#   - Coverage tracks:  bam/SRR*_human.bw

## Output Files

Week10/
├── reads/
│   ├── SRR3194430.fastq
│   ├── fastqc_reports/
│   │   └── SRR3194430_fastqc.html
├── ref/
│   ├── human_hg38.fa
│   ├── human_hg38.gtf.gz
│   └── ZikaVirus_genome.fa
├── bam/
│   ├── SRR3194430_human.bam
│   ├── SRR3194430_human.bam.bai
│   ├── SRR3194430_human.bw
│   └── SRR3194430_human_align_stats.txt
├── vcf/
│   ├── SRR3194430_human.vcf.gz
│   ├── SRR3194430_human.vcf.gz.tbi
│   └── multisample_merged.vcf.gz
└── design.csv


---

#Screenshot from IGV visualizing three SRR bam files and their bw files
<img src="img1.png" alt="image" width="800">

**Note:**
- The Makefile automatically handles both single-end and paired-end data.
- Adjust the number of jobs in GNU Parallel as needed for your system.
- All output directories are created automatically if they do not exist.
