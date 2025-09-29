# Week 05 Assignment :

## This assignment follows the workflow for analyzing the **Zika virus MR766 strain** used in  Tang et al. *Cell Stem Cell* (2016). We explore genome download, annotation, sequencing data retrieval,  basic statistics, and quality control. Zika virus, MR766 strain, from the Cell Stem Cell 2016 paper (Tang et al. Cell Stem Cell 2016).

## Identify the BioProject and SRR accession numbers for the sequencing data associated with the publication.
### Bioproject - PRJNA313294
### SRR accession numbers - SRR3194431 (Single-end Illumina dataset, NextSeq platform)
### Genome accession:  - GenBank: GCA_000882815.1  - RefSeq: GCF_000882815.3 

## Step 1: Write a Bash script:

### Make a new file
```
nano zika_pipeline.sh
```
### activation requirements before running the script
module load anaconda
conda activate bioinfo
conda install -c conda-forge jq
conda install -c bioconda seqtk cutadapt fastqc ncbi-datasets-cli sra-tools

### Content of the script is given below. 
### Paste from the next line
#!/usr/bin/env bash

#Requirements
#module load anaconda
#conda activate bioinfo
#conda install -c conda-forge jq
#conda install -c bioconda seqkit


# ---- 1. Download genome + annotation ----
echo ">> Downloading genome and annotation..."
datasets summary genome accession GCA_000882815.1 | jq
datasets download genome accession GCF_000882815.3 --include genome,gff3,gtf
unzip ncbi_dataset.zip -d ncbi_dataset -x README.md
# Rename for convenience
mv ncbi_dataset/ncbi_dataset/data/GCF_000882815.3/GCF_000882815.3_ViralProj36615_genomic.fna Zikagenome.fa
mv ncbi_dataset/ncbi_dataset/data/GCF_000882815.3/genomic.gff Zikagenome.gff

# ---- 2. Genome size ----
GENOME_SIZE=$(grep -v ">" Zikagenome.fa | tr -d '\n' | wc -c)
echo "Genome size (bp): ${GENOME_SIZE}" > genome_summary.txt

# ---- 3. Count features in GFF ----
echo "Genomic details"
cat Zikagenome.gff | cut -f 1 | sort | uniq

echo "Genome features are:"
cat Zikagenome.gff | grep -v '#' | cut -f 3 | sort | uniq -c | sort -nr

echo "The gene being coded is:"
grep -v '#' Zikagenome.gff | awk '$3=="gene" {print $9}'

# ---- 4.What other accessions are there ----
echo "Some other accessions include:" 
datasets summary genome taxon "Zika virus" \
  | grep -oE 'GC[FA]_[0-9]+\.[0-9]' \
  | sort -u

# ---- 5. Download sequence dataset (Bioproject PRJNA313294) ----
SRR_ID="SRR3194430" # paired-end Illumina

echo ">> Prefetching ${SRR_ID}..."
prefetch ${SRR_ID}


# ---- 6. Subsample reads (~10x coverage) ----
NREADS=1000
echo ">> Subsampling ${NREADS} reads ..."
mkdir -p reads

fasterq-dump --split-files ${SRR_ID} -O reads/
echo "Single-end detected."
seqkit stats reads/${SRR_ID}_subsample.fastq

### Bash script ends here
### Save and exit
### Make it executable 
```
chmod +x zika_pipeline.sh
```
## Run the script 
```
./zika_pipeline.sh
```
### The outputs should be - Zikagenome.fa, Zikagenome.gff, genome_summary.txt, read_stats.txt, trimmed/*.fastq

# ----7. Statistics ----
cd reads
echo ">> Basic statistics on ${SRR_ID}:" 
seqkit stats ${SRR_ID}_subsample.fastq
cd ..

# ---- 8. Quality Control: Evaluating FASTQC report ----
echo ">> Running FastQC..."
mkdir -p fastqc_reports
fastqc -t 4 -o fastqc_reports reads/${SRR_ID}_subsample.fastq

# ----9. Adapter trimming ----
echo ">> Running cutadapt..."
mkdir -p trimmed 
cutadapt -a AGATCGGAAGAGC \
  -o trimmed/${SRR_ID}.trimmed.fastq \
  reads/${SRR_ID}_subsample.fastq > trimmed/cutadapt_report.txt

  # Re-run FastQC on trimmed reads
fastqc -t 4 -o fastqc_reports trimmed/${SRR_ID}.trimmed.fastq 

echo "Pipeline Complete"


# Searched for another dataset SRA for the same genome but uses a different sequencing platform
## SRA1: SRR3194431 (single-end Illumina nextSeq dataset)
## SRA : SRR3191544 (paired-end Illumina MiSeq dataset)

SRR2="SRR3191544"
prefetch ${SRR2}
fasterq-dump --split-files ${SRR2} -O SRR2reads/
### Basic statistics of SRR2
seqkit stats SRR2reads/${SRR2}_1.fastq SRR2reads/${SRR2}_2.fastq > SRR2reads/${SRR2}_stats.txt
fastqc -o fastqc_reports SRR2reads/${SRR2}_1.fastq SRR2reads/${SRR2}_2.fastq

### Comparison notes:
NextSeq single-end (SRR3194431): shorter reads (1×150 bp), uniform quality, smaller files.

MiSeq paired-end (SRR3191544): longer reads (2×250 bp), more overlap, but quality drops at 3′ ends.

Both adequately cover the small Zika genome (~10.8 kb). MiSeq may require more trimming.
