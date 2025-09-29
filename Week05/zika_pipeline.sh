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
SRR_ID= "SRR3194430" # paired-end Illumina

echo ">> Prefetching ${SRR_ID}..."
prefetch ${SRR_ID}


# ---- 6. Subsample reads (~10x coverage) ----
NREADS=1000
echo ">> Subsampling ${NREADS} reads ..."
mkdir -p reads

fasterq-dump --split-files ${SRR_ID} -O reads/
echo "Single-end detected."
seqkit sample -n ${NREADS} reads/${SRR_ID}_1.fastq > reads/${SRR_ID}_subsample.fastq