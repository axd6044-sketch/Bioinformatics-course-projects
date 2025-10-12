#This is a makefile for Wee07 assignment
#The purpose of this makefile is : 


#The outputs the makefile will produce:

SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

# ------ No changes in code beyond this point -----


all: get_genome	get_fastq get_pairedfastq index alignsingle alignpaired bigwig
	echo ">> All steps completed successfully."

# Download reference genome
# Usage: make get_genome genome=GCF_000882815.3
get_genome:
	mkdir -p ref
	datasets download genome accession $(genome) --include genome,gff3 --filename ${genome}_genome.zip
	unzip -n -o ${genome}_genome.zip -d ref
	rm ${genome}_genome.zip
	mv ref/ncbi_dataset/data/$(genome)/*.fna ref/genome.fa
	mv ref/ncbi_dataset/data/$(genome)/*.gff ref/genome.gff
	echo "Genome saved to ref/genome.fa | GFF saved to ref/genome.gff"

# Download FULL FASTQ 
# Usage: make get_fastq fastq=SRR3194430
get_fastq:
	mkdir -p reads
	fasterq-dump $(fastq) -O reads/
	seqkit stats reads/$(fastq).fastq > reads/$(fastq)_stats.txt
	@echo ">> Running FastQC..."
	mkdir -p reads/fastqc_reports
	fastqc -o reads reads/$(fastq).fastq
	echo "FASTQ downloaded + QC run"

# Download PAIRED FASTQ
# Usage: make get_pairedfastq fastq=SRR3191544
get_pairedfastq:
	mkdir -p reads
	fasterq-dump --split-files $(fastq) -O reads/
	seqkit stats reads/$(fastq)_1.fastq reads/$(fastq)_2.fastq > reads/$(fastq)_stats.txt
	@echo ">> Running FastQC..."
	mkdir -p reads/fastqc_reports
	fastqc -o reads/fastqc_reports reads/$(fastq)_1.fastq reads/$(fastq)_2.fastq
	echo "Paired FASTQ downloaded + QC run"

# Index genome for BWA and samtools
# Usage: make index genome_fa=ref/genome.fa
index:
	@echo ">> Indexing genome with BWA..."
	bwa index -p ref/genome_index $(genome_fa)
	@echo " Genome indexing complete."

# Align FULL FASTQ to genome
# Usage: make align fastq=SRR3194430

alignsingle:
	@echo ">> Aligning reads to genome..."
	mkdir -p alignments
	bwa mem ref/genome_index reads/$(fastq).fastq | samtools sort -o alignments/$(fastq)_sorted.bam
	samtools index alignments/$(fastq)_sorted.bam
	samtools flagstat alignments/$(fastq)_sorted.bam > alignments/$(fastq)_align_stats.txt
	echo "BAM + index + stats ready in alignments/"

# Usage: make align fastq=SRR3191544
alignpaired:
	@echo ">> Aligning reads to genome..."
	mkdir -p alignments
	bwa mem ref/genome_index reads/$(fastq)_1.fastq reads/$(fastq)_2.fastq | samtools sort -o alignments/$(fastq)_sorted.bam
	samtools index alignments/$(fastq)_sorted.bam
	samtools flagstat alignments/$(fastq)_sorted.bam > alignments/$(fastq)_align_stats.txt
	echo "BAM + index + stats ready in alignments/"


# Convert BAM to BedGraph + BigWig for IGV/Genome Browser
# Usage: make bigwig fastq=SRR3194430 genome_fa=ref/genome.fa
# Usage: make bigwigfastq=SRR3191544 genome_fa=ref/genome.fa
bigwig:
	# Ensure fai index exists
	test -f $(genome_fa).fai || samtools faidx $(genome_fa)
	echo ">> Generating BedGraph..."
	bedtools genomecov -ibam alignments/$(fastq)_sorted.bam -bg | sort -k1,1 -k2,2n > alignments/$(fastq).bedgraph
	echo ">> Converting to BigWig..."
	bedGraphToBigWig alignments/$(fastq).bedgraph $(genome_fa).fai alignments/$(fastq).bw
	echo "BigWig generated: alignments/$(fastq).bw | Load into IGV!"

.PHONY: get_genome get_fastq get_pairedfastq index alignsingle alignpaired bigwig 
