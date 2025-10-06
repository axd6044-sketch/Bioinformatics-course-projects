# This week we will make a make file for running script and align sequence data to genome and finally visualize the BAM file in IGV

## First:  Making a makefile 
### the makefile will do the following :
downloads the reference genome, indexes it, downloads raw reads from the NCBI SRA database, performs read subsampling, aligns reads to the genome, and computes coverage statistics.

### Things to do before calling makefile 
```
micromamba activate bioinfo
```
Prepare the makefile and run all 
```
make all
```

## Stats 
```
samtools flagstat SRR3194430.bam 
```

76300000 + 0 in total (QC-passed reads + QC-failed reads)
76299868 + 0 primary
0 + 0 secondary
132 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
202364 + 0 mapped (0.27% : N/A)
202232 + 0 primary mapped (0.27% : N/A)
0 + 0 paired in sequencing
0 + 0 read1
0 + 0 read2
0 + 0 properly paired (N/A : N/A)
0 + 0 with itself and mate mapped
0 + 0 singletons (N/A : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

### What percentage of reads aligned to the genome?
0.27% which is extremely less 
###What was the expected average coverage?
N = number of reads (expected 1500 subsampled reads, or 76,300,000 full reads)
L = average read length (~75 bp for SRR3194430)
G = genome length (~10,794 bp for Zika MR766)

### What is the observed average coverage?
X+ observed expected​=0.0027×530,000≈1,431× . So the expected average coverage for the mapped fraction is about 1,400×.
How much does the coverage vary across the genome? (Provide a visual estimate
```
samtools depth SRR3194430.bam | awk '{sum+=$3} END { print "Average coverage:", sum/NR }'
```
Average coverage: 1409.39


