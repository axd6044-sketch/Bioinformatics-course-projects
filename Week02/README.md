# Week 02: Genome Annotation Analysis (Anabas testudineus)

## 1.Go to Ensemble 
### http://ftp.ensembl.org/pub/current_gff3/

## 2. Organism chosen -
### anabas_testudineus

## 3.Download GFF file using command
```
wget http://ftp.ensembl.org/pub/current_gff3/anabas_testudineus/Anabas_testudineus.fAnaTes1.3.115.gff3.gz
```
## 4.About the organism
```
echo "Anabas testudineus is the scientific name for the climbing perch,
a freshwater fish belonging to the family Anabantidae.
Found widely in South and Southeast Asia "
```
## 5.How many sequence regions (chromosomes)? Does that match with the expectation for this organism?  
### Input
```
grep '^##sequence-region' Anabas_testudineus.fAnaTes1.3.115.gff3 | grep -v OOHO | cut -f 1 | wc -l
grep '^##sequence-region' Anabas_testudineus.fAnaTes1.3.115.gff3 | grep -v OOHO | cut -f 1 | #this will show the details of the chromosome and we can see chromosome number 20 is missing. And the total becomes 23.
```
### Output
```
23

# Interpretation :It matches with the expectation of the organism whose diploid chromosome number is 46. 
# While parsing the annotation file Anabas_testudineus.fAnaTes1.3.115.gff3, I noticed that the ##sequence-region metadata declares only 23 chromosomes:1–19, 21–24 There is no line for chromosome 20.
```
## 6.How many features does the file contain?
### Input
```
cat Anabas_testudineus.fAnaTes1.3.115.gff3 | grep -v '#' > Anabas.gff3
ls -lh
cat Anabas.gff3 | cut -f 3 | sort | uniq -c 
cat Anabas.gff3 | cut -f 3 | sort | uniq | wc -l
```
### Output 
```
total 647M
-rw-rw-r-- 1 axd6044 axd6044_collab 323M Sep  7 14:27 Anabas.gff3
-rw-rw-r-- 1 axd6044 axd6044_collab 324M Jul 10 19:02 Anabas_testudineus.fAnaTes1.3.115.gff3
-rw-rw-r-- 1 axd6044 axd6044_collab 1.2K Sep  7 09:52 README.md

#The features of this file
 914601 CDS 
 938673 exon
  51602 five_prime_UTR
  23996 gene
   1832 lnc_RNA
     97 miRNA
  64517 mRNA
   1904 ncRNA_gene
    137 pseudogene
    137 pseudogenic_transcript
     50 region
    138 rRNA
      7 scRNA
    203 snoNA
    207 snRNA
  43764 three_prime_UTR
     10 transcript
      4 V_gene_segment
      1 Y_RNogenic_transcript

#Counting the number of features
19

```
## 7.How many genes are listed for this organism?
### Input
```
grep -v '^#' Anabas.gff3 | awk '$3=="gene"' | wc -l
```
### Output 
```
23996
```
## 8. Is there a feature type that you may have not heard about before? What is the feature and how is it defined?  
```
#Feature- Pseudogene
“A sequence that closely resembles a known functional gene, but is generally non-functional as a consequence of mutations that prevent its expression or translation.”
How they arise:
Duplication — a gene is copied, but one copy accumulates mutations that disable it.
Retrotransposition — an mRNA is reverse-transcribed back into DNA and inserted into the genome, but lacks proper regulatory elements.
Unitary loss — the only copy of a gene in the genome becomes nonfunctional through mutations.
In annotations (like this  GFF3):
A pseudogene feature marks the genomic coordinates of such a gene. Often, you’ll also see pseudogenic_transcript entries beneath them, representing the “ghost” transcripts they would have produced.
```
## 9. What are the top-ten most annotated feature types (column 3) across the genome?
### Input
```
cat Anabas.gff3 | cut -f 3 | sort | uniq -c | sort -rn
cat Anabas.gff3 | cut -f3 | sort | uniq -c | sort -nr | head -10
```
### Output 
```
#Order of annotated features
938673 exon
914601 CDS
  64517 mRNA
  51602 five_prime_UTR
  43764 three_prime_UTR
  23996 gene
   1904 ncRNA_gene
   1832 lnc_RNA
    207 snRNA
    203 snoRNA
    138 rRNA
    137 pseudogenic_transcript
    137 pseudogene
     97 miRNA
     50 region
     10 transcript
      7 scRNA
      4 V_gene_segment
      1 Y_RNA

#TOP-10 annotated features
938673 exon
 914601 CDS
  64517 mRNA
  51602 five_prime_UTR
  43764 three_prime_UTR
  23996 gene
   1904 ncRNA_gene
   1832 lnc_RNA
    207 snRNA
    203 snoRNA
```
## 10. Having analyzed this GFF file, does it seem like a complete and well-annotated organism?
```
The GFF3 file appears reasonably comprehensive, with ~24,000 genes, ~64,000 mRNAs, and nearly one million exons and CDS entries. The presence of noncoding RNA categories (lnc_RNA, miRNA, snoRNA, snRNA, rRNA) and pseudogenes indicates that the annotation pipeline attempted to capture a broad range of genomic features, not just protein-coding genes.  It does seem like a complete and well-annotated organism. 
```

## 11.Other insights
```
Only a handful of miRNAs (97), snoRNAs (203), and snRNAs (207) are annotated, while most vertebrates harbor hundreds to thousands.
The missing chromosome 20 is a red flag; downstream tools may assume consecutive numbering, so analyses should be cautious and perhaps use sequence IDs directly rather than numeric assumptions.
Not all transcripts have both 5′ and 3′ UTRs, which may affect transcript structure analysis and regulatory motif searches.
```

