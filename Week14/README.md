## Week 14 assignment 

Exploratory and Functional Analysis of RNA-seq Differential Expression Results
Abstract

Exploratory visualization and functional enrichment are essential for interpreting RNA-seq differential expression results. Here, we analyzed significantly differentially expressed genes derived from an edgeR output to assess sample relationships and biological function. Principal component analysis (PCA) and heatmap visualization were used to characterize transcriptional structure, followed by functional enrichment analysis using g:Profiler and Enrichr. All analyses were automated using a Makefile to ensure reproducibility.

Introduction

RNA-seq experiments generate high-dimensional datasets that require dimensionality reduction and functional interpretation to uncover biologically meaningful patterns. PCA and heatmaps provide complementary views of global expression variation and gene-level structure, while enrichment analyses link differentially expressed genes to known biological pathways. In this study, we applied a reproducible Biostar-based workflow to explore differential expression results from a human RNA-seq dataset.

Methods
Input Data

edger.csv: Differential expression results containing Ensembl gene identifiers and FDR values

design.csv: Experimental design file specifying sample names and group labels (HBR and UHR)

Genes with FDR < 0.05 were retained for downstream analyses.

Software Environment

All analyses were performed in the Biostar bioinfo environment:

micromamba activate bioinfo

Workflow Automation

All steps were executed using a Makefile to ensure reproducibility:

make all


This command generates:

PCA plot

Heatmap

g:Profiler enrichment results

Enrichr enrichment results

Principal Component Analysis

PCA was performed using the top 500 most variable genes.

Command:
```bash
src/r/plot_pca.r -c edger.csv -d design.csv -o pca.pdf

```
This script uses DESeq2-based variance estimation and ggplot2 for visualization.

Heatmap Visualization

A heatmap of significantly differentially expressed genes was generated using scaled expression values.

Command:
```bash
src/r/plot_heatmap.r -c edger.csv -d design.csv -o heatmap.pdf
```

Hierarchical clustering was applied to both genes and samples.

Functional Enrichment Analysis
g:Profiler

Significantly differentially expressed genes were submitted to g:Profiler for functional enrichment analysis.
```bash
bio gprofiler -c edger.csv -d hsapiens
```

Output:

gprofiler.csv

Enrichr

A complementary enrichment analysis was performed using Enrichr.
```
bio enrichr -c edger.csv
```

Output:

enrichr.csv

Results
Principal Component Analysis

PCA revealed strong separation between experimental groups (Figure 1). The first principal component (PC1) explained 97% of the total variance, clearly distinguishing HBR and UHR samples. PC2 explained 1% of the variance, capturing minor within-group variation. Tight clustering of biological replicates indicates high reproducibility.

Figure 1. Principal component analysis of RNA-seq samples showing clear separation between HBR and UHR groups.

Heatmap Visualization

Heatmap analysis revealed structured gene expression patterns that segregated by experimental condition (Figure 2). Distinct gene clusters showed reciprocal expression between HBR and UHR samples, consistent with the PCA results.

Figure 2. Heatmap of significantly differentially expressed genes (FDR < 0.05). Red indicates higher relative expression and green indicates lower relative expression.

Functional Enrichment

A total of 288 genes met the FDR < 0.05 threshold.

g:Profiler identified 352 enriched functional terms

Enrichr identified 95 enriched categories

Enriched pathways were consistent with known biological processes relevant to the experimental conditions, providing functional context for the observed transcriptional differences.

Discussion

The agreement between PCA clustering, heatmap structure, and enrichment results demonstrates the robustness of the differential expression analysis. The dominance of PC1 indicates that experimental condition is the primary driver of transcriptional variation. Using two independent enrichment tools increased confidence in biological interpretation, while Makefile automation ensured reproducibility.

Conclusion

This study presents a reproducible workflow for exploratory and functional analysis of RNA-seq differential expression results. PCA and heatmap analyses revealed clear group-specific transcriptional patterns, and enrichment analyses connected these patterns to biologically meaningful pathways. The Makefile-based framework provides a scalable approach for RNA-seq data interpretation.
