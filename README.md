# RNA-seq analysis pipeline, featureCounts input, edgeR, functional enrichment, visualization

- `Figure_clusterProfiler_nes.Rmd` - Takes the results of edgeR analysis from an Excel file, performs GO and KEGG GSEA and plots the results as horizontal barplots, sorted by normalized enrichment score (NES). [Example](Figure_clusterProfiler.pdf)

- `Figure_clusterProfiler_asis.Rmd` - Takes the results of edgeR analysis from an Excel file, performs GO and KEGG GSEA and plots the results as horizontal barplots, sorted by p-value, as they come out of the enrichment analysis.


## Misc

- DESeq results to pathways in 60 Seconds with the fgsea package, https://stephenturner.github.io/deseq-to-fgsea/