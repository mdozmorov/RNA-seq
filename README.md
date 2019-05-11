# RNA-seq analysis pipeline, featureCounts input, edgeR, functional enrichment, visualization

- `Analysis.Rmd` - main RNA-seq analysis pipeline, for `featureCount` counts. Prerequisites:
    - A path to data folder. This folder should have 3 subfolders:
        - `03_featureCount` - gzipped count files outputted by `featureCount`
        - `results` - folder where the results will be stored
        - `data` - Must have `sample_annotation.csv` file. Annotation file should have "Sample" column with sample names, and any other annotation columns. Include "Group" column containing covariate of interest. Example:
```
# Sample,Group
VLI10_AA_S61_L006_R1_001.txt.gz,AA
VLI10_AA_S61_L007_R1_001.txt.gz,AA
VLI10_AA_S61_L008_R1_001.txt.gz,AA
VLI11_C_S62_L006_R1_001.txt.gz,C
VLI11_C_S62_L007_R1_001.txt.gz,C
VLI11_C_S62_L008_R1_001.txt.gz,C
```

- `calcTPM.R` - a function to calculate TPMs from gene counts

- `Figure_clusterProfiler_nes.Rmd` - Takes the results of edgeR analysis from an Excel file, performs GO and KEGG GSEA and plots the results as horizontal barplots, sorted by normalized enrichment score (NES). [Example](Figure_clusterProfiler.pdf)

- `Figure_clusterProfiler_asis.Rmd` - Takes the results of edgeR analysis from an Excel file, performs GO and KEGG GSEA and plots the results as horizontal barplots, sorted by p-value, as they come out of the enrichment analysis.

- `utils.R` - helper functions

## Misc

- DESeq results to pathways in 60 Seconds with the fgsea package, https://stephenturner.github.io/deseq-to-fgsea/