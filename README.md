# RNA-seq analysis pipeline, STAR input, edgeR, functional enrichment, visualization

- [RNAseq-workflow](https://github.com/twbattaglia/RNAseq-workflow) - A repository for setting up a RNAseq workflow. Detailed instructions and code for each analysis and visualization step.

- [Analysis_STAR.Rmd](Analysis_STAR.Rmd) - RNA-seq analysis pipeline for `STAR` counts. Prerequisites:
    - A path to data folder. This folder should have 3 subfolders:
        - `02_STAR-align` - gzipped count files with `.tab` extension outputted by `STAR` aligner
        - `results` - folder where the results will be stored
        - `data` - Must have `sample_annotation.csv` file, example below

- [GSEA.Rmd](GSEA.Rmd) - EnrichR (non-directional) and GSEA (directional) analysis using KEGG, GO, MSigDb.

- [GSEA_figures.Rmd](GSEA_figures.Rmd) - Visualization of GSEA enrichment results as horizontal barplots.

- [Figure_heatmap.Rmd](Figure_heatmap.Rmd) - make heatmap of top 50 differentially expressed genes. Uses `TMP.xlsx` produced by `Analysis*.Rmd`. May use a custom signature of genes. Includes [EnhancedVolcano](https://bioconductor.org/packages/EnhancedVolcano/) and boxplots of selected genes.

- [oncoEnrichR.Rmd](oncoEnrichR.Rmd) - Cancer-dedicated gene set interpretation using the [oncoEnrichR](https://github.com/sigven/oncoEnrichR) R package

- [Pathview.Rmd](Pathview.Rmd) - visualization of top KEGG pathways. Uses `DEGs.xlsx` produced by `Analysis*.Rmd`. [Example](examples/pathview_example.pdf)

- [calcTPM.R](calcTPM.R) - a function to calculate TPMs from gene counts

- [utils.R](utils.R) - helper functions

### [misc](misc) - Outdated scripts

- `Analysis_featurecounts.Rmd` - RNA-seq analysis pipeline for `featureCount` counts. Prerequisites:
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

- `Figure_clusterProfiler_nes.Rmd` - Takes the results of edgeR analysis from an Excel file, performs GO and KEGG GSEA and plots the results as horizontal barplots, sorted by normalized enrichment score (NES). [Example](examples/Figure_clusterProfiler.pdf)

- `Figure_clusterProfiler_asis.Rmd` - Takes the results of edgeR analysis from an Excel file, performs GO and KEGG GSEA and plots the results as horizontal barplots, sorted by p-value, as they come out of the enrichment analysis.

- [enrichR_analysis.Rmd](enrichR_analysis.Rmd) - Analyze gene lists using [enrichR](https://cran.r-project.org/package=enrichR). Analyze all genes, and up- and downregulated genes separately. Uses `DEGs.xlsx` produced by `Analysis*.Rmd`.

- [enrichR_plot.Rmd](enrichR_plot.Rmd) - barplot of selected enrichment results, similar to [Example](examples/Figure_clusterProfiler.pdf). WIP

### [scripts](scripts)

Scripts for running RNA-seq preprocessing steps on a cluster using PBS job submission system. `subread-featurecounts` scripts are in the [dcaf/ngs.rna-seq](https://github.com/mdozmorov/dcaf/tree/master/ngs.rna-seq) repository

- [submit00_fastqc.sh](scripts/submit00_fastqc.sh) - FASTQC on raw FASTQ files
- [MultiQC](https://multiqc.info/) commands to summarize QC reports generated by TrimGalore and STAR
```
multiqc --filename multiqc_01_trimmed.html --outdir multiqc_01_trimmed 01_trimmed/
multiqc --filename multiqc_02_STAR-align.html --outdir multiqc_02_STAR-align 02_STAR-align/
```
- [submit01_trimgalore.sh](scripts/submit01_trimgalore.sh) - Adapter trimming using [TrimGalore](https://github.com/FelixKrueger/TrimGalore)
- [submit02_STAR-index.sh](scripts/submit02_STAR-index.sh) - Index the genome for the [STAR](https://github.com/alexdobin/STAR) aligner
- [submit02_STAR.sh](scripts/submit02_STAR.sh) - Align samples using [STAR](https://github.com/alexdobin/STAR). Requires `input01_toStarAlign.list` text file with the list of input files, each string contains (comma-separated) file name(s), space separates first and second read pairs

### [CaSpER](https://github.com/akdess/CaSpER) pipeline detecting CNVs from RNA-seq data

Dedicated repository with detailed instructions: [mdozmorov/CaSpER_pipeline](https://github.com/mdozmorov/CaSpER_pipeline)

- [submit05_BAFExtract-index.sh](scripts/submit05_BAFExtract-index.sh) - indexing the genome for [BAFExtract](https://github.com/akdess/BAFExtract) 
- [submit05_BAFExtract.sh](scripts/submit05_BAFExtract.sh) - [BAFExtract](https://github.com/akdess/BAFExtract) run

## Misc

- DESeq results to pathways in 60 Seconds with the fgsea package, https://stephenturner.github.io/deseq-to-fgsea/

- [A Shiny app for visualizing DESeq2 results](https://jokergoo.github.io/InteractiveComplexHeatmap/articles/deseq2_app.html) by Zuguang Gu. [Tweet](https://twitter.com/jokergoo_gu/status/1373013068821725189?s=20)
