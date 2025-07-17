# Helper functions


# Function to compare groups
#' @param object glmLTR result
#' @param p_adj_cutoff FDR adjusted p-value cutoff
#' @param subset_by logical vectore defined earlier and used to subset the whole dataset
#' @param comparison text defining what comparison has been made
#' @param nplot how many top differential genes to plot
p.vals <- function(object, p_adj_cutoff = 0.1, subset_by = subset_by, comparison = "DEGs", nplot = 50, print_to_file = TRUE, merge_by = "symbol") {
  
  genes.ABC <- topTags(object, n = Inf, p.value = p_adj_cutoff) %>% as.data.frame # Get all significant genes
  # Proceed only if some results returned
  if (nrow(genes.ABC) > 0) {
    if (merge_by == "symbol") {
      genes.ABC <- left_join(data.frame(genes.ABC, stringsAsFactors = FALSE), 
                             gene_annotations, by = c("genes" = "symbol")) # Attach annotations
    }
    if (merge_by == "ensgene") {
      genes.ABC <- left_join(data.frame(genes.ABC, stringsAsFactors = FALSE), 
                             gene_annotations, by = c("genes" = "ensgene")) # Attach annotations
    }
    # Sort by LR
    genes.ABC <- genes.ABC[order(genes.ABC$LR, decreasing = TRUE), ]
    
    if (print_to_file) {
      # Save most significant genes
      save_res(genes.ABC, fileName_rna, wb = wb, sheetName = comparison) # Save gene list
    }
  }
  # Return actual gene list
  return(genes.ABC)
}


# Function to plot and save selected gene lists
genes_to_heatmap <- function(object, edgeR.dgelist = edgeR.dgelist, nplot = 50, subset_by = c("A", "B"), comparison = "C231_noCRISPR_vs_Negative", clustmethod = "ward", width = 6, height = 8, print_to_file = TRUE, merge_by = "symbol") {
  # Select columns to subset
  if (any(!is.na(subset_by))){
    subset_index <- grepl(paste(subset_by, collapse = "|"), sample_annotation_subset$Group) 
  } else {
    subset_index <- 1:ncol(cpm(edgeR.dgelist, log = TRUE))
  }
  # # Order object by the most on-average significant genes
  # object <- object[ order( (rank(object$FDR.x) + rank(object$FDR.y)) / 2), ] 
  # # Save most significant genes
  # fileName <- paste0("results/DEGs_", comparison,  "_hg19.xlsx") # Construct filename
  # unlink(fileName) # Delete file, if exists
  # write.xlsx2(object, fileName, sheetName = comparison, row.names = FALSE) # Save gene list
  
  # Plot top genes
  genes.to.plot <- object[1:nplot, c("genes"), drop = FALSE] 
  if (merge_by == "symbol") {
    genes.to.plot <- left_join(genes.to.plot, gene_annotations, by = c("genes" = "symbol"))
  }
  if (merge_by == "ensgene") {
    genes.to.plot <- left_join(genes.to.plot, gene_annotations, by = c("genes" = "ensgene"))
  }
  genes.to.plot <- genes.to.plot[match(edgeR.dgelist@.Data[[3]]$genes[edgeR.dgelist@.Data[[3]]$genes %in% genes.to.plot$genes], genes.to.plot$genes), , drop = FALSE] # Make the same order as the genes subsetted from the whole dataset
  matrix.to.plot <- cpm(edgeR.dgelist, log = TRUE, normalized.lib.sizes=TRUE)[edgeR.dgelist@.Data[[3]]$genes %in% genes.to.plot$genes, subset_index] # * edgeR.dgelist@.Data[[2]]$norm.factors[subset_index] # Subset matrix to plot
  if (merge_by == "symbol") {
    rownames(matrix.to.plot) <- genes.to.plot$genes
  }
  if (merge_by == "ensgene") {
    rownames(matrix.to.plot) <- genes.to.plot$symbol
  }
  
  annotation_col <- data.frame(Group = sample_annotation_subset$Group[subset_index])
  rownames(annotation_col) <- colnames(cpm(edgeR.dgelist)[, subset_index])
  
  # If too many genes, do not plot gene names
  if (nplot > 50) {
    genes.to.plot$symbol <- ""
  }
  
  # Save to PDF
  if (print_to_file) {
    pdf(paste0("results/Figure_Heatmap_", comparison, "_", nplot, ".pdf"), width = width, height = height)
    pheatmap(matrix.to.plot, color = col3, clustering_method = "ward", treeheight_row = FALSE, treeheight_col = FALSE, annotation_col = annotation_col, labels_row = genes.to.plot$symbol, scale = "row")
    dev.off()
  } else {
    pheatmap(matrix.to.plot, color = col3, clustering_method = "ward", treeheight_row = FALSE, treeheight_col = FALSE, annotation_col = annotation_col, labels_row = genes.to.plot$symbol, scale = "row")
  }
}



# A function to plot a barplot of selected genes in selected conditions
genes_to_boxplot <- function(selected_genes = c("NF1", "SND1", "MTDH"), subset_by = c("C", "AA"), merge_by = "symbol") {
  # Select columns to subset
  # subset_index <- grepl(paste(subset_by, collapse = "|"), sample_annotation$Group) 
  subset_index <- sample_annotation_subset$Group == subset_by[1] | sample_annotation_subset$Group == subset_by[2]
  # Select EntrezIDs to plot
  if (merge_by == "symbol") {
    genes.to.plot <- gene_annotations[gene_annotations$symbol %in% selected_genes, ]
    genes.to.plot <- genes.to.plot[match(edgeR.dgelist@.Data[[3]]$genes[edgeR.dgelist@.Data[[3]]$genes %in% genes.to.plot$symbol], genes.to.plot$symbol), ] # Make the same order as the genes subsetted from the whole dataset
    index_to_keep  <- edgeR.dgelist@.Data[[3]]$genes %in% genes.to.plot$symbol
  }
  if (merge_by == "ensgene") {
    genes.to.plot <- gene_annotations[gene_annotations$ensgene %in% selected_genes, ]
    genes.to.plot <- genes.to.plot[match(edgeR.dgelist@.Data[[3]]$genes[edgeR.dgelist@.Data[[3]]$genes %in% genes.to.plot$ensgene], genes.to.plot$ensgene), ] # Make the same order as the genes subsetted from the whole dataset
    index_to_keep  <- edgeR.dgelist@.Data[[3]]$genes %in% genes.to.plot$ensgene
  }
  # Expression of genes to plot
  matrix.to.plot <- cpm(edgeR.dgelist, log = TRUE)
  matrix.to.plot <- matrix.to.plot[index_to_keep, subset_index, drop = FALSE] 
  rownames(matrix.to.plot) <- edgeR.dgelist@.Data[[3]]$genes[index_to_keep]
  colnames(matrix.to.plot) <- sample_annotation_subset$Group[subset_index]
  # Add genes
  matrix.to.plot <- data.frame(Gene = rownames(matrix.to.plot), matrix.to.plot, check.names = FALSE)
 # Make long format
  matrix.to.plot_melted <- tidyr::pivot_longer(matrix.to.plot, cols = colnames(matrix.to.plot)[colnames(matrix.to.plot) != "Gene"])
  matrix.to.plot_melted$Gene <- as.factor(matrix.to.plot_melted$Gene) # factor(matrix.to.plot_melted$Gene, levels = unique(matrix.to.plot_melted$Gene))
  matrix.to.plot_melted$name <- as.factor(matrix.to.plot_melted$name) # factor(matrix.to.plot_melted$variable %>% as.character(), levels = unique(matrix.to.plot_melted$variable %>% as.character()))
  # Plot
  ggplot(matrix.to.plot_melted, aes(x = name, y = value, group = name)) +
    geom_boxplot(aes(fill = name)) +
      facet_wrap(~ Gene, ncol = 3, scales = "free_y")
}


## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

