
######## THE FOLLOWING FUNCTIONS ALLOW USER TO #######
# * Create a DESeq data object given a gse and filter out low gene counts
# * Create a variance stabilized representation of a given count matrix
# * Create a results object given a contrast, plot volcano, and save to file
# * Order significant DE gene results only by lfc or pval and save to file
######################################################

#############################################
##### CREATE DESEQ DATASET AND FILTER  ######
#############################################
# desig should be c() vector with list of design elements (could be just one), eg c("macaque", "condition")
# hist is whether or not to plot a histogram of the frequency of count totals for each gene
create_dds <- function(gse, desig, hist=FALSE) { 
  str <- paste('~', desig[1], sep = ' ')
  if (length(desig)!=1) {
    for (element in desig) {
      str <- paste(str, element, sep=" + ")
    }
  } # note: "~ var1 + var2" measures the effect of var2 while controlling for var1
  dds <- DESeqDataSet(gse, design = as.formula(str))
  counts_per_gene <- rowSums(counts(dds)) # for every gene, get the sum of counts for all samples
  # note: counts(dds) is the same thing as using assays for gse
  if (hist==TRUE) {
    # Plot histogram to visualize counts per gene
    hist(log(counts_per_gene,base=10), xlab='Log counts per gene', 
         ylab='Frequency', main='Count Frequency Before Filtering')
  }
  n_gene_initial <- nrow(dds) # starting number of total genes
  # get the count that represents the bottom 5% of all count sums
  count_threshold <- quantile(counts_per_gene, 0.05) # quantile in distribution of log counts per genes
  if (count_threshold==0) {
    count_threshold==1
  }
  # Filter on number of genes counted more than the quantile above
  keep <- rowSums(counts(dds)) > count_threshold
  dds <- dds[keep, ]
  n_gene_final <- nrow(dds) # final number of genes after being filtered
  cat("Filtered from initial number of genes:", n_gene_initial, "to", n_gene_final, 'genes.\n')
  return(dds)
}

#############################################
# VARIANCE STABILIZING TRANSFORMATION (vst) #
#############################################
## - Calculated from fitted dispersion mean relation (in DESeq object I believe)
## - Transforms normalized count data to yeild homoskedastic values
## - (homoskedastic refers to having the same variance at different mean values)
## - (in RNA seq data, variance grows with the mean which must be corrected for)
var_transform <- function(dds, plotFit=FALSE) {
  vsd <- DESeq2::vst(dds, blind = TRUE) # easier to compute but less sensitive than rlog to high count outliers
  # blind determines whether to blind the transformation to the experimental design
  # TRUE used for comparing samples in unbiased manner
  # blind=FALSE used for transforming data for downstream analysis where full design info is made
  # (skips re-estimation of dispersion trend if already calculated)
  # (do this if many genes have large differences in counts due to experimental design)
  if (plotFit==TRUE) {
    meanSdPlot(vsd@assays@data[[1]], ranks=FALSE) # plot a mean vs variance graph of genes
  }
  # Extract dictionary between ENSEMBL IDs and HGNC symbols
  ensembl_ids <- vsd@rowRanges@elementMetadata@listData$Ensembl
  vsd_symbols <- vsd@rowRanges@elementMetadata@listData$Symbol
  vsd_symbols[is.na(vsd_symbols)] <- ensembl_ids[is.na(vsd_symbols)]
  rownames(vsd) <- vsd_symbols # rename rows to HGNC symbols
  return(vsd) 
}

#############################################
############ COMPUTE DGE RESULTS ############
#############################################
compute_dge <- function(dds, contr=c("genotype", "HT", "WT"), a=0.1, 
                        titl="DGE", co=2,dir=file.path(results_dir, 'DGE'),  
                        savetoFile=FALSE, showSummary=FALSE) { 
  # contr: a c() beginning with contrast condition name, then independent var, then control
  # alpha: false discovery rate threshold (set to default of 0.1)
  # co: # of measurements in your design (eg ~genotype is co=2 while ~batch+genotype is co=3)
  
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast=contr, alpha=a)
  res <- lfcShrink(dds, coef=co, res=res, type = 'apeglm')
  # steps for above deseq2 function:
  # (1) estimates size factors (controlling for differences in sequencing depth of samples)
  # (2) estimation of dispersion values for each gene
  # (3) fitting a generalized linear model
  
  #mcols command below gives definition of each column in res...
  #mcols(res, use.names = TRUE)
  ensembl_IDs <- rownames(res) # Create dictionary from ENSEMBL genes to gene symbols
  if(all(grepl('ENSMMUG', ensembl_IDs))) { # only correct the gene IDs if they are ENSEMBL IDs
    symbols <- dds@rowRanges@elementMetadata$Symbol
    res$ensembl_ID <- ensembl_IDs
    res$symbols <- symbols
    rownames(res) <- symbols
  } else {
    cat('Gene names are already gene symbols.  Nothing to do!')
  }
  # The annotation code below doesn't work (still trying to debug)
  # Query ENSEMBL database for gene annotations
  #ann_df <- add_anns(ensembl_ids = ensembl_IDs)
  #anns <- ann_df$description; names(anns) <- ann_df$ensembl_gene_id
  res$descriptions <- dds@rowRanges@elementMetadata$Desc
  
  p <- EnhancedVolcano(res, lab = res$symbols, x = 'log2FoldChange', y = 'pvalue')
  plot(p)
  if (showSummary==TRUE) {
    summary(res)
  }
  if (savetoFile==TRUE) {
    if(!dir.exists(dir)) {
      dir.create(dir)
    }
    titl <- paste(titl, ".txt", sep="")
    write.table(res, file = file.path(dir, titl),append = FALSE, 
                quote = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE)
  }
  cat("Considering a faction of 10% False Positives as acceptable, number of genes with an Adjusted P-Value below 10% is: ")
  cat(sum(res$padj < 0.1, na.rm=TRUE))
  cat("\n")
  return(res)
}

#############################################
######### SAVE SIGNIFICANT DGE ONLY #########
#############################################
compute_sig <- function(res, cutof=0.1, order='log2FoldChange', 
                        titl="significant", savetoFile=FALSE, 
                        dir=file.path(results_dir, "DGE")) {
  # order MUST be either log2FoldChange or padj!!!!
  sig <- subset(res, padj<cutof)
  if (order=='log2FoldChange') {
    sig <- sig[order(sig$log2FoldChange), ]
  } else {
    sig <- sig[order(sig$padj), ]
  }
  rownames(sig) <- sig$symbols
  
  sig <- subset(sig, !grepl("^ENSMMUG0", sig$symbols))
  sig <- subset(sig, !grepl("^LOC", sig$symbols))
  
  if (savetoFile==TRUE) {
    if(!dir.exists(dir)) {
      dir.create(dir)
    }
    titl <- paste(titl, ".txt", sep="")
    write.table(sig, file = file.path(dir, titl),append = FALSE, 
                quote = FALSE, sep = '\t', col.names = TRUE, row.names = TRUE)
  }
  return(sig)
}