
######## THE FOLLOWING FUNCTIONS ALLOW USER TO #######
# * Create a DESeq data object given a gse and filter out low gene counts
# * Create a variance stabilized representation of a given count matrix
# * Create a results object given a contrast, plot volcano, and save to file
# * Order significant DE gene results only by lfc or pval and save to file
######################################################

#############################
##### CREATE DGE LIST  ######
#############################
# NOTE: This code is not utilized in markdown!!!!
# (instead a function already exists to make DGE directly from gse)
create_dgl <- function(counts,metadata,condition) {
  ensembl_ids <- rownames(counts) # create gene annotations to add to list
  entrez_ids <- as.vector(create_ENS_Entrez_dict(ensembl_IDs = ensembl_ids))
  symbols <- as.vector(create_ENS_HGNC_dict(ensembl_IDs = ensembl_ids))
  genes <- data.frame(ENSEMBL = ensembl_ids, 
                      EntrezGene = entrez_ids, 
                      Symbol = symbols)
  dgl <- DGEList(counts = counts, 
                 genes = genes,
                 group = condition)
  cat(paste("Number of starting genes is: ", nrow(dgl), "\n", sep=""))
  return(dgl)
}

###########################################
##### FILTER DUPLICATE GENES IN DGL  ######
###########################################
# different refseq transcripts for same gene symbol predominantly count the same reads!!
# as a result, we choose to just keep one transcript for each symbol 
# (choosing the transcript with the highest overall count to represent)
filter_duplicates <- function(dgl) {
  num1 <- nrow(dgl) # original number of genes
  order_dgl <- order(rowSums(dgl$counts), decreasing=TRUE)
  dgl <- dgl[order_dgl,] # order dgl by highest counts
  
  dup <- duplicated(dgl$genes$Symbol)
  dgl <- dgl[!dup,] # filtered from 22520 to 22479 genes
  num2 <- nrow(dgl) # filtered number
  cat(paste("Filtered out ", (num1-num2), 
            " duplicated genes out for a new total of: ", num2, "\n", sep=""))
  return(dgl)
}

###########################################
##### FILTER OUT LOW COUNTS IN DGL  #######
###########################################
filter_lowcts <- function(dgl) {
  num1 <- nrow(dgl) # original number of genes
  keep <- filterByExpr(dgl) # filter out genes expressed at low levels!
  dgl <- dgl[keep, , keep.lib.sizes=FALSE] # apply filter to dgl
  dgl$samples$lib.size <- colSums(dgl$counts) # recompute library sizes after filtering
  rownames(dgl$counts) <- rownames(dgl$genes) <- dgl$genes$EntrezGene # use Entrez Gene IDs as row names
  #dgl$genes$EntrezGene <- NULL
  num2 <- nrow(dgl)
  cat(paste("Filtered out ", (num1-num2), " low count genes out for a new total of: ", num2, "\n", sep=""))
  return(dgl)
}

