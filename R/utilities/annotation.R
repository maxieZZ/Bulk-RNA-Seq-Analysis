
######## THE FOLLOWING FUNCTIONS ALLOW USER TO #######
# 1) Get associated gene descriptions given a list of ENSEMBL IDs 
# (saves output to an RDS file to be accessed in later attempts)
# 2) Convert to gene symbols given a list of ENSEMBL IDs
# 3) Convert to ENTREZ IDs given a list of ENSEMBL IDs
######################################################

######################################
# Add Gene Descriptions from ENSEMBL #
######################################
add_anns <- function(ensembl_ids, force_remote_query = FALSE) {
  gene_annot_file <- file.path(rds_dir, 'gene_annotations.RDS')
  if(!file.exists(gene_annot_file) || force_remote_query) { # determine whether necessary to perform remote query
    ensembl_dataset <- useDataset(dataset = 'mmulatta_gene_ensembl', # specify ENSEMBL data-set (assumes macaque is species)
                                  mart = useMart('ENSEMBL_MART_ENSEMBL', host = 'www.ensembl.org'))
    anns <- getBM(
      attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'),
      filters = "ensembl_gene_id",
      values = ensembl_ids,
      mart = ensembl_dataset) # query to get descriptions
    clean.descriptions <- as.vector(sapply(anns$description, function(d){ # clean up description output
      strsplit(d, " \\[Source:")[[1]][1]
    }))
    anns$description <- clean.descriptions
    anns <- anns[match(ensembl_ids, anns[, 1]), ] # match description annotation with results
    rownames(anns) <- ensembl_ids # set the row names to the ENSEMBL ID
    save(anns, file = gene_annot_file) # save annotations in a file in the rds directory to be accessed in previous tries
  } else {
    load(file = gene_annot_file) # if file was already saved, we don't need to query again!
  }
  invisible(anns)
}

#########################################################
##### Convert given ENSEMBL list to HGNC symbol list ####
#########################################################
create_ENS_HGNC_dict <- function(ensembl_IDs,manual_dict_fn = NULL) {
  symbols <- mapIds(org.Mmu.eg.db,
                    keys = ensembl_IDs,
                    column = 'SYMBOL',
                    keytype = 'ENSEMBL',
                    multiVals = 'first')
  missing <- which(is.na(symbols)) # set missing symbols to corresponding ENSEMBL ID
  symbols[missing] <- ensembl_IDs[missing]
  if(!is.null(manual_dict_fn)) { # check if there is dictionary of manual translations
    # can be added in the future
  }
  invisible(symbols)
}

#############################################
#### Convert ENSEMBL IDs to Entrez IDs ######
#############################################
create_ENS_Entrez_dict <- function(ensembl_IDs) {
  entrez_ids <- mapIds(org.Mmu.eg.db,
                       keys = ensembl_IDs,
                       column = 'ENTREZID',
                       keytype = 'ENSEMBL',
                       multiVals = 'first')
  # set missing symbols to corresponding ENSEMBL ID
  missing <- which(is.na(entrez_ids))
  entrez_ids[missing] <- ensembl_IDs[missing]
  invisible(entrez_ids)
}
  