##################
# LOAD PACKAGES ##
##################

library('gplots')
library('ggplot2')
library('DESeq2')
library('limma')
library('reshape2')
library('RColorBrewer')
library('WGCNA')
library('edgeR')
library("vsn")
library("fgsea")
library("msigdbr")
library("org.Mmu.eg.db")
library("pathview")
library("grid")
library("png")

##################
### USER PATHS ###
##################
# NOTE: Please change the paths below to match the given users local environment 
# base directory
local_dir <- '/Users/gricemz/Desktop/RNA_CFU/FINAL_VERSION' 
# folder for json, cdna, index, gtf and other transcript info (copied from biowulf or found online)
transcript_dir <- '/Users/gricemz/Documents/transcripts' 
# directory containing raw quant folders (output from salmon)
quant_dir <- '/Users/gricemz/Desktop/RNA_CFU/quant' 
results_dir <- file.path(local_dir, 'results') 
if(!dir.exists(results_dir)) {
  dir.create(results_dir)
}

rds_dir <- file.path(local_dir, 'rds')
#meta_file <- file.path(local_dir, 'metadata/metadata.txt')
broad_geneset_dir <- '/Users/gricemz/Documents/NHLBI/reference_data/broad_genesets' # broad gene sets
species_dataset <- "mmulatta_gene_ensembl" # ENSEMBL dataset for this species




