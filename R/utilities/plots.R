
######## THE FOLLOWING FUNCTIONS ALLOW USER TO #######
# * Plot a euclidean distance map given variance stabilized count data
# * Plot a heat map given a list of custom genes and a count matrix
# * Plot top two principle components (PCA) with designated metadata grouping
# * Plot a volcano plot and save to file given a results object and annotation
# * Compare two gene lists, then print and save output of intersection
# * Plot the top results of undirected go analysis on barchart
######################################################

#############################################
## PLOT EUCLIDEAN DISTANCE AND SAVE TO PDF ##
#############################################
# dir=the directory the pdf is saved to (IF savetoPDF=TRUE)
# vsd=variance stabilized count matrix & titl=name of pdf file
# rownames=the labels you would like in the plot you create
euclidean_plot <- function(vsd, dir=file.path(local_dir, 'results'), 
                          titl, rownames=vsd@colData@listData$names, # normally where metadata is stored
                          savetoPDF=FALSE) { # titl is what the pdf file should be called
  
  sampleDists <- dist(t(assay(vsd))) # Compute Euclidean Distances
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- rownames
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  p <- pheatmap(sampleDistMatrix,
                clustering_distance_rows = sampleDists,
                clustering_distance_cols = sampleDists,
                col = colors)
  if(savetoPDF==TRUE) {
    if(!dir.exists(dir)) {
      dir.create(dir)
    }
    fname <- paste(dir,"/",titl, ".pdf", sep="")
    pdf(fname, width=7, height=7)
    grid::grid.newpage()
    grid::grid.draw(p$gtable)
    dev.off()
  }
}

#############################################
### DISPLAY AND SAVE PRINCIPLE COMPONENTS ###
#############################################
# num_catagories tells the function how many metadata variables its showing on the plot
# eg when =2 it will use shapes for one and colors for the other but if =1 will just use color
# catogories tells funciton which catagories it is labeling through shapes/colors on plot
plot_pca <- function(vsd, num_catagories=2, catagories=c("genotype","batch"), 
                     dir=file.path(local_dir, 'results/plots'), titl="ddsPCA", savetoPDF=TRUE) {
  pcaData <- plotPCA(vsd, intgroup = catagories, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaData[,4] <- factor(pcaData[,4])
  if (num_catagories==2) {
    pcaData[,5] <- factor(pcaData[,5])
  } 
  if (num_catagories==2) {
    ggplot(pcaData, aes(x = PC1, y = PC2, color =pcaData[,4], shape = pcaData[,5])) +
      geom_point(size =3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed() +
      ggtitle("Principle Component Analysis")
  } else {
    ggplot(pcaData, aes(x = PC1, y = PC2, color = pcaData[,4])) +
      geom_point(size =3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed() +
      ggtitle("Principle Component Analysis")
  }
  if (savetoPDF==TRUE) {
    if(!dir.exists(dir)) {
      dir.create(dir)
    }
    titl <- paste(titl, ".pdf", sep="")
    ggsave(filename=titl, device="pdf",path=dir)
  }
}

#############################################
# PLOT CUSTOM HEATMAP GIVEN A LIST OF GENES #
#############################################
custom_heatmap <- function(vsd, lst, dir=file.path(local_dir, 'results/plots/heatmaps'), 
                           titl, anno, group_num=2) { 
  # lst=list of interesting genes to plot
  # group_num= what number of biggest clusters you want to split the heatmap by
  # norm_counts=matrix of normalized counts
  # anno=the metadata features you want to label on the heatmap 
  # (eg label by batch and genotype so we need a table with one column 
  # of batch values and one for genotype values with sample names as rownames)
  .fs <- function(n) { # shrink font size depending on number of genes 
    if(n < 100) { fs <- 6 } else if(n < 200) {
      fs <- 5 } else if(n < 500) { fs <- 4 } 
    else if(n < 1000) { fs <- 3 } else { fs <- 3 
    }
    return(fs)
  }
  norm_counts <- as.matrix(assay(vsd))
  colnames(norm_counts) <- colData(vsd)$names
  
  custom_genes <- intersect(lst,rownames(vsd))
  # make sure these genes are also in the vst (should be!)
  submatrix <- subset(norm_counts,rownames(norm_counts)%in%custom_genes) 
  mat <- submatrix - rowMeans(submatrix)
  
  if(!dir.exists(dir)) {
    dir.create(dir)
  }
  titl <- paste(titl,".pdf",sep="")
  pheatmap(mat, annotation_col = anno, cutree_cols=group_num,
           fontsize_row = .fs(length(custom_genes)), #cutree_rows=2,
           filename=paste(dir,titl,sep="/"),cluster_rows = TRUE)#,
           #cluster_cols=FALSE)
}

##################################################
############## PLOT VOLCANO PLOT #################
##################################################
# dir is the directory to save the PDF to (if TRUE) and titl is name of file
plot_volcano <- function(data,gene_anno=data$symbols,dir=file.path(local_dir, 'results'),
                         xlab='log2FoldChange',ylab='padj', savetoPDF=TRUE, titl) {
  p <- EnhancedVolcano(data, lab=gene_anno, x=xlab, y=ylab)
  if(savetoPDF==TRUE) {
    if(!dir.exists(dir)) {
      dir.create(dir)
    }
    fname <- paste(dir,"/",titl,".pdf",sep="")
    pdf(fname, width=11, height=7)
    plot(p)
    dev.off() 
  }
}

###############################################
############## PLOT GO TERMS ##################
###############################################
# go is the list of terms obtained by running undirected gsea
# n is the number of terms you want to plot
# col is the color you want the bars to be
plotGSEA_undirected <- function(go, n1, n2, col, titl, saveToPDF=TRUE, 
                                dir=file.path(local_dir, 'results'), 
                                label="Term", direction="up") {
  if (direction=="up") {
    goSub <- go[order(go$P.Up, decreasing=FALSE, na.last=TRUE), ]
  } else {
    goSub <- go[order(go$P.Down, decreasing=FALSE, na.last=TRUE), ]
  }
  goSub <- go[n1:n2,]
  pathways <- c()
  for (i in 1:length(rownames(goSub))) {
    name <- paste(rownames(goSub)[i],goSub[i,label],sep=": ")
    pathways <- c(pathways, name)
  }
  #pathways <- paste(rownames(goSub),goSub[label],sep=": ")
  prob <- goSub$P.Up
  plotGo <- data.frame(pathways=pathways,prob=prob)
  p <- ggplot(data=plotGo, aes(x=pathways, y=-log10(prob)))  + 
    geom_bar(stat="identity", color=col, fill=col, width=0.5) + 
    coord_flip() + xlab("GO Pathways") + ylab("-log10 P Value") 
  if (saveToPDF==TRUE) {
    filename <- file.path(dir,paste(titl,".pdf",sep=""))
    pdf(filename, width=7, height=7)
    print(p)
    dev.off()
  } else {
    plot(p)
  }
}
###############################################
############## COMPARE GENE LISTS #############
###############################################
# Given a interesting subset list (compare_lst), see how many genes are also in the parent list
# Returns the final shared list of genes between the two sets and prints results
compare_lists <- function(parent_lst, compare_lst, 
                          compare_titl="the", parent_titl="significant DGE") {
  final <- parent_lst[which(parent_lst%in%compare_lst==TRUE)]
  cat(paste("Out of the ", length(compare_lst), " genes from ", compare_titl, "gene set, ", 
            length(final), " of these genes were found in the ", parent_titl, "!\n", sep=""))
  return(final)
}

###############################################
########## MAKE PATHMAP WITH PATHVIEW #########
###############################################
path_map <- function(data, id, dir=file.path(results_dir, "plots")) {
  if(!dir.exists(dir)) {
    dir.create(dir)
  }
  setwd(dir) # so that the map saves to this location
  logFC <- data$table$logFC
  names(logFC) <- data$table$Entrez
  pathview(gene.data=logFC,pathway.id=id,
           species="mcc", kegg.native=FALSE)
}
