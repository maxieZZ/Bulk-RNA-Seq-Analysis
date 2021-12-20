# Bulk-RNA-Seq-Analysis

> This code was made to analyze the Bulk Expression RNA data from two data-sets sequenced seperately via Novaseq.
> Previous to this analysis, both data-sets had been analyzed seperately via DESeq2 and a tab deliminated text file 
> of significantly differentially expresssed gene symbols shared between the two results was saved (labeled in the 
> results folder as "sharedGenes.txt." Because the data is still in the process of publishing, the original expression 
> datasets themselves have not been uploaded and data that is visible has been altered. 
>
> The code begins by reading in the two expression datasets (saved as rangeSummerizedExperiment RDS files as output
> from tximeta), creating shared columns between the two datasets, and then finally combining these data sets based
> on shared columns into one single rangeSummerizedExperiment (saved to a file). Dataset1 has 7 samples, all cells from
> the same macaque, 4/7 of which have a Mutant (Mut) genotype and 3/7 of which have a Wild type (HT) genotype. The seond 
> The second data set has 6 samples, 3/6 of which have a Mutant (Mut) genotype and 3/6 of which have a Wild type (HT) genotype.
> The second differs slightly in its labeling of "genotype" as "condition" (this is corrected to keep consistency
> between the labeling of the two datasets.
>
> After combining the expression datasets and adding annotation, exploratory analysis is preformed to visualize 
> top principle components and to observe trends in euclidean distance between genotype. After intial visualizations
> (using variance stabilized ie normalized expression data) methods DESeq2 and EdgeR are both used seperately to 
> perform differential gene expression analysis, the results of which may be saved in heatmaps, volcanoplots, and files.
> After performing each method seperately, significant and differentially expressed genes common between the two methods
> are found and methods are compared. 
>
> Undirected gene set analysis was also performed on the differential expression results of EdgeR via goanna and kegga 
> functions and can be visualized via custom built barchart functions. Pathway tool was also used breifly to visualize 
> the differential expression of genes in pathway maps of interest obtained from kegg.
>
> A final list of significant differentially expressed genes shared between DESeq2, EdgeR, and also the old separated 
> method of analysis was then used as the starting point for module analysis. 11 gene modules from this shared list where
> created using a custom function analogus to a pearson correlation blend. Subsets of the similarity and adjacency matrixes 
> were plotted for better visualization and the final module results were analyzed via heatmap plots. Heatmap plots and data
> were saved to files and three modules were selected based on geotype clustering for further analysis and subsetting.

>
![Similarity Matrix](https://github.com/maxieZZ/Bulk-RNA-Seq-Analysis/tree/master/results/images/similarityMatrixSubset.png)
![Adjacency Matrix](https://github.com/maxieZZ/Bulk-RNA-Seq-Analysis/tree/master/results/images/adjacencyMatrixSubset.png)


