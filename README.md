# Overview
This is a wrapper package to normalize and integrate common omics data types collected during toxicological exposure studies. This package contains two functions: 
	- import_MO(), see R/import_MO(), for importing data with an option to normalize
	- integrate_MO(), see R/integrate_MO(), to integrate omics layers

"wrapper package" means that algoriths used to integrate data employed within integrate_MO() come from other packages. These statistical methods were developed by other groups and deployed within R packages/scripts. integrateMO is a package developed to streamline the use of these disperate methods for omics data collected on toxicological exposure data. Using integrateMO means you're actually using one of the "wrapped" packages and citations should be made accordingly, see "References & links for information on integration methods".

# Install importMO & associated packages
Most packages will be installed within the installation of integrateMO. Several packages must be installed seperately which can be accomplished by running install_run.R. The following code can be run to install integrateMO from github.
```
# Install the package from GitHub & add to library
devtools::install_github("omtorano/integrateMO")
library(integrateMO)
```

# Step 1: Import data with import_MO()
Step one of the integration workflow is to import omics data & associated metadata with the importMO() function. Within this function data are reformatted and an optional normalization step is executed, see Normalization section for further details. importMO() must be run before integrateMO().

import_MO() has eight parameter inputs:

- rnaseq_counts: Data matrix of RNAseq counts 
- metab_peaks: Data matrix of metabolite peaks
- rrbs_mvals: Data matrix of RRBS m-values
- mirna_counts: Data matrix of microRNA counts
- pirna_counts: Data matrix of piRNA counts
- meta: Metadata, must have sample ids in first column matching those in omics data and treatment col called TRT, if batch correction is required must also have unique batch column for each omics layer with batches 
- norm: Normalization, true or false, default true
- batch: Omics layers with batches corresponding to "batch_omic" columns in metadata, column labels as follows: batch_rna, batch_metab, batch_rrbs

Example usage
```
import_MO(rnaseq_counts = counts, rrbs_mvals = mvals, meta = meta) #minimum of 2 omics required, default normalization = TRUE
import_MO(rnaseq_counts = counts, rrbs_mvals = mvals, metab_peaks = metab, meta = meta, batch = “batch_rna”) #will normalize all three data layers and correct for batch in RNAseq data
```
The orientation of the omics data matricies does not matter. Samples can be either row names or column names and features can be row names or column names. Sample names must be either row or column names (i.e. not the first row or first column) and must match the row names of the metadata.
## Function output
Running import_MO() will save a data list named “data_list” to the global environment. The elements of this list are composed of the omics layers provided and will be used automatically as the input for integrateMO(). If the normalization option is set to “TRUE” this function will also output a MOnorm folder to the current working directory. This folder will be labeled with the current date and time, so rerunning will not overwrite previous results. The MOnorm folder will contain visualizations of the omics layers provided including boxplots, PCA, scree, and MDS plots. 

## Normalization 
The following normalization steps are carried out for the given omics layers when norm = TRUE.
### RNAseq counts
Features with low counts are removed with the filterByExpr() function in edgeR. Log2 counts per million adjusted for library size are calculated and saved for integration. Library size is normalized using the trimmed means of m-values method via calcNormFactors() in edgeR. Counts are extracted using cpm(x, log = TRUE, normalized.lib.sizes = TRUE).
If batch correction is indicated with batch = “rnaseq_counts” correction is carried out with ComBat() from the sva package.

Links & references
- Robinson MD, McCarthy DJ, Smyth GK (2010). “edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.” Bioinformatics, 26(1), 139-140. doi:10.1093/bioinformatics/btp616.
- Leek, J. T., Johnson, W. E., Parker, H. S., Fertig, E. J., Jaffe, A. E., Storey, J. D., ... & Torres, L. C. (2019). sva: Surrogate variable analysis. R package version, 3(0), 882-883.

### Metabolite sum peak area
Metabolite feature processing was chosen to match standard protocol already in use within EPA. Feature values are first multiplied by 1000, this is done so results match those of metaboanalyst, which is commonly used within EPA. Then rows are mean centered and values are pareto scaled using pareto_scale() from the IMIFA package.

Links & references
- Murphy, K., Viroli, C., & Gormley, I. C. (2020). Infinite mixtures of infinite factor analysers.

### miRNA and piRNA
There is currently no normalization procedure carried out by import_MO() for these data types. MIB has yet to measure these data types in a multi-omics study and they may be added or removed in future versions of this package.

### RRBS M-values
These features do not require normalization. To generate M-values from RRBS Bismark files see the formatRRBS package https://github.com/omtorano/formatRRBS.

# Step 2: Integrate data with integrate_MO()
After running import_MO() the imported omics data sets will automatically be detected by integrate_MO(). integrate_MO() therefore has only two parameter inputs: int_method (integration method) allows users to select the data integration method to be used, see "integration method" section, and the optional RRBS_feature_map parameter which allows users to specify the path to RRBS feature map input.

- int_method: sPLS-DA, MOFA2, WGCNA, SNF, or iPCA
- RRBS_feature_map: dataframe of RRBS features to genes, must match format of Unique_Features_to_Genes.csv from formatRRBS output

Example usage
```
integrate_MO(int_method = "WGCNA") #default RRBS_feature_map = NA
integrate_MO(int_method = "sPLS-DA", RRBS_feature_map = ML-0_1000-1000_0.8_Unique_Features_to_Genes) 
```

## Function output
Running integrate_MO() will generate a folder in the current working directory with the following naming convention "integrateMO_*integration method*_*year-month-day time*". The contents of the folder depend on the integration method chosen, details below.

## Integration method
Integration method specifies which package & method will be employed to integrate data. Note that integrate_MO() is a wrapper for the packages described below and appropriate citations to the original packages should be used.
### sPLS-DA  
Multiblock sparse partial least squares discriminant anslysis from the mixOmics package (DIABLO N-integration). The steps performed in the integrate_MO() function closely follow those described in the mixOmics vignette linked below.
Important decision points and outputs:
- Features from each omic layer are limited to top 10,000 features ordered by median absolute deviation.
- The design matrix is determined using mixOmics' "data driven" approach. For each omics group cross comparison, PLS with one component is calculated, the
cross-correlations between components are averaged and used as the off diagonal values of the design matrix.  
- The number of components is first set to # of treatments + 1, performance evaluation is done with perf(), folds = number of samples in treatment level
with minimum replicates, nrepeat = 10, final number of components and distance metric selected will be saved in the figure title of the classification_error_rate plot.
- The number of variables to select is determined with tune.block.splsda(), folds = 5 and nrepeat = 10. The sequence of numbers tested is c(5:9, seq(10, 50, 5)). The final number
of variable selected on each component will match the number of rows in the splsda_variables_comp.csv files.
- The design, number of components, and number of variables to select are saved and input as parameters in the final block.splsda() model. Visualization outputs from final model are weighted average plotIndiv() plots, correlation circle
plots plotVar(), and clustered image maps cimDiablo() of each component. In tests with tox data sets, the number of components in the final model has been <= 4. Therefore weighted average and correlation circle plots are created for pairs of
components up to ncomp = 4.  

Estimated run time and resourse usage - 
On a 11th Gen Intel(R) Core(TM) i7-1185G7 PC with 4 cores and 8 threads a test run took 3915.49 seconds, 148.90 seconds user CPU time, and 6.11 system CPU time. Saved output from tests is <2MB.

Links & references
- https://mixomicsteam.github.io/mixOmics-Vignette/id_06.html#id_06:diablo  
- Rohart, F., Gautier, B., Singh, A., & Lê Cao, K. A. (2017). mixOmics: An R package for ‘omics feature selection and multiple data integration. PLoS computational biology, 13(11), e1005752.  
- Singh, A., Shannon, C. P., Gautier, B., Rohart, F., Vacher, M., Tebbutt, S. J., & Lê Cao, K. A. (2019). DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. Bioinformatics, 35(17), 3055-3062.  

### MOFA2
multi-omic factor analysis from MOFA2
Links & references
- https://biofam.github.io/MOFA2/index.html
- https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html  
- https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/downstream_analysis.html  
- Argelaguet, R., Velten, B., Arnol, D., Dietrich, S., Zenz, T., Marioni, J. C., ... & Stegle, O. (2018). Multi‐Omics Factor Analysis—a framework for unsupervised integration of multi‐omics data sets. Molecular systems biology, 14(6), e8124.  
- Argelaguet, R., Arnol, D., Bredikhin, D., Deloro, Y., Velten, B., Marioni, J. C., & Stegle, O. (2020). MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. Genome biology, 21(1), 1-17.  


### WGCNA
Weighted gene correlation network analysis from the WGCNA package. The steps performed in the integrate_MO() function are based on the WGCNA tutorial in WGCNAtutorials.
Note that this method is currently set up for pairwise comparisons between omics layers.    

Steps, important decision points, and output:  
- If present, RRBS features are filtered to top 25% most variable features ordered by median absolute deviation. 
- Omics features are separately clustered for viz purposes with 7 different clustering algorithms within the hclust() function, output is saved in hclust_sampleTree file. Information below gathered from https://cicerocq.wordpress.com/wp-content/uploads/2019/05/cluster-analysis_5ed_everitt.pdf Chapter 4 Hierarchical Clustering, Statistics and Data Analysis in Geology Chapter 6 and Unsupervised Methods Allen (slides in resources).
	- Single aka nearest neighbor: distance between groups is defined as that of the closest pair of individuals. Sensitive to outliers so can be useful in identifying potential outliers.
	- Complete aka furthest neighbor: distance between two groups is defined as most distant pair of individuals. 
	- Average - UPGMA unweighted pair-group method using the average approach: distance between two clusters is average of the distance.  between all pairs of individuals that are made up of one individual from each group, in short the average of all pairwise distances.  
	- Centroid UPGMC unweighted pair-group method using centroid approach: merges clusters with most similar mean vectors.
	- Median WPGMC weighted pair-group methods using the centroid approach: weights centroids equally to produce new centroid of merged cluster. This avoids clusters with numerous objects dominating cluster with few objects when merged.
	- Mcquitty WPGMA weighted average linkage: also weights clusters, points in small clusters weighted more highly than points in large clusters.
	- Ward.D2 fusion of clusters is based on the size of error sum-of-square criterion. Similar to centroids method but centroids are weights. Sensitive to outliers.
- For omic layers with between 100 and 5,000 features, modules are determined with moduleEigengenes() and saved to a lowcount_omics_MEs object. For omic layers with >5,000 features, modules are created with blockwiseModules() and saved to a
block_MEs object.
	- lowcount_omics_MEs: soft thresholding powers tested are c(1:10, seq(from = 12, to = 40, by = 2)), power chosen with pickSoftThreshold(). The power that is chosen will be red in resulting soft_thresholding.svg figure. The
	minimum module size is set to 30 features. Eigengenes are clustered using method = "average". The dissimilarity threshold is set to 0.25, which is used for tree cut height. Uses default pearson unsigned network construction.
	- block_MEs: same as above, additional options include maxBlockSize = 5000.
- Output includes hclust_sampleTree.svg showing hierarchical clustering of each omic lyer,
Module_Eigengenes_heatplot.svg for each omic layer showing the relationship between module eigengenes and samples, Module_Eigengenes.csv contains the eigengene values,
Module_membership.csv has correlation values of a gene to a module eigengene, Module-module.csv of correlation and p-values show these values for the relationship between
modules belonging to each omic layer. The dendro_*omic_layer*.pdf files show module clustering for each omic layer.  

Estimated run time and resource usage  
On a 11th Gen Intel(R) Core(TM) i7-1185G7 PC with 4 cores and 8 threads a test run took 1624.73 seconds, 1535.94 seconds user CPU time, and 66.34 system CPU time. Saved output from tests is <2MB.


Links & references
- https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/index.html  -- it looks like this is offline as of 11/7/2023, there are other tutorials but no permanant replacement yet. The contents of this tutorial are available in this repo, see WGCNAtutorials (see https://bioinformatics.stackexchange.com/questions/21885/where-to-access-the-wgcna-tutorial-documents-horvath-lab-site-down). Other tutorials are also available including https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0, https://edu.sib.swiss/pluginfile.php/158/course/section/65/_01_SIB2016_wgcna.pdf, https://pages.stat.wisc.edu/~yandell/statgen/ucla/WGCNA/wgcna.pdf.
- Zhang B and Horvath S (2005) A General Framework for Weighted Gene Co-Expression Network Analysis, Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 PMID: 16646834  
- Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)

### SNF
Similarity fusion network from SNFtool package.

Steps, important decision points, and output:
- If present, RRBS features are filtered to top 25% most variable features ordered by median absolute deviation. 
- Parameter settings: Number of nearest neighbors is set to 20 or the number of samples - 1, if fumber of samples is less than 20. The hyperparameter alpha is kept at default recommended 0.5, the iteration parameter T is set to 20. The number of clusters is set to the number of unique treatments in metadata.
- Pair-wise distances and similarity graphs are calculated for each omic layer with the dist2() and affinityMatirx() functions. These are then fused with SNF() into a unified graph. 
- Spectral clustering is then performed with spectralClustering() and normalized mutual information (NMI) is calculated to measure clustering results, a value close to 1 indicates good clustering.
- The values from the unified graph and a heatmap of these values are saved as SNF_cluster.csv and SNF_cluster_heatmap.svg respectively. In addition, a network plot of a subset of the highest values from the unified graph is saved to network_map.svg. NMI information is included in this plot.

Links & references
- https://www.rdocumentation.org/packages/SNFtool/versions/2.3.1/topics/SNF  
- https://github.com/cran/SNFtool  
- Wang, B., Mezlini, A. M., Demir, F., Fiume, M., Tu, Z., Brudno, M., ... & Goldenberg, A. (2014). Similarity network fusion for aggregating data types on a genomic scale. Nature methods, 11(3), 333-337.

### iPCA
or integrated principal components analysis from iPCA.
Links & references
- https://github.com/DataSlingers/iPCA  
- Tang, T. M., & Allen, G. I. (2021). Integrated principal components analysis. The Journal of Machine Learning Research, 22(1), 8953-9023.


## RRBS_feature_map
If RRBS methylation loci features are being input, the resulting output will have the "chromosome-location" format for any output with feature names. 
Including the RRBS_feature_map option will include associated genes in .csv output with feature names. If RRBS M-values were formatted with the formatRRBS 
package, the path to the Unique_Features_to_Genes.csv file should be the input for this parameter. If RRBS M-values were not formatted with the formatRRBS
package, the input for this parameter must be a .csv of methylation feature loci and associated genes. The format of this spreadsheet must have chromosome locations
in a column labeled "chrom", location along the chromosome in a column labeled "loc", and a column of associated genes as the fourth column of the spreadsheet.

If formatRRBS was used, "association" of methylation loci and genes is defined as a methylation location falling within the defined base pair window of a given gene. Default parameters for formatRRBS
define this window as 1000bp, so a methylation location is associted with a gene if its location falls within 1000bp upstream, within the gene body, or within 1000bp downstream of the gene. 
