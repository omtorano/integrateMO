# integrateMO 
This is a wrapper package to normalize and integrate common omics data types collected during toxicological exposure studies. This package contains two functions import_MO(), see R/import_MO(), for importing data with an option to normalize, and integrate_MO(), see R/integrate_MO(), to integrate omics layers.

What does "wrapper package" mean here? It means that algoriths used to integrate data employed within integrate_MO() come from other packages. These statistical methods were developed by other groups and deployed within R packages/scripts. integrateMO is a package developed to streamline the use of these disperate methods for omics data collected on toxicological exposure data. Using integrateMO means you're actually using one of the "wrapped" packages and citations should be made accordingly, see "References & links for information on integration methods".

# Import
import_MO() has eight parameter inputs

- rnaseq_counts: Data matrix of RNAseq counts 
- metab_peaks: Data matrix of metabolite peaks
- rrbs_mvals: Data matrix of RRBS m-values
- mirna_counts: Data matrix of microRNA counts
- pirna_counts: Data matrix of piRNA counts
- meta: Metadata, must have sample ids in first column matching those in omics data and treatment col called TRT, if batch correction is required must also have unique batch column for each omics layer with batches. 
- norm: Normalization, true or false, default true
- batch: Omics layers with batches corresponding to "batch_omic" columns in metadata, column labels as follows: batch_rna, batch_metab, batch_rrbs

Example usage
```
import_MO(rnaseq_counts = counts, rrbs_mvals = mvals, meta = meta) #minimum of 2 omics required, default normalization = TRUE
import_MO(rnaseq_counts = counts, rrbs_mvals = mvals, metab_peaks = metab, meta = meta, batch = “batch_rna”) #will normalize all three data layers and correct for batch in RNAseq data
```
The orientation of the omics data matricies does not matter. Samples can be either row names or column names and features can be row names or column names. Sample names must be either row or column names (i.e. not the first row or first column) and must match the row names of the metadata.
## Output
This function will save a data list “data_list” to the global environment. The elements of this list are composed of the omics layers provided. If the normalization option is set to “true” this function will also output a MOnorm folder to the current working directory. This folder will be labeled with the current date and time, so rerunning will not overwrite previous results. This folder will contain visualizations of the omics layers provided including boxplots, PCA, scree, and MDS plots. 

## Normalization 
The following steps are executed if norm = TRUE
### RNAseq counts
Features with low counts are removed with the filterByExpr() function in edgeR. Log2 counts per million adjusted for library size are calculated and saved for integration. Library size is normalized using the trimmed means of m-values method via calcNormFactors() in edgeR. Counts are extracted using cpm(x, log = TRUE, normalized.lib.sizes = TRUE).
If batch correction is indicated with batch = “rnaseq_counts” correction is carried out with ComBat() from the sva package.
### Metabolite sum peak area
Metabolite feature processing was chosen to match standard protocol already in use within EPA. Feature values are first multiplied by 1000, this is done so results match those of metaboanalyst, which is commonly used within EPA. Then rows are mean centered and values are pareto scaled using pareto_scale() from the IMIFA package.

### miRNA and piRNA
There is currently no normalization procedure carried out by import_MO() for these data types. MIB has yet to measure these data types in a multi-omics study and they may be removed in future versions of this package.

### RRBS M-values
These features do not require normalization. To generate m-values from RRBS Bismark files see the formatRRBS package.

# Integrate
integrate_MO() has two parameter inputs

- int_method: sPLS-DA, MOFA2, WGCNA, SNF, or iPCA
- RRBS_feature_map: dataframe of RRBS features to genes, must match format of Unique_Features_to_Genes.csv from formatRRBS output


Example usage
```
integrate_MO(int_method = "SNF") #default RRBS_feature_map = NA
integrate_MO(int_method = "sPLS-DA", RRBS_feature_map = FL-F_1000-1000_0.8_Unique_Features_to_Genes) 
```

Integration method specifies which package/method will be employed to integrate data. sPLS-DA is sparse partial least squares discriminant anslysis from the mixOmics package, multi-omic factor analysis from MOFA2, weighted gene correlation network analysis from WGCNA, similarity fusion network from SNF, or integrated principal components analysis from iPCA.
output 

## References & links for information on integration methods
Note that integrate_MO() is a wrapper for the below packages and appropriate citations to the original packages should be used.

mixOmics sPLS-DA (DIABLO N-integration)
https://mixomicsteam.github.io/mixOmics-Vignette/id_06.html#id_06:diablo


WGCNA
https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/index.html
Zhang B and Horvath S (2005) A General Framework for Weighted Gene Co-Expression Network Analysis, Statistical Applications in Genetics and Molecular Biology: Vol. 4: No. 1, Article 17 PMID: 16646834
Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559)

MOFA2

SNF

iPCA
