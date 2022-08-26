### Raw Reads
Paired-end 150 bp mRNA reads from *Poecilia mexicana* gill tissues can be downloaded using the NCBI accession PRJNA290391.

The variable `$sample` should be replaced with the full sample name.

### Trim
Trim FASTQ files containing the raw reads using TrimGalore! (v.0.4.2) with FastQC (v.0.11.4) and Cutadapt (v.1.9).

Adapter trim.
```{bash}
trim_galore \
  --output_dir /path/to/outdir \
  --quality 0 \
  --fastqc \
  --fastqc_args "--nogroup --noextract --outdir /path/to/outdir" \
  -path_to_cutadapt /path/to/cutadapt \
  --illumina \
  --stringency 6 \
  --clip_R1 11 \
  --clip_R2 11 \
  --paired \
  --gzip ${sample}_1.fq.gz ${sample}_2.fq.gz
```

Quality trim.
```{bash}
trim_galore \
  --output_dir /path/to/outdir \
  --quality 24 \
  --fastqc \
  --fastqc_args "--nogroup --noextract --outdir /path/to/outdir" \
  -path_to_cutadapt /path/to/cutadapt \
  --illumina \
  --stringency 6 \
  --length 50 \
  --paired \
  --gzip ${sample}_1_val_1.fq.gz ${sample}_2_val_2.fq.gz
```

### Map
Index the *Poecilia mexicana* reference genome using HISAT2 (v.2.1.0). Reference genome and annotations in GFF format were downloaded from NCBI (Accession: GCF_001443325.1) with mitochondrial sequence appended (Accession: KC992995.1).
```{bash}
hisat2-build \
  -f GCF_001443325.1_P_mexicana-1.0_genomic_with_mito_sequence.fa \
  GCF_001443325.1_P_mexicana-1.0_genomic_with_mito_sequence
```

Align reads to the reference genome.
```{bash}
hisat2 \
  -q \
  --threads 4 \
  -k 5 \
  --downstream-transcriptome-assembly \
  --fr \
  --rg SM:${sample} \
  --rg ID:${rgid} \
  --rg PL:ILLUMINA \
  -x GCF_001443325.1_P_mexicana-1.0_genomic_with_mito_sequence \
  -1 ${sample}_1_val_1_val_1.fq.gz \
  -2 ${sample}_2_val_2_val_2.fq.gz \
  -S ${sample}_HISat_out.sam \
  --summary-file ${sample}_statistics.txt
```

### Modify SAM files
Performed using samtools (v.1.9).<br>
Convert to BAM.
```{bash}
samtools view \
  -bSh ${sample}_HISat_out.sam > \
  ${sample}_HISat_out.bam
```

Sort by coordinate.
```{bash}
samtools sort \
  -@ 8 \
  ${sample}_HISat_out.bam \
  -o ${sample}_sorted.bam
```

Merge BAM files from the same individuals.
```{bash}
samtools merge \
  -r \
  ${sample}_merged_sorted.bam \
  ${sample}_sorted.bam \
  ${sample}_sorted.bam
```

### Create gene and transcript counts matrices
Generate GTF files for each sample using StringTie (v.2.0.3) and associated Python script prepDE.py.
```{bash}
mkdir -p ballgown

stringtie \
  ${sample}_sorted.bam \
  -o ballgown/${sample}/${sample}.gtf \
  -p 8 \
  -G GCF_001443325.1_P_mexicana-1.0_genomic_with_mito_sequence-edited.gff \
  -B \
  -e
```
Generate gene and transcript counts matrices. `sample_list.txt` contains sample prefixes in column 1 and paths to GTF files in column 2.
```{bash}
module load python/2.7.10

python prepDE.py \
  -i sample_list.txt \
  -g gene_count_matrix.csv \
  -t transcript_count_matrix.csv
```

Delete the STRG genes (mitochondrial) from the gene counts matrix.
```{bash}
# Gene counts matrix
grep -v "STRG" gene_count_matrix.csv > gene_count_matrix_no_STRG.csv

# Transcript counts matrix
grep -v "STRG" transcript_count_matrix.csv > transcript_count_matrix_no_STRG.csv
```

All subsequent analyses were run using R (v.3.6.0 for WGCNA, v.4.0.3 for EdgeR).

Install packages.
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("data.table")
install.packages("dplyr")
BiocManager::install("edgeR")
install.packages("ggplot2")
BiocManager::install("limma")
install.packages("scales")
install.packages("splitstackshape")
install.packages("tidyr")
BiocManager::install("WGCNA") 
```

Load packages.
```{r}
library(data.table)
library(dplyr)
library(edgeR)
library(ggplot2)
library(limma)
library(scales)
library(splitstackshape)
library(tidyr)
library(WGCNA)
```

### Network Analysis (Weighted Gene Co-expression Network Analysis, WGCNA)
#### WGCNA Part 1

Read in gene counts matrix.
```{r}
wild_data <- read.csv(file = 'gene_count_matrix_no_STRG.csv', row.names = 1)
```

Create grouping factors to be incorporated later into the DGELists.<br>
W = wild-caught samples, S = sulfidic habitat, NS = non-sulfidic habitat.<br>
Rosita and Banos = Pichucalco Drainage, LaLluvia and Trib1 = Puyacatengo Drainage, Bonita and ElAzufreI = Tacotalpa Drainage
```{r}
WildGroup = c(rep("W_Rosita_NS", 6), rep("W_Banos_S", 6), rep("W_LaLluvia_S", 5), rep("W_Trib1_NS", 6), rep("W_Bonita_NS", 6), rep("W_ElAzufreI_S", 6))
```

Create DGEList object.
```{r}
y_wild <- DGEList(counts=wild_data, group=WildGroup)
```

Filtering calculations.
```{r}
# Sort samples by library size
y_wild_samples_sorted <- y_wild$samples[order(y_wild$samples$lib.size) , ]

# Select lowest library size
wild_lowest_library_size <- y_wild_samples_sorted$lib.size[1]

# Calculate minimum counts per million (cpm) to have ~5 counts in the smallest library, as per EdgeR manual recommendations
x_wild <- (5*1000000)/wild_lowest_library_size

# Samples must have a cpm of at least 1.144 (~5 counts in smallest library) and be present in at least 5 samples
keep_wild <- rowSums(cpm(y_wild)>=1.144) >= 5
filtered_y_wild <- y_wild[keep_wild, , keep.lib.sizes=FALSE]
```

Calculates normalization factors to scale raw library sizes, TMM is default method.
```{r}
filtered_y_wild <- calcNormFactors(filtered_y_wild)
```

Multidimensional scaling plot (MDS) of the top 500 genes (the default).
```{r}
pdf("0_MDSplot_500genes_wild_symbols_ecotype.pdf")
par(cex.axis=1, cex=1)
# NS plotted in blue, S plotted in gold
plotMDS(filtered_y_wild, gene.selection='common', top=500, method='logFC', cex=2, 
                  col=alpha(c(rep("blue", 6), rep("gold", 6), rep("gold", 5), rep("blue", 6), rep("blue", 6), rep("gold", 6)), 0.7),
                  pch=c(rep(19, 6), rep(19, 6), rep(17, 5),rep(17, 6),rep(15, 6),rep(15, 6)),
                  xlim=c(-2.5, 2.5),
                  ylim=c(-2.5, 2.5))
legend("bottomright", 
       inset = 0.02, 
       legend = c("Pichucalco", "Puyacatengo", "Tacotalpa", "Non-sulfidic", "Sulfidic"), 
       pch = c(19, 17, 15, NA, NA), 
       lty = c(NA, NA, NA, 1, 1),
       lwd = c(NA, NA, NA, 4, 4),
       col = c(rep("black", 3), "blue", "gold"), 
       cex = 1,
       pt.cex = 1.2,
       bg = "white")
dev.off()
```

Note, a majority of the WGCNA code and annotation below was copied from https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/, and slightly modified to fit the design of this study.

#### Data input, cleaning, and pre-processing.
The following setting is important, do not omit
```{r}
options(stringsAsFactors = FALSE)
```
Make a data frame of log counts per million from EdgeR.
```{r}
counts_per_mil_wild <- cpm(filtered_y_wild, normalized.lib.sizes=TRUE, log=TRUE, prior.count=1)
counts_per_mil_wild <- as.data.frame(counts_per_mil_wild[,])
```

Transpose the data frame.
```{r}
datExprWild0 = as.data.frame(t(counts_per_mil_wild[,]))
```

#### Check data for excessive missing values and identify outlier samples
Check for genes and samples with too many missing values.
```{r}
gsg_wild = goodSamplesGenes(datExprWild0, verbose = 3)
cat("Do all wild-caught genes pass cuts?", gsg_wild$allOK)
```

If above is TRUE, all genes passed cuts; if not, bad genes and samples are removed.
```{r}
if (!gsg_wild$allOK)
{
  # Optionally, print the gene and sample names that were removed
  if (sum(!gsg_wild$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExprWild0)[!gsg_wild$goodGenes], collapse = ", ")))
  if (sum(!gsg_wild$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExprWild0)[!gsg_wild$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data
  datExprWild0 = datExprWild0[gsg_wild$goodSamples, gsg_wild$goodGenes]
}
```

Cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers.
```{r}
sampleTreeWild = hclust(dist(datExprWild0), method = "average")
```

Plot the sample tree.
```{r}
pdf(file = "1_wild_sample_clustering.pdf", width = 12, height = 9)
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTreeWild, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Choose a height cut that will remove any offending samples (in this case none) and use a branch cut at that height
# Plot a line to show the cut
abline(h = 220, col = "red")
dev.off()
```

Determine cluster under the line.
```{r}
clust_wild = cutreeStatic(sampleTreeWild, cutHeight = 220, minSize = 10)
print(table(clust_wild), row.names=getOption("datatable.print.rownames"))
```

Cluster 1 contains the samples we want to keep.
datExprWild contains the expression data that can be used for network analysis.
```{r}
keepSamplesWild = (clust_wild==1)
datExprWild = datExprWild0[keepSamplesWild, ]
nGenesWild = ncol(datExprWild)
nSamplesWild = nrow(datExprWild)
```

#### Input trait data
Make a data frame of traits.
```{r}
# Habitat: non-sulfidic (NS) = 0, sulfidic (S) = 1
NS_vs_S_wild <- c(rep(0,6), rep(1,6), rep(1,5), rep(0,6), rep(0,6), rep(1, 6))

# Pichucalco (Pich) drainage = 1, all other drainages = 0
Pich_vs_all_wild <- c(rep(1,12), rep(0,11), rep(0,12))

# Puyacatengo (Puya) drainage = 1, all other drainages = 0
Puya_vs_all_wild <- c(rep(0,12), rep(1,11), rep(0,12))

# Tacotalpa (Taco) drainage = 1, all other drainages = 0
Taco_vs_all_wild <- c(rep(0,12), rep(0,11), rep(1,12))

datTraitsWild <- data.frame(NS_vs_S_wild, Pich_vs_all_wild, Puya_vs_all_wild, Taco_vs_all_wild)
```

Visualize how the clinical trails relate to the dendrogram.
```{r}
# Re-cluster samples
sampleTree2Wild = hclust(dist(datExprWild), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColorsWild = numbers2colors(datTraitsWild, signed = FALSE)

# Plot the sample dendrogram and the colors underneath
pdf(file = "2_wild_dendrogram_and_trait_heatmap.pdf")
plotDendroAndColors(sampleTree2Wild, traitColorsWild,
                    groupLabels = names(datTraitsWild),
                    main = "Sample dendrogram and trait heatmap"
dev.off()
```

Save as RData file.
```{r}
save(datExprWild, datTraitsWild, file = "wild-01-dataInput.RData")
```

#### Step-by-step network construction and module detection

#### Preliminaries: setting up the R session
Allow multi-threading within WGCNA; at present this call is necessary. Caution: skip this line if you run RStudio or other third-party R environments.
```{r}
enableWGCNAThreads()
```

Load the data saved in the first part.
```{r}
lnames_wild = load(file = "wild-01-dataInput.RData")
```

#### Step-by-step construction of the gene network and identification of modules

#### Choosing the soft-thresholding power: analysis of network topology
Choose a set of soft-thresholding powers.
```{r}
powers_wild = c(c(1:10), seq(from = 12, to=20, by=2))
```

Call the network topology analysis function.
```{r}
# Signed network
sftWild = pickSoftThreshold(datExprWild, powerVector = powers_wild, networkType = "signed", verbose = 5)
```

Plot the results.
```{r}
pdf(file = "3_wild_scale_independence_and_connectivity.pdf", width = 9, height = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sftWild$fitIndices[,1], -sign(sftWild$fitIndices[,3])*sftWild$fitIndices[,2],
  xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",type="n",
  main = paste("Scale independence"))
text(sftWild$fitIndices[,1], -sign(sftWild$fitIndices[,3])*sftWild$fitIndices[,2],
  labels=powers_wild,cex=cex1,col="red")

# R^2 cutoff = 0.8460, which corresponds to a power (beta) of 14
abline(h=0.8460,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sftWild$fitIndices[,1], sftWild$fitIndices[,5],
  xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
  main = paste("Mean connectivity"))
text(sftWild$fitIndices[,1], sftWild$fitIndices[,5], labels=powers_wild, cex=cex1, col="red")
dev.off()
```

#### Co-expression similarity and adjacency
Calculate the adjacencies, using the soft thresholding power of 14.
```{r}
softPower_wild = 14
adjacency_wild = adjacency(datExprWild, power = softPower_wild, type = "signed")
```

#### Topological Overlap Matrix (TOM)
Turn adjacency into topological overlap.
```{r}
TOM_wild = TOMsimilarity(adjacency_wild, TOMType = "unsigned")
```

Calculate the corresponding dissimilarity.
```{r}
dissTOM_wild = 1-TOM_wild
```

#### Clustering using TOM
Call the hierarchical clustering function.
```{r}
geneTree_wild = hclust(as.dist(dissTOM_wild), method = "average")
```

Plot the resulting clustering tree (dendrogram).
```{r}
pdf(file = "4_wild_gene_clustering_on_TOM-based_dissimilarity.pdf", width = 12, height = 9)
sizeGrWindow(12,9)
plot(geneTree_wild, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()
```

Set the minimum module size.
```{r}
minModuleSize_wild = 30
```

Module identification using dynamic tree cut.
```{r}
dynamicMods_wild = cutreeDynamic(dendro = geneTree_wild, distM = dissTOM_wild,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize_wild)
table(dynamicMods_wild)
```

Convert numeric labels into colors.
```{r}
dynamicColors_wild = labels2colors(dynamicMods_wild)
table(dynamicColors_wild)
```

Plot the dendrogram and colors underneath.
```{r}
pdf(file = "5_wild_gene_dendrogram_and_module_colors.pdf", width = 8, height = 6)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree_wild, dynamicColors_wild, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
```

#### Merging of modules whose expression profiles are very similar
To quantify co-expression similarity of entire modules, calculate their eigengenes and cluster them based on their correlation.

Calculate eigengenes.
```{r}
MEList_wild = moduleEigengenes(datExprWild, colors = dynamicColors_wild)
MEs_wild = MEList_wild$eigengenes
```

Calculate dissimilarity of module eigengenes.
```{r}
MEDiss_wild = 1-cor(MEs_wild)
```

Cluster module eigengenes.
```{r}
METree_wild = hclust(as.dist(MEDiss_wild), method = "average")
```

Plot the result.
```{r}
pdf(file = "6_wild_clustering_of_module_eigengenes.pdf", width = 7, height = 6)
sizeGrWindow(7, 6)
plot(METree_wild, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# I chose a height cut of 0.25, corresponding to a correlation of 0.75, to merge
MEDissThres_wild = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres_wild, col = "red")
dev.off()
```

Call an automatic merging function.
```{r}
merge_wild = mergeCloseModules(datExprWild, dynamicColors_wild, cutHeight = MEDissThres_wild, verbose = 3)
```

The merged module colors.
```{r}
mergedColors_wild = merge_wild$colors
```

Eigengenes of the new merged modules.
```{r}
mergedMEs_wild = merge_wild$newMEs
```

To see what the merging did to the module colors, plot the gene dendrogram again, with the original and merged module colors underneath.
```{r}
pdf(file = "7_wild_gene_dendgrogram_original_and_merged.pdf", width = 12, height = 9)
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree_wild, cbind(dynamicColors_wild, mergedColors_wild),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
```

From now on, this analysis will use the merged module colors in mergedColors_wild.

Rename to moduleColors_wild.
```{r}
moduleColors_wild = mergedColors_wild
```

Construct numerical labels corresponding to the colors.
```{r}
colorOrder_wild = c("grey", standardColors(50))
moduleLabels_wild = match(moduleColors_wild, colorOrder_wild)-1
MEs_wild = mergedMEs_wild
write.csv(x = MEs_wild, file = "7_2_module_eigengenes_wild.csv")
```

Save module colors and labels for use in subsequent parts.
```{r}
save(MEs_wild, moduleLabels_wild, moduleColors_wild, geneTree_wild, file = "wild-02-networkConstruction-stepByStep.RData")
```

#### WGCNA Part 2

#### Preliminaries: set up the R session and load results of previous parts
The following setting is important, do not omit.
```{r}
options(stringsAsFactors = FALSE)
```

Load the expression and trait data saved in the first part.
```{r}
lnames_wild = load(file = "wild-01-dataInput.RData")
```

Load network data saved in the second part.
```{r}
lnames_wild = load(file = "wild-02-networkConstruction-stepByStep.RData")
```

#### Relate modules to external traits

#### Quantify module-trait associations
Identify modules that are significantly associated with the measured traits (sulfidic vs. non-sulfidic, drainage). To do this, correlate eigengenes with external traits.

Define numbers of genes and samples.
```{r}
nGenes_wild = ncol(datExprWild)
nSamples_wild = nrow(datExprWild)

# Use the MEs calculated after merging that were saved in the previous script
MEs_wild = orderMEs(MEs_wild)
moduleTraitCor_wild = cor(MEs_wild, datTraitsWild, use = "p")
moduleTraitPvalue_wild = corPvalueStudent(moduleTraitCor_wild, nSamples_wild)
```

Color code associations between modules and traits using their correlation values.
```{r}
pdf(file = "8_wild_correlations_and_p_values.pdf", width = 10, height = 16)
sizeGrWindow(10,16)
# Will display correlations and their p-values
textMatrix_wild = paste(signif(moduleTraitCor_wild, 2), "\n(",
signif(moduleTraitPvalue_wild, 1), ")", sep = "")
dim(textMatrix_wild) = dim(moduleTraitCor_wild)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor_wild,
xLabels = names(datTraitsWild),
yLabels = names(MEs_wild),
ySymbols = names(MEs_wild),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix_wild,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()
```

#### Gene relationship to trait and important modules: Gene Significance and Module Membership

#### (1) Comparing ecotypes (sulfidic habitat to non-sulfidic habitat)
```{r}
# Define variable habitat_wild containing the NS_vs_S_wild column of datTraitsWild
habitat_wild = as.data.frame(datTraitsWild$NS_vs_S_wild)
names(habitat_wild) = "habitat_wild"

# Names (colors) of the modules
modNames_wild = substring(names(MEs_wild), 3)

geneModuleMembership_wild = as.data.frame(cor(datExprWild, MEs_wild, use = "p"))
MMPvalue_wild = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_wild), nSamples_wild))

names(geneModuleMembership_wild) = paste("MM", modNames_wild, sep="")
names(MMPvalue_wild) = paste("p.MM", modNames_wild, sep="")

geneTraitSignificance_habitat_wild = as.data.frame(cor(datExprWild, habitat_wild, use = "p"))
GSPvalue_habitat_wild = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_habitat_wild), nSamples_wild))

names(geneTraitSignificance_habitat_wild) = paste("GS.", names(habitat_wild), sep="")
names(GSPvalue_habitat_wild) = paste("p.GS.", names(habitat_wild), sep="")
```

#### Intramodular analysis: identifying genes with high GS and MM
Identify genes that have a high significance for habitat as well as high module membership in interesting modules.

Most positively correlated module with habitat.
```{r}
moduleTraitCor_wild2 <- as.data.frame(moduleTraitCor_wild)

# Sort by NS_vs_S_wild column in descending order
moduleTraitCor_wild2_desc_habitat <- moduleTraitCor_wild2[order(-moduleTraitCor_wild2$NS_vs_S_wild),]

# Print ordered modules
cat("Most positively correlated module with habitat:\n")
moduleTraitCor_wild2_desc_habitat

# Move row names to their own column called "Module"
moduleTraitCor_wild2_desc_habitat <- tibble::rownames_to_column(moduleTraitCor_wild2_desc_habitat, "Module")

# Pull out module color of most positively correlated module with habitat
positive_wild_habitat <- moduleTraitCor_wild2_desc_habitat$Module[1]

# Remove ME from string, left with only the name of the color
modulePositiveWildHabitat <- gsub("ME", "", positive_wild_habitat)

columnPositiveWildHabitat = match(modulePositiveWildHabitat, modNames_wild)
moduleGenesPositiveWildHabitat = moduleColors_wild==modulePositiveWildHabitat
```

Plot.
```{r}
pdf(file = "9_wild_habitat_positive_corr_module.pdf", width = 7, height = 7)
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership_wild[moduleGenesPositiveWildHabitat, columnPositiveWildHabitat]),
abs(geneTraitSignificance_habitat_wild[moduleGenesPositiveWildHabitat, 1]),
xlab = paste("Module Membership in", modulePositiveWildHabitat, "module"),
ylab = "Gene significance for habitat (NS vs S)",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = modulePositiveWildHabitat)
dev.off()
```

Most negatively correlated module with habitat.
```{r}
moduleTraitCor_wild2 <- as.data.frame(moduleTraitCor_wild)

# Sort by NS_vs_S_wild column in ascending order
moduleTraitCor_wild2_asc_habitat <- moduleTraitCor_wild2[order(moduleTraitCor_wild2$NS_vs_S_wild),]

# Print ordered modules
cat("Most negatively correlated module with habitat:\n")
moduleTraitCor_wild2_asc_habitat

# Move row names to their own column called "Module"
moduleTraitCor_wild2_asc_habitat <- tibble::rownames_to_column(moduleTraitCor_wild2_asc_habitat, "Module")

# Pull out module color of most negatively correlated module with habitat
negative_wild_habitat <- moduleTraitCor_wild2_asc_habitat$Module[1]

# Remove ME from string, left with only the name of the color
moduleNegativeWildHabitat <- gsub("ME", "", negative_wild_habitat)

columnNegativeWildHabitat = match(moduleNegativeWildHabitat, modNames_wild)
moduleGenesNegativeWildHabitat = moduleColors_wild==moduleNegativeWildHabitat
```

Plot.
```{r}
pdf(file = "9_wild_habitat_negative_corr_module.pdf", width = 7, height = 7)
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership_wild[moduleGenesNegativeWildHabitat, columnNegativeWildHabitat]),
abs(geneTraitSignificance_habitat_wild[moduleGenesNegativeWildHabitat, 1]),
xlab = paste("Module Membership in", moduleNegativeWildHabitat, "module"),
ylab = "Gene significance for habitat (NS vs S)",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleNegativeWildHabitat)
dev.off()
```

#### Summary output of network analysis results
Import annotation file. See `Table S2` in the `README.md`.
```{r}
annot = read.csv(file = "Table\ S2")
```

Calculate the number of genes without annotation--in this case, the number of genes without a BLAST hit in `Table S2`.
```{r}
genesWild = names(datExprWild)
genes2annotWild = match(genesWild, annot$geneID)

cat("This is the number of genes without annotation:\n")
sum(is.na(genes2annotWild))
```

```{r}
# Create the starting data frame.
geneInfoHabitatWild0 = data.frame(geneID = genesWild,
geneName = annot$geneName[genes2annotWild],
subjectSequenceID = annot$SubjectSequenceID[genes2annotWild],
proteinAnnotations = annot$ProteinAnnotations[genes2annotWild],
moduleColors = moduleColors_wild,
geneTraitSignificance_habitat_wild,
GSPvalue_habitat_wild)

# Order modules by their significance for habitat
modOrder_habitat_wild = order(-abs(cor(MEs_wild, habitat_wild, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership_wild))
{
oldNames_habitat_wild = names(geneInfoHabitatWild0)
geneInfoHabitatWild0 = data.frame(geneInfoHabitatWild0, geneModuleMembership_wild[, modOrder_habitat_wild[mod]],
MMPvalue_wild[, modOrder_habitat_wild[mod]])
names(geneInfoHabitatWild0) = c(oldNames_habitat_wild, paste("MM.", modNames_wild[modOrder_habitat_wild[mod]], sep=""),
paste("p.MM.", modNames_wild[modOrder_habitat_wild[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder_habitat_wild = order(geneInfoHabitatWild0$moduleColors, -abs(geneInfoHabitatWild0$GS.habitat_wild))
geneInfoHabitatWild = geneInfoHabitatWild0[geneOrder_habitat_wild, ]
```

Write data frame into a text-format spreadsheet.
```{r}
write.csv(geneInfoHabitatWild, file = "10_geneInfoHabitatWild.csv")
```

#### Replace NAs in 10_geneInfoHabitatWild.csv with their gene name/locus IDs from the *P. mexicana* GFF file
NAs are from genes that do not have a BLAST hit.

Read in annotations pulled from the *P. mexicana* GFF file. See `PmexGeneNameMatching.csv` in the `README.md`.
```{r}
PmexGeneNameMatching <- read.csv(file = "PmexGeneNameMatching.csv")
```

Make row 1 the headers for the columns.
```{r}
names(PmexGeneNameMatching) <- as.matrix(PmexGeneNameMatching[1, ])
PmexGeneNameMatching <- PmexGeneNameMatching[-1, ]
PmexGeneNameMatching[] <- lapply(PmexGeneNameMatching, function(x) type.convert(as.character(x)))
```

Merge `PmexGeneNameMatching` and `geneInfoHabitatWild` by gene ID.
```{r}
mergedHabitatWild <- merge(x = geneInfoHabitatWild, y = PmexGeneNameMatching, by.x = "geneID", by.y = "gene.ID", all.x = TRUE)
```

Select columns in a new order to replace NAs in the geneName column of mergedHabitatWild with gene.name.
```{r}
mergedHabitatWild2 <- select(mergedHabitatWild, geneID,gene.name,subjectSequenceID,proteinAnnotations,moduleColors,GS.habitat_wild,p.GS.habitat_wild,matches("MM.*"),matches("p.MM.*"))
```

Save as a CSV file.
```{r}
write.csv(mergedHabitatWild2, file = "11_geneInfoHabitatWild_NAs_replaced_with_LOCs.csv")
```

#### (2) Comparing drainages (Pichucalco, Puyacatengo, Tacotalpa)
***Pichucalco Drainage***<br>
```{r}
# Define variable pich_wild containing the Pich_vs_all_wild column of datTraitsWild
pich_wild = as.data.frame(datTraitsWild$Pich_vs_all_wild)
names(pich_wild) = "pich_wild"

# Names (colors) of the modules
modNames_wild = substring(names(MEs_wild), 3)

geneModuleMembership_wild = as.data.frame(cor(datExprWild, MEs_wild, use = "p"))
MMPvalue_wild = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_wild), nSamples_wild))

names(geneModuleMembership_wild) = paste("MM", modNames_wild, sep="")
names(MMPvalue_wild) = paste("p.MM", modNames_wild, sep="")

geneTraitSignificance_pich_wild = as.data.frame(cor(datExprWild, pich_wild, use = "p"))
GSPvalue_pich_wild = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_pich_wild), nSamples_wild))

names(geneTraitSignificance_pich_wild) = paste("GS.", names(pich_wild), sep="")
names(GSPvalue_pich_wild) = paste("p.GS.", names(pich_wild), sep="")
```

***Puyacatengo Drainage***<br>
```{r}
# Define variable puya_wild containing the Puya_vs_all_wild column of datTraitsWild
puya_wild = as.data.frame(datTraitsWild$Puya_vs_all_wild)
names(puya_wild) = "puya_wild"

# Names (colors) of the modules
modNames_wild = substring(names(MEs_wild), 3)

geneModuleMembership_wild = as.data.frame(cor(datExprWild, MEs_wild, use = "p"))
MMPvalue_wild = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_wild), nSamples_wild))

names(geneModuleMembership_wild) = paste("MM", modNames_wild, sep="")
names(MMPvalue_wild) = paste("p.MM", modNames_wild, sep="")

geneTraitSignificance_puya_wild = as.data.frame(cor(datExprWild, puya_wild, use = "p"))
GSPvalue_puya_wild = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_puya_wild), nSamples_wild))

names(geneTraitSignificance_puya_wild) = paste("GS.", names(puya_wild), sep="")
names(GSPvalue_puya_wild) = paste("p.GS.", names(puya_wild), sep="")
```

***Tacotalpa Drainage***<br>
```{r}
# Define variable taco_wild containing the Taco_vs_all_wild column of datTraitsWild
taco_wild = as.data.frame(datTraitsWild$Taco_vs_all_wild)
names(taco_wild) = "taco_wild"
# Names (colors) of the modules
modNames_wild = substring(names(MEs_wild), 3)

geneModuleMembership_wild = as.data.frame(cor(datExprWild, MEs_wild, use = "p"))
MMPvalue_wild = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_wild), nSamples_wild))

names(geneModuleMembership_wild) = paste("MM", modNames_wild, sep="")
names(MMPvalue_wild) = paste("p.MM", modNames_wild, sep="")

geneTraitSignificance_taco_wild = as.data.frame(cor(datExprWild, taco_wild, use = "p"))
GSPvalue_taco_wild = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_taco_wild), nSamples_wild))

names(geneTraitSignificance_taco_wild) = paste("GS.", names(taco_wild), sep="")
names(GSPvalue_taco_wild) = paste("p.GS.", names(taco_wild), sep="")
```

#### Intramodular analysis: identifying genes with high GS and MM
Identify genes that have a high significance for drainage as well as high module membership in interesting modules.

***Pichucalco Drainage***<br>
Most positively correlated module with Pichucalco drainage.
```{r}
moduleTraitCor_wild2 <- as.data.frame(moduleTraitCor_wild)

# Sort by Pich_vs_all_wild column in descending order
moduleTraitCor_wild2_desc_pich <- moduleTraitCor_wild2[order(-moduleTraitCor_wild2$Pich_vs_all_wild),]

# Print ordered modules
cat("Most positively correlated module with Pichucalco drainage:\n")
moduleTraitCor_wild2_desc_pich

# Move row names to their own column called "Module"
moduleTraitCor_wild2_desc_pich <- tibble::rownames_to_column(moduleTraitCor_wild2_desc_pich, "Module")

# Pull out module color of most positively correlated module with Pichucalco drainage
positive_wild_pich <- moduleTraitCor_wild2_desc_pich$Module[1]

# Remove ME from string, left with only the name of the color
modulePositiveWildPich <- gsub("ME", "", positive_wild_pich)

columnPositiveWildPich = match(modulePositiveWildPich, modNames_wild)
moduleGenesPositiveWildPich = moduleColors_wild==modulePositiveWildPich
```

Plot.
```{r}
pdf(file = "9_wild_pich_positive_corr_module.pdf", width = 7, height = 7)
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership_wild[moduleGenesPositiveWildPich, columnPositiveWildPich]),
abs(geneTraitSignificance_pich_wild[moduleGenesPositiveWildPich, 1]),
xlab = paste("Module Membership in", modulePositiveWildPich, "module"),
ylab = "Gene significance for habitat (NS vs S)",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = modulePositiveWildPich)
dev.off()
```

Most negatively correlated module with Pichucalco drainage.
```{r}
moduleTraitCor_wild2 <- as.data.frame(moduleTraitCor_wild)

# Sort by Pich_vs_all_wild column in ascending order
moduleTraitCor_wild2_asc_pich <- moduleTraitCor_wild2[order(moduleTraitCor_wild2$Pich_vs_all_wild),]

# Print ordered modules
cat("Most negatively correlated module with Pichucalco drainage:\n")
moduleTraitCor_wild2_asc_pich

# Move row names to their own column called "Module"
moduleTraitCor_wild2_asc_pich <- tibble::rownames_to_column(moduleTraitCor_wild2_asc_pich, "Module")

# Pull out module color of most negatively correlated module with Pichucalco drainage
negative_wild_pich <- moduleTraitCor_wild2_asc_pich$Module[1]

# Remove ME from string, left with only the name of the color
moduleNegativeWildPich <- gsub("ME", "", negative_wild_pich)

columnNegativeWildPich = match(moduleNegativeWildPich, modNames_wild)
moduleGenesNegativeWildPich = moduleColors_wild==moduleNegativeWildPich
```

Plot.
```{r}
pdf(file = "9_wild_pich_negative_corr_module.pdf", width = 7, height = 7)
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership_wild[moduleGenesNegativeWildPich, columnNegativeWildPich]),
abs(geneTraitSignificance_pich_wild[moduleGenesNegativeWildPich, 1]),
xlab = paste("Module Membership in", moduleNegativeWildPich, "module"),
ylab = "Gene significance for habitat (NS vs S)",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleNegativeWildPich)
dev.off()
```

***Puyacatengo Drainage***<br>
Most positively correlated module with Puyacatengo drainage.
```{r}
moduleTraitCor_wild2 <- as.data.frame(moduleTraitCor_wild)

# Sort by Puya_vs_all_wild column in descending order
moduleTraitCor_wild2_desc_puya <- moduleTraitCor_wild2[order(-moduleTraitCor_wild2$Puya_vs_all_wild),]

# Print ordered modules
cat("Most positively correlated module with Puyacatengo drainage:\n")
moduleTraitCor_wild2_desc_puya

# Move row names to their own column called "Module"
moduleTraitCor_wild2_desc_puya <- tibble::rownames_to_column(moduleTraitCor_wild2_desc_puya, "Module")

# Pull out module color of most positively correlated module with Puyacatengo drainage
positive_wild_puya <- moduleTraitCor_wild2_desc_puya$Module[1]

# Remove ME from string, left with only the name of the color
modulePositiveWildPuya <- gsub("ME", "", positive_wild_puya)

columnPositiveWildPuya = match(modulePositiveWildPuya, modNames_wild)
moduleGenesPositiveWildPuya = moduleColors_wild==modulePositiveWildPuya
```

Plot.
```{r}
pdf(file = "9_wild_puya_positive_corr_module.pdf", width = 7, height = 7)
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership_wild[moduleGenesPositiveWildPuya, columnPositiveWildPuya]),
abs(geneTraitSignificance_puya_wild[moduleGenesPositiveWildPuya, 1]),
xlab = paste("Module Membership in", modulePositiveWildPuya, "module"),
ylab = "Gene significance for habitat (NS vs S)",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = modulePositiveWildPuya)
dev.off()
```

Most negatively correlated module with Puyacatengo drainage.
```{r}
moduleTraitCor_wild2 <- as.data.frame(moduleTraitCor_wild)
# Sort by Puya_vs_all_wild column in ascending order
moduleTraitCor_wild2_asc_puya <- moduleTraitCor_wild2[order(moduleTraitCor_wild2$Puya_vs_all_wild),]
# Print ordered modules
cat("\n")
cat("Most negatively correlated module with Puyacatengo drainage:\n")
moduleTraitCor_wild2_asc_puya
cat("\n")
# Move row names to their own column called "Module"
moduleTraitCor_wild2_asc_puya <- tibble::rownames_to_column(moduleTraitCor_wild2_asc_puya, "Module")
#colnames(moduleTraitCor_wild2_asc_puya)
# Pull out module color of most negatively correlated module with Puyacatengo drainage
negative_wild_puya <- moduleTraitCor_wild2_asc_puya$Module[1]
#negative_wild_puya
# Remove ME from string, left with only the name of the color
moduleNegativeWildPuya <- gsub("ME", "", negative_wild_puya)
#moduleNegativeWildPuya
columnNegativeWildPuya = match(moduleNegativeWildPuya, modNames_wild)
moduleGenesNegativeWildPuya = moduleColors_wild==moduleNegativeWildPuya
```

Plot.
```{r}
pdf(file = "9_wild_puya_negative_corr_module.pdf", width = 7, height = 7)
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership_wild[moduleGenesNegativeWildPuya, columnNegativeWildPuya]),
abs(geneTraitSignificance_puya_wild[moduleGenesNegativeWildPuya, 1]),
xlab = paste("Module Membership in", moduleNegativeWildPuya, "module"),
ylab = "Gene significance for habitat (NS vs S)",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleNegativeWildPuya)
dev.off()
```

***Tacotalpa Drainage***<br>
Most positively correlated module with Tacotalpa drainage.
```{r}
moduleTraitCor_wild2 <- as.data.frame(moduleTraitCor_wild)

# Sort by Taco_vs_all_wild column in descending order
moduleTraitCor_wild2_desc_taco <- moduleTraitCor_wild2[order(-moduleTraitCor_wild2$Taco_vs_all_wild),]

# Print ordered modules
cat("Most positively correlated module with Tacotalpa drainage:\n")
moduleTraitCor_wild2_desc_taco

# Move row names to their own column called "Module"
moduleTraitCor_wild2_desc_taco <- tibble::rownames_to_column(moduleTraitCor_wild2_desc_taco, "Module")

# Pull out module color of most positively correlated module with Tacotalpa drainage
positive_wild_taco <- moduleTraitCor_wild2_desc_taco$Module[1]

# Remove ME from string, left with only the name of the color
modulePositiveWildTaco <- gsub("ME", "", positive_wild_taco)

columnPositiveWildTaco = match(modulePositiveWildTaco, modNames_wild)
moduleGenesPositiveWildTaco = moduleColors_wild==modulePositiveWildTaco
```

Plot.
```{r}
pdf(file = "9_wild_taco_positive_corr_module.pdf", width = 7, height = 7)
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership_wild[moduleGenesPositiveWildTaco, columnPositiveWildTaco]),
abs(geneTraitSignificance_taco_wild[moduleGenesPositiveWildTaco, 1]),
xlab = paste("Module Membership in", modulePositiveWildTaco, "module"),
ylab = "Gene significance for habitat (NS vs S)",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = modulePositiveWildTaco)
dev.off()
```

Most negatively correlated module with Tacotalpa drainage.
```{r}
moduleTraitCor_wild2 <- as.data.frame(moduleTraitCor_wild)

# Sort by Taco_vs_all_wild column in ascending order
moduleTraitCor_wild2_asc_taco <- moduleTraitCor_wild2[order(moduleTraitCor_wild2$Taco_vs_all_wild),]

# Print ordered modules
cat("Most negatively correlated module with Tacotalpa drainage:\n")
moduleTraitCor_wild2_asc_taco

# Move row names to their own column called "Module"
moduleTraitCor_wild2_asc_taco <- tibble::rownames_to_column(moduleTraitCor_wild2_asc_taco, "Module")

# Pull out module color of most negatively correlated module with Tacotalpa drainage
negative_wild_taco <- moduleTraitCor_wild2_asc_taco$Module[1]

# Remove ME from string, left with only the name of the color
moduleNegativeWildTaco <- gsub("ME", "", negative_wild_taco)

columnNegativeWildTaco = match(moduleNegativeWildTaco, modNames_wild)
moduleGenesNegativeWildTaco = moduleColors_wild==moduleNegativeWildTaco
```

Plot.
```{r}
pdf(file = "9_wild_taco_negative_corr_module.pdf", width = 7, height = 7)
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership_wild[moduleGenesNegativeWildTaco, columnNegativeWildTaco]),
abs(geneTraitSignificance_taco_wild[moduleGenesNegativeWildTaco, 1]),
xlab = paste("Module Membership in", moduleNegativeWildTaco, "module"),
ylab = "Gene significance for habitat (NS vs S)",
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = moduleNegativeWildTaco)
dev.off()
```

#### Summary output of network analysis results
***Pichucalco Drainage***<br>
```{r}
# Create the starting data frame
geneInfoPichWild0 = data.frame(geneID = genesWild,
geneName = annot$geneName[genes2annotWild],
subjectSequenceID = annot$SubjectSequenceID[genes2annotWild],
proteinAnnotations = annot$ProteinAnnotations[genes2annotWild],
moduleColors = moduleColors_wild,
geneTraitSignificance_pich_wild,
GSPvalue_pich_wild)

# Order modules by their significance for the Pichucalco drainage
modOrder_pich_wild = order(-abs(cor(MEs_wild, pich_wild, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership_wild))
{
oldNames_pich_wild = names(geneInfoPichWild0)
geneInfoPichWild0 = data.frame(geneInfoPichWild0, geneModuleMembership_wild[, modOrder_pich_wild[mod]],
MMPvalue_wild[, modOrder_pich_wild[mod]])
names(geneInfoPichWild0) = c(oldNames_pich_wild, paste("MM.", modNames_wild[modOrder_pich_wild[mod]], sep=""),
paste("p.MM.", modNames_wild[modOrder_pich_wild[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder_pich_wild = order(geneInfoPichWild0$moduleColors, -abs(geneInfoPichWild0$GS.pich_wild))
geneInfoPichWild = geneInfoPichWild0[geneOrder_pich_wild, ]
```

Write data frame into a text-format spreadsheet.
```{r}
write.csv(geneInfoPichWild, file = "10_geneInfoPichWild.csv")
```

***Puyacatengo Drainage***<br>
```{r}
# Create the starting data frame
geneInfoPuyaWild0 = data.frame(geneID = genesWild,
geneName = annot$geneName[genes2annotWild],
subjectSequenceID = annot$SubjectSequenceID[genes2annotWild],
proteinAnnotations = annot$ProteinAnnotations[genes2annotWild],
moduleColors = moduleColors_wild,
geneTraitSignificance_puya_wild,
GSPvalue_puya_wild)

# Order modules by their significance for the Puyacatengo drainage
modOrder_puya_wild = order(-abs(cor(MEs_wild, puya_wild, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership_wild))
{
oldNames_puya_wild = names(geneInfoPuyaWild0)
geneInfoPuyaWild0 = data.frame(geneInfoPuyaWild0, geneModuleMembership_wild[, modOrder_puya_wild[mod]],
MMPvalue_wild[, modOrder_puya_wild[mod]])
names(geneInfoPuyaWild0) = c(oldNames_puya_wild, paste("MM.", modNames_wild[modOrder_puya_wild[mod]], sep=""),
paste("p.MM.", modNames_wild[modOrder_puya_wild[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder_puya_wild = order(geneInfoPuyaWild0$moduleColors, -abs(geneInfoPuyaWild0$GS.puya_wild))
geneInfoPuyaWild = geneInfoPuyaWild0[geneOrder_puya_wild, ]
```

Write data frame into a text-format spreadsheet.
```{r}
write.csv(geneInfoPuyaWild, file = "10_geneInfoPuyaWild.csv")
```

***Tacotalpa Drainage***<br>
```{r}
# Create the starting data frame
geneInfoTacoWild0 = data.frame(geneID = genesWild,
geneName = annot$geneName[genes2annotWild],
subjectSequenceID = annot$SubjectSequenceID[genes2annotWild],
proteinAnnotations = annot$ProteinAnnotations[genes2annotWild],
moduleColors = moduleColors_wild,
geneTraitSignificance_taco_wild,
GSPvalue_taco_wild)

# Order modules by their significance for the Tacotalpa drainage
modOrder_taco_wild = order(-abs(cor(MEs_wild, taco_wild, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership_wild))
{
oldNames_taco_wild = names(geneInfoTacoWild0)
geneInfoTacoWild0 = data.frame(geneInfoTacoWild0, geneModuleMembership_wild[, modOrder_taco_wild[mod]],
MMPvalue_wild[, modOrder_taco_wild[mod]])
names(geneInfoTacoWild0) = c(oldNames_taco_wild, paste("MM.", modNames_wild[modOrder_taco_wild[mod]], sep=""),
paste("p.MM.", modNames_wild[modOrder_taco_wild[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder_taco_wild = order(geneInfoTacoWild0$moduleColors, -abs(geneInfoTacoWild0$GS.taco_wild))
geneInfoTacoWild = geneInfoTacoWild0[geneOrder_taco_wild, ]
```

Write data frame into a text-format spreadsheet.
```{r}
write.csv(geneInfoTacoWild, file = "10_geneInfoTacoWild.csv")
```

#### Replace NAs in 10_geneInfo\*Wild.csv with their gene name/locus IDs from the P. mexicana GFF file
NAs are genes that do not have a BLAST hit.

***Pichucalco Drainage***<br>
Merge `PmexGeneNameMatching` and `geneInfoPichWild` by gene ID.
```{r}
mergedPichWild <- merge(x = geneInfoPichWild, y = PmexGeneNameMatching, by.x = "geneID", by.y = "gene.ID", all.x = TRUE)
```

Select columns in a new order to replace NAs in the geneName column of mergedPichWild with gene.name.
```{r}
mergedPichWild2 <- select(mergedPichWild, geneID,gene.name,subjectSequenceID,proteinAnnotations,moduleColors,GS.pich_wild,p.GS.pich_wild,matches("MM.*"),matches("p.MM.*"))
```

Save as a CSV file.
```{r}
write.csv(mergedPichWild2, file = "11_geneInfoPichWild_NAs_replaced_with_LOCs.csv")
```

***Puyacatengo Drainage***<br>
Merge `PmexGeneNameMatching` and `geneInfoPuyaWild` by gene ID.
```{r}
mergedPuyaWild <- merge(x = geneInfoPuyaWild, y = PmexGeneNameMatching, by.x = "geneID", by.y = "gene.ID", all.x = TRUE)
```

Select columns in a new order to replace NAs in the geneName column of mergedPuyaWild with gene.name.
```{r}
mergedPuyaWild2 <- select(mergedPuyaWild, geneID,gene.name,subjectSequenceID,proteinAnnotations,moduleColors,GS.puya_wild,p.GS.puya_wild,matches("MM.*"),matches("p.MM.*"))
```

Save as a CSV file.
```{r}
write.csv(mergedPuyaWild2, file = "11_geneInfoPuyaWild_NAs_replaced_with_LOCs.csv")
```

***Tacotalpa Drainage***<br>
Merge `PmexGeneNameMatching` and `geneInfoTacoWild` by gene ID.
```{r}
mergedTacoWild <- merge(x = geneInfoTacoWild, y = PmexGeneNameMatching, by.x = "geneID", by.y = "gene.ID", all.x = TRUE)
```

Select columns in a new order to replace NAs in the geneName column of mergedTacoWild with gene.name.
```{r}
mergedTacoWild2 <- select(mergedTacoWild, geneID,gene.name,subjectSequenceID,proteinAnnotations,moduleColors,GS.taco_wild,p.GS.taco_wild,matches("MM.*"),matches("p.MM.*"))
```

Save as a CSV file.
```{r}
write.csv(mergedTacoWild2, file = "11_geneInfoTacoWild_NAs_replaced_with_LOCs.csv")
```

#### Merge WGCNA output with TRRUST v2 human dataset to identify genes in significant modules that are transcription factors (TFs)
Read in TRRUST v2 human dataset in TSV format (from https://www.grnpedia.org/trrust/downloadnetwork.php), downloaded 6/3/2020 (Release note (2018.04.16)).
```{r}
trrust <- read.csv("trrust_rawdata.human.tsv", sep = "\t", header = FALSE)
```

Add column names.
```{r}
colnames(trrust) <- c("transcription_factor", "target_gene", "relationship", "pubmed_interaction_IDs")
```

#### (1) Comparing ecotypes (sulfidic habitat to non-sulfidic habitat)
Load CSV file with module information generated via WGCNA.
```{r}
geneInfoHabitatWild <- read.csv("11_geneInfoHabitatWild_NAs_replaced_with_LOCs.csv")
```

Pull out NCBI gene IDs from WGCNA dataset.
```{r}
# Split SubjectSequenceID by |
geneInfoHabitatWild_split <- cSplit(geneInfoHabitatWild, 'subjectSequenceID', sep = "|", type.convert = FALSE)

# Split SubjectSequenceID_3 by _
geneInfoHabitatWild_split2 <- cSplit(geneInfoHabitatWild_split, 'subjectSequenceID_3', sep = "_", type.convert = FALSE)

# Subset columns to keep subjectSequenceID_3_1 (the NCBI gene IDs)
geneInfoHabitatWild_subset <- select(geneInfoHabitatWild_split2, geneID,gene.name,subjectSequenceID_3_1,proteinAnnotations,moduleColors,GS.habitat_wild,p.GS.habitat_wild,matches("MM.*"),matches("p.MM.*"))

# Rename subjectSequenceID_3_1
colnames(geneInfoHabitatWild_subset)[3] <- "subjectSequenceID"
```

Merge WGCNA output and TRRUST dataframe, retain all lines in geneInfoHabitatWild_subset.
```{r}
geneInfoHabitatWild_merged <- merge(x = geneInfoHabitatWild_subset, y = trrust, by.x = "subjectSequenceID", by.y = "transcription_factor", all.x = TRUE, allow.cartesian = TRUE)
```

#### Pull out TFs in modules significantly correlated to habitat
```{r}
# Pull out Pearson correlation p-values for all modules, rounded to 1 significant figure
moduleTraitPvalue_wild2 <- as.data.frame(signif(moduleTraitPvalue_wild, 1))

# Sort by NS_vs_S_wild column in ascending order
moduleTraitPvalue_wild2_asc_habitat <- moduleTraitPvalue_wild2[order(moduleTraitPvalue_wild2$NS_vs_S_wild),]

# Print ordered modules
cat("Habitat modules ordered by p-value:\n")
moduleTraitPvalue_wild2_asc_habitat

moduleTraitPvalue_wild2_asc_habitat <- tibble::rownames_to_column(moduleTraitPvalue_wild2_asc_habitat, "Module")

# Pull out colors of modules with significant correlations (p < 0.05) to habitat
moduleTraitPvalue_wild2_habitat_sig_only <- moduleTraitPvalue_wild2_asc_habitat[moduleTraitPvalue_wild2_asc_habitat$NS_vs_S_wild < 0.05,]
sig_wild_habitat <- moduleTraitPvalue_wild2_habitat_sig_only$Module

# Remove ME from string, left with only the name of the color
sigModulesWildHabitat <- gsub("ME", "", sig_wild_habitat)

# Pull out TFs from significantly correlated modules
# Subset merged data frame by columns
# Add prefixes to modules significantly correlated with habitat
MM.sigModulesWildHabitat <- gsub("^", "MM.", sigModulesWildHabitat)
p.MM.sigModulesWildHabitat <- gsub("^", "p.MM.", sigModulesWildHabitat)

# Pull out only MM and p-values for significant columns
geneInfoHabitatWild_merged_sig <- select(geneInfoHabitatWild_merged, subjectSequenceID,geneID,gene.name,proteinAnnotations,moduleColors,all_of(MM.sigModulesWildHabitat),all_of(p.MM.sigModulesWildHabitat),target_gene,relationship,pubmed_interaction_IDs)

# Subset to only the genes in significant modules using the column moduleColors
geneInfoHabitatWild_merged_sig2 <- geneInfoHabitatWild_merged_sig %>% filter(moduleColors %in% (!!sigModulesWildHabitat)) %>% group_by(moduleColors)

# Keep only rows with target genes (a way to get only the TFs)
geneInfoHabitatWild_merged_sig3 <- geneInfoHabitatWild_merged_sig2[complete.cases(geneInfoHabitatWild_merged_sig2$target_gene),]

# Unique TFs
cat("Unique TFs for WildHabitat from significant modules:\n")
unique(geneInfoHabitatWild_merged_sig3$subjectSequenceID)

# Number of unique TFs
cat("Number of unique TFs for WildHabitat from significant modules:\n")
length(unique(geneInfoHabitatWild_merged_sig3$subjectSequenceID))

# Save subset merged dataframe as a CSV file (TFs are in the subjectSequenceID column)
write.csv(x = geneInfoHabitatWild_merged_sig3, file = "14_wild_habitat_all_sig_corr_modules_TFs.csv")
```

#### Make table of significant module correlations and p-values
Correlations.
```{r}
moduleTraitCor_wild2 <- as.data.frame(moduleTraitCor_wild)
moduleTraitCor_wild2 <- tibble::rownames_to_column(moduleTraitCor_wild2, "Module")
HabitatWild_sig_corr <- moduleTraitCor_wild2 %>% filter(Module %in% (!!sig_wild_habitat))

# Keep only modules and NS_vs_S_wild column
HabitatWild_sig_corr2 <- HabitatWild_sig_corr[,c(1,2)]
```

p-values.
```{r}
moduleTraitPvalue_wild3 <- tibble::rownames_to_column(moduleTraitPvalue_wild2, "Module")
HabitatWild_sig_p <- moduleTraitPvalue_wild3 %>% filter(Module %in% (!!sig_wild_habitat))

# Keep only modules and NS_vs_S_wild column
HabitatWild_sig_p2 <- HabitatWild_sig_p[,c(1,2)]
```

Combine correlations and p-values.
```{r}
HabitatWild_corr_p <- merge(x = HabitatWild_sig_corr2, y = HabitatWild_sig_p2, by.x = "Module", by.y = "Module", all = TRUE)

# Rename columns
names(HabitatWild_corr_p)[2] <- "Correlation"
names(HabitatWild_corr_p)[3] <- "p-value"

# Remove "ME" prefix in Module column
HabitatWild_corr_p[,1] <- sub("ME", "", HabitatWild_corr_p[,1])

# Add counts of unique TFs
HabitatWild_uniqueTFs_by_module <- geneInfoHabitatWild_merged_sig3 %>% group_by(moduleColors) %>% summarise(Unique_TFs = n_distinct(subjectSequenceID))
HabitatWild_corr_p_uniq <- merge(x = HabitatWild_corr_p, y = HabitatWild_uniqueTFs_by_module, by.x = "Module", by.y = "moduleColors", all.x = TRUE)

# Replace any NAs in Unique_TFs column with 0s
HabitatWild_corr_p_uniq$Unique_TFs[is.na(HabitatWild_corr_p_uniq$Unique_TFs)] <- 0

# Order by correlation in descending order
HabitatWild_corr_p_uniq_desc <- HabitatWild_corr_p_uniq[order(-HabitatWild_corr_p_uniq$Correlation),]

# Save as a CSV file
write.csv(x = HabitatWild_corr_p_uniq_desc, file = "15_wild_habitat_all_sig_modules_corr_pvalues_tfs.csv")
```

#### (2) Comparing drainages (Pichucalco, Puyacatengo, Tacotalpa)
***Pichucalco Drainage***<br>
Load CSV file with module information generated via WGCNA.
```{r}
geneInfoPichWild <- read.csv("11_geneInfoPichWild_NAs_replaced_with_LOCs.csv")
```

Pull out NCBI gene IDs from WGCNA dataset.
```{r}
# Split SubjectSequenceID by |
geneInfoPichWild_split <- cSplit(geneInfoPichWild, 'subjectSequenceID', sep = "|", type.convert = FALSE)

# Split SubjectSequenceID_3 by _
geneInfoPichWild_split2 <- cSplit(geneInfoPichWild_split, 'subjectSequenceID_3', sep = "_", type.convert = FALSE)

# Subset columns to keep subjectSequenceID_3_1 (the NCBI gene IDs)
geneInfoPichWild_subset <- select(geneInfoPichWild_split2, geneID,gene.name,subjectSequenceID_3_1,proteinAnnotations,moduleColors,GS.pich_wild,p.GS.pich_wild,matches("MM.*"),matches("p.MM.*"))

# Rename subjectSequenceID_3_1
colnames(geneInfoPichWild_subset)[3] <- "subjectSequenceID"
```

Merge WGCNA output and TRRUST dataframe, retain all lines in geneInfoPichWild_subset.
```{r}
geneInfoPichWild_merged <- merge(x = geneInfoPichWild_subset, y = trrust, by.x = "subjectSequenceID", by.y = "transcription_factor", all.x = TRUE, allow.cartesian = TRUE)
```

#### Pull out TFs in modules significantly correlated to the Pichucalco drainage
```{r}
# Pull out Pearson correlation p-values for all modules, rounded to 1 significant figure
moduleTraitPvalue_wild2 <- as.data.frame(signif(moduleTraitPvalue_wild, 1))

# Sort by Pich_vs_all_wild column in ascending order
moduleTraitPvalue_wild2_asc_pich <- moduleTraitPvalue_wild2[order(moduleTraitPvalue_wild2$Pich_vs_all_wild),]

# Print ordered modules
cat("Pich modules ordered by p-value:\n")
moduleTraitPvalue_wild2_asc_pich

moduleTraitPvalue_wild2_asc_pich <- tibble::rownames_to_column(moduleTraitPvalue_wild2_asc_pich, "Module")

# Pull out colors of modules with significant correlations (p < 0.05) to the Pichucalco drainage
moduleTraitPvalue_wild2_pich_sig_only <- moduleTraitPvalue_wild2_asc_pich[moduleTraitPvalue_wild2_asc_pich$Pich_vs_all_wild < 0.05,]
sig_wild_pich <- moduleTraitPvalue_wild2_pich_sig_only$Module

# Remove ME from string, left with only the name of the color
sigModulesWildPich <- gsub("ME", "", sig_wild_pich)

# Pull out TFs from significantly correlated modules
# Subset merged data frame by columns
# Add prefixes to modules significantly correlated with the Pichucalco drainage
MM.sigModulesWildPich <- gsub("^", "MM.", sigModulesWildPich)
p.MM.sigModulesWildPich <- gsub("^", "p.MM.", sigModulesWildPich)

# Pull out only MM and p-values for significant columns
# select(data.frame, desired columns)
geneInfoPichWild_merged_sig <- select(geneInfoPichWild_merged, subjectSequenceID,geneID,gene.name,proteinAnnotations,moduleColors,all_of(MM.sigModulesWildPich),all_of(p.MM.sigModulesWildPich),target_gene,relationship,pubmed_interaction_IDs)

# Subset to only the genes in significant modules using the column moduleColors
geneInfoPichWild_merged_sig2 <- geneInfoPichWild_merged_sig %>% filter(moduleColors %in% (!!sigModulesWildPich)) %>% group_by(moduleColors)

# Keep only rows with target genes (a way to get only the TFs)
geneInfoPichWild_merged_sig3 <- geneInfoPichWild_merged_sig2[complete.cases(geneInfoPichWild_merged_sig2$target_gene),]

# Unique TFs
cat("Unique TFs for WildPich from significant modules:\n")
unique(geneInfoPichWild_merged_sig3$subjectSequenceID)

# Number of unique TFs

cat("Number of unique TFs for WildPich from significant modules:\n")
length(unique(geneInfoPichWild_merged_sig3$subjectSequenceID))

# Save subset merged dataframe as a CSV file (TFs are in the subjectSequenceID column)
write.csv(x = geneInfoPichWild_merged_sig3, file = "14_wild_pich_all_sig_corr_modules_TFs.csv")
```

#### Make table of significant module correlations and p-values
Correlations.
```{r}
moduleTraitCor_wild2 <- as.data.frame(moduleTraitCor_wild)
moduleTraitCor_wild2 <- tibble::rownames_to_column(moduleTraitCor_wild2, "Module")
PichWild_sig_corr <- moduleTraitCor_wild2 %>% filter(Module %in% (!!sig_wild_pich))

# Keep only modules and Pich_vs_all_wild column
PichWild_sig_corr2 <- PichWild_sig_corr[,c(1,3)]
```

p-values.
```{r}
moduleTraitPvalue_wild3 <- tibble::rownames_to_column(moduleTraitPvalue_wild2, "Module")
PichWild_sig_p <- moduleTraitPvalue_wild3 %>% filter(Module %in% (!!sig_wild_pich))

# Keep only modules and Pich_vs_all_wild column
PichWild_sig_p2 <- PichWild_sig_p[,c(1,3)]
```

Combine correlations and p-values.
```{r}
PichWild_corr_p <- merge(x = PichWild_sig_corr2, y = PichWild_sig_p2, by.x = "Module", by.y = "Module", all = TRUE)

# Rename columns
names(PichWild_corr_p)[2] <- "Correlation"
names(PichWild_corr_p)[3] <- "p-value"

# Remove "ME" prefix in Module column
PichWild_corr_p[,1] <- sub("ME", "", PichWild_corr_p[,1])

# Add counts of unique TFs
PichWild_uniqueTFs_by_module <- geneInfoPichWild_merged_sig3 %>% group_by(moduleColors) %>% summarise(Unique_TFs = n_distinct(subjectSequenceID))
PichWild_corr_p_uniq <- merge(x = PichWild_corr_p, y = PichWild_uniqueTFs_by_module, by.x = "Module", by.y = "moduleColors", all.x = TRUE)

# Replace any NAs in Unique_TFs column with 0s
PichWild_corr_p_uniq$Unique_TFs[is.na(PichWild_corr_p_uniq$Unique_TFs)] <- 0

# Order by correlation in descending order
PichWild_corr_p_uniq_desc <- PichWild_corr_p_uniq[order(-PichWild_corr_p_uniq$Correlation),]

# Save as a CSV file
write.csv(x = PichWild_corr_p_uniq_desc, file = "15_wild_pich_all_sig_modules_corr_pvalues_tfs.csv")
```

***Puyacatengo Drainage***<br>
Load CSV file with module information generated via WGCNA.
```{r}
geneInfoPuyaWild <- read.csv("11_geneInfoPuyaWild_NAs_replaced_with_LOCs.csv")
```

Pull out NCBI gene IDs from WGCNA dataset.
```{r}
# Split SubjectSequenceID by |
geneInfoPuyaWild_split <- cSplit(geneInfoPuyaWild, 'subjectSequenceID', sep = "|", type.convert = FALSE)

# Split SubjectSequenceID_3 by _
geneInfoPuyaWild_split2 <- cSplit(geneInfoPuyaWild_split, 'subjectSequenceID_3', sep = "_", type.convert = FALSE)

# Subset columns to keep subjectSequenceID_3_1 (the NCBI gene IDs)
geneInfoPuyaWild_subset <- select(geneInfoPuyaWild_split2, geneID,gene.name,subjectSequenceID_3_1,proteinAnnotations,moduleColors,GS.puya_wild,p.GS.puya_wild,matches("MM.*"),matches("p.MM.*"))

# Rename subjectSequenceID_3_1
colnames(geneInfoPuyaWild_subset)[3] <- "subjectSequenceID"
```

Merge WGCNA output and TRRUST dataframe, retain all lines in geneInfoPuyaWild_subset.
```{r}
geneInfoPuyaWild_merged <- merge(x = geneInfoPuyaWild_subset, y = trrust, by.x = "subjectSequenceID", by.y = "transcription_factor", all.x = TRUE, allow.cartesian = TRUE)
```

#### Pull out TFs in modules significantly correlated to the Puyacatengo drainage
```{r}
# Pull out Pearson correlation p-values for all modules, rounded to 1 significant figure
moduleTraitPvalue_wild2 <- as.data.frame(signif(moduleTraitPvalue_wild, 1))

# Sort by Puya_vs_all_wild column in ascending order
moduleTraitPvalue_wild2_asc_puya <- moduleTraitPvalue_wild2[order(moduleTraitPvalue_wild2$Puya_vs_all_wild),]

# Print ordered modules
cat("Puya modules ordered by p-value:\n")
moduleTraitPvalue_wild2_asc_puya

moduleTraitPvalue_wild2_asc_puya <- tibble::rownames_to_column(moduleTraitPvalue_wild2_asc_puya, "Module")

# Pull out colors of modules with significant correlations (p < 0.05) to the Puyacatengo drainage
moduleTraitPvalue_wild2_puya_sig_only <- moduleTraitPvalue_wild2_asc_puya[moduleTraitPvalue_wild2_asc_puya$Puya_vs_all_wild < 0.05,]
sig_wild_puya <- moduleTraitPvalue_wild2_puya_sig_only$Module

# Remove ME from string, left with only the name of the color
sigModulesWildPuya <- gsub("ME", "", sig_wild_puya)

# Pull out TFs from significantly correlated modules
# Subset merged data frame by columns
# Add prefixes to modules significantly correlated with the Puyacatengo drainage
MM.sigModulesWildPuya <- gsub("^", "MM.", sigModulesWildPuya)
p.MM.sigModulesWildPuya <- gsub("^", "p.MM.", sigModulesWildPuya)

# Pull out only MM and p-values for significant columns
geneInfoPuyaWild_merged_sig <- select(geneInfoPuyaWild_merged, subjectSequenceID,geneID,gene.name,proteinAnnotations,moduleColors,all_of(MM.sigModulesWildPuya),all_of(p.MM.sigModulesWildPuya),target_gene,relationship,pubmed_interaction_IDs)

# Subset to only the genes in significant modules using the column moduleColors
geneInfoPuyaWild_merged_sig2 <- geneInfoPuyaWild_merged_sig %>% filter(moduleColors %in% (!!sigModulesWildPuya)) %>% group_by(moduleColors)

# Keep only rows with target genes (a way to get only the TFs)
geneInfoPuyaWild_merged_sig3 <- geneInfoPuyaWild_merged_sig2[complete.cases(geneInfoPuyaWild_merged_sig2$target_gene),]

# Unique TFs
cat("Unique TFs for WildPuya from significant modules:\n")
unique(geneInfoPuyaWild_merged_sig3$subjectSequenceID)

# Number of unique TFs
cat("Number of unique TFs for WildPuya from significant modules:\n")
length(unique(geneInfoPuyaWild_merged_sig3$subjectSequenceID))

# Save subset merged dataframe as a CSV file (TFs are in the subjectSequenceID column)
write.csv(x = geneInfoPuyaWild_merged_sig3, file = "14_wild_puya_all_sig_corr_modules_TFs.csv")
```

#### Make table of significant module correlations and p-values
Correlations.
```{r}
moduleTraitCor_wild2 <- as.data.frame(moduleTraitCor_wild)
moduleTraitCor_wild2 <- tibble::rownames_to_column(moduleTraitCor_wild2, "Module")
PuyaWild_sig_corr <- moduleTraitCor_wild2 %>% filter(Module %in% (!!sig_wild_puya))

# Keep only modules and Puya_vs_all_wild column
PuyaWild_sig_corr2 <- PuyaWild_sig_corr[,c(1,4)]
```

p-values.
```{r}
moduleTraitPvalue_wild3 <- tibble::rownames_to_column(moduleTraitPvalue_wild2, "Module")
PuyaWild_sig_p <- moduleTraitPvalue_wild3 %>% filter(Module %in% (!!sig_wild_puya))

# Keep only modules and Puya_vs_all_wild column
PuyaWild_sig_p2 <- PuyaWild_sig_p[,c(1,4)]
```

Combine correlations and p-values.
```{r}
PuyaWild_corr_p <- merge(x = PuyaWild_sig_corr2, y = PuyaWild_sig_p2, by.x = "Module", by.y = "Module", all = TRUE)

# Rename columns
names(PuyaWild_corr_p)[2] <- "Correlation"
names(PuyaWild_corr_p)[3] <- "p-value"

# Remove "ME" prefix in Module column
PuyaWild_corr_p[,1] <- sub("ME", "", PuyaWild_corr_p[,1])

# Add counts of unique TFs
PuyaWild_uniqueTFs_by_module <- geneInfoPuyaWild_merged_sig3 %>% group_by(moduleColors) %>% summarise(Unique_TFs = n_distinct(subjectSequenceID))
PuyaWild_corr_p_uniq <- merge(x = PuyaWild_corr_p, y = PuyaWild_uniqueTFs_by_module, by.x = "Module", by.y = "moduleColors", all.x = TRUE)

# Replace any NAs in Unique_TFs column with 0s
PuyaWild_corr_p_uniq$Unique_TFs[is.na(PuyaWild_corr_p_uniq$Unique_TFs)] <- 0

# Order by correlation in descending order
PuyaWild_corr_p_uniq_desc <- PuyaWild_corr_p_uniq[order(-PuyaWild_corr_p_uniq$Correlation),]

# Save as a CSV file
write.csv(x = PuyaWild_corr_p_uniq_desc, file = "15_wild_puya_all_sig_modules_corr_pvalues_tfs.csv")
```

***Tacotalpa Drainage***<br>
Load CSV file with module information generated via WGCNA.
```{r}
geneInfoTacoWild <- read.csv("11_geneInfoTacoWild_NAs_replaced_with_LOCs.csv")
```

Pull out NCBI gene IDs from WGCNA dataset.
```{r}
# Split SubjectSequenceID by |
geneInfoTacoWild_split <- cSplit(geneInfoTacoWild, 'subjectSequenceID', sep = "|", type.convert = FALSE)

# Split SubjectSequenceID_3 by _
geneInfoTacoWild_split2 <- cSplit(geneInfoTacoWild_split, 'subjectSequenceID_3', sep = "_", type.convert = FALSE)

# Subset columns to keep subjectSequenceID_3_1 (the NCBI gene IDs)
geneInfoTacoWild_subset <- select(geneInfoTacoWild_split2, geneID,gene.name,subjectSequenceID_3_1,proteinAnnotations,moduleColors,GS.taco_wild,p.GS.taco_wild,matches("MM.*"),matches("p.MM.*"))

# Rename subjectSequenceID_3_1
colnames(geneInfoTacoWild_subset)[3] <- "subjectSequenceID"
```

Merge WGCNA output and TRRUST dataframe, retain all lines in geneInfoTacoWild_subset.
```{r}
geneInfoTacoWild_merged <- merge(x = geneInfoTacoWild_subset, y = trrust, by.x = "subjectSequenceID", by.y = "transcription_factor", all.x = TRUE, allow.cartesian = TRUE)
```

#### Pull out TFs in modules significantly correlated to the Tacotalpa drainage
```{r}
# Pull out Pearson correlation p-values for all modules, rounded to 1 significant figure
moduleTraitPvalue_wild2 <- as.data.frame(signif(moduleTraitPvalue_wild, 1))

# Sort by Puya_vs_all_wild column in ascending order
moduleTraitPvalue_wild2_asc_taco <- moduleTraitPvalue_wild2[order(moduleTraitPvalue_wild2$Taco_vs_all_wild),]

# Print ordered modules
cat("Taco modules ordered by p-value:\n")
moduleTraitPvalue_wild2_asc_taco

moduleTraitPvalue_wild2_asc_taco <- tibble::rownames_to_column(moduleTraitPvalue_wild2_asc_taco, "Module")

# Pull out colors of modules with significant correlations (p < 0.05) to the Tacotalpa drainage
moduleTraitPvalue_wild2_taco_sig_only <- moduleTraitPvalue_wild2_asc_taco[moduleTraitPvalue_wild2_asc_taco$Taco_vs_all_wild < 0.05,]
sig_wild_taco <- moduleTraitPvalue_wild2_taco_sig_only$Module

# Remove ME from string, left with only the name of the color
sigModulesWildTaco <- gsub("ME", "", sig_wild_taco)

# Pull out TFs from significantly correlated modules
# Subset merged data frame by columns
# Add prefixes to modules significantly correlated with the Tacotalpa drainage
MM.sigModulesWildTaco <- gsub("^", "MM.", sigModulesWildTaco)
p.MM.sigModulesWildTaco <- gsub("^", "p.MM.", sigModulesWildTaco)

# Pull out only MM and p-values for significant columns
geneInfoTacoWild_merged_sig <- select(geneInfoTacoWild_merged, subjectSequenceID,geneID,gene.name,proteinAnnotations,moduleColors,all_of(MM.sigModulesWildTaco),all_of(p.MM.sigModulesWildTaco),target_gene,relationship,pubmed_interaction_IDs)

# Subset to only the genes in significant modules using the column moduleColors
geneInfoTacoWild_merged_sig2 <- geneInfoTacoWild_merged_sig %>% filter(moduleColors %in% (!!sigModulesWildTaco)) %>% group_by(moduleColors)

# Keep only rows with target genes (a way to get only the TFs)
geneInfoTacoWild_merged_sig3 <- geneInfoTacoWild_merged_sig2[complete.cases(geneInfoTacoWild_merged_sig2$target_gene),]

# Unique TFs
cat("Unique TFs for WildTaco from significant modules:\n")
unique(geneInfoTacoWild_merged_sig3$subjectSequenceID)

# Number of unique TFs
cat("Number of unique TFs for WildTaco from significant modules:\n")
length(unique(geneInfoTacoWild_merged_sig3$subjectSequenceID))

# Save subset merged dataframe as a CSV file (TFs are in the subjectSequenceID column)
write.csv(x = geneInfoTacoWild_merged_sig3, file = "14_wild_taco_all_sig_corr_modules_TFs.csv")
```

#### Make table of significant module correlations and p-values
Correlations.
```{r}
moduleTraitCor_wild2 <- as.data.frame(moduleTraitCor_wild)
moduleTraitCor_wild2 <- tibble::rownames_to_column(moduleTraitCor_wild2, "Module")
TacoWild_sig_corr <- moduleTraitCor_wild2 %>% filter(Module %in% (!!sig_wild_taco))

# Keep only modules and Taco_vs_all_wild column
TacoWild_sig_corr2 <- TacoWild_sig_corr[,c(1,5)]
```

p-values.
```{r}
moduleTraitPvalue_wild3 <- tibble::rownames_to_column(moduleTraitPvalue_wild2, "Module")
TacoWild_sig_p <- moduleTraitPvalue_wild3 %>% filter(Module %in% (!!sig_wild_taco))

# Keep only modules and Taco_vs_all_wild column
TacoWild_sig_p2 <- TacoWild_sig_p[,c(1,5)]
```

Combine correlations and p-values.
```{r}
TacoWild_corr_p <- merge(x = TacoWild_sig_corr2, y = TacoWild_sig_p2, by.x = "Module", by.y = "Module", all = TRUE)

# Rename columns
names(TacoWild_corr_p)[2] <- "Correlation"
names(TacoWild_corr_p)[3] <- "p-value"

# Remove "ME" prefix in Module column
TacoWild_corr_p[,1] <- sub("ME", "", TacoWild_corr_p[,1])

# Add counts of unique TFs
TacoWild_uniqueTFs_by_module <- geneInfoTacoWild_merged_sig3 %>% group_by(moduleColors) %>% summarise(Unique_TFs = n_distinct(subjectSequenceID))
TacoWild_corr_p_uniq <- merge(x = TacoWild_corr_p, y = TacoWild_uniqueTFs_by_module, by.x = "Module", by.y = "moduleColors", all.x = TRUE)

# Replace any NAs in Unique_TFs column with 0s
TacoWild_corr_p_uniq$Unique_TFs[is.na(TacoWild_corr_p_uniq$Unique_TFs)] <- 0

# Order by correlation in descending order
TacoWild_corr_p_uniq_desc <- TacoWild_corr_p_uniq[order(-TacoWild_corr_p_uniq$Correlation),]

# Save as a CSV file
write.csv(x = TacoWild_corr_p_uniq_desc, file = "15_wild_taco_all_sig_modules_corr_pvalues_tfs.csv")
```




























### Differential Gene Expression in EdgeR
Read in gene counts matrix from StringTie, above.
```{r}
wild_data <- read.csv(file = 'gene_count_matrix_no_STRG.csv', row.names = 1)
```

Create grouping factors to be incorporated later into the DGELists
```{r}
# S = sulfidic, NS = non-sulfidic
group = c(rep("NS", 6), rep("S", 6), rep("S", 5), rep("NS", 6), rep("NS", 6), rep("S", 6))
```

Create DGEList object.
```{r}
y <- DGEList(counts=wild_data, group=group)
```

Calculate dimensions of dataset before filtering.
```{r}
cat("y dimensions pre-filtering:", dim(y$counts), "genes and samples, respectively")
```

Filter data. 
```{r}
# Sort samples by library size
y_sorted <- y$samples[order(y$samples$lib.size) , ]
y_sorted
# Select lowest library size
lowest_library_size <- y_sorted$lib.size[1]
# Print lowest library size
cat("Smallest wild library size:", lowest_library_size)
# Calculate minimum counts per million (cpm) to have ~5 counts in the smallest library, as per edgeR manual recommendations
  # (5)/(lowest_library_size)=(x_wild)/(1,000,000)
  # Solve for x_wild
x <- (5*1000000)/lowest_library_size
# Print minimum cpm
cat("\n")
cat("Minimum wild cpm:", x)
# Samples must have a cpm of at least 1.144 (~5 counts in smallest library) and be present in at least 5 samples
keep <- rowSums(cpm(y) >= 1.144) >= 5
filtered_y <- y[keep, , keep.lib.sizes=FALSE]
# Library sizes post-filtering
filtered_y$samples
```

Dimensions of dataset after filtering.
```{r}
cat("filtered_y dimensions post-filtering:", dim(filtered_y$counts), "genes and samples, respectively")
```

Normalize for RNA composition (i.e., prevents genes from appearing falsely downregulated). Calculates normalization factors to scale raw library sizes, TMM is default method, to create effective library sizes.
```{r}
filtered_y <- calcNormFactors(filtered_y)
```

Create design matrix.
```{r}
design <- model.matrix(~0+group, data=y$samples) # 0 is the intercept
colnames(design) <- levels(y$samples$group)
design
```
Estimating Dispersions - CR method. Uses GLM instead of qCML method because testing multiple factors here. Also uses Cox-Reid profile-adjusted likelihood.
```{r}
# Estimate common dispersion and tagwise dispersions in one run
filtered_y <- estimateDisp(filtered_y, design)
# Common dispersions
filtered_y$common.dispersion
# Tagwise (gene-specific) dispersions
summary(filtered_y$tagwise.dispersion)
# Estimated prior degrees of freedom
filtered_y$prior.df
```

BCV plot.
```{r}
BCV <- plotBCV(filtered_y)
```

Save counts per million (CPM).
```{r}
cpm <- cpm(filtered_y)
write.csv(x=cpm, file = "3_wild_cpm.csv")
```

Measure differential gene expression using quasi-likelihood F-tests. Quasi-likelihood F-tests are better for bulk RNA-seq data because they have stricter error rate control, account for uncertainty in dispersion estimation.
```{r}
# Habitat 
# Make QL model representing the study design fitted to the data
fit <- glmQLFit(filtered_y, design)

# Make a contrast vector to compare SULFIDIC samples to NON-SULFIDIC samples
S_vs_NS_contrast <- makeContrasts(S-NS, levels = design)
S_vs_NS_contrast
qlf_habitat <- glmQLFTest(fit, contrast=S_vs_NS_contrast)

# n = Inf means report all genes
habitat_S_vs_NS <- topTags(qlf_habitat, n = Inf)
habitat_S_vs_NS <- data.frame(habitat_S_vs_NS)

# Write CSV of QL F-test results
write.csv(x=habitat_S_vs_NS, file = "~/Documents/WSU/RESEARCH/pmex_tf_biomed/8_edgeR_for_DGE/4_habitat_wild_qlf_S_vs_NS_RAW.csv")

# Move gene IDs to column 1
habitat_S_vs_NS <- setDT(habitat_S_vs_NS, keep.rownames = TRUE)[]

# Rename first column
names(habitat_S_vs_NS)[1] <- "geneID"
head(habitat_S_vs_NS)

# Read in PmexGeneNameMatching from Courtney Passow
PmexGeneNameMatching <- read.csv(file = "~/Documents/WSU/RESEARCH/pmex_tf_biomed/8_edgeR_for_DGE/scripts/PmexGeneNameMatching.csv")

# Make row 1 the headers for the columns
names(PmexGeneNameMatching) <- as.matrix(PmexGeneNameMatching[1, ])
PmexGeneNameMatching <- PmexGeneNameMatching[-1, ]
PmexGeneNameMatching[] <- lapply(PmexGeneNameMatching, function(x) type.convert(as.character(x)))
head(PmexGeneNameMatching)

# Merge habitat_wild_S_vs_NS and PmexGeneNameMatching by gene ID
merged_habitat_S_vs_NS <- merge(x = habitat_S_vs_NS, y = PmexGeneNameMatching, by.x = "geneID", by.y = "gene.ID", all.x = TRUE)

# Read in BLAST2GO file from Courtney Passow
anno <- read.csv(file = "~/Documents/WSU/RESEARCH/pmex_tf_biomed/8_edgeR_for_DGE/scripts/annotations_from_cpassow.csv", row.names = 1)

# Merge merged_habitat_wild_S_vs_NS and anno by gene ID
merged_x2_habitat_S_vs_NS <- merge(x = merged_habitat_S_vs_NS, y = anno, by.x = "geneID", by.y = "geneID", all.x = TRUE)

# Subset data frame to only contain relevant columns in an order that makes sense
merged_x2_habitat_S_vs_NS_subset <- merged_x2_habitat_S_vs_NS[,c(1,7,11,2:6,22)]

# Re-order by FDR
# FDR < 0.05 considered significant
sorted_habitat_S_vs_NS <- merged_x2_habitat_S_vs_NS_subset[order(FDR),]

# Write CSV of QL F-test results with annotations
write.csv(x=sorted_habitat_S_vs_NS, file = "~/Documents/WSU/RESEARCH/pmex_tf_biomed/8_edgeR_for_DGE/5_habitat_qlf_S_vs_NS_WITH_ANNOTATIONS.csv")
```


