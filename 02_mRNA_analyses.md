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

***FIX*** rgid !!!!!!!!!!!!!!!!
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
  {sample}_merged_sorted.bam \
  {sample}_sorted.bam \
  {sample}_sorted.bam
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

### Network Analysis (Weighted Gene Co-expression Network Analysis, WGCNA)
Subsequent analyses in this section were run using R (v.3.6.0).

Install packages.
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("WGCNA") 

install.packages("scales")
install.packages("tidyr")
install.packages("ggplot2")
```

Load packages.
```{r}
library(edgeR)
library(scales)
library(tidyr)
library(ggplot2)
library(WGCNA)
```

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

# Show library sizes and normalization factors
filtered_y_wild$samples
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

WGCNA. Note, a majority of the code and annotation was copied from https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/, and slightly modified to fit the design of this study.

#### Data input, cleaning, and pre-processing.
The following setting is important, do not omit
```{r}
options(stringsAsFactors = FALSE)
```
Make a data frame of log counts per million from EdgeR.
```{r}
counts_per_mil_wild <- cpm(filtered_y_wild, normalized.lib.sizes=TRUE, log=TRUE, prior.count=1)
counts_per_mil_wild <- as.data.frame(counts_per_mil_wild[,])
write.csv(counts_per_mil_wild, file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/log_cpm_wild.csv")
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
                    main = "Sample dendrogram and trait heatmap")
dev.off()
```

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
pdf(file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/4_wild_gene_clustering_on_TOM-based_dissimilarity.pdf", width = 12, height = 9)
sizeGrWindow(12,9)
plot(geneTree_wild, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()
```

# Set the minimum module size
minModuleSize_wild = 30

# Module identification using dynamic tree cut
dynamicMods_wild = cutreeDynamic(dendro = geneTree_wild, distM = dissTOM_wild,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize_wild)
table(dynamicMods_wild)

# Convert numeric labels into colors
dynamicColors_wild = labels2colors(dynamicMods_wild)
table(dynamicColors_wild)

# Plot the dendrogram and colors underneath
pdf(file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/5_wild_gene_dendrogram_and_module_colors.pdf", width = 8, height = 6)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree_wild, dynamicColors_wild, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

#### Merging of modules whose expression profiles are very similar -----------------------------------------------------------------------------------------------------------------

# To quantify co-expression similarity of entire modules, calculate their eigengenes and cluster them based on their correlation

# Calculate eigengenes
MEList_wild = moduleEigengenes(datExprWild, colors = dynamicColors_wild)
MEs_wild = MEList_wild$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss_wild = 1-cor(MEs_wild)

# Cluster module eigengenes
METree_wild = hclust(as.dist(MEDiss_wild), method = "average")

# Plot the result
pdf(file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/6_wild_clustering_of_module_eigengenes.pdf", width = 7, height = 6)
sizeGrWindow(7, 6)
plot(METree_wild, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# I chose a height cut of 0.25, corresponding to a correlation of 0.75, to merge
MEDissThres_wild = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres_wild, col = "red")
dev.off()

# Call an automatic merging function
merge_wild = mergeCloseModules(datExprWild, dynamicColors_wild, cutHeight = MEDissThres_wild, verbose = 3)

# The merged module colors
mergedColors_wild = merge_wild$colors

# Eigengenes of the new merged modules
mergedMEs_wild = merge_wild$newMEs

# To see what the merging did to the module colors, plot the gene dendrogram again, with the original and merged module colors underneath
pdf(file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/7_wild_gene_dendgrogram_original_and_merged.pdf", width = 12, height = 9)
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree_wild, cbind(dynamicColors_wild, mergedColors_wild),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# From now on, this analysis will use the merged module colors in mergedColors_wild

# Rename to moduleColors_wild
moduleColors_wild = mergedColors_wild

# Construct numerical labels corresponding to the colors
colorOrder_wild = c("grey", standardColors(50))
moduleLabels_wild = match(moduleColors_wild, colorOrder_wild)-1
MEs_wild = mergedMEs_wild

write.csv(x = MEs_wild, file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/7_2_module_eigengenes_wild.csv")

# Save module colors and labels for use in subsequent parts
save(MEs_wild, moduleLabels_wild, moduleColors_wild, geneTree_wild, file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/wild-02-networkConstruction-stepByStep.RData")







#### WGCNA Lab Exposure Samples -----------------------------------------------------------------------------------------------------------------

# Note, a majority of the code and annotation was copied from https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/, and slightly modified to fit the design of this study

# Install required packages
#BiocManager::install("WGCNA") 
#install.packages("flashClust", "cluster", lib="~/kerry.mcgowan/lib/R_libs/")

# Load required packages
library(WGCNA)
#library(flashClust)
#library(cluster)

#### Data input, cleaning, and pre-processing -----------------------------------------------------------------------------------------------------------------

#### Load expression data -----------------------------------------------------------------------------------------------------------------

# The following setting is important, do not omit
options(stringsAsFactors = FALSE)

# Make a data frame of log counts per million from edgeR
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html recommends "start with normalized counts (or RPKM/FPKM data) and log-transform them using log2(x+1)" for RNA-seq data
counts_per_mil_lab <- cpm(filtered_y_lab, normalized.lib.sizes=TRUE, log=TRUE, prior.count=1)
counts_per_mil_lab <- as.data.frame(counts_per_mil_lab[,])
write.csv(counts_per_mil_lab, file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/log_cpm_lab.csv")
#is.data.frame(counts_per_mil_lab)

# Remove samples BONGM2_ACTCAG_merged, BONGM5_ACAGTC_L002, LLUGL1_ACTGAC_AAGTCC_merged, LLUGH2_CTTGTA_L004, and LLUGH5_TGACAC_L006 from counts_per_mil_lab because they were shown to be mis-labeled samples in an admixture analysis done by J. Landers
cat("List sample names in counts_per_mil_lab to remove outliers:\n")
colnames(counts_per_mil_lab)
cat("\n")
counts_per_mil_lab2 <- counts_per_mil_lab[,-c(19,22,53,50,51)]
cat("List sample names in counts_per_mil_lab2 to check that the outliers were successfully removed:\n")
colnames(counts_per_mil_lab2)
cat("\n")
write.csv(counts_per_mil_lab2, file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/log_cpm_lab2.csv")

# Take a quick look at what is in the data set
cat("\n")
cat("counts_per_mil_lab2 dimensions:", dim(counts_per_mil_lab2), "genes and samples, respectively")
cat("\n")
cat("sample names in counts_per_mil_lab2:", names(counts_per_mil_lab2))
cat("\n")

# Transpose the data frame
datExprLab0 = as.data.frame(t(counts_per_mil_lab2[,]))

#### Check data for excessive missing values and identify outlier samples -----------------------------------------------------------------------------------------------------------------

# Check for genes and samples with too many missing values
gsg_lab = goodSamplesGenes(datExprLab0, verbose = 3)
cat("\n")
cat("Do all lab exposure genes pass cuts?", gsg_lab$allOK)
cat("\n")

# If above is TRUE, all genes passed cuts; if not, bad genes and samples are removed
if (!gsg_lab$allOK)
{
  # Optionally, print the gene and sample names that were removed
  if (sum(!gsg_lab$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExprLab0)[!gsg_lab$goodGenes], collapse = ", ")))
  if (sum(!gsg_lab$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExprLab0)[!gsg_lab$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data
  datExprLab0 = datExprLab0[gsg_lab$goodSamples, gsg_lab$goodGenes]
}

# Cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers
sampleTreeLab = hclust(dist(datExprLab0), method = "average")

# Plot the sample tree: open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small
pdf(file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/1_lab_sample_clustering.pdf", width = 12, height = 9)
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTreeLab, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Sample PSOGC1_TAGTGC_L001 looks like an outlier, so it is removed
# Choose a height cut that will remove the offending sample and use a branch cut at that height
# Plot a line to show the cut
abline(h = 300, col = "red")
dev.off()

# Determine cluster under the line
clust_lab = cutreeStatic(sampleTreeLab, cutHeight = 300, minSize = 10)
print(table(clust_lab), row.names=getOption("datatable.print.rownames"))

# Cluster 1 contains the samples we want to keep
# datExprLab contains the expression data that can be used for network analysis
keepSamplesLab = (clust_lab==1)
datExprLab = datExprLab0[keepSamplesLab, ]
nGenesLab = ncol(datExprLab)
nSamplesLab = nrow(datExprLab)

#### Input trait data -----------------------------------------------------------------------------------------------------------------

# Make a data frame of traits
# Remember, BONGM2_ACTCAG_merged, BONGM5_ACAGTC_L002, LLUGL1_ACTGAC_AAGTCC_merged, LLUGH2_CTTGTA_L004, LLUGH5_TGACAC_L006, and PSOGC1_TAGTGC_L001 were removed
cat("\n")
cat("Here are the row names of datExprLab to make a data frame of traits:\n")
rownames(datExprLab)
cat("\n")
# Habitat (meaning treatment--C, L, M, H):
  # non-sulfidic (NS) = 0, sulfidic (S) = 1
NS_vs_S_lab <- c(rep(0,6), rep(1,15), rep(0,6), rep(1,14), rep(0,6), rep(1,12), rep(0,5), rep(1,17))
# Drainage:
  # Ixtapangajoya (Ixta) drainage = 0, Puyacatengo (Puya) drainage = 1, Tacotalpa (Taco) drainage = 2 <--- OLD! JLK says WGCNA doesn't handle 0,1,2 coding well (only continuous or binary variables work from her experience)
#drainage_lab <- c(rep(2,23), rep(0,20), rep(1,21), rep(2,22)) <--- OLD!
  # Ixtapangajoya (Ixta) drainage = 1, all other drainages = 0
Ixta_vs_all_lab <- c(rep(0,21), rep(1,20), rep(0,18), rep(0,22))
  # Puyacatengo (Puya) drainage = 1, all other drainages = 0
Puya_vs_all_lab <- c(rep(0,21), rep(0,20), rep(1,18), rep(0,22))
  # Tacotalpa (Taco) drainage = 1, all other drainages = 0
Taco_vs_all_lab <- c(rep(1,21), rep(0,20), rep(0,18), rep(1,22))
# Ancestry:
  # non-sulfidic (NS) = 0, sulfidic (S) = 1
ancestry_lab <- c(rep(0,21), rep(0,20), rep(1,18), rep(1,22))
# datTraitsLab <- data.frame(NS_vs_S_lab, drainage_lab) <--- OLD!
datTraitsLab <- data.frame(NS_vs_S_lab, Ixta_vs_all_lab, Puya_vs_all_lab, Taco_vs_all_lab, ancestry_lab)
#print(datTraitsLab)
#cat("\n")

# Visualize how the clinical trails relate to the dendrogram
# Re-cluster samples
sampleTree2Lab = hclust(dist(datExprLab), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColorsLab = numbers2colors(datTraitsLab, signed = FALSE)

# Plot the sample dendrogram and the colors underneath
pdf(file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/2_lab_dendrogram_and_trait_heatmap.pdf")
plotDendroAndColors(sampleTree2Lab, traitColorsLab,
                    groupLabels = names(datTraitsLab),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

save(datExprLab, datTraitsLab, file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/lab-01-dataInput.RData")

#### Step-by-step network construction and module detection -----------------------------------------------------------------------------------------------------------------

#### Preliminaries: setting up the R session -----------------------------------------------------------------------------------------------------------------

# Allow multi-threading within WGCNA; at present this call is necessary
# Any error here may be ignored but you may want to update WGCNA if you see one
# Caution: skip this line if you run RStudio or other third-party R environments
enableWGCNAThreads()

# Load the data saved in the first part
lnames_lab = load(file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/lab-01-dataInput.RData")

# The variable lnames_lab contains the names of loaded variables
cat("\n")
cat("print lnames_lab:", lnames_lab)
cat("\n")

#### Step-by-step construction of the gene network and identification of modules -----------------------------------------------------------------------------------------------------------------

#### Choosing the soft-thresholding power: analysis of network topology -----------------------------------------------------------------------------------------------------------------

# Choose a set of soft-thresholding powers
powers_lab = c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
# Unsigned network
#sftLab = pickSoftThreshold(datExprLab, powerVector = powers_lab, verbose = 5)
# Signed network
sftLab = pickSoftThreshold(datExprLab, powerVector = powers_lab, networkType = "signed", verbose = 5)

# Plot the results
pdf(file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/3_lab_scale_independence_and_connectivity.pdf", width = 9, height = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sftLab$fitIndices[,1], -sign(sftLab$fitIndices[,3])*sftLab$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",type="n",
     main = paste("Scale independence"))
text(sftLab$fitIndices[,1], -sign(sftLab$fitIndices[,3])*sftLab$fitIndices[,2],
     labels=powers_lab,cex=cex1,col="red")
# This line corresponds to using an R^2 cut-off of h
# In Zhang and Horvath: Weighted Gene Co-Expression Network Analysis, they state "In our applications, we use the first parameter [power] value where saturation is reached as long as it is above [an R^2 value of] 0.8."
# Therefore, I chose a R^2 cutoff = 0.896, which corresponds to a power (beta) of 14
# See https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html #6 for how I chose the power for a signed network
abline(h=0.896,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sftLab$fitIndices[,1], sftLab$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sftLab$fitIndices[,1], sftLab$fitIndices[,5], labels=powers_lab, cex=cex1, col="red")
dev.off()

#### Co-expression similarity and adjacency -----------------------------------------------------------------------------------------------------------------

# Calculate the adjacencies, using the soft thresholding power of 14
softPower_lab = 14
# Unsigned network
#adjacency_lab = adjacency(datExprLab, power = softPower_lab)
# Signed network
adjacency_lab = adjacency(datExprLab, power = softPower_lab, type = "signed")

#### Topological Overlap Matrix (TOM) -----------------------------------------------------------------------------------------------------------------

# The TOM minimizes the effects of noise and spurious associations

# Turn adjacency into topological overlap
# Unsigned network
TOM_lab = TOMsimilarity(adjacency_lab, TOMType = "unsigned")
# Signed network, TOM keeps track of the sign of the adjacency between neighbors
# In signed networks, the signed and unsigned TOMs are practically identical (from "Signed vs. Unsigned Topological Overlap Matrix Technical Report")
#TOM_lab = TOMsimilarity(adjacency_lab, TOMType = "signed")

# Calculate the corresponding dissimilarity
dissTOM_lab = 1-TOM_lab

#### Clustering using TOM -----------------------------------------------------------------------------------------------------------------

# Use hierarchical clustering to produce a hierarchical clustering tree (dendrogram) of genes

# Call the hierarchical clustering function
geneTree_lab = hclust(as.dist(dissTOM_lab), method = "average")

# Plot the resulting clustering tree (dendrogram)
pdf(file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/4_lab_gene_clustering_on_TOM-based_dissimilarity.pdf", width = 12, height = 9)
sizeGrWindow(12,9)
plot(geneTree_lab, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

# Set the minimum module size
minModuleSize_lab = 30

# Module identification using dynamic tree cut
dynamicMods_lab = cutreeDynamic(dendro = geneTree_lab, distM = dissTOM_lab,
                                 deepSplit = 2, pamRespectsDendro = FALSE,
                                 minClusterSize = minModuleSize_lab)
table(dynamicMods_lab)

# Convert numeric labels into colors
dynamicColors_lab = labels2colors(dynamicMods_lab)
table(dynamicColors_lab)

# Plot the dendrogram and colors underneath
pdf(file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/5_lab_gene_dendrogram_and_module_colors.pdf", width = 8, height = 6)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree_lab, dynamicColors_lab, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

#### Merging of modules whose expression profiles are very similar -----------------------------------------------------------------------------------------------------------------

# To quantify co-expression similarity of entire modules, calculate their eigengenes and cluster them based on their correlation

# Calculate eigengenes
MEList_lab = moduleEigengenes(datExprLab, colors = dynamicColors_lab)
MEs_lab = MEList_lab$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss_lab = 1-cor(MEs_lab)

# Cluster module eigengenes
METree_lab = hclust(as.dist(MEDiss_lab), method = "average")

# Plot the result
pdf(file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/6_lab_clustering_of_module_eigengenes.pdf", width = 7, height = 6)
sizeGrWindow(7, 6)
plot(METree_lab, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
# I chose a height cut of 0.25, corresponding to a correlation of 0.75, to merge
MEDissThres_lab = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres_lab, col = "red")
dev.off()

# Call an automatic merging function
merge_lab = mergeCloseModules(datExprLab, dynamicColors_lab, cutHeight = MEDissThres_lab, verbose = 3)

# The merged module colors
mergedColors_lab = merge_lab$colors

# Eigengenes of the new merged modules
mergedMEs_lab = merge_lab$newMEs

# To see what the merging did to the module colors, plot the gene dendrogram again, with the original and merged module colors underneath
pdf(file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/7_lab_gene_dendgrogram_original_and_merged.pdf", width = 12, height = 9)
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree_lab, cbind(dynamicColors_lab, mergedColors_lab),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# From now on, this analysis will use the merged module colors in mergedColors_lab

# Rename to moduleColors_lab
moduleColors_lab = mergedColors_lab

# Construct numerical labels corresponding to the colors
colorOrder_lab = c("grey", standardColors(50))
moduleLabels_lab = match(moduleColors_lab, colorOrder_lab)-1
MEs_lab = mergedMEs_lab

write.csv(x = MEs_lab, file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/7_2_module_eigengenes_lab.csv")

# Save module colors and labels for use in subsequent parts
save(MEs_lab, moduleLabels_lab, moduleColors_lab, geneTree_lab, file = "/data/kelley/projects/kerry/pmex_tf_biomed/7_edgeR_WGCNA/lab-02-networkConstruction-stepByStep.RData")
