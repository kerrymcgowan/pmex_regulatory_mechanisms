### Raw Reads
Raw single-end 75 bp input and csRNA reads from gill tissues can be downloaded using the NCBI accession PRJNA743555. mRNA reads from the same gill tissues were used as part of the "full approach" for peak finding, described below, can be downloaded using the NCBI accession PRJNA290391.

All analyses use HOMER Tools (v.4.11) unless otherwise specified.

The variable `$sample` should be replaced with the full sample name.

### Trim
Trim FASTQ files containing raw reads.
```{bash}
homerTools trim \
  -3 AGATCGGAAGAGCACACGTCT \
  -mis 2 \
  -minMatchLength 4 \
  -min 20 \
  -q 20 \
  $sample.fastq.gz
```

### Map
Index the *Poecilia mexicana* reference genome using STAR (v.2.7.6). Reference genome and annotations in GFF format were downloaded from NCBI (Accession: GCF_001443325.1).
```{bash}
STAR \
  --runThreadN 7 \
  --runMode genomeGenerate \
  --genomeDir genomeDir \
  --genomeFastaFiles GCF_001443325.1_P_mexicana-1.0_genomic.fna \
  --sjdbGTFfile GCF_001443325.1_P_mexicana-1.0_genomic.gff \
  --sjdbOverhang 83 \
  --genomeSAindexNbases 13
```
Align trimmed input, csRNA, and mRNA FASTA files to the reference genome.
```{bash}
STAR \
	--genomeDir genomeDir \
	--runThreadN 24 \
	--readFilesIn $sample.fastq.gz.trimmed \
	--outFileNamePrefix $sample. \
	--outSAMstrandField intronMotif \
	--outMultimapperOrder Random \
	--outSAMmultNmax 1 \
	--outFilterMultimapNmax 10000 \
	--limitOutSAMoneReadBytes 10000000
```

### Make tag directories
Make an input, csRNA, and mRNA tag directories for each sample using the aligned reads.
```{bash}
makeTagDirectory \
  $sample\_tagDir \
  $sample.Aligned.out.sam \
  -genome $reference_genome \
  -checkGC \
  -unique \
  -fragLength 30 \
  -single
```

### Identify peaks
Convert the *Poecilia mexicana* GFF file to GTF format using Gffread (v.0.9.9).
```{bash}
gffread \
	-T GCF_001443325.1_P_mexicana-1.0_genomic.gff \
	-o GCF_001443325.1_P_mexicana-1.0_genomic.gtf
```

Peaks are clusters of transcription initiation, called using the "full approach" of the HOMER pipeline. This incorporates mRNA tag directories into the csRNA peak calling to help eliminate false positives like exomic contaminants.
```{bash}
findcsRNATSS.pl \
  $sample\_csRNA_tagDir  \
  -o $sample \
  -i $sample\_input_tagDir \
  -rna $sample\_rna_tagDir \
  -gtf $annotations \
  -genome GCF_001443325.1_P_mexicana-1.0_genomic.fna \
  -ntagThreshold 7
```

### Merge peaks
Produce a non-redundant list of peaks from all samples; 12 samples in this study.
```{bash}
mergePeaks \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	$sample.tss.txt \
	-strand \
	> merged.tss.txt
```

### Annotate peaks
Quantify csRNA read counts and annotate peaks.
```{bash}
annotatePeaks.pl \
	merged.tss.txt \
	GCF_001443325.1_P_mexicana-1.0_genomic.fna \
	-gtf GCF_001443325.1_P_mexicana-1.0_genomic.gtf \
	-strand + \
	-fragLength 1 \
	-raw \
	-d $sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	$sample\_csRNA_tagDir \
	> rawcounts.txt
```

### Calculate differential initiation
Find differentially regulated features using R (v.3.6.0) EdgeR (v.3.28.1). "H2S" and "noH2S" indicate if the sample came from a sulfidic or non-sulfidic environment, respectively.
```{bash}
getDiffExpression.pl \
	rawcounts.txt \
	noH2S noH2S H2S H2S H2S H2S noH2S noH2S noH2S noH2S H2S H2S \
	-edgeR \
	> ns_vs_s_diffOutput_edgeR.txt
```

Subsequent analyses in this section were run in RStudio (v.2022.02.1) using R (v.4.0.3).
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Install required packages.
```{r, message=FALSE, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("limma")

install.packages("xlsx")
install.packages("RColorBrewer")
install.packages("pheatmap")
install.packages("splitstackshape")
install.packages("factoextra")
install.packages("dplyr")
install.packages("tidyr")
```

Load required packages.
```{r, message = FALSE}
library(edgeR)
library(xlsx)
library(RColorBrewer)
library(pheatmap)
library(limma)
library(splitstackshape)
library(factoextra)
library(dplyr)
library(tidyr)
```

Read in EdgeR output from HOMER's getDiffExpression.pl.
**NOTE**: The output from getDiffExpression compared sulfidic to non-sulfidic samples, so for example, a peak with a positive log-fold change (logFC) is upregulated in sulfidic samples compared to non-sulfidic samples.
```{r}
diff_exp_eco <- read.delim("~/Documents/WSU/RESEARCH/pmex_gills_csrna_2021/09_homer_getdiffexpression/ns_vs_s_diffOutput_edgeR.txt", sep = "\t")
```

Replace column names with readable labels.
```{r}
colnames(diff_exp_eco) <- c("PeakID", "Chr", "Start", "End", "Strand", "Peak_Score", "Focus_Ratio/Region_Size", "Annotation", "Detailed_Annotation", "Distance_to_TSS", "Nearest_PromoterID", "Entrez_ID", "Nearest_Unigene", "Nearest_Refseq", "Nearest_Ensembl", "Gene_Name", "Gene_Alias", "Gene_Description", "Gene_Type", "MX04", "MX05", "MX31", "MX32", "MX46", "MX48", "MX52", "MX53", "MX60", "MX62", "MX76", "MX77", "noH2S.vs..H2S.Log2.Fold.Change", "noH2S.vs..H2S.p.value", "noH2S.vs..H2S.adj..p.value")
names(diff_exp_eco)
```

Make the first column the row names.
```{r}
rownames(diff_exp_eco) <- diff_exp_eco[,1]
diff_exp_eco <- diff_exp_eco[,-1]
```

Sort by adjusted p-value.
```{r}
sorted_eco <- diff_exp_eco[order(diff_exp_eco$`noH2S.vs..H2S.adj..p.value`),]

# Subset to significant peaks only (FDR < 0.05)
sig_eco <- subset(x = sorted_eco, subset = `noH2S.vs..H2S.adj..p.value` < 0.05)
#names(sig_eco)

# Subset to significantly upregulated peaks only (positive logFC, FDR < 0.05)
up_sig_eco <- sig_eco[sig_eco$noH2S.vs..H2S.Log2.Fold.Change > 0, ]

# Subset to significantly downregulated peaks only (negative logFC, FDR < 0.05)
down_sig_eco <- sig_eco[sig_eco$noH2S.vs..H2S.Log2.Fold.Change < 0, ]
```

Subset to only normalized read counts for all peaks and significant peaks only for plotting.
```{r}
# All peaks
subset_all_eco <- as.matrix(diff_exp_eco[,19:30])

# Significant peaks
subset_sig_eco <- as.matrix(sig_eco[,19:30])
```

Make a heatmap of all peaks and significant peaks only showing differential transcription initiation. Scale bars are z-scores.
```{r}
heat.colors <- rev(brewer.pal(6, "RdBu"))

# All peaks
# Plot top 10,000 most expressed genes
subset_all_eco2 <- subset_all_eco[order(-rowSums(subset_all_eco)),]
subset_all_eco3 <- subset_all_eco2[1:10000,]

eco_all_pheatmap <- pheatmap(subset_all_eco3, color = heat.colors, cluster_rows = T,
                             show_rownames = F, border_color=NA, fontsize = 10, 
                             scale="row", fontsize_row = 10, height=20)
eco_all_pheatmap

pdf("01_heatmap_all_ecotype.pdf")
print(eco_all_pheatmap)
dev.off()

# Significant peaks
eco_sig_pheatmap <- pheatmap(subset_sig_eco, color = heat.colors, cluster_rows = T, show_rownames = F,
         border_color=NA, fontsize = 10, scale="row",
         fontsize_row = 10, height=20)
eco_sig_pheatmap

pdf("01_heatmap_sig_only_ecotype.pdf")
print(eco_sig_pheatmap)
dev.off()
```

# Add PCA!!!

Merge significant peaks with H<sub>2</sub>S-related candidate list. The Sulfide Detox/Response Gene Set is relevant for this study. The NuclearRef and Broughton Gene Sets are reference sets.
```{r}
# Read in candidate list
cand_list <- read.csv(OXPHOS_Reference_and_Detox_Gene_IDS.csv', header = 1)

# Subset to only include Sulfide Detox/Response Gene Set
subset_cand <- subset(cand_list, Gene.Set == "Sulfide Detox/Response")

# Change name of sulfide:quinone oxidoreductase from "sqor" to "sqrdl" in Gene.name.from.accession column
subset_cand[subset_cand == "sqor"] <- "sqrdl"

# Merge H2S candidate list with significant csRNA peaks
merged_eco <- merge(x = sig_eco, y = subset_cand, by.x = "Gene_Name", by.y = "Gene.name.from.accession")
```

Merge with Blast annotations from Passow *et al.*
```{r}
pmex_anno <- read.csv('Table S2 - Poecilia mexicana Annotations.csv')

anno_eco <- merge(x = merged_eco, y = pmex_anno, by.x = "Gene_Name", by.y = "gene.name")

# Sort by logFC
anno_eco_sorted <- anno_eco[order(anno_eco$noH2S.vs..H2S.Log2.Fold.Change, decreasing = TRUE),]

write.csv(04_h2s_candidate_sig_peaks_intersection.csv')
```

Determine if the number of differentially expressed peaks in the Sulfide Detox/Response Gene Set is significant using a Fisher's Exact Test.

2 x 2 contingency table (Set = Sulfide Detox/Response Gene Set):

|                | Differential Expression  | No Differential Expression |
| -------------- | ------------------------ | -------------------------- |
| In Set         | 18                       | 346                        |
| Not in Set     | 1,488                    | 65,049                     |

```{r}
# Intersect non-significant peaks with the Sulfide Detox/Response Gene Set to complete the 2 x 2 contingency table above
# Subset to non-significant peaks only (FDR >= 0.05)
non_sig_eco <- subset(sorted_eco, `noH2S.vs..H2S.adj..p.value` >= 0.05)

# Merge H2S candidate list with non-significant csRNA peaks
merged_non_sig_eco <- merge(x = non_sig_eco, y = subset_cand, by.x = "Gene_Name", by.y = "Gene.name.from.accession")

# Make the 2 x 2 contingency table as a matrix
table <- matrix(c(18, 1488, 346, 65049),
                nrow = 2,
                dimnames = list(c('in_set', 'not_in_set'),
                                c('diff_exp', 'no_diff_exp')))
# Fisher's Exact Test for count data
fisher.test(table, alternative = "two.sided")
```

Intersect the list of significant peaks with Blast annotations.
```{r}
anno_sig_eco <- merge(x = sig_eco, y = pmex_anno, by.x = "Gene_Name", by.y = "gene.name", all.x = TRUE)
# Split Subject.sequence.ID to only contain the name of the human ortholog
anno_sig_eco <- cSplit(anno_sig_eco, "Subject.sequence.ID", sep = "|", type.convert = FALSE)
anno_sig_eco <- cSplit(anno_sig_eco, 'Subject.sequence.ID_3', sep = "_", type.convert = FALSE)
# Subset to desired columns
#names(anno_sig_eco)
anno_sig_eco <- anno_sig_eco[,c(1:33,48,45)]
write.csv(anno_sig_eco, file = '05_annotated_significant_peaks.csv')
```

Get a unique list of XM (messenger RNA) and XR (non-coding RNA) accessions for all peaks to be used in OmicsBox functional enrichment analyses.
```{r}
# Isolate XM and XR accessions
sorted_eco2 <- separate(sorted_eco, Nearest_Unigene, into = paste0("Nearest_Unigene", 1:2), sep = '-')

# Pull unique list of accessions
unigene_list_all <- as.data.frame(unique(sorted_eco2$Nearest_Unigene2))

# Remove NA
#sum(is.na(unigene_list_all$`unique(sorted_eco2$Nearest_Unigene2)`))
unigene_list_all <- unigene_list_all %>% drop_na()

# Remove pesudogenes (locus, or LOC, identifiers)
unigene_list_all <- unigene_list_all %>% filter(!grepl("LOC", `unique(sorted_eco2$Nearest_Unigene2)`))

# Remove tRNAs
unigene_list_all <- unigene_list_all %>% filter(!grepl("trna", `unique(sorted_eco2$Nearest_Unigene2)`))

# Save all XM peaks
unigene_list_all_XM <- as.data.frame(unigene_list_all[grepl("XM", unigene_list_all$`unique(sorted_eco2$Nearest_Unigene2)`),])
write.table(x = unigene_list_all_XM, file = '06_unigene_list_for_blast2go_all_XM_peaks.csv', sep="\t", row.names = FALSE, col.names = FALSE)

# Save all XR peaks
unigene_list_all_XR <- as.data.frame(unigene_list_all[grepl("XR", unigene_list_all$`unique(sorted_eco2$Nearest_Unigene2)`),])
write.table(x = unigene_list_all_XR, file = '06_unigene_list_for_blast2go_all_XR_peaks.csv', sep="\t", row.names = FALSE, col.names = FALSE)
```

### OmicsBox Functional Enrichment Analysis
Formerly Blast2GO. Pull FASTA sequences for all XM and XR accessions from above. `sequence.fasta` contains all of the RefSeq FASTA sequences for *Poecilia mexicana* downloaded from NCBI's Nucleotide database.
```{blast}
# Turn multi-line FASTA to single line FASTA
awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < sequence.fasta > sequence_single_line.fasta

# Subset to mRNAs corresponding to csRNA-seq peaks
while read line || [ -n "$line" ];
do
        gene=$(echo $line | awk '{print $1}')

        grep $gene -A 1 sequence_single_line.fasta >> nearest_unigenes_to_all_peaks_XM.fa

done < 06_unigene_list_for_blast2go_all_XM_peaks.txt

# Subset to ncRNAs corresponding to csRNA-seq peaks
while read line || [ -n "$line" ];
do
        gene=$(echo $line | awk '{print $1}')

        grep $gene -A 1 sequence_single_line.fasta >> nearest_unigenes_to_all_peaks_XR.fa

done < 06_unigene_list_for_blast2go_all_XR_peaks.txt
```

The resulting subset FASTA files were used as input for a functional enrichment analysis in OmicsBox.











