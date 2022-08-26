### Intersect Differentially Initiated Transcripts and Differentially Expressed Genes
Subsequent analyses in this section were run using R (v.4.0.3).

Install packages.
```{r}
install.packages("DescTools")
install.packages("ggrepel")
install.packages("ggpubr")
install.packages("mblm")
install.packages("RColorBrewer")
install.packages("tidyverse")
install.packages("viridis")
```

Load packages.
```{r}
library(DescTools)
library(ggrepel)
library(ggpubr)
library(mblm)
library(RColorBrewer)
library(tidyverse)
library(viridis)
```

#### Differentially Initiated (DI) peaks from csRNA-seq
Read in DI peaks (csRNA-seq) analyzed using edgeR.
```{r}
DI <- read.delim("ns_vs_s_diffOutput_edgeR.txt", sep = "\t")

# Rename columns
colnames(DI) <- c("PeakID", "Chr", "Start", "End", "Strand", "Peak_Score", "Focus_Ratio/Region_Size", "Annotation", "Detailed_Annotation", "Distance_to_TSS", "Nearest_PromoterID", "Entrez_ID", "Nearest_Unigene", "Nearest_Refseq", "Nearest_Ensembl", "Gene_Name", "Gene_Alias", "Gene_Description", "Gene_Type", "MX04", "MX05", "MX31", "MX32", "MX46", "MX48", "MX52", "MX53", "MX60", "MX62", "MX76", "MX77", "noH2S.vs..H2S.Log2.Fold.Change", "noH2S.vs..H2S.p.value", "noH2S.vs..H2S.adj..p.value")

# Split Annotation column, this is needed to later summarize the annotations of the peaks
DI_split <- DI %>% separate(col = Annotation, into = c("Annotation", "Anno_RefSeq_Accession"), sep = " ")

# Capitalize the first letter of all annotations
str_sub(DI_split$Annotation, 1, 1) <- str_sub(DI_split$Annotation, 1, 1) %>% str_to_upper()

#Sort by FDR (adjusted p-value)
DI_sorted <- DI_split[order(DI_split$`noH2S.vs..H2S.adj..p.value`),]

# Subset to significantly DI peaks only (FDR < 0.05)
DI_sig <- subset(x = DI_sorted, subset = `noH2S.vs..H2S.adj..p.value` < 0.05)
print(paste0("Number of significantly DI peaks: ", nrow(DI_sig)))
```

Summarize the annotations of all significantly DI peaks (i.e., what region of the genome the peaks map to).
```{r}
# Count unnanotated DI peaks (NAs)
print(paste0("Number of unannotated peaks: ", sum(is.na(DI_sig$Annotation))))

# Sum occurrences by genomic peak location
DI_sig_gen_loc <- as.data.frame(table(DI_sig$Annotation))
# Change column header
names(DI_sig_gen_loc)[1] <- "Annotation"
DI_sig_gen_loc
```

#### Differentially Expressed (DE) genes from RNA-seq
Read in DE genes (RNA-seq) analyzed using edgeR.
```{r}
DE <- read.csv("5_habitat_qlf_S_vs_NS_WITH_ANNOTATIONS.csv")

# Split SubjectSequenceID column to only include the gene name.
DE_split <- DE %>% 
                  separate(SubjectSequenceID, into = paste0("SubjectSequenceID", 1:4), sep = "[|_]") %>%
                  dplyr::select(-c(SubjectSequenceID1, SubjectSequenceID2, SubjectSequenceID4)) %>%
                  rename(SubjectSequenceID = SubjectSequenceID3)


# Subset to significantly DE genes only (FDR < 0.05)
DE_sig <- subset(x = DE_split, subset = `FDR` < 0.05)
print(paste0("Number of significantly DE genes: ", nrow(DE_sig)))
```

#### Intersect significantly DI peaks with DE genes
Note, we start with only significantly DI peaks but all DE genes, regardless of their significance.
```{r}
all_DI_DE <- merge(x = DI_sig, y = DE_split, by.x = "Gene_Name", by.y = "gene.name", all = TRUE)

# Subset columns in an order that makes sense
all_DI_DE_subset <- all_DI_DE %>% dplyr::select(geneID, PeakID, Chr, Start, End, Strand, Peak_Score, Annotation, Anno_RefSeq_Accession, Distance_to_TSS, Entrez_ID, SubjectSequenceID, Nearest_PromoterID, Nearest_Unigene, logFC, PValue, FDR, noH2S.vs..H2S.Log2.Fold.Change, noH2S.vs..H2S.p.value, noH2S.vs..H2S.adj..p.value, ProteinAnnotations)

# Rename columns
colnames(all_DI_DE_subset) <- c("RNA_GeneID", "csRNA_PeakID", "csRNA_Chr", "csRNA_Start", "csRNA_End", "csRNA_Strand", "csRNA_Peak_Score", "csRNA_Annotation", "csRNA_Anno_RefSeq_Accession", "csRNA_Distance_to_TSS", "Gene_Name_aka_Entrez_ID", "RNA_Subject_Sequence_ID", "csRNA_Nearest_Promoter_ID", "csRNA_Nearest_Unigene", "RNA_LogFC", "RNA_Pvalue", "RNA_FDR", "csRNA_LogFC", "csRNA_Pvalue", "csRNA_FDR", "RNA_Protein_Annotations")
```

Subset significantly DI peaks to significantly DE genes. 
```{r}
DI_sig_DE_sig <- subset(x = all_DI_DE_subset, `csRNA_FDR` < 0.05 & `RNA_FDR` < 0.05 )
print(paste0("Number of significantly DI peaks with significant DE: ", nrow(DI_sig_DE_sig)))
write.csv(x = DI_sig_DE_sig, file = "01_DI_sig_DE_sig.csv")
```

Summarize the annotations of significantly DI peaks that have significant DE (i.e., what region of the genome the peaks map to).
```{r}
# Count unoannotated DI peaks (NAs)
print(paste0("Number of unannotated peaks: ", sum(is.na(DI_sig_DE_sig$csRNA_Annotation))))

# Sum occurrences by genomic peak location
DI_sig_DE_sig_gen_loc <- as.data.frame(table(DI_sig_DE_sig$csRNA_Annotation))
# Change column header
names(DI_sig_DE_sig_gen_loc)[1] <- "Annotation"
# Add column of percentages
DI_sig_DE_sig_gen_loc_w_percents <- DI_sig_DE_sig_gen_loc %>% 
                                    mutate(Perc = `Freq` / sum(`Freq`)) %>% 
                                    arrange(Perc) %>% 
                                    mutate(Labels = scales::percent(Perc))
# Add column with data frame name
DI_sig_DE_sig_gen_loc_w_percents$df <- "DI_sig_DE_sig"
DI_sig_DE_sig_gen_loc_w_percents
```

G-test (log likelihood ratio) of goodness of fit for significantly DI peaks that have significant DE. Do the annotations fit the distribution of all DI peaks?
```{r}
# Total number of annotated peaks
#259 + 251 + 499 + 387 + 61

# Expected frequencies
# In the order: Exon, Intergenic, Intron, Promoter-TSS, TTS
expected_all_DI_peaks <- c((259/1457), (251/1457), (499/1457), (387/1457), (61/1457))
expected_all_DI_peaks
# Observed frequencies
observed_DI_sig_DE_sig <- c(111, 78, 190, 215, 24)
GTest(x = observed_DI_sig_DE_sig,
      p = expected_all_DI_peaks,
      correct = "none")
```

Subset significantly DI peaks to non-significantly DE genes.
```{r}
DI_sig_DE_nonsig <- subset(x = all_DI_DE_subset, `csRNA_FDR` < 0.05 & `RNA_FDR` > 0.05)
print(paste0("Number of significantly DI peaks without significant DE: ", nrow(DI_sig_DE_nonsig)))
write.csv(x = DI_sig_DE_nonsig, file = "01_DI_sig_DE_nonsig.csv")
```

Summarize the annotations of significantly DI peaks that have non-significant DE (i.e., what region of the genome the peaks map to).
```{r}
# Count unoannotated DI peaks (NAs)
print(paste0("Number of unannotated peaks: ", sum(is.na(DI_sig_DE_nonsig$csRNA_Annotation))))

# Sum occurrences by genomic peak location
DI_sig_DE_nonsig_gen_loc <- as.data.frame(table(DI_sig_DE_nonsig$csRNA_Annotation))
# Change column header
names(DI_sig_DE_nonsig_gen_loc)[1] <- "Annotation"
# Add column of percentages
DI_sig_DE_nonsig_gen_loc_w_percents <- DI_sig_DE_nonsig_gen_loc %>% 
                                    mutate(Perc = `Freq` / sum(`Freq`)) %>% 
                                    arrange(Perc) %>% 
                                    mutate(Labels = scales::percent(Perc))
# Add column with data frame name
DI_sig_DE_nonsig_gen_loc_w_percents$df <- "DI_sig_DE_nonsig"
DI_sig_DE_nonsig_gen_loc_w_percents
```

G-test (log likelihood ratio) of goodness of fit for significantly DI peaks that have non-significant DE. Do the annotations fit the distribution of all DI peaks?
```{r}
# Total number of annotated peaks
#259 + 251 + 499 + 387 + 61

# Expected frequencies
# In the order: Exon, Intergenic, Intron, Promoter-TSS, TTS
expected_all_DI_peaks <- c((259/1457), (251/1457), (499/1457), (387/1457), (61/1457))
expected_all_DI_peaks
# Observed frequencies
observed_DI_sig_DE_nonsig <- c(105, 98, 222, 139, 26)
GTest(x = observed_DI_sig_DE_nonsig,
      p = expected_all_DI_peaks,
      correct = "none")
```

Subset significantly DI peaks that did not match a DE gene.
```{r}
DI_sig_DE_nomatch <- subset(x = all_DI_DE_subset, `csRNA_FDR` < 0.05 & is.na(`RNA_FDR`))
print(paste0("Number of significantly DI peaks without a DE match: ", nrow(DI_sig_DE_nomatch)))
write.csv(x = DI_sig_DE_nomatch, file = "01_DI_sig_DE_nomatch.csv")
```

Summarize the annotations of significantly DI peaks that did not match a DE gene (i.e., what region of the genome the peaks map to).
```{r}
# Count unannotated DI peaks (NAs)
print(paste0("Number of unannotated peaks: ", sum(is.na(DI_sig_DE_nomatch$csRNA_Annotation))))

# Sum occurrences by genomic peak location
DI_sig_DE_nomatch_gen_loc <- as.data.frame(table(DI_sig_DE_nomatch$csRNA_Annotation))
# Change column header
names(DI_sig_DE_nomatch_gen_loc)[1] <- "Annotation"
# Manually add in row for unannotated peaks
unanno_row <- data.frame(Annotation = "Unannotated", Freq = 31)
DI_sig_DE_nomatch_gen_loc <- rbind(DI_sig_DE_nomatch_gen_loc, unanno_row)
# Add column of percentages
DI_sig_DE_nomatch_gen_loc_w_percents <- DI_sig_DE_nomatch_gen_loc %>% 
                                    mutate(Perc = `Freq` / sum(`Freq`)) %>% 
                                    arrange(Perc) %>% 
                                    mutate(Labels = scales::percent(Perc))
# Add column with data frame name
DI_sig_DE_nomatch_gen_loc_w_percents$df <- "DI_sig_DE_nomatch"
DI_sig_DE_nomatch_gen_loc_w_percents
```

G-test (log likelihood ratio) of goodness of fit for significantly DI peaks that did not match a DE gene. Do the annotations fit the distribution of all DI peaks?
```{r}
# Total number of annotated peaks
#259 + 251 + 499 + 387 + 61

# Expected frequencies
# In the order: Exon, Intergenic, Intron, Promoter-TSS, TTS
expected_all_DI_peaks <- c((259/1457), (251/1457), (499/1457), (387/1457), (61/1457))
expected_all_DI_peaks
# Observed frequencies
observed_DI_sig_DE_nomatch <- c(43, 75, 87, 33, 11)
GTest(x = observed_DI_sig_DE_nomatch,
      p = expected_all_DI_peaks,
      correct = "none")
```

Generate a stacked bar chart of genomic locations. Colorblind friendly palette from: (https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible).
```{r}
# Combine data frames
gen_loc_w_percents <- rbind(DI_sig_DE_sig_gen_loc_w_percents, DI_sig_DE_nonsig_gen_loc_w_percents, DI_sig_DE_nomatch_gen_loc_w_percents)

# Color palette
colorBlind   <- c("#E69F00", "#56B4E9", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#F0E442")

# Plot
stacked_bar <- ggplot(data = gen_loc_w_percents, aes(fill = Annotation, x = df, y = Perc)) +
  # Stacked bar chart, use position = "stack" for frequency and position = "fill" for percentages
  geom_bar(position = "fill", stat = "identity", width = 0.3) +
  #scale_fill_viridis(discrete = TRUE, option = "cividis") +
  scale_fill_manual(values = colorBlind) +
  xlab("Intersection") +
  ylab("Frequency") +
  theme(aspect.ratio = .1) +
  theme_classic() +
  # Flip x and y axes
  coord_flip()
stacked_bar

ggsave(file = "02_genomic_locations_stacked_bar_DE_or_not.pdf", stacked_bar, width = 9, height = 5)
```

#### Subset significantly DI peaks with significantly DE genes by the directionality of their log2-fold change (logFC)

Both DI and DE upregulated (+logFC) in sulfidic populations.
```{r}
DI_sig_DE_sig_upreg <- subset(x = DI_sig_DE_sig, `csRNA_LogFC` > 0 & `RNA_LogFC` > 0)
print(paste0("Upregulated DI and DE: ", nrow(DI_sig_DE_sig_upreg)))
write.csv(x = DI_sig_DE_sig_upreg, file = "03_DI_sig_DE_sig_upreg.csv")
```

Summarize the annotations of both DI and DE upregulated (i.e., what region of the genome the peaks map to).
```{r}
# Sum occurrences by genomic peak location
DI_sig_DE_sig_upreg_gen_loc <- as.data.frame(table(DI_sig_DE_sig_upreg$csRNA_Annotation))

# Change column header
names(DI_sig_DE_sig_upreg_gen_loc)[1] <- "Annotation"
DI_sig_DE_sig_upreg_gen_loc

# Add additional columns
DI_sig_DE_sig_upreg_gen_loc$DE <- "1_DE_genes"
DI_sig_DE_sig_upreg_gen_loc$Direction <- "1_Upregulated"
DI_sig_DE_sig_upreg_gen_loc

# Add column of percentages
DI_sig_DE_sig_upreg_gen_loc_w_percents <- DI_sig_DE_sig_upreg_gen_loc %>% 
                                    mutate(Perc = `Freq` / sum(`Freq`)) %>% 
                                    arrange(Perc) %>% 
                                    mutate(Labels = scales::percent(Perc))
DI_sig_DE_sig_upreg_gen_loc_w_percents
```

G-test (log likelihood ratio) of goodness of fit for both DI and DE upregulated. Do the annotations fit the distribution of all significantly DI peaks that have significant DE?
```{r}
# Total number of annotated peaks
#111 + 78 + 190 + 215 + 24

# Expected frequencies
# In the order: Exon, Intergenic, Intron, Promoter-TSS, TTS
expected_DI_sig_DE_sig_peaks <- c((111/618), (78/618), (190/618), (215/618), (24/618))
expected_DI_sig_DE_sig_peaks
# Observed frequencies
observed_DI_sig_DE_sig_upreg <- c(64, 42, 101, 140, 13)
GTest(x = observed_DI_sig_DE_sig_upreg,
      p = expected_DI_sig_DE_sig_peaks,
      correct = "none")
```

Both DI and DE downregulated (-logFC) in sulfidic populations.
```{r}
DI_sig_DE_sig_downreg <- subset(x = DI_sig_DE_sig, `csRNA_LogFC` < 0 & `RNA_LogFC` < 0)
print(paste0("Downregulated DE and DI: ", nrow(DI_sig_DE_sig_downreg)))
write.csv(x = DI_sig_DE_sig_downreg, file = "03_DI_sig_DE_sig_downreg.csv")
```

Summarize the annotations of both DI and DE downregulated (i.e., what region of the genome the peaks map to).
```{r}
# Sum occurrences by genomic peak location
DI_sig_DE_sig_downreg_gen_loc <- as.data.frame(table(DI_sig_DE_sig_downreg$csRNA_Annotation))
# Change column header
names(DI_sig_DE_sig_downreg_gen_loc)[1] <- "Annotation"
DI_sig_DE_sig_downreg_gen_loc

# Add additional columns
DI_sig_DE_sig_downreg_gen_loc$DE <- "1_DE_genes"
DI_sig_DE_sig_downreg_gen_loc$Direction <- "2_Downregulated"
DI_sig_DE_sig_downreg_gen_loc

# Add column of percentages
DI_sig_DE_sig_downreg_gen_loc_w_percents <- DI_sig_DE_sig_downreg_gen_loc %>% 
                                    mutate(Perc = `Freq` / sum(`Freq`)) %>% 
                                    arrange(Perc) %>% 
                                    mutate(Labels = scales::percent(Perc))
DI_sig_DE_sig_downreg_gen_loc_w_percents
```

G-test (log likelihood ratio) of goodness of fit for both DI and DE downregulated. Do the annotations fit the distribution of all significantly DI peaks that have significant DE?
```{r}
# Total number of annotated peaks
#111 + 78 + 190 + 215 + 24

# Expected frequencies
# In the order: Exon, Intergenic, Intron, Promoter-TSS, TTS
expected_DI_sig_DE_sig_peaks <- c((111/618), (78/618), (190/618), (215/618), (24/618))
expected_DI_sig_DE_sig_peaks
# Observed frequencies
observed_DI_sig_DE_sig_downreg <- c(34, 24, 62, 68, 6)
GTest(x = observed_DI_sig_DE_sig_downreg,
      p = expected_DI_sig_DE_sig_peaks,
      correct = "none")
```

DI upregulated (+logFC) and DE downregulated (-logFC) in sulfidic populations.
```{r}
DI_sig_upreg_DE_sig_downreg <- subset(x = DI_sig_DE_sig, `csRNA_LogFC` > 0 & `RNA_LogFC` < 0)
print(paste0("Upregulated DI, downregulated DE: ", nrow(DI_sig_upreg_DE_sig_downreg)))
write.csv(x = DI_sig_upreg_DE_sig_downreg, file = "03_DI_sig_upreg_DE_sig_downreg.csv")
```

DI downregulated (-logFC) and DE upregulated (+logFC) in sulfidic populations.
```{r}
DI_sig_downreg_DE_sig_upreg <- subset(x = DI_sig_DE_sig, `csRNA_LogFC` < 0 & `RNA_LogFC` > 0)
print(paste0("Downregulated DI, upregulated DE: ", nrow(DI_sig_downreg_DE_sig_upreg)))
write.csv(x = DI_sig_downreg_DE_sig_upreg, file = "03_DI_sig_downreg_DE_sig_upreg.csv")
```

Combine divergent DI and DE in sulfidic populations.
```{r}
DI_sig_DE_sig_divergent <- rbind(DI_sig_upreg_DE_sig_downreg, DI_sig_downreg_DE_sig_upreg)
```

Summarize the annotations of DI and DE divergent (i.e., what region of the genome the peaks map to).
```{r}
# Sum occurrences by genomic peak location
DI_sig_DE_sig_divergent_gen_loc <- as.data.frame(table(DI_sig_DE_sig_divergent$csRNA_Annotation))
# Change column header
names(DI_sig_DE_sig_divergent_gen_loc)[1] <- "Annotation"
DI_sig_DE_sig_divergent_gen_loc

# Add additional columns
DI_sig_DE_sig_divergent_gen_loc$DE <- "1_DE_genes"
DI_sig_DE_sig_divergent_gen_loc$Direction <- "3_Divergent"
DI_sig_DE_sig_divergent_gen_loc

# Add column of percentages
DI_sig_DE_sig_divergent_gen_loc_w_percents <- DI_sig_DE_sig_divergent_gen_loc %>% 
                                    mutate(Perc = `Freq` / sum(`Freq`)) %>% 
                                    arrange(Perc) %>% 
                                    mutate(Labels = scales::percent(Perc))
DI_sig_DE_sig_divergent_gen_loc_w_percents
```

G-test (log likelihood ratio) of goodness of fit for  DI and DE divergent. Do the annotations fit the distribution of all significantly DI peaks that have significant DE?
```{r}
# Total number of annotated peaks
#111 + 78 + 190 + 215 + 24

# Expected frequencies
# In the order: Exon, Intergenic, Intron, Promoter-TSS, TTS
expected_DI_sig_DE_sig_peaks <- c((111/618), (78/618), (190/618), (215/618), (24/618))
expected_DI_sig_DE_sig_peaks
# Observed frequencies
observed_DI_sig_DE_sig_divergent <- c(13, 12, 27, 7, 5)
GTest(x = observed_DI_sig_DE_sig_divergent,
      p = expected_DI_sig_DE_sig_peaks,
      correct = "none")
```

Generate a stacked bar chart of genomic locations using a colorblind friendly palette.
```{r}
# Combine data frames
gen_loc_w_percents2 <- rbind(DI_sig_DE_sig_upreg_gen_loc_w_percents, DI_sig_DE_sig_downreg_gen_loc_w_percents, DI_sig_DE_sig_divergent_gen_loc_w_percents)

# Color palette
colorBlind2   <- c("#E69F00", "#56B4E9", "#D55E00", "#009E73", "#CC79A7", "#F0E442")

# Plot
stacked_bar2 <- ggplot(data = gen_loc_w_percents2, aes(fill = Annotation, x = reorder(Direction, desc(Direction)), y = Perc)) +
  # Stacked bar chart, use position = "stack" for frequency and position = "fill" for percentages
  geom_bar(position = "fill", stat = "identity", width = 0.3) +
  #scale_fill_viridis(discrete = TRUE, option = "cividis") +
  scale_fill_manual(values = colorBlind) +
  xlab("Intersection") +
  ylab("Frequency") +
  theme(aspect.ratio = .1) +
  theme_classic() +
  coord_flip()
stacked_bar2

ggsave(file = "02_genomic_locations_stacked_bar_upreg_downreg_div.pdf", stacked_bar2, width = 9, height = 5)
```

Nodes and edges for significant DI but not significant DE.
```{r}
# Add additional columns
DE <- rep("2_Not_DE", 5)
DI_sig_DE_nonsig_gen_loc2 <- cbind(DI_sig_DE_nonsig_gen_loc, DE)
DI_sig_DE_nonsig_gen_loc2$Direction <- "4_None"
DI_sig_DE_nonsig_gen_loc2
```

Nodes and edges for significant DI but no match to DE.
```{r}
# Add additional columns
DE <- rep("3_No_match", 6)
DI_sig_DE_nomatch_gen_loc2 <- cbind(DI_sig_DE_nomatch_gen_loc, DE)
DI_sig_DE_nomatch_gen_loc2$Direction <- "4_None"
DI_sig_DE_nomatch_gen_loc2
```

#### Combine nodes and edges to make an alluvial plot (below)

```{r}
# Bind nodes and edges
alluvial_df <- rbind(DI_sig_DE_sig_upreg_gen_loc, DI_sig_DE_sig_downreg_gen_loc, DI_sig_DE_sig_divergent_gen_loc, DI_sig_DE_nonsig_gen_loc2, DI_sig_DE_nomatch_gen_loc2)

# Re-order columns
alluvial_df <- alluvial_df[, c(1, 3, 4, 2)]
alluvial_df

# Save as CSV
write.csv(x = alluvial_df, file = "03_alluvial_nodes_and_edges.csv", row.names = FALSE)
```

### CiiiDER Transcription Factor Binding Site (TFBS) Enrichment Analysis
Subsequent analyses in this section were run using R (v.4.0.3).<br>

CiiiDER was run using a GUI. Read in CiiiDER analysis (permutation 1). The foreground for this analysis were regions -150, +50 bp around the TSS from significantly *upregulated* DI peaks with significantly *upregulated* DE genes.
```{r}
ciiider_upreg <- read.csv(file = "Enrichment: background_per1_MostSigDeficit.csv")
```

Separate into significantly enriched TFBSs (Gene.P.Value < 0.05) and non-significantly enriched TFBSs (Gene.P.Value > 0.05).
```{r}
# Significant
ciiider_upreg_sig <- subset(x = ciiider_upreg, subset = `Gene.P.Value` < 0.05)

# Sort by enrichment
ciiider_upreg_sig_sorted <- ciiider_upreg_sig[order(-ciiider_upreg_sig$Log2.Enrichment),]
write.csv(x = ciiider_upreg_sig_sorted, file = "04_ciiider_upreg_significantly_enriched_tfbs_by_enrichment.csv")

# Non-significant
ciiider_upreg_nonsig <- subset(x = ciiider_upreg, subset = `Gene.P.Value` > 0.05)
```

Sort significantly enriched TFBSs by their enrichment.
```{r}
# Positively enriched
ciiider_upreg_sig_pos_enrich <- subset(x = ciiider_upreg_sig, `Log2.Enrichment` > 0)

# Zero enrichment
ciiider_upreg_sig_zero_enrich <- subset(x = ciiider_upreg_sig, `Log2.Enrichment` == 0)

# Negatively enriched
ciiider_upreg_sig_neg_enrich <- subset(x = ciiider_upreg_sig, `Log2.Enrichment` < 0)
```

#### Intersect CiiiDER results with WGCNA results

CiiiDER results.
```{r}
# Remove (var.#) at the end of some TFs and capitalize all
ciiider_upreg_sig_pos_enrich$Transcription.Factor.Simple <- toupper(gsub("\\(.*", "", ciiider_upreg_sig_pos_enrich$Transcription.Factor.Name))

# Split co-binding TFs by '::' into 2 rows
ciiider_upreg_sig_pos_enrich_split <- separate_rows(ciiider_upreg_sig_pos_enrich, Transcription.Factor.Simple, sep = '::')
```

Read in WGCNA.
```{r}
wgcna <- read.csv(file = '~/Documents/WSU/RESEARCH/pmex_tf_biomed/7_edgeR_WGCNA/11_geneInfoHabitatWild_NAs_replaced_with_LOCs.csv')

# Split SubjectSequenceID column to only include the gene name.
wgcna_split <- wgcna %>% 
                  separate(subjectSequenceID, into = paste0("subjectSequenceID", 1:4), sep = "[|_]") %>% 
                  dplyr::select(-subjectSequenceID1, -subjectSequenceID2, -subjectSequenceID4) %>% 
                  rename(subjectSequenceID = subjectSequenceID3)
```

Pull module membership (MM) and p-values (p.MM) for each gene.
```{r}
# Add MODULE_MEMBERSHIP column that copies the module membership from the column that matches the gene's module color
wgcna_split_MM <- data.frame(wgcna_split,
        MODULE_MEMBERSHIP = wgcna_split[cbind(seq_len(nrow(wgcna_split)),
        match(as.character(paste("MM.", wgcna_split$moduleColors, sep = "")), colnames(wgcna_split)))])

# Add P_VALUE column that copies the p-values from the column that matches the gene's module color
wgcna_split_MM_pval <- data.frame(wgcna_split_MM,
        P_VALUE = wgcna_split_MM[cbind(seq_len(nrow(wgcna_split_MM)),
        match(as.character(paste("p.MM.", wgcna_split$moduleColors, sep = "")), colnames(wgcna_split_MM)))])
wgcna_split_MM_pval <- transform(wgcna_split_MM_pval, 
                                 MODULE_MEMBERSHIP = as.numeric(MODULE_MEMBERSHIP),
                                 P_VALUE = as.numeric(P_VALUE))

wgcna_split_MM_pval_sig <- subset(x = wgcna_split_MM_pval, subset = `P_VALUE` < 0.05)
```

Merge CiiiDER results and WGCNA.
```{r}
ciiider_upreg_wgcna <- merge(x = ciiider_upreg_sig_pos_enrich_split, y = wgcna_split_MM_pval_sig, by.x = "Transcription.Factor.Simple", by.y = "subjectSequenceID", all = FALSE)
```

Plot a histogram of WGCNA module membership to see distribution.
```{r}
ggplot(data = ciiider_upreg_wgcna, aes(MODULE_MEMBERSHIP)) +
  geom_histogram()
```


Plot CiiiDER enrichment versus WGCNA module membership.
```{r}
# Add column noting if modules were positively, negatively, or not correlated with habitat
ciiider_upreg_wgcna <- mutate(ciiider_upreg_wgcna,
                            significance = case_when(moduleColors == "black" ~ "sig_pos_corr",
                                                     moduleColors == "salmon4" ~ "sig_pos_corr",
                                                     moduleColors == "lavenderblush3" ~ "sig_pos_corr",
                                                     moduleColors == "brown" ~ "sig_pos_corr",
                                                     moduleColors == "darkturquoise" ~ "sig_pos_corr",
                                                     moduleColors == "yellow" ~ "sig_pos_corr",
                                                     moduleColors == "skyblue" ~ "sig_pos_corr",
                                                     moduleColors == "ivory" ~ "sig_neg_corr",
                                                     moduleColors == "lightsteelblue1" ~ "sig_neg_corr",
                                                     moduleColors == "bisque4" ~ "sig_neg_corr",
                                                     moduleColors == "darkmagenta" ~ "sig_neg_corr",
                                                     moduleColors == "blue" ~ "sig_neg_corr",
                                                     moduleColors == "darkorange2" ~ "sig_neg_corr",
                                                     TRUE ~ "non-sig"))
# Check module colors for plotting
#unique(ciiider_upreg_wgcna$moduleColors)

# Change the module colors "white" to "gray95" which is more visible
ciiider_upreg_wgcna_col_changed <- mutate(ciiider_upreg_wgcna, 
                          moduleColors = replace(moduleColors, moduleColors == "white", "gray95"))

# Plot
p <- ggplot(data = ciiider_upreg_wgcna_col_changed,
            aes(x = Log2.Enrichment, y = MODULE_MEMBERSHIP, 
                shape = significance, label = Transcription.Factor.Name)) +
  # Scatterplot, color points by module color
  geom_point(color = ciiider_upreg_wgcna_col_changed$moduleColors, size = 3) +
  # Regression line (linear model) without confidence interval
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, aes(color = significance)) +
  # Colors of regression lines, legend format
  scale_color_manual(values = c("gray75", "gray40", "black"), 
                     name = "Network Correlation \nto Habitat", 
                     labels = c("p > 0.05", "Negative, p < 0.05", "Positive, p < 0.05")) +
  # Shapes of points, legend format
  scale_shape_manual(values = c(1, 17, 15), 
                     name = "Network Correlation \nto Habitat", 
                     labels = c("p > 0.05", "Negative, p < 0.05", "Positive, p < 0.05")) +
  # Add TF labels
  geom_text_repel(data = subset(x = ciiider_upreg_wgcna, subset = significance == "sig_pos_corr"),
                  size = 3, hjust = 0, vjust = 0, segment.color = NA) +
  # Add R^2 values
  #stat_regline_equation(aes(label = ..rr.label..)) +
  # x-axis title
  xlab(bquote(~ Log[2] ~ "Enrichment")) +
  # y-axis title
  ylab("Module Membership") +
  # Theme
  theme_classic()
p
ggsave(filename = "~/Documents/WSU/RESEARCH/pmex_gills_csrna_2021/19_intersect_DI_and_DE_genes/05_module_membership_vs_log2_enrichment.pdf", plot = p, width = 10, height = 7)
```

Linear regressions.
```{r}
lm_sig_pos_corr_wgcna <- lm(data = subset(ciiider_upreg_wgcna, significance == "sig_pos_corr"), MODULE_MEMBERSHIP ~ Log2.Enrichment)
summary(lm_sig_pos_corr_wgcna)
```
```{r}
lm_sig_neg_corr_wgcna <- lm(data = subset(ciiider_upreg_wgcna, significance == "sig_neg_corr"), MODULE_MEMBERSHIP ~ Log2.Enrichment)
summary(lm_sig_neg_corr_wgcna)
```
```{r}
lm_nonsig_corr_wgcna <- lm(data = subset(ciiider_upreg_wgcna, significance == "non-sig"), MODULE_MEMBERSHIP ~ Log2.Enrichment)
summary(lm_nonsig_corr_wgcna)
```

Plot same as above but split by network correlation to habitat.
```{r}
# Re-title plots
plot_names <- c(`non-sig` = "Not correlated, p > 0.05",
                `sig_neg_corr` = "Negatively correlated, p < 0.05",
                `sig_pos_corr` = "Positively correlated, p < 0.05")

# Make a data frame of slopes and intercepts from the non-parametric linear regressions.
#ablines <- data.frame(significance = c("non-sig", "sig_neg_corr", "sig_pos_corr"),
#                    slope = c(as.numeric(summary(lm_nonpar_nonsig_corr_wgcna)$coefficients[2,1]),
#                            as.numeric(summary(lm_nonpar_sig_neg_corr_wgcna)$coefficients[2,1]),
#                            as.numeric(summary(lm_nonpar_sig_pos_corr_wgcna)$coefficients[2,1])),
#                    int = c(as.numeric(summary(lm_nonpar_nonsig_corr_wgcna)$coefficients[1,1]),
#                            as.numeric(summary(lm_nonpar_sig_neg_corr_wgcna)$coefficients[1,1]),
#                            as.numeric(summary(lm_nonpar_sig_pos_corr_wgcna)$coefficients[1,1])))

# Plot
# Removed <label = Transcription.Factor.Name> to get rid of TF labels
p_fac_wrap <- ggplot(data = ciiider_upreg_wgcna_col_changed, aes(x = Log2.Enrichment, y = MODULE_MEMBERSHIP)) +
  # Scatterplot, color points by module color
  geom_point(color = ciiider_upreg_wgcna_col_changed$moduleColors, size = 2.2) +
  # Facet wrap
  facet_wrap(~significance, labeller = as_labeller(plot_names)) +
  # Regression line (linear model) without confidence interval
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, color = "black") +
  # Colors of regression lines, legend format
  # Add non-parametric regression lines <- OLD
  #geom_abline(data = ablines, aes(slope = slope, intercept = int)) +
  # Add TF labels
  #geom_text_repel(data = subset(x = ciiider_upreg_wgcna),
  #                size = 3, hjust = 0, vjust = 0, segment.color = NA) +
  # x-axis title
  xlab(bquote(~ Log[2] ~ "Enrichment")) +
  # y-axis title
  ylab("Module Membership") + 
  # Theme
  theme_classic()
p_fac_wrap
ggsave(filename = "~/Documents/WSU/RESEARCH/pmex_gills_csrna_2021/19_intersect_DI_and_DE_genes/05_module_membership_vs_log2_enrichment_facet_wrap.pdf", plot = p_fac_wrap, width = 10, height = 7)
```
Kendall–Theil Sen Siegel nonparametric linear regressions. See [this website](https://rcompanion.org/handbook/F_12.html) for help.
```{r}
#lm_nonpar_sig_pos_corr_wgcna = mblm(data = subset(ciiider_upreg_wgcna, significance == "sig_pos_corr"), MODULE_MEMBERSHIP ~ Log2.Enrichment)
#summary(lm_nonpar_sig_pos_corr_wgcna)
```

```{r}
#lm_nonpar_sig_neg_corr_wgcna = mblm(data = subset(ciiider_upreg_wgcna, significance == "sig_neg_corr"), MODULE_MEMBERSHIP ~ Log2.Enrichment)
#summary(lm_nonpar_sig_neg_corr_wgcna)
```

```{r}
#lm_nonpar_nonsig_corr_wgcna <- mblm(data = subset(ciiider_upreg_wgcna, significance == "non-sig"), MODULE_MEMBERSHIP ~ Log2.Enrichment)
#summary(lm_nonpar_nonsig_corr_wgcna)
```

## Plot module membership per module for TF and non-TF genes

Read in TRRUST v2 human dataset in TSV format (from https://www.grnpedia.org/trrust/downloadnetwork.php), downloaded 6/3/2020 (Release note (2018.04.16))
```{r}
trrust <- read.csv("~/Documents/WSU/RESEARCH/pmex_tf_biomed/7_edgeR_WGCNA/scripts/trrust_rawdata.human.tsv", sep = "\t", header = FALSE)

# Add column names
colnames(trrust) <- c("transcription_factor", "target_gene", "relationship", "pubmed_interaction_IDs")

# Generate a list of all TFs
tfs <- as.data.frame(unique(trrust$transcription_factor))

# Rename tfs column
names(tfs)[1] <- "subjectSequenceID"
```

Label all genes from the WGCNA with TF or not_TF.
```{r}
# All TFs
wgcna_tfs <- merge(x = wgcna_split_MM_pval_sig, y = tfs, by.x = "subjectSequenceID", by.y = "subjectSequenceID")
# Add TF column 
wgcna_tfs$TF <- "TF"

# All non-TFs 
wgcna_not_tfs <- anti_join(x = wgcna_split_MM_pval_sig, y = tfs, by = "subjectSequenceID")
# Add TF column 
wgcna_not_tfs$TF <- "not_TF"

# Merge TF and not_TF dataframes
wgcna_tfs_labeled <- rbind(wgcna_tfs, wgcna_not_tfs)
```

Plot gene significance and module membership across all network modules.
```{r}
module_colors <- sort(unique(wgcna_split_MM_pval_sig$moduleColors))
module_colors


ggplot(data = wgcna_tfs_labeled, aes(y = GS.habitat_wild, group = TF)) +
  geom_boxplot() +
  ylab("Gene Significance") +
  facet_wrap(~ moduleColors, nrow = 1) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
```
Add column with correlation to habitat.
```{r}
wgcna_tfs_sig_labeled <- mutate(wgcna_tfs_labeled,
                            significance = case_when(moduleColors == "black" ~ "sig_pos_corr",
                                                     moduleColors == "salmon4" ~ "sig_pos_corr",
                                                     moduleColors == "lavenderblush3" ~ "sig_pos_corr",
                                                     moduleColors == "brown" ~ "sig_pos_corr",
                                                     moduleColors == "darkturquoise" ~ "sig_pos_corr",
                                                     moduleColors == "yellow" ~ "sig_pos_corr",
                                                     moduleColors == "skyblue" ~ "sig_pos_corr",
                                                     moduleColors == "ivory" ~ "sig_neg_corr",
                                                     moduleColors == "lightsteelblue1" ~ "sig_neg_corr",
                                                     moduleColors == "bisque4" ~ "sig_neg_corr",
                                                     moduleColors == "darkmagenta" ~ "sig_neg_corr",
                                                     moduleColors == "blue" ~ "sig_neg_corr",
                                                     moduleColors == "darkorange2" ~ "sig_neg_corr",
                                                     TRUE ~ "non-sig"))
```

Group by correlation to habitat.
```{r}
# Gene significance
ggplot(data = wgcna_tfs_sig_labeled, aes(y = GS.habitat_wild), fill = moduleColors) +
  geom_boxplot(aes(fill = moduleColors)) +
  facet_grid(~ significance)

# Module membership
ggplot(data = wgcna_tfs_sig_labeled, aes(y = MODULE_MEMBERSHIP), fill = moduleColors) +
  geom_boxplot(aes(fill = moduleColors)) +
  facet_grid(~ significance)

plot_names <- c(`non-sig` = "Not correlated, p > 0.05",
                `sig_neg_corr` = "Negatively correlated, p < 0.05",
                `sig_pos_corr` = "Positively correlated, p < 0.05")

# Gene significance
p_gs <- ggplot(data = wgcna_tfs_sig_labeled, aes(y = MODULE_MEMBERSHIP, fill = TF)) +
  geom_boxplot() +
  facet_wrap(~ significance, nrow = 1, labeller = as_labeller(plot_names)) +
  scale_fill_manual(values = c("#404788FF", "#FDE725FF"), name = "Gene Type", labels = c("Non-TF", "TF")) +
  xlab("Modules") +
  ylab("Module Membership") +
  theme_classic()
p_gs + theme(axis.text.x=element_blank(),
             axis.ticks.x=element_blank())
ggsave(filename = "~/Documents/WSU/RESEARCH/pmex_gills_csrna_2021/19_intersect_DI_and_DE_genes/06_tf_vs_non-tf_gene_significance_facet_wrap.pdf", plot = p_gs + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()))

# Module membership
p_mm <- ggplot(data = wgcna_tfs_sig_labeled, aes(y = GS.habitat_wild, fill = TF)) +
  geom_boxplot() +
  facet_wrap(~ significance, nrow = 1, labeller = as_labeller(plot_names)) +
  scale_fill_manual(values = c("#404788FF", "#FDE725FF"), name = "Gene Type", labels = c("Non-TF", "TF")) +
    xlab("Modules") +
  ylab("Gene Significance") +
  theme_classic()
p_mm + theme(axis.text.x=element_blank(),
             axis.ticks.x=element_blank())
ggsave(filename = "~/Documents/WSU/RESEARCH/pmex_gills_csrna_2021/19_intersect_DI_and_DE_genes/06_tf_vs_non-tf_module_membership_facet_wrap.pdf", plot = p_mm + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()))
```

### Alluvial Diagram
Subsequent analyses in this section were run using R (v.4.0.3).

Install packages.
```{r}
install.packages("ggforce")
```
Load packages.
```
library(ggforce)
```

Read in data frame to be plotted. Variables are numbered so they order properly in alluvial diagram.
```{r}
alluvial_df <- read.csv("03_alluvial_nodes_and_edges_edited.csv", header = TRUE)
```

Add required rows and columns for ggforce.
```{r}
# Copy rows of the dataset 3 times
alluvial_df_x3 <- rbind(alluvial_df, alluvial_df, alluvial_df)

# Add ID column
alluvial_df_x3$id <- rep(c(1:30), times = 3)

# Add x column
alluvial_df_x3$x <- c(rep("Annotation", 30), rep("DE", 30), rep("Direction", 30))

# Add y column
alluvial_df_x3$y <- c(rep(c("Exon", "Intergenic", "Intron", "Promoter-TSS", "TTS", "Unannotated"), 5),
                      rep("1_DE_genes", 18), rep("2_Not_DE", 6), rep("3_No_match", 6),
                      rep("1_Upregulated", 6), rep("2_Downregulated", 6), rep("3_Divergent", 6), rep("4_None", 12))              
```

Plot.
```{r}
p <- ggplot(alluvial_df_x3, aes(x, id = id, split = y, value = Freq)) +
  geom_parallel_sets(aes(fill = Annotation), alpha = 0.3, axis.width = 0.1) +
  geom_parallel_sets_axes(aes(fill=y), axis.width = 0.1) +
  geom_parallel_sets_labels(colour = 'black', size = 4) + 
  theme_void()
p2 <- p + theme(legend.position="none")
ggsave(filename = "~/Documents/WSU/RESEARCH/pmex_gills_csrna_2021/23_alluvial_diagram/01_alluvial_plot.pdf", plot = p2, width = 15, height = 10)
```
