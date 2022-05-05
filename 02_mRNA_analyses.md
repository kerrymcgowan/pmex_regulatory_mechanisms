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

### Create gene counts matrix
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





