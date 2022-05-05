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
