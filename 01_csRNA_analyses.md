### Raw Reads
Raw single-end 75 bp reads from gill tissues can be downloaded using the NCBI accession PRJNA743555.

All analyses use HOMER (v.4.11) unless otherwise specified.

### Trimming
Trim raw reads.
```{bash}
homerTools trim \
  -3 AGATCGGAAGAGCACACGTCT \
  -mis 2 \
  -minMatchLength 4 \
  -min 20 \
  -q 20 \
  $sample.fastq.gz
```

### Mapping
Index the *Poecilia mexicana* reference genome using STAR (v.2.7.6).
```{bash}
STAR \
        --runThreadN 7 \
        --runMode genomeGenerate \
        --genomeDir /data/kelley/projects/kerry/pmex_gills_csrna_2021/03_star/genomeDir \
        --genomeFastaFiles /data/kelley/projects/kerry/pmex_reference/GCF_001443325.1_P_mexicana-1.0_genomic.fna \
        --sjdbGTFfile /data/kelley/projects/kerry/pmex_reference/GCF_001443325.1_P_mexicana-1.0_genomic.gff \
        --sjdbOverhang 83 \
        --genomeSAindexNbases 13
```
