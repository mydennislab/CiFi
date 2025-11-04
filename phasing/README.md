# Phasing CiFi segments using Whatshap

Whatshap documentation [here](https://whatshap.readthedocs.io/en/latest/)

Whatshap haplotag documentation [here](https://whatshap.readthedocs.io/en/latest/guide.html#whatshap-haplotag)

### 1. Whatshap haplotag for phasing mapped segment bam

Phasing with whatshap:

```
whatshap haplotag -o ./GM12878_DpnII_T2Tref_haplotagged.bam \
--reference /quobyte/mydennisgrp/grp/datasets/t2t/assembly/v2.0/chm13v2.0.fa \
--output-threads 60 \
--ignore-read-groups \
--skip-missing-contigs  \
--output-haplotag-list GM12878_DpnII_T2Tref_haplotagged_IDs.txt \
NA12878.dip.vcf.gz #Diploid VCF \
GM12878_merged_all_DpnII_T2Tref.cs.bam #Input bam\
```

### 2. Convert segment bams to bed

Expand each phased segments to +- 30 Mb from start and end:

```
samtools view -@ 60 GM12878_DpnII_T2Tref_haplotagged.bam | awk 'BEGIN{OFS="\t"}{
    # get HP tag
    hp="NA";
    for(i=12;i<=NF;i++){
        if($i ~ /^HP:i:/){ hp=$i; break }
    }
    if(hp=="NA") next;

    # parse read length from read name suffix ...:start:end
    # works for IDs like m84050_.../ccs:2178:2470
    n = split($1, parts, ":");
    readlen = 0;
    if(n >= 2){
        start_off = parts[n-1]+0;
        end_off   = parts[n]+0;
        readlen   = end_off - start_off;
    }
    if(readlen <= 0) next;   # skip if we can’t parse a sensible length

    # BED start is 0-based; SAM pos ($4) is 1-based
    bedStart = $4 - 1;
    bedEnd   = bedStart + readlen;   # half-open end

    # expand ±30mb
    bedStart -= 30000000;
    if (bedStart < 0) bedStart = 0;
    bedEnd   += 30000000;

    # chr start stop hp(read as HP:i:1/2) read
    print $3, bedStart, bedEnd, hp, $1;
}' > GM12878_DpnII_T2Tref_haplotagged_expanded30mb.bed
```

Split hap1 and hap2:

```
awk '$4 == "HP:i:1" {print > "GM12878_DpnII_T2Tref_haplotagged_expanded30mb_hap1.bed"}
     $4 == "HP:i:2" {print > "GM12878_DpnII_T2Tref_haplotagged_expanded30mb_hap2.bed"}' \
     GM12878_DpnII_T2Tref_haplotagged_expanded30mb.bed
```

Filter to extract unphased segments:

```
cut -f 5 GM12878_DpnII_T2Tref_haplotagged_expanded30mb.bed > phased_segment_IDs.txt
samtools view -@ 60 -N ^phased_segment_IDs.txt -b -o unphased_GM12878_haplotagged.bam GM12878_DpnII_T2Tref_haplotagged.bam
```

Create unique chr_Read.ID to simplify overlap:

Hap1:

```
awk 'BEGIN{OFS="\t"} {
    split($5, a, "/");            # Split on "/"
    readID = a[1] "_" a[2];       # Keep only first two parts
    $1 = $1"_"readID;             # Append to column 1
    print $0
}' GM12878_DpnII_T2Tref_haplotagged_expanded30mb_hap1.bed > GM12878_merged_all_DpnII_T2Tref_haplotagged_expanded30mb_hap1_read_chr.bed

sort -k1,1 -k2,2n GM12878_merged_all_DpnII_T2Tref_haplotagged_expanded30mb_hap1_read_chr.bed > GM12878_DpnII_T2Tref_haplotagged_expanded30mb_hap1_read_chr.sorted.bed
```

Hap2:

```
awk 'BEGIN{OFS="\t"} {
    split($5, a, "/");            # Split on "/"
    readID = a[1] "_" a[2];       # Keep only first two parts
    $1 = $1"_"readID;             # Append to column 1
    print $0
}' GM12878_DpnII_T2Tref_haplotagged_expanded30mb_hap2.bed > GM12878_merged_all_DpnII_T2Tref_haplotagged_expanded30mb_hap2_read_chr.bed

sort -k1,1 -k2,2n GM12878_merged_all_DpnII_T2Tref_haplotagged_expanded30mb_hap2_read_chr.bed > GM12878_DpnII_T2Tref_haplotagged_expanded30mb_hap2_read_chr.sorted.bed
```

Repeat chr_Read.ID for unphased segments:

```
bedtools bamtobed -i unphased_GM12878_haplotagged.bam \
| awk 'BEGIN{OFS="\t"}{
    split($4, a, "/");                 # e.g. m84050_..._s4 / 194315137 / ccs:...
    rid = (length(a)>=2 ? a[1] "_" a[2] : $4);  # m84050_..._s4_194315137
    $1 = $1 "_" rid;                   # chr1 -> chr1_m84050_..._194315137
    print
}' > unphased_DpnII_GM12878_haplotagged.chrRID.bed

sort -k1,1 -k2,2n unphased_DpnII_GM12878_haplotagged.chrRID.bed > unphased_DpnII_GM12878_haplotagged.chrRID.sorted.bed
```

### 3. Impute phasing based on 30Mbp overlap

Hap 1:
```
bedtools intersect -a unphased_DpnII_GM12878_haplotagged.chrRID.sorted.bed -b GM12878_DpnII_T2Tref_haplotagged_expanded30mb_hap1_read_chr.sorted.bed -wa -wb | cut -f 4,10,11 > DpnII_GM12878_haplotagged_imputed_phase_30mb_overlap_updated_chr_read.ID.hap1.txt
```

Hap 2:
```
bedtools intersect -a unphased_DpnII_GM12878_haplotagged.chrRID.sorted.bed -b GM12878_DpnII_T2Tref_haplotagged_expanded30mb_hap2_read_chr.sorted.bed -wa -wb | cut -f 4,10,11 > DpnII_GM12878_haplotagged_imputed_phase_30mb_overlap_updated_chr_read.ID.hap2.txt
```

Count total number of segments with imputed phase. Only keep segments where every phased segment within 30 Mbp has the same phase:

```{r}
Hap1_imputed_phase = fread("DpnII_GM12878_haplotagged_imputed_phase_30mb_overlap_updated_chr_read.ID.hap1.txt", header = F)
Hap2_imputed_phase = fread("DpnII_GM12878_haplotagged_imputed_phase_30mb_overlap_updated_chr_read.ID.hap2.txt", header = F)

CiFi_imputed_phasing = rbind(Hap1_imputed_phase,Hap2_imputed_phase)

CiFi_imputed_phasing_unique_ID = CiFi_imputed_phasing %>%
  group_by(V1) %>%
  summarise(n_phase_tags = n_distinct(V2))

CiFi_imputed_phasing_unique_ID %>%
 filter(n_phase_tags == 1) %>%
 nrow()
```
