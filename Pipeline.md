New CiFi Pipeline

Processing Workflow

This workflow covers CiFi digestion, FASTQ conversion, fragment mapping, filtering and conversion to pairs, and conversion to Hi-C format.

> **Note:** For Mamba environment YAML files, use:
>
> ``` bash
> mamba env create -f environment.yml
> ```
>
> Then activate the environment using the value of the `name:` field inside the YAML file. The examples below assume the environment name matches the YAML filename without `.yml`.

## 1. Create and Activate the Pore-C Tools Environment

```bash
mamba env create -f /quobyte/mydennisgrp/users/smcginty/pore-c-tools_env.yml
mamba activate pore-c-tools_env
```

## 2. Digest CiFi Reads

``` bash
pore-c-py digest \
  --output OUTPUT_name_segments.bam \
  --threads THREADS \
  INPUT.fastq \
  RESTRICTION_ENZYME
```

The input may be either a FASTQ or BAM file.

## 3. Convert BAM to FASTQ

Load Samtools:

``` bash
module load samtools
```

Convert the digested BAM file to FASTQ:

``` bash
samtools fastq \
  -@ THREADS \
  OUTPUT_name_segments.bam \
  > OUTPUT_name_segments.fastq
```

## 4. Create and Activate the Minimap2 Environment

``` bash
mamba env create -f /quobyte/mydennisgrp/users/smcginty/minimap2_env.yml
mamba activate minimap2_env
```

## 5. Map Fragments

Map the fragments with Minimap2:

``` bash
minimap2 \
  -t THREADS \
  -ax map-pb \
  REF.fa \
  OUTPUT_name_segments.fastq \
  > OUTPUT_name_segments.sam
```

Sort the alignments by genomic coordinate:

``` bash
samtools sort \
  -@ THREADS \
  -o OUTPUT_name_segments.sorted.bam \
  OUTPUT_name_segments.sam
```

Then sort by read name:

``` bash
samtools sort \
  -n \
  -@ THREADS \
  -o OUTPUT_name_segments.name.sorted.bam \
  OUTPUT_name_segments.sorted.bam
```

## 6. Filter Alignments and Assign Segment Order

``` bash
samtools view \
  -@ THREADS \
  -F 2304 \
  OUTPUT_name_segments.name.sorted.bam \
  | awk -v OFS='\t' '{
      read = $1
      sub(/:[^:]+:[^:]+$/, "", read)

      if (read == previous_read) {
          segment++
      } else {
          segment = 1
          previous_read = read
      }

      print $1, $2, $3, $4, $5, segment
  }' \
  > OUTPUT.name.sorted.segment.order.txt
```

## 7. Create and Activate the Pairtools Environment

``` bash
mamba env create -f /quobyte/mydennisgrp/users/smcginty/pairtools_env.yml
mamba activate pairtools_env
```

## 8. Convert Segments to Pairs

``` bash
/quobyte/mydennisgrp/grp/datasets/CiFi/Pore-C-Nextflow/results/GM12878_DpnII_merged_5_runs/neighbor_order_phased_bedtools/segments_to_all_pairs.sh \
  -i OUTPUT.name.sorted.segment.order.txt \
  -o OUTPUT.name.sorted.segment.order.pairs \
  -q 1
```

Compress the pairs file:

``` bash
bgzip OUTPUT.name.sorted.segment.order.pairs
```

Sort the pairs file:

``` bash
pairtools sort \
  --nproc THREADS \
  --output OUTPUT.name.sorted.segment.order.sorted.pairs.gz \
  OUTPUT.name.sorted.segment.order.pairs.gz
```

Index the sorted pairs file:

``` bash
pairix OUTPUT.name.sorted.segment.order.sorted.pairs.gz
```

## 9. Create and Activate the Java Environment

``` bash
mamba env create -f /quobyte/mydennisgrp/users/smcginty/java_env.yml
mamba activate java_env
```

## 10. Convert Pairs to Hi-C

``` bash
java -jar juicer_tools.1.9.9_jcuda.0.8.jar pre \
  -t /scratch/ \
  OUTPUT.name.sorted.segment.order.sorted.pairs.gz \
  OUTPUT.name.sorted.segment.order.pairs.hic \
  /quobyte/mydennisgrp/grp/datasets/CiFi/Pore-C-Nextflow/results/GM12878_DpnII_merged_5_runs/neighbor_distance/segment_pairs_by_order_all_MAPQ1_max_250/sizes.genome
```

## Environment Command Reference

For any Conda/Mamba YAML environment file:

``` bash
mamba env create -f ENVIRONMENT.yml
mamba activate ENVIRONMENT_NAME
```

If the environment already exists and you want to update it from the YAML file:

``` bash
mamba env update -f ENVIRONMENT.yml --prune
```
