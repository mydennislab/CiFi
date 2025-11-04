# Summary of the distance and number of intra vs interchromosomal pairs

### 1. Extract distance between intrachromosomal pairs

The input is a mapq filtered paired end bam.

Mapq filtered paired-end bam is here: `/RESULTS/FOLDER_NAME/paired_end/output.mapq_filtered_paired_end.bam`

Filter paired end bam with -f 67:

| Bit Position | Decimal Value | Meaning                         |
|--------------|---------------|---------------------------------|
| 1 (2\^0)     | 1             | Read is paired                  |
| 2 (2\^1)     | 2             | Read is mapped in a proper pair |
| 7 (2\^6)     | 64            | Read is the first in a pair     |

Then we calculate take the difference between the coordinate of the read and its pair:

```         
MAPQ_filt_PE_bam="results/Merged_GM12878_CiFi/paired_end/GM12878_CiFi_mapq1_paired_end.bam"
samtools view -@ 150 -f 67 $MAPQ_filt_PE_bam | awk '{print $8 - $4}' > GM12878_CiFi_mapq1_paired_end_intra_pairs_distance.txt
```

### 2. Read and plot distance between intrachromosomal pairs

Read in file of distance between intrachromosomal pairs:

```{r}
library(data.table)
library(tidyverse)
intra_pairs_distance = fread("GM12878_CiFi_mapq1_paired_end_intra_pairs_distance.txt.txt",header = F)
```

Plot a histogram of these distances with a bin width of 5e6:

```{r}
# Create a histogram to visualize intra-pair distances
intra_distance_ggplot = intra_pairs_distance %>%
  
  # Take absolute value of distances
  mutate(V1 = abs(V1)) %>%
  
  ggplot(aes(x = V1, fill = "Intra Distance")) +
  geom_histogram(
    col = I("white"),
    alpha = 0.8,
    position = "identity",  # Overlay histogram without stacking
    binwidth = 5e6     # Set bin width to 5 million
  ) +
  
  # Customize the x-axis scale
  scale_x_continuous(
    name = "Intra-Distance",   
    breaks = seq(0, 12e7, by = 2e7)  # Define tick marks every 20 million
    # No explicit limits set to allow ggplot to determine range dynamically
  ) +
  scale_fill_brewer(palette = "Set2") +
  ggtitle("Histogram of Intra-Distance (Counts)") +
  
  # Customize plot appearance with a clean theme
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
    axis.text.x = element_text(size = 13),   # Adjust x-axis text size
    axis.text.y = element_text(size = 13),   # Adjust y-axis text size
    panel.background = element_rect(fill = "white", color = NA), # White panel background
    plot.background = element_rect(fill = "white", color = NA),  # White plot background
    panel.grid.major = element_line(color = NA), # Remove major grid lines
    panel.grid.minor = element_line(color = NA), # Remove minor grid lines
    legend.text = element_text(size = 13)  # Adjust legend text size
  ) +
  
  # Update y-axis label to indicate the count of observations
  ylab("Counts")

```

### 3. Compare number of intrachromosomal contacts to number of total contacts

Dump contacts at lowest resolution (2.5 Mbp):

```
JAR="/PATH/TO/juicer_tools.jar"
HIC="/PATH/TO/INPUT.hic"
NORM="NONE"
MAP="observed"
UNIT="BP"
BIN=2500000                 # 2.5 Mb
LABEL="2.5_mbp"             # filename label

# Chromosomes 1–22 only
CHRS=( {1..22} )

# Loop over all chromosome pairs in upper triangle (including self)
for (( i=0; i<${#CHRS[@]}; i++ )); do
  for (( j=i; j<${#CHRS[@]}; j++ )); do
    a="${CHRS[$i]}"
    b="${CHRS[$j]}"
    out="./chr${a}_chr${b}_${LABEL}_counts.txt"
    
    echo "Starting chr${a} vs chr${b}..."
    java -jar "$JAR" dump "$MAP" "$NORM" "$HIC" "chr${a}" "chr${b}" "$UNIT" "$BIN" > "$out"
    echo "Finished chr${a} vs chr${b} → $out"
  done
done
```

Summarize into a single file:

```
#!/usr/bin/env bash
set -euo pipefail
OUTPUT_NAME="" #Fill in output file name
SUMMARY=${OUTPUT_NAME}"_chr_pair_sums.tsv"

# Header
echo -e "chrA\tchrB\tsum" > "$SUMMARY"

shopt -s nullglob
for f in chr*_chr*_*_counts.txt; do
  # Extract chr indices from the filename: chrA_chrB_...
  if [[ "$f" =~ ^chr([0-9]+)_chr([0-9]+)_ ]]; then
    a="${BASH_REMATCH[1]}"
    b="${BASH_REMATCH[2]}"
  else
    echo "Skipping unrecognized file name: $f" >&2
    continue
  fi

  # Sum column 3; ignores non-numeric values automatically
  sum=$(awk '{ if ($3+0 == $3) s += $3 } END { print s+0 }' "$f")

  echo -e "chr${a}\tchr${b}\t${sum}" >> "$SUMMARY"
  echo "Summed column 3 for $f → chr${a} chr${b} = ${sum}"
done

# Optional: sort by chrA, then chrB (version sort handles chr10 after chr9)
# keep header at top
{ head -n1 "$SUMMARY"; tail -n +2 "$SUMMARY" | sort -t$'\t' -k1,1V -k2,2V; } > "${SUMMARY}.sorted"
mv "${SUMMARY}.sorted" "$SUMMARY"

echo "Summary written to: $SUMMARY"
```

Read and count intra vs inter CiFi contacts:

```
chrom_contacts = fread("/PATH/TO/inter_intra_dump/OUTPUT_NAME_chr_pair_sums.tsv")

Total_contacts = chrom_contacts %>% pull(sum) %>% sum()

Total_contacts

intra_contacts = chrom_contacts %>%
	filter(chrA == chrB) %>%
	pull(sum) %>%
	sum()

intra_CiFi_contacts

inter_contacts = chrom_contacts %>%
	filter(chrA != chrB) %>%
	pull(sum) %>%
	sum()

inter_contacts

#intra contact percentage
intra_contacts/Total_contacts
```
