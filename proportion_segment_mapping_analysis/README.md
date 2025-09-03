# Identify percent of segments within a read mapping to the same chromosome across different genomic bins

### 1. Download bedfiles and define bins (bed files) to loop through:

Download bedfiles

```
wget -r --level=10 -nH -nc --cut-dirs=5 --no-parent --reject "wget_index.html" --no-check-certificate --header "Cookie: sessionid=l0xwzgn8fnp9ke2fo1tyn9kule0jwwnj;" https://bioshare.bioinformatics.ucdavis.edu/bioshare/wget/cpqqdfge5lfvovq/Sean/t2t_coords/wget_index.html
```

Define genomic bins, input bam, and bedfile locations

```
# Define the list of spaces (files) to loop through
spaces=("LINE.t2tv2.0.merge.bed" "SINE.t2tv2.0.merge.bed" "SD.t2tv2.0.merge.bed" "cenSat.t2tv2.0.merge.bed")

# Define the input BAM file location
bam_file="/PATH/TO/RESULTS/FOLDER_NAME/bams/<sample>.cs.bam"

# Define the input bed directory
bed_dir="/share/dennislab/users/smcginty/t2t_coords/"
```

### 2. Run intersection across each genomic bin

Run loop to intersect input bam with each genomic bin extracting segment ID, chromosome, and MAPQ
```
# Loop through each space and process the corresponding bed file
for space in "${spaces[@]}"
do
    # Perform the intersection of BAM file with BED file to only keep alignments intersecting each genomic bin
    bedtools intersect -a "$bam_file" -b "${bed_dir}${space}" | \
    
    # Use samtools view on the intersection result and filter for mapped reads (-F 4 excludes unmapped)
    samtools view -F 4 - | \
    
    # Extract segment ID, chromosome, and MAPQ
    cut -f 1,3,5 > "mapped_segment_IDs_chr_${space}.ns.txt"
done

samtools view $bam_file | cut -f 1,3,5 > mapped_segment_IDs_chr_all.ns.txt
```

### 3. Load and process segment mapping data

Load required packages:

```         
library(tidyverse)
library(data.table)
```


Read and process data in R:

(Repeat for each restriction enzyme)

```{r}
SINE_space_mapping = fread("mapped_segment_IDs_chr_SINE.t2tv2.0.merge.bed.ns.txt") %>%
  mutate(
    # Extract read ID 
    read_ID = str_extract(V1, "/([:digit:]+)/", group = 1),
    
    # Extract the segment start position
    segment_start = str_extract(V1, ":([:digit:]+):", group = 1) %>% as.numeric(),
    
    # Extract the segment end position 
    segment_end = str_extract(V1, ":([:digit:]+)$", group = 1) %>% as.numeric(),
    
    segment_ID = paste(read_ID,"ccs",segment_start,segment_end,sep = ":"),
    read_ID_seq_run = str_extract(V1, "^[^/]+/[^/]+") %>% gsub("/", ".", .),
    region = "SINE")

LINE_space_mapping = fread("mapped_segment_IDs_chr_LINE.t2tv2.0.merge.bed.ns.txt") %>%
  mutate(
    # Extract read ID 
    read_ID = str_extract(V1, "/([:digit:]+)/", group = 1),
    
    # Extract the segment start position
    segment_start = str_extract(V1, ":([:digit:]+):", group = 1) %>% as.numeric(),
    
    # Extract the segment end position 
    segment_end = str_extract(V1, ":([:digit:]+)$", group = 1) %>% as.numeric(),
    
    segment_ID = paste(read_ID,"ccs",segment_start,segment_end,sep = ":"),
    read_ID_seq_run = str_extract(V1, "^[^/]+/[^/]+") %>% gsub("/", ".", .),
    region = "LINE")


SD_space_mapping = fread("mapped_segment_IDs_chr_SD.t2tv2.0.merge.bed.ns.txt") %>%
  mutate(
    # Extract read ID 
    read_ID = str_extract(V1, "/([:digit:]+)/", group = 1),
    
    # Extract the segment start position
    segment_start = str_extract(V1, ":([:digit:]+):", group = 1) %>% as.numeric(),
    
    # Extract the segment end position 
    segment_end = str_extract(V1, ":([:digit:]+)$", group = 1) %>% as.numeric(),
    
    segment_ID = paste(read_ID,"ccs",segment_start,segment_end,sep = ":"),
    read_ID_seq_run = str_extract(V1, "^[^/]+/[^/]+") %>% gsub("/", ".", .),
    region = "SD")

cenSat_space_mapping = fread("mapped_segment_IDs_chr_cenSat.t2tv2.0.merge.bed.ns.txt") %>%
  mutate(
    # Extract read ID 
    read_ID = str_extract(V1, "/([:digit:]+)/", group = 1),
    
    # Extract the segment start position
    segment_start = str_extract(V1, ":([:digit:]+):", group = 1) %>% as.numeric(),
    
    # Extract the segment end position 
    segment_end = str_extract(V1, ":([:digit:]+)$", group = 1) %>% as.numeric(),
    
    segment_ID = paste(read_ID,"ccs",segment_start,segment_end,sep = ":"),
    read_ID_seq_run = str_extract(V1, "^[^/]+/[^/]+") %>% gsub("/", ".", .),
    region = "cenSat")

Unique_space_mapping = fread("All_mapped_segment_IDs.chr.cs.txt") %>%
  mutate(
    # Extract read ID 
    read_ID = str_extract(V1, "/([:digit:]+)/", group = 1),
    
    # Extract the segment start position
    segment_start = str_extract(V1, ":([:digit:]+):", group = 1) %>% as.numeric(),
    
    # Extract the segment end position 
    segment_end = str_extract(V1, ":([:digit:]+)$", group = 1) %>% as.numeric(),
    
    segment_ID = paste(read_ID,"ccs",segment_start,segment_end,sep = ":"),
    read_ID_seq_run = str_extract(V1, "^[^/]+/[^/]+") %>% gsub("/", ".", .),
    region = "unique")

#Merge all regions but have each segment ID appear once. Any segment not in SINE, LINE, SD, or cenSat is considered unique.

merged_all_regions_with_unique.df = rbind(SINE_space_mapping,LINE_space_mapping,SD_space_mapping,cenSat_space_mapping,Unique_space_mapping) %>%
	distinct(V1, .keep_all = TRUE)
```

### 4. Count the percentage of segments that share a chromosome with each other segment within the same read

```{r}
DpnII_Unique_segment_mapping_chr_MAPQ1.df = merged_all_regions_with_unique.df %>%

  # Keep only rows where MAPQ >= 1
  filter(V3 >= 1) %>%
  
  # Group data by sequencing run ID
  group_by(read_ID_seq_run) %>%
  
  # Summarise within each group
  summarise(
    read_ID_seq_run = first(read_ID_seq_run),    # Preserve run ID
    segment_ID = segment_ID,                     # Keep segment IDs
    chrom = V2,                                  # Chromosome info
    concat_chr = paste(chrom, collapse = "_"),   # Concatenate chromosomes
    n = length(segment_ID),                      # Count of segments
    region = region                              # Region type
  ) %>%
  
  mutate(
    same_chr_count = (str_count(concat_chr, chrom) - 1),     # How often same chromosome repeats
    same_chr_pct = same_chr_count / (n - 1) * 100            # % of segments on the same chromosome
  ) %>%
  
  # Drop the concatenated chromosome column (no longer needed)
  select(-concat_chr) %>%

  mutate(region = factor(region, levels = c("unique", "SINE", "LINE", "SD", "cenSat")),RE = "DpnII")

HindIII_Unique_segment_mapping_chr_MAPQ1.df = ...

Merged_DpnII_HindIII_segment_mapping = rbind(DpnII_Unique_segment_mapping_chr_MAPQ1.df, HindIII_Unique_segment_mapping_chr_MAPQ1.df)

```


### 5. Make violin plot

Plot DpnII and HindIII data and save

```{r}
Merged_segment_mapping_mean_lines.ggplot <- ggplot(
  Merged_DpnII_HindIII_segment_mapping,
  aes(x = region, y = same_chr_pct, fill = RE)
) +
  geom_violin(position = position_dodge(width = 0.9)) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.25,                            # short horizontal mean line
    position = position_dodge(width = 0.9),
    colour = "black",
    fatten = 2
  ) +
  scale_fill_brewer(palette = "Paired") +
  labs(
    title = "Same Chromosome Percentage by Region for DpnII and HindIII CiFi",
    x = "Region",
    y = "Same Chromosome Percentage",
    fill = "Restriction Enzyme"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(
  filename = "Merged_DpnII_HindIII_CiFi_segment_mapping_ALL_chr.mean_lines.pdf",
  plot = Merged_segment_mapping_mean_lines.ggplot,
  width = 12, height = 12, dpi = 300, units = "in"
)
```
