# CiFi Coverage Analysis

Aligned coordinate sorted segment bam is here: `/PATH/TO/RESULTS/FOLDER_NAME/bams/<sample>.cs.bam`

`mosdepth` will be used to assess coverage.

### 1. Generate bedfile of coverage windows

Read in chromosome sizes file:

```{r}
# Load necessary libraries for data manipulation
library(tidyverse) 
library(data.table)

# Read the chromosome size data from a file named "sizes.genome"
Chrom_sizes.df = fread("sizes.genome") %>%
  # Rename the first and second columns
  rename(CHROM = V1, End_coord = V2) %>%
  # Set the "CHROM" column as row names for indexing
  column_to_rownames("CHROM")

```

Specify a window size and an overlap size. Use a window size of 5Kbp with 1Kbp overlap:

```{r}
# Define the window size and overlap for creating genomic windows
window_size = 5000  # Size of each window
overlap = 1000  # Overlap between consecutive windows

# Extract the list of chromosome names from the row names of Chrom_sizes.df
chrom_list = rownames(Chrom_sizes.df)

# Generate a list of tibbles, each containing window coordinates for a chromosome
mosdepth_windows_list = chrom_list %>%
  lapply(function(x) {  
    tibble(  
      chr = x,  # Chromosome name  
      start = seq(1, Chrom_sizes.df[[x, "End_coord"]] - window_size, overlap),  # Start positions  
      end = start + window_size  # End positions  
    )  
  })

```

Write to bed:

```{r}
# Combine all chromosome-specific window tibbles into a single data frame
mosdepth_windows_list %>%
  bind_rows() %>%  # Merge the list of tibbles into a single tibble
  # Write the resulting tibble to a BED file
  write_tsv(file = "~/DennisLab/Cifi/coverage/all_chr_coverage_windows.5000.bed", 
            col_names = FALSE, quote = "none")
```

### 2. Download and activate conda environment for coverage analysis

mosdepth Github page: [mosdepth](https://github.com/brentp/mosdepth)

Bedtools man page: [BEDTools suite](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)

Download `CiFi_coverage_analysis.yml` and activate:

```         
wget https://raw.githubusercontent.com/mydennislab/CiFi/refs/heads/main/coverage_analysis/CiFi_coverage_analysis.yml
conda env create -f 
conda activate CiFi_coverage_analysis
```

### 3. Filter bed windows and run `mosdepth`

For the analysis in this study we are only interested in studying "unique" space windows or windows that do not share overlap with segmental duplications (SDs) or centromere satellites (cenSat). This filtering is only necessary if you wish to retain specific windows.

Remove SD and cenSat overlapping windows with `bedtools`:

```         
bedtools intersect -a all_chr_coverage_windows.5000.bed -b merge_bed_SD_SD98_cenSat_cenSat.noCT.sorted.merged.bed -v > unique_all_chr_coverage_windows.5000.bed
```

Run `mosdepth` on Standard data without whole-genome amplification: `GM12878_standard_SequelII`

```         
mosdepth -t 150 -b unique_all_chr_coverage_windows.5000.bed  -n GM12878_standard_SequelII /results/GM12878_standard/bams/GM12878_standard.cs.bam
```

Run `mosdepth` on CiFi data with whole-genome amplification: `GM12878_ultralow_SequelII`

```         
mosdepth -t 150 -b unique_all_chr_coverage_windows.5000.bed  -n GM12878_ultralow_SequelII /results/GM12878_ultralow/bams/GM12878_ultralow.cs.bam
```

Mosdepth output:

-   <name>.mosdepth.global.dist.txt

-   <name>.mosdepth.region.dist.txt

-   <name>.mosdepth.summary.txt

-   <name>.regions.bed.gz

-   <name>.regions.bed.gz.csi

### 4. Read `mosdepth` output into R

Load in necessary packages:

```{r}
# Load necessary libraries for data manipulation
library(tidyverse) 
library(data.table)
```

Read in Standard CiFi data:

```{r}
# Read the BED file containing standard coverage data for GM12878 cells
GM12878_standard_coverage_unique = 
  fread("GM12878_standard_SequelII.regions.bed") %>%
  # Rename columns
  dplyr::rename(Chrom = V1, Start = V2, End = V3, Standard_Depth = V4) %>%
  # Remove entries corresponding to chromosome X,Y, and mitochondrial chromosome (chrM)
  filter(!(Chrom %in% c("chrX",chrY", "chrM")))

# Normalize the coverage depth values and compute Z-scores
GM12878_standard_coverage_unique = GM12878_standard_coverage_unique %>%
  mutate(
    # Normalize coverage depth by dividing by the mean coverage depth
    norm_Standard_Depth = Standard_Depth / mean(Standard_Depth))
  )
```

Read in Ultralow CiFi data:

```{r}
# Read the BED file containing standard coverage data for GM12878 cells
GM12878_ultralow_coverage_unique = 
  fread("GM12878_ultralow_SequelII.regions.bed") %>%
  # Rename columns
  dplyr::rename(Chrom = V1, Start = V2, End = V3, Ultralow_Depth = V4) %>%
  # Remove entries corresponding to chromosome X,Y, and mitochondrial chromosome (chrM)
  filter(!(Chrom %in% c("chrX",chrY", "chrM")))

# Normalize the coverage depth values and compute Z-scores
GM12878_ultralow_coverage_unique = GM12878_ultralow_coverage_unique %>%
  mutate(
    # Normalize coverage depth by dividing by the mean coverage depth
    norm_Ultralow_Depth = Ultralow_Depth / mean(Ultralow_Depth))
  )
```

Merge the 2 dataframes:

```{r}
merged_coverage = GM12878_standard_coverage_unique %>%
	left_join(GM12878_ultralow_sequel_coverage_unique) %>%
	mutate(Middle = (Start + End)/2,.before = End) %>%
	mutate(Coverage_difference = norm_Standard_Depth - norm_Ultralow_Depth,
				 Coverage_ratio_ult_sta = norm_Ultralow_Depth/norm_Standard_Depth,
				 Coverage_ratio_sta_ult = norm_Standard_Depth/norm_Ultralow_Depth) %>%
	mutate(UCSC_coord = paste(Chrom,":",Start,"-",End,sep = ""),.before = Standard_Depth)
```

### 5. Plot mosdepth data

Plot one or more chromosomes:

```{r}
chromosome = 1 #plot just chromosome 1
chromosome = c(1,2,3,4,5,6) # plot chromosomes 1-5

CiFi_coverage_plot = merged_coverage %>%
  filter(Chrom %in% paste("chr", chromosome, sep = "")) %>%
  mutate(Chrom_num = str_extract(Chrom, "[:digit:]+") %>% as.numeric()) %>%
  # Remove raw depth values and coverage difference
  dplyr::select(-c(Standard_Depth, Ultralow_Depth, Coverage_difference)) %>%
  # Reshape data from wide to long format for plotting
  pivot_longer(c(norm_Standard_Depth, norm_Ultralow_Depth), 
               names_to = "Sample", values_to = "norm_depth") %>%
  # Clean up sample names by removing "norm_" and "_Depth" substrings
  mutate(Sample = str_remove_all(Sample, "norm_|_Depth")) %>%
  # Create a ggplot object
  ggplot(aes(x = Middle, y = norm_depth)) +
  geom_line(aes(col = Sample), alpha = 0.3, linewidth = 0.3) +
  geom_smooth(aes(col = Sample), alpha = 0.7, linewidth = 0.35) +
  xlab("Window middle (bp)") +
  ylab("Normalized Depth") +
  # Customize plot theme for better readability
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 19),  # Centered, bold title
    axis.text.x = element_text(angle = 50, hjust = 1, size = 13),  # Tilted x-axis labels
    axis.text.y = element_text(size = 13),  # Y-axis text size
    axis.title.x = element_text(size = 16),  # X-axis title size
    axis.title.y = element_text(size = 16),  # Y-axis title size
    legend.title = element_text(size = 16),  # Legend title size
    legend.text = element_text(size = 16),  # Legend text size
    panel.background = element_rect(fill = "white", color = NA),  # White panel background
    plot.background = element_rect(fill = "white", color = NA),  # White plot background
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),  # Remove minor gridlines
    panel.border = element_blank()  # Remove panel border
  ) +
  # Set plot title
  ggtitle("GM12878 CiFi normalized coverage (No SDs or Censat)") +
  # Manually set colors for different samples
  scale_color_manual(values = c("blue", "orange")) +
  # Create facets by Chrom_num, allowing free scaling along the x-axis
  facet_grid(. ~ Chrom_num, scales = "free_x", space = "free_x")
```
