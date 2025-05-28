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

### 3. Compare number of intrachromosoal contacts to number of total contacts

Total number of intrachromosomal contacts:

```         
samtools view -@ 150 -f 67 -c $MAPQ_filt_PE_bam
```

Total number of contacts:

```         
samtools view -@ 150 -f 65 -c $MAPQ_filt_PE_bam
```
