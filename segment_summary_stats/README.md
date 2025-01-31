# CiFi Segment Summary Analysis

Aligned coordinate sorted segment bam is here: `/PATH/TO/RESULTS/FOLDER_NAME/bams/<sample>.cs.bam`

`samtools` will be used to manipulating bams

### 1. Download and activate conda environment for segment summary analysis

Samtools man page: [samtools suite](https://www.htslib.org/doc/samtools.html)

Download `CiFi_analysis.yml` and activate:

```         
wget https://raw.githubusercontent.com/mydennislab/CiFi/refs/heads/main/CiFi_Nextflow_env.yml
conda env create -f CiFi_analysis.yml
conda activate CiFi_analysis
```

### 2. Extract segment IDs from Aligned coordinate sorted segment bam

Example segment IDs:\
\

|                   Segment ID                    |  Read ID  | Segment Start | Segment End | Segment Length |
|:-------------------------------------:|:---------:|:---------:|:---------:|:---------:|
| m84036_241120_004303_s3/185995108/ccs:0000:6821 | 185995108 |     0000      |    6821     |      6821      |
| m84036_241120_004303_s3/105581362/ccs:2391:6812 | 105581362 |     2391      |    6812     |      4421      |
| m84036_241120_004303_s3/76220098/ccs:1685:6136  | 76220098  |     1685      |    6136     |      4451      |
| m84036_241120_004303_s3/255530899/ccs:0000:3508 | 255530899 |     0000      |    3508     |      3508      |
| m84036_241120_004303_s3/147129652/ccs:4631:7706 | 147129652 |     4631      |    7706     |      3075      |

Filter out unmapped reads and extract segment ID:

```         
samtools view -F4 /results/GM12878_ultralow_Revio_DpnII/bams/GM12878_ultralow_Revio_DpnII.cs.bam | \
cut -f 1  > GM12878_ultralow_Revio_DpnII.mapped_segment_IDs.txt
samtools view -F4 /results/GM12878_ultralow_Revio_HindIII/bams/GM12878_ultralow_Revio_HindIII.cs.bam | \
cut -f 1  > GM12878_ultralow_Revio_HindIII.mapped_segment_IDs.txt
```

### 3. Read in mapped_segment_IDs.txt

Extract segment end and segment start to get read length. Also extract read ID

```{r}
# Read the file containing mapped segment IDs for DpnII as a cutter
DpnII_segment_lengths.df = fread("GM12878_ultralow_Revio_DpnII.mapped_segment_IDs.txt", header = FALSE) %>%
  mutate(
    # Extract read ID 
    read_ID = str_extract(V1, "/([:digit:]+)/", group = 1),
    
    # Extract the segment start position
    segment_start = str_extract(V1, ":([:digit:]+):", group = 1) %>% as.numeric(),
    
    # Extract the segment end position 
    segment_end = str_extract(V1, ":([:digit:]+)$", group = 1) %>% as.numeric(),
    
    # Calculate segment length
    segment_length = segment_end - segment_start,
    sequencer = "Revio",
    RE = "DpnII"
  )

# Read the file containing mapped segment IDs for HindIII as a cutter
HindIII_segment_lengths.df = fread("GM12878_ultralow_Revio_HindIII.mapped_segment_IDs.txt", header = FALSE) %>%
  mutate(
    # Extract read ID 
    read_ID = str_extract(V1, "/([:digit:]+)/", group = 1),
    
    # Extract the segment start position
    segment_start = str_extract(V1, ":([:digit:]+):", group = 1) %>% as.numeric(),
    
    # Extract the segment end position 
    segment_end = str_extract(V1, ":([:digit:]+)$", group = 1) %>% as.numeric(),
    
    # Calculate segment length
    segment_length = segment_end - segment_start,
    sequencer = "Revio",
    RE = "HindIII"
  )
```

### 4. Count the number of segments per read


Group by read ID to count the number of segments per read

```{r}
DpnII_segments_per_read.df = DpnII_segment_lengths.df %>%
	group_by(read_ID,sequencer,RE) %>%
	tally()
	
HindIII_segments_per_read.df = HindIII_segment_lengths.df %>%
	group_by(read_ID,sequencer,RE) %>%
	tally()	
```

### 5. Generate summary statistics

Summary statistics function:

```{r}
summary_stats = function(x){
	list(
  Mean = mean(x),
  Std_Dev = sd(x),
  Median = median(x),
  Min = min(x),
  Max = max(x),
  P90 = quantile(x, 0.9),
  P95 = quantile(x, 0.95)
)
}
```

Segment length statistics:

```{r}
DpnII_segment_length_stats = DpnII_segment_lengths.df %>% pull(segment_length)

HindIII_segment_length_stats = HindIII_segment_lengths.df %>% pull(segment_length)

summary_stats(DpnII_segment_length_stats)
summary_stats(HindIII_segment_length_stats)
```

Number of segments per read statistics:

```{r}
DpnII_segments_per_read_stats = DpnII_segments_per_read.df %>% pull(n)

HindIII_segments_per_read_stats = HindIII_segments_per_read.df %>% pull(n)

summary_stats(DpnII_segments_per_read_stats)
summary_stats(HindIII_segments_per_read_stats)
```

### 6. Plotting data


Merge the HindIII and DpnII datasets so they can be plotted together:

```{r}
Merged_DpnII_HindIII_segment_lens.df = rbind(DpnII_segment_lengths.df,HindIII_segment_lengths.df)

Merged_DpnII_HindIII_num_segments_per_read.df = rbind(DpnII_segments_per_read.df,HindIII_segments_per_read.df)
```


Create a histogram of the segment length comparing DpnII and HindIII

```{r}

# Compute the median segment length for each restriction enzyme
median_segment_lengths = Merged_DpnII_HindIII_segment_lens.df %>%
  group_by(RE) %>%  # Group data by restriction enzyme
  summarize(median_segment_len = median(segment_len))

CiFi_segment_length_ggplot = Merged_DpnII_HindIII_segment_lens.df %>%
  ggplot(aes(x = segment_len, fill = RE)) +  # Set x-axis as segment length and fill by RE
  
  # Create a histogram with density normalization
  geom_histogram(
    aes(y = ..density..),
    binwidth = 50,
    position = "identity",
    alpha = 0.8           
  ) +
  
  # Add vertical lines representing median segment length for each RE
  geom_vline(
    data = median_segment_lengths,
    aes(xintercept = median_segment_len, color = RE),
    size = 1  # Set line thickness
  ) +
  
  # Customize x-axis: Set name, breaks at every 200 bp, and limit range to 0-10,000 bp
  scale_x_continuous(
    name = "Segment Length (bp)",
    breaks = seq(0, 10000, 200),
    limits = c(0, 10000)
  ) +
  
  # Use color palettes for fill and line colors
  scale_fill_brewer(palette = "Set2") +  
  scale_color_brewer(palette = "Set2") + 
  
  # Set plot title
  ggtitle("Segment length of mapped segments for GM12878") +
  
  # Customize theme for better visualization
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),  
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),   
    panel.grid.major = element_line(color = NA),  
    panel.grid.minor = element_line(color = NA),  
    legend.text = element_text(size = 13) 
  ) +
  
  # Set y-axis label
  ylab("Density")
```


Create a histogram of the number of segments per read comparing HindIII and DpnII:

```{r}

medians_segment_count <- Merged_DpnII_HindIII_num_segments_per_read.df %>%
	group_by(RE) %>%
	summarize(median_num_segments = median(n))

CiFi_num_segments_per_read_ggplot = Merged_DpnII_HindIII_num_segments_per_read.df %>%
  ggplot(aes(x = n, fill = RE)) +
  
  # Create a histogram with density normalization
  geom_histogram(
    aes(y = ..density..),  # Convert y-axis to density scale
    binwidth = 1,   
    col = I("white"),     
    alpha = 0.8,          
    position = "identity"  
  ) +
  
  # Add vertical lines representing median segment count per RE
  geom_vline(
    data = medians_segment_count, 
    aes(xintercept = median_num_segments, color = RE),
    size = 1  # Set line thickness
  ) +
  
  # Customize x-axis: Set name, breaks at every 2 segments, and limit range to 0-44
  scale_x_continuous(
    name = "Number of segments per read",
    breaks = seq(0, 44, 2),
    limits = c(0, 44)
  ) +
  
  # Use color palettes for fill and line colors
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  
  # Set plot title
  ggtitle("Number of mapped segments per read for GM12878 CiFi") +
  
  # Customize theme for better visualization
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),  
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),  
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),   
    panel.grid.major = element_line(color = NA),  
    panel.grid.minor = element_line(color = NA),  
    legend.text = element_text(size = 13) 
  ) +
  
  # Set y-axis label
  ylab("Count")
```
