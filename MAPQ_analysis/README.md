# Comparing MAPQ of segments across different genomic regions

for CiFi:

### 1. Save merged bedfiles of different repetitive genomic regions:

Segmental Duplications (SDs) and Centromeric Satellites (cenSats)
```         
SD="SD.t2tv2.0.merge.bed"
SD98="SD98.t2tv2.0.merge.bed"
cenSat="cenSat.t2tv2.0.merge.bed"
cenSat_noCT="cenSat_noCT.t2tv2.0.merge.bed"
```

Short interspersed nuclear elements (SINEs)
```
Alu="Alu.t2tv2.0.merge.bed"
tRNA_RTE="tRNA_RTE.t2tv2.0.merge.bed"
SINE_5S_Deu_L2="SINE_5S_Deu_L2.t2tv2.0.merge.bed"
tRNA="tRNA.t2tv2.0.merge.bed"
tRNA_Deu="tRNA_Deu.t2tv2.0.merge.bed"
```

Long interspersed nuclear elements (LINEs)
```
L1="L1.t2tv2.0.merge.bed"
L2="L2.t2tv2.0.merge.bed"
CR1="CR1.t2tv2.0.merge.bed"
RTE_X="RTE_X.t2tv2.0.merge.bed"
RTE_BovB="RTE_BovB.t2tv2.0.merge.bed"
Penelope="Penelope.t2tv2.0.merge.bed"
Dong_R4="Dong_R4.t2tv2.0.merge.bed"
L1_Tx1="L1_Tx1.t2tv2.0.merge.bed"
I_Jockey="I_Jockey.t2tv2.0.merge.bed"
```

Unique space represents any region not including the one listed above
```
unique_space="unique_space.t2tv2.0.bed"
```

### 2. Extract the MAPQ of segments that overlap specific regions

Use `bedtools` and `samtools` to extract the mapq for each segments depending on which region it shares overlap with

MAPQ.sh:

```         
# Usage: sh MAPQ.sh input output

# Check if correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: sh MAPQ.sh <input> <output>"
    exit 1
fi

input=$1
output=$2


#SDs and cenSats
for space in $SD $SD98 $cenSat $cenSat_noCT;
    do 
    echo $space;
    name=$(basename $space | cut -d'.' -f1);
    bedtools intersect -a ${input} -b ${space} | samtools view -@ 200 | cut -f5 | awk -v name="$name" '{print $1, name, "SD_cenSat"}' > ${output}.${name}.txt;
done

#SINEs
for space in $Alu $MIR $tRNA_RTE $SINE_5S_Deu_L2 $tRNA $tRNA_Deu;
    do 
    echo $space;
    name=$(basename $space | cut -d'.' -f1);
    bedtools intersect -a ${input} -b ${space} | samtools view -@ 200 | cut -f5 | awk -v name="$name" '{print $1, name, "SINE"}' > ${output}.${name}.txt;
done

#LINEs
for space in $L1 $L2 $CR1 $RTE_X $RTE_BovB $Penelope $Dong_R4 $L1_Tx1 $I_Jockey;
    do 
    echo $space;
    name=$(basename $space | cut -d'.' -f1);
    bedtools intersect -a ${input} -b ${space} | samtools view -@ 200 | cut -f5 | awk -v name="$name" '{print $1, name, "LINE"}' > ${output}.${name}.txt;
done

#Unique Space
for space in $unique_space;
    do 
    echo $space;
    name=$(basename $space | cut -d'.' -f1);
    bedtools intersect -a ${input} -b ${space} | samtools view | cut -f5 | awk -v name="$name" '{print $1, name, "unique_space"}' > ${output}.${name}.txt
```

Example run for DpnII CiFi, HindIII CiFi, and DpnII Hi-C:

```         
sh MAPQ.sh GM12878_CiFi_ultralow_DpnII.cs.bam   GM12878_CiFi_DpnII.mapq
sh MAPQ.sh GM12878_CiFi_ultralow_HindIII.cs.bam GM12878_CiFi_HindIII.mapq
sh MAPQ.sh GM12878_Illumina_HiC_DpnII.PE.bam    GM12878_Illumina_HiC_DpnII.mapq
```

### 3. Load MAPQ data

Load required packages:

```         
library(tidyverse)
library(data.table)
```

Read in files
```{r}
# List all files in the current directory and store them as a character vector
all_files = system("ls", intern = T)

# Convert the list of files into a data frame with a single column named 'file'
files_df = data.frame(file = all_files)

# Filter files that contain "CiFi_DpnII" in their names
DpnII_CiFi_mapq_files = files_df %>%
  filter(str_detect(file, "mapq"),      
         str_detect(file, "CiFi_DpnII")) 
  pull(file)  # Extract the matching filenames as a character vector

# Filter files that contain "Illumina_HiC_DpnII" in their names
HiC_Illumina_mapq_files = files_df %>%
  filter(str_detect(file, "mapq"),            
         str_detect(file, "Illumina_HiC_DpnII")) 
  pull(file)  # Extract the matching filenames as a character vector

# Filter files that contain "CiFi_HindIII" in their names
HindIII_CiFi_mapq_files = files_df %>%
  filter(str_detect(file, "mapq"),        
         str_detect(file, "CiFi_HindIII")) 
  pull(file)  # Extract the matching filenames as a character vector
```

Read and combine data from all DpnII_CiFi_mapq files:

```{r}
DpnII_CiFi_mapq.df <- lapply(DpnII_CiFi_mapq_files, function(x) { 
  message("Processing file: ", x) # Print the filename being processed
  fread(paste(x)) # Read the file using fread from data.table package
}) %>%
  bind_rows() %>%
  mutate(seq = "DpnII CiFi") # Add sequence type column

# Read and combine data from all HiC_Illumina_mapq files
HiC_Illumina_mapq.df <- lapply(HiC_Illumina_mapq_files, function(x) { 
  message("Processing file: ", x)
  fread(paste(x))
}) %>%
  bind_rows() %>%
  mutate(seq = "Illumina Hi-C")

# Read and combine data from all HindIII_CiFi_mapq files
HindIII_CiFi_mapq.df <- lapply(HindIII_CiFi_mapq_files, function(x) { 
  message("Processing file: ", x)
  fread(paste(x))
}) %>%
  bind_rows() %>%
  mutate(seq = "HindIII CiFi")
```

Merge all MAPQ files together:

```{r}
merged_mapq = rbind(DpnII_CiFi_mapq.df,HiC_Illumina_mapq.df,HindIII_CiFi_mapq.df) %>%
	rename(MAPQ = V1,
				 Region = V2,
				 Type = V3)
```

Calculate the percenbtage of reads in a region that a variety of MAPQ thresholds:

```{r}
mapq_pct_by_region.df = merged_mapq %>%
  # Filter for specific regions of interest
  filter(Region %in% c("unique_space","SD","SD98","Alu","MIR","cenSat","cenSat_noCT","L1","L2")) %>%
  
  # Rename "unique_space" to "Unique Space"
  mutate(Region = case_when(Region == "unique_space" ~ "Unique Space",
                            .default = Region),
         # Convert 'seq' into a factor with specific levels for ordering
         seq = factor(seq, levels = c("Illumina Hi-C","DpnII CiFi","HindIII CiFi"))) %>%

  # Define the factor levels for 'Region' to ensure a consistent order
  mutate(Region = factor(Region, levels = c("Unique Space",
                                            "Alu","MIR","L1",
                                            "L2","SD","SD98",
                                            "cenSat","cenSat_noCT"))) %>%
  
  # Group data by sequencing method and region
  group_by(seq, Region) %>%
  
  # Compute summary statistics for MAPQ values
  summarise(n = length(MAPQ),
            num_above_60 = sum(MAPQ >= 60),
            pct_mapq_60 = (num_above_60 / n) * 100,
            num_above_50 = sum(MAPQ >= 50),
            pct_mapq_50 = (num_above_50 / n) * 100,
            num_at_zero = sum(MAPQ == 0),
            pct_mapq_0 = (num_at_zero / n) * 100,
            num_above_20 = sum(MAPQ >= 20),
            pct_mapq_above_20 = (num_above_20 / n) * 100,
      

```

### 4. Plot with `ggplot2`

Percentage of reads with a MAPQ ≥ 1:

```{r}
mapq_1.ggplot <- mapq_pct_by_region.df %>% 
  mutate(MAPQ_pct = round(MAPQ_pct, 1)) %>%  # Round percentage values to 1 decimal place
  filter(MAPQ_measure == "MAPQ > 1") %>%  # Filter only for MAPQ > 1 category
  ggplot(aes(y = MAPQ_pct, x = Region, fill = seq)) +  # Set up ggplot with axes and fill color
  geom_bar(stat = "identity", position = "dodge") +  # Create grouped bar plot
  
  # Add labels and legend title
  labs(fill = "Sequencing method") +  
  
  # Customize theme for better readability
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 20),
    axis.text.y = element_text(size = 20),  
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white")
  ) +  

  # Add title and axis labels
  ggtitle("Comparison of MAPQ across different genomic regions (MAPQ > 0)") +  
  ylab("Percentage of reads") +  
  xlab("MAPQ Category") +  

  # Use a color palette for fill aesthetics
  scale_fill_brewer(palette = "Paired")

```

Percentage of reads with a MAPQ ≥ 10:

```{r}
mapq_10.ggplot <- mapq_pct_by_region.df %>% 
  mutate(MAPQ_pct = round(MAPQ_pct, 1)) %>%
  filter(MAPQ_measure == "MAPQ > 10") %>% 
  ggplot(aes(y = MAPQ_pct, x = Region, fill = seq)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  #geom_text(aes(label = paste(MAPQ_pct, "%", sep = "")), 
            #position = position_dodge(width = 0.9), vjust = 2, size = 4) +  # Much smaller label size
  labs(fill = "Sequencing method") + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
    axis.text.x = element_text(angle = 60, hjust = 1, size = 20), 
    axis.text.y = element_text(size = 20), 
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white")
  ) + 
  ggtitle("Comparison of MAPQ across different genomic regions (MAPQ > 9)") +  # Use expression for ≥ symbol
  ylab("Percentage of reads") + 
  xlab("MAPQ Category") + 
  scale_fill_brewer(palette = "Paired")
```

Percentage of reads with a MAPQ ≥ 20:

```{r}
mapq_20.ggplot <- mapq_pct_by_region.df %>% 
  mutate(MAPQ_pct = round(MAPQ_pct, 1)) %>%
  filter(MAPQ_measure == "MAPQ > 20") %>% 
  ggplot(aes(y = MAPQ_pct, x = Region, fill = seq)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  #geom_text(aes(label = paste(MAPQ_pct, "%", sep = "")), 
            #position = position_dodge(width = 0.9), vjust = 2, size = 4) +  # Much smaller label size
  labs(fill = "Sequencing method") + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
    axis.text.x = element_text(angle = 60, hjust = 1, size = 20), 
    axis.text.y = element_text(size = 20), 
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white")
  ) + 
  ggtitle("Comparison of MAPQ across different genomic regions (MAPQ > 19)") +  # Use expression for ≥ symbol
  ylab("Percentage of reads") + 
  xlab("MAPQ Category") + 
  scale_fill_brewer(palette = "Paired")
```

Percentage of reads with a MAPQ ≥ 30:

```{r}
mapq_30.ggplot <- mapq_pct_by_region.df %>% 
  mutate(MAPQ_pct = round(MAPQ_pct, 1)) %>%
  filter(MAPQ_measure == "MAPQ > 30") %>% 
  ggplot(aes(y = MAPQ_pct, x = Region, fill = seq)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  #geom_text(aes(label = paste(MAPQ_pct, "%", sep = "")), 
            #position = position_dodge(width = 0.9), vjust = 2, size = 4) +  # Much smaller label size
  labs(fill = "Sequencing method") + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
    axis.text.x = element_text(angle = 60, hjust = 1, size = 20), 
    axis.text.y = element_text(size = 20), 
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white")
  ) + 
  ggtitle("Comparison of MAPQ across different genomic regions (MAPQ > 29)") +  # Use expression for ≥ symbol
  ylab("Percentage of reads") + 
  xlab("MAPQ Category") + 
  scale_fill_brewer(palette = "Paired")
```

Percentage of reads with a MAPQ = 60:

```{r}
mapq_60.ggplot <- mapq_pct_by_region.df %>% 
  mutate(MAPQ_pct = round(MAPQ_pct, 1)) %>%
  filter(MAPQ_measure == "MAPQ = 60") %>% 
  ggplot(aes(y = MAPQ_pct, x = Region, fill = seq)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  #geom_text(aes(label = paste(MAPQ_pct, "%", sep = "")), 
            #position = position_dodge(width = 0.9), vjust = 2, size = 4) +  # Much smaller label size
  labs(fill = "Sequencing method") + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
    axis.text.x = element_text(angle = 60, hjust = 1, size = 20), 
    axis.text.y = element_text(size = 20), 
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white")
  ) + 
  ggtitle("Comparison of MAPQ across different genomic regions (MAPQ > 59)") +  # Use expression for ≥ symbol
  ylab("Percentage of reads") + 
  xlab("MAPQ Category") + 
  scale_fill_brewer(palette = "Paired")
```
