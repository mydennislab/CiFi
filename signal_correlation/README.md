# Comparing chromatin contact signal between 2 .hic files

### 1. Dump intra-chromosome contact numbers for each resolution in the .hic file

Define variables:

```
JUICER_TOOLS="/PATH/TO/juicer_tools.jar"
HIC_FILE="/PATH/TO/input.hic"            # Path to your .hic file
OUTPUT_DIR="/PATH/TO/OUTPUT_juicer_dump"

#For Human:
CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X

#For Mosquito:
CHROMOSOMES=(OX030893 OX030891 OX030892)

# List of resolutions in base pairs
RESOLUTIONS=(2500000 1000000 500000 250000 100000 50000 25000 10000 5000)
```

Loop through resolutions:

```
for RES in "${RESOLUTIONS[@]}"; do
    # Create a subdirectory for each resolution
    RESOLUTION_DIR="${OUTPUT_DIR}/res_${RES}"
    mkdir -p "$RESOLUTION_DIR"
    echo "Processing resolution: $RES"
    
    # Loop through chromosomes
    for CHR in "${CHROMOSOMES[@]}"; do
        OUTPUT_FILE="${RESOLUTION_DIR}/chr${CHR}_res${RES}.txt"
        echo "Processing chromosome: $CHR at resolution: $RES"
        
        # Execute Juicer Tools dump
        java -jar "$JUICER_TOOLS" dump observed KR "$HIC_FILE" "$CHR" "$CHR" BP "$RES" "$OUTPUT_FILE"
        
        # Check exit status and handle errors
        if [ $? -ne 0 ]; then
            echo "Error processing chromosome $CHR at resolution $RES"
        else
            echo "Extracted intrachromosomal contacts for chromosome $CHR at resolution $RES into $OUTPUT_FILE"
        fi
    done
done
```

### 2. Merge intra-chromosome contact numbers for each resolution

Set variables:

```
BASE_DIR="/PATH/TO/OUTPUT_juicer_dump"
OUTPUT_DIR="${BASE_DIR}/merged_contacts" 
mkdir -p "$OUTPUT_DIR"
RESOLUTION_FOLDERS=(res_10000 res_100000 res_1000000 res_25000 res_250000 res_2500000 res_5000 res_50000 res_500000)
```

Loop through each resolution folder:

```
for RES_FOLDER in "${RESOLUTION_FOLDERS[@]}"; do
    RES_DIR="${BASE_DIR}/${RES_FOLDER}"  # Full path to the resolution folder
    MERGED_FILE="${OUTPUT_DIR}/${RES_FOLDER}_merged.txt"  # Output file for this resolution
    
    echo "Merging files in $RES_DIR into $MERGED_FILE"
    
    # Create/overwrite the merged file
    > "$MERGED_FILE"  # Empty the file if it exists
    
    # Find and sort files by chromosome numerically (handles chr1, chr2, ..., chr22, chrX)
    FILES=$(ls "$RES_DIR"/chr*_res*.txt | sort -V)

    # Loop through sorted chromosome files
    for FILE in $FILES; do
        CHR=$(basename "$FILE" | cut -d'_' -f1)  # Extract chromosome name (e.g., "chr10")
        
        # Append data with an additional column for the chromosome
        awk -v chr="$CHR" '{print $1, $2, $3, chr}' "$FILE" >> "$MERGED_FILE"
    done
    
    echo "Finished merging for $RES_FOLDER. Output: $MERGED_FILE"
done
```

### 3. Read in contact number data

Example comparing HiC and CiFi:

```{r}

library(data.table)
library(tidyverse)

options(scipen = 999)
resolutions=c(2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000)
CiFi_dump_res.list = lapply(resolutions, function(x) {
  print(paste("Reading resolution:", x))  # Message to indicate which resolution is being read
  fread(paste("res_", x, "_merged.txt", sep = "")) %>% 
    mutate(res = x) %>%
  	rename(Win1 = V1, Win2 = V2, `DpnII CiFi` = V3, chrom = V4)
})
CiFi_dump_res.df = bind_rows(CiFi_dump_res.list)

Illumina_HiC_dump_res.list = lapply(resolutions, function(x) {
  print(paste("Reading resolution:", x))  # Message to indicate which resolution is being read
  fread(paste("res_", x, "_merged.txt", sep = "")) %>% 
    mutate(res = x) %>%
  	rename(Win1 = V1, Win2 = V2, `DpnII Hi-C` = V3, chrom = V4)
})
Illumina_HiC_dump_res.df = bind_rows(Illumina_HiC_dump_res.list)

merged_juicer_dump = left_join(CiFi_dump_res.df,Illumina_HiC_dump_res.df, by = c("Win1","Win2","chrom","res")) %>%
	replace(is.na(.), 0)
```

### 4. Plot and calculate correlation coefficient

Example comparing HiC and CiFi plot in ggplot2:

```
library(ggpmisc) # Ensure ggpmisc is installed and loaded

merged_juicer_dump %>%
filter(res == 2500000) %>%
  ggplot(aes(x = `DpnII CiFi`, y = `DpnII Hi-C`)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = FALSE) + # Adds a blue linear trend line
  theme(
    panel.background = element_rect(fill = "white"), # White background
    panel.grid = element_blank(), # Removes gridlines
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13), # Centers and bolds the title
    axis.title = element_text(size = 14), # Customizes axis titles
    axis.text = element_text(size = 12) # Customizes axis text
  ) +
  stat_poly_eq(
    aes(label = after_stat(paste("R^2: ", ..rr.label..))),
    formula = y ~ x,
    parse = TRUE,
    label.x = 0.9, # Adjusts the x position of the label
    label.y = 0.15, # Adjusts the y position of the label
    size = 4) +
  labs(
    title = "DpnII CiFi vs DpnII Hi-C at 2.5 Mb resolution",
    x = "DpnII CiFi contact number",
    y = "DpnII Hi-C contact number"
  )
```


Example comparing HiC and CiFi calculate correlation coefficients:


```
with(merged_juicer_dump %>% filter(res == 2500000), cor(`DpnII CiFi`,`DpnII Hi-C`))^2
with(merged_juicer_dump %>% filter(res == 1000000), cor(`DpnII CiFi`,`DpnII Hi-C`))^2
with(merged_juicer_dump %>% filter(res == 500000), cor(`DpnII CiFi`,`DpnII Hi-C`))^2
with(merged_juicer_dump %>% filter(res == 250000), cor(`DpnII CiFi`,`DpnII Hi-C`))^2
with(merged_juicer_dump %>% filter(res == 100000), cor(`DpnII CiFi`,`DpnII Hi-C`))^2
with(merged_juicer_dump %>% filter(res == 50000), cor(`DpnII CiFi`,`DpnII Hi-C`))^2
with(merged_juicer_dump %>% filter(res == 25000), cor(`DpnII CiFi`,`DpnII Hi-C`))^2
with(merged_juicer_dump %>% filter(res == 10000), cor(`DpnII CiFi`,`DpnII Hi-C`))^2
with(merged_juicer_dump %>% filter(res == 5000), cor(`DpnII CiFi`,`DpnII Hi-C`))^2
```
