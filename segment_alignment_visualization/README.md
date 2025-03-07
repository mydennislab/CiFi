# Visualizing a CiFi read using SVbyEye

SVbyEye described [here](https://htmlpreview.github.io/?https://github.com/daewoooo/SVbyEye/blob/master/man/doc/SVbyEye.html)

### 1. Extract the aligned segments for a specific read

Aligned coordinate sorted segment bam is here: `/PATH/TO/RESULTS/FOLDER_NAME/bams/<sample>.cs.bam`

`samtools` will be used to manipulating bams

Select a read of interest. Extract all segments with that read ID from the coordinate sorted bam of segments. Then create a sam file with only those segments. 
```
Read_ID="m84050_240627_172653_s4/99878549/"
Input_CS_segment_bam="/results/GM12878_ultralow_Revio_DpnII/bams/GM12878_ultralow_Revio_DpnII.cs.bam"
samtools view $Input_CS_segment_bam chr2 | cut -f 1 | grep $Read_ID > segment_IDs.txt 
samtools view -h -N segment_IDs.txt $Input_CS_segment_bam chr2 > CiFi_one_read_segments.sam
```

### 2. Convert sam to paf

Use `paftools` to convert to paf
```
paftools.js sam2paf CiFi_one_read_segments.sam > CiFi_one_read_segments.paf
```

### 3. Install SVbyEye and load

SVbyEye described [here](https://htmlpreview.github.io/?https://github.com/daewoooo/SVbyEye/blob/master/man/doc/SVbyEye.html)

Install `SVbyEye`

```{r}
library(devtools)

## Install from GitHub repository
devtools::install_github("daewoooo/SVbyEye", branch="master")

library(SVbyEye)
```


### 4. Read in paf file and modify

Read in Paf file
```{r}
paf.table_one_CiFi_read <- readPaf(
    paf.file = "CiFi_one_read_segments.paf",
    include.paf.tags = TRUE, restrict.paf.tags = "cg"
)
```

Clean paf file and modify query start, end, and name:

```{r}
# Process and filter PAF (Pairwise Alignment Format)
filt_paf.table_one_CiFi_read = paf.table_one_CiFi_read %>%
  
  # Arrange the table in descending order
  arrange(desc(q.name)) %>%
  
  # Extract and mutate new columns from q.name
  mutate(
    # Extract the query start position from q.name
    q.start = str_extract(q.name, "ccs:([:digit:]+)", group = 1) %>% as.numeric(),
    
    # Extract the query end position from q.name
    q.end = str_extract(q.name, "([:digit:]+)$", group = 1) %>% as.numeric(),
    
    # Extract the base query name before "/ccs"
    q.name = str_extract(q.name, "(.+)/ccs", group = 1)
  )
```

In order to show chromosome context as well as the CiFi read on the same plot we need to adjust the scale. Set scale and adjust the sizes of the segments by this scale:

```{r}
# Calculate the total read size
read_size = max(filt_paf.table_one_CiFi_read$end) - min(filt_paf.table_one_CiFi_read$start)

# Calculate the region size
region_size = max(filt_paf.table_one_CiFi_read$end) - min(filt_paf.table_one_CiFi_read$start)

# Compute the scaling factor by dividing region size by read size
scale = round(region_size / read_size, 0)

# Apply scaling to query start, end, and length values in the filtered PAF table
scale_paf.table_one_CiFi_read = filt_paf.table_one_CiFi_read %>%
  mutate(
    q.start = q.start * scale,  # Scale query start position
    q.end = q.end * scale,      # Scale query end position
    q.len = q.len * scale       # Scale query length
  )

```

### 5. Plot the segment alignments

Plot with `plotMiro`:

```{r}
plotMiro(paf.table = scale_paf.table_one_CiFi_read, color.by = "direction") +
	ggtitle("DpnII GM12878 CiFi read in SVbyEYE") +
	theme(plot.title = element_text(hjust = 0.5,face = "bold"),
				axis.text.x = element_blank()
```				
				
*Note* the final read size scale will need to be manually added back to the plot on the X axis
