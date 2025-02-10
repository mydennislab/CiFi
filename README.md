# CiFi Computational Pipeline

### 1. Create nextflow Conda env and activate:

`Conda` will be used to install software. The pipeline relies on `singularity` and `nextflow` . Installation for both of these as well as any necessary modules are included in the `conda` environment.

```         
wget https://raw.githubusercontent.com/mydennislab/CiFi/refs/heads/main/CiFi_Nextflow_env.yml
conda env create -f CiFi_Nextflow_env.yml
conda activate CiFi_Nextflow_env
```

### 2. Create paths to tmp files

\
By default the pipeline will store tmp files in `/tmp`. If you wish to write to new temporary files for `SINGULARITY_TMPDIR` and `SINGULARITY_CACHEDIR` save them here:

```         
export SINGULARITY_TMPDIR=/PATH/TO/TMPDIR
export SINGULARITY_CACHEDIR=/PATH/TO/CACHE
```

### 3. Create path to results folder

```         
mkdir /PATH/TO/RESULTS/FOLDER_NAME
```

When complete this folder will contain these folders:

-   bams

-   execution

-   filtered_out

-   ingress_results

-   paired_end

As well as this report: `wf-pore-c-report.html`

### 4. Run `nextflow` command and check flags


The initial processing will use the wf-pore-c pipeline. See pipeline [here](https://github.com/epi2me-labs/wf-pore-c) \


nf-core framework :https://doi.org/10.1038/s41587-020-0439-x

Check flags:

```         
nextflow run epi2me-labs/wf-pore-c \--help
```

Run Command:

```         
nextflow run epi2me-labs/wf-pore-c -work-dir ./workflow/  -c  nextflow.config \
--bam <INPUT_BAM> \ 
--ref <PATH_TO_REF_GENOME> \
--cutter <RESTRICTION ENZYME NAME> \
--out_dir <./results/FOLDER_NAME> \
--threads 120 \
--minimap2_settings '-x map-pb' \
--paired_end \ #generates paired-end file that will be converted to Hi-C
--sample <SAMPLE_NAME> # can match folder
```

### 5. Filter the paired-end bam by mapping quality (MAPQ)

This software needs to be downloaded separated. See github [here](https://github.com/mydennislab/2024-sep-mapqfilter/tree/main)

Mapq filter only retains reads where BOTH reads in the pair cross the mapq threshold as well as removing unpaired reads

Download software\
\
Environment Setup

```         
mamba create -n bamfilter bioconda::samtools bioconda::htslib zlib gcc 
```

#### Compilation

```         
g++ -O2 -o filter_bam bam_filter.cpp -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib -lhts -lz -lpthread -Wl,-rpath,$CONDA_PREFIX/lib
```

#### Usage

Run the program with the following command:

```         
./filter_bam input.bam output.bam MAPQ_threshold [num_threads]
```

-   `input.bam`: Path to the input BAM file
-   `output.bam`: Path for the output BAM file
-   `MAPQ_threshold`: Minimum MAPQ value to keep a read pair
-   `num_threads` (optional): Number of threads to use (default is 1)

Run within the paired_end folder within results with a mapq filter of 1:

```         
./filter_bam </RESULTS/FOLDER_NAME/paired_end/input.paired_end.bam> </RESULTS/FOLDER_NAME/paired_end/output.mapq_filtered_paired_end.bam> <1> <num_threads>
```

### 6. Convert filtered paired-end bam to pairs

Reactivate `CiFi_Nextflow_env` conda environment:

```         
conda activate CiFi_Nextflow_env
```

Generate chromosome sizes file from index fasta reference

```         
cut -f1,2 fasta.fai > sizes.genome
```

Run `bam2pairs` to convert filtered paired-end file into `.pairs.gz` file

```         
bam2pairs -l -c <sizes.genome> <output.mapq_filtered_paired_end.bam> <sample.bam2pairs>
```

### 7. Convert pairs file to HiC file

Install [java](https://www.java.com/en/download/help/linux_x64_install.html#download) and load

Download [juicertools](https://github.com/aidenlab/juicer/wiki/Download)\
\
`Juicertools` takes in pairs input to output a `.hic` file. This `.hic` file can be loaded into `Juicebox` or `UCSC` for visualization.\
\
Run `juicertools`:

```
java -jar juicer_tools.VERSION_jcuda.0.8.jar pre -t </PATH/TO/WORKING/DIR> <sample.bam2pairsbsorted.pairs.gz> <output_sample_name.hic> <sizes.genome>
```


### Downloading CiFi Data

##### `GM12878` CiFi HiFi v2.0 with DpnII as a restriction enzyme

```
GM12878_v2_path="https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/rlc692m7tk5cibb/ciFi/GM12878_CiFi_v2/"

wget ${GM12878_v2_path}GM12878_DpnII_CiFi.HiFi_v2.0_hifi_reads.bam #Raw CiFi HiFi reads
wget ${GM12878_v2_path}GM12878_DpnII_CiFi.HiFi_v2.0_hifi_reads.bam.pbi #Raw CiFi HiFi reads index
wget ${GM12878_v2_path}GM12878_DpnII_CiFi.HiFi_v2.0_aligned_segments.cs.bam #Aligned (coordinate sorted) segments
wget ${GM12878_v2_path}GM12878_DpnII_CiFi.HiFi_v2.0_aligned_segments.cs.bam.csi #Aligned (coordinate sorted) segments index
wget ${GM12878_v2_path}GM12878_DpnII_CiFi.HiFi_v2.0_paired_end.ns.bam #Mock paired end bam
```

##### `HG002` CiFi HiFi v2.0 with DpnII as a restriction enzyme

```
HG002_v2_path="https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/rlc692m7tk5cibb/ciFi/HG002_CiFi_v2/"

wget ${HG002_v2_path}HG002_DpnII_CiFi.HiFi_v2.0_hifi_reads.bam #Raw CiFi HiFi reads
wget ${HG002_v2_path}HG002_DpnII_CiFi.HiFi_v2.0_hifi_reads.bam.pbi #Raw CiFi HiFi reads index
wget ${HG002_v2_path}HG002_DpnII_CiFi.HiFi_v2.0_aligned_segments.cs.bam #Aligned (coordinate sorted) segments
wget ${HG002_v2_path}HG002_DpnII_CiFi.HiFi_v2.0_aligned_segments.cs.bam.csi #Aligned (coordinate sorted) segments index
wget ${HG002_v2_path}HG002_DpnII_CiFi.HiFi_v2.0_paired_end.ns.bam #Mock paired end bam
```


##### `GM12878` CiFi HiFi v3.0 with DpnII as a restriction enzyme

```
GM12878_v3_path="https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/rlc692m7tk5cibb/ciFi/GM12878_CiFi_v3/"

wget ${GM12878_v3_path}GM12878_DpnII_CiFi.HiFi_v3.0_hifi_reads.bam #Raw CiFi HiFi reads
wget ${GM12878_v3_path}GM12878_DpnII_CiFi.HiFi_v3.0_hifi_reads.bam.pbi #Raw CiFi HiFi reads index
wget ${GM12878_v3_path}GM12878_DpnII_CiFi.HiFi_v3.0_aligned_segments.cs.bam #Aligned (coordinate sorted) segments
wget ${GM12878_v3_path}GM12878_DpnII_CiFi.HiFi_v3.0_aligned_segments.cs.bam.csi #Aligned (coordinate sorted) segments index
wget ${GM12878_v3_path}GM12878_DpnII_CiFi.HiFi_v3.0_paired_end.ns.bam #Mock paired end bam
```

##### `GM12878` CiFi HiFi v3.0 with HindIII as a restriction enzyme

```
GM12878_v3_path="https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/rlc692m7tk5cibb/ciFi/GM12878_CiFi_v3/"

wget ${GM12878_v3_path}GM12878_HindIII_CiFi.HiFi_v3.0_hifi_reads.bam #Raw CiFi HiFi reads
wget ${GM12878_v3_path}GM12878_HindIII_CiFi.HiFi_v3.0_hifi_reads.bam.pbi #Raw CiFi HiFi reads index
wget ${GM12878_v3_path}GM12878_HindIII_CiFi.HiFi_v3.0_aligned_segments.cs.bam #Aligned (coordinate sorted) segments
wget ${GM12878_v3_path}GM12878_HindIII_CiFi.HiFi_v3.0_aligned_segments.cs.bam.csi #Aligned (coordinate sorted) segments index
wget ${GM12878_v3_path}GM12878_HindIII_CiFi.HiFi_v3.0_paired_end.ns.bam #Mock paired end bam
```


##### `HG002` CiFi HiFi v3.0 with DpnII as a restriction enzyme

```
HG002_v3_path="https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/rlc692m7tk5cibb/ciFi/HG002_CiFi_v3/"

wget ${HG002_v3_path}HG002_DpnII_CiFi.HiFi_v3.0_hifi_reads.bam #Raw CiFi HiFi reads
wget ${HG002_v3_path}HG002_DpnII_CiFi.HiFi_v3.0_hifi_reads.bam.pbi #Raw CiFi HiFi reads index
wget ${HG002_v3_path}HG002_DpnII_CiFi.HiFi_v3.0_aligned_segments.cs.bam #Aligned (coordinate sorted) segments
wget ${HG002_v3_path}HG002_DpnII_CiFi.HiFi_v3.0_aligned_segments.cs.bam.csi #Aligned (coordinate sorted) segments index
wget ${HG002_v3_path}HG002_DpnII_CiFi.HiFi_v3.0_paired_end.ns.bam #Mock paired end bam
```

##### `HG002` CiFi HiFi v3.0 with HindIII as a restriction enzyme

```
HG002_v3_path="https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/rlc692m7tk5cibb/ciFi/HG002_CiFi_v3/"

wget ${HG002_v3_path}HG002_HindIII_CiFi.HiFi_v3.0_hifi_reads.bam #Raw CiFi HiFi reads
wget ${HG002_v3_path}HG002_HindIII_CiFi.HiFi_v3.0_hifi_reads.bam.pbi #Raw CiFi HiFi reads index
wget ${HG002_v3_path}HG002_HindIII_CiFi.HiFi_v3.0_aligned_segments.cs.bam #Aligned (coordinate sorted) segments
wget ${HG002_v3_path}HG002_HindIII_CiFi.HiFi_v3.0_aligned_segments.cs.bam.csi #Aligned (coordinate sorted) segments index
wget ${HG002_v3_path}HG002_HindIII_CiFi.HiFi_v3.0_paired_end.ns.bam #Mock paired end bam
```


`GM12878` data is included in this study and to please cite upon use:
https://www.biorxiv.org/content/10.1101/2025.01.31.635566v1.article-info

`HG002` data remains unpublished. If planning to use in publication please reach out to the Dennis lab
