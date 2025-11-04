# Call TADs, measure concordance, and compare gaps between TADs

### 1. Call TADs using TopDom

[TopDom](https://github.com/HenrikBengtsson/TopDom)

Call domains on .hic matrix (Knight-Ruiz normalization) for all chromosomes at give resolution:

```
bash hic2topdom.sh sample1 /path/to/sample1.hic sample2 /path/to/sample2.hic [resolution in bp]
```

### 2. Calculate Measure of Concordance MoC

Calculate measure of concordance (MoC) for all chromosomes, as implemented by Zufferey et al.

R script is available [here](https://github.com/CSOgroup/TAD-benchmarking-scripts/blob/master/Figure2/fig2_fig3_fig4_fig5_moc_calc.R).


Define get_MoC function:

```{r}
library(TopDom)
library(foreach)
library(doMC)

get_MoC <- function(file1, file2, chrSize, nCpu = 1, fillInter=TRUE, meanWithInter = FALSE, correctClust = FALSE, binSize = NA, noMinusOne = FALSE, fillInterGapZero = FALSE) {

  library(foreach)
  library(doMC)

  if(correctClust)
    if(is.na(binSize))
      stop("should provide binSize for correctClust!\n")

  if(fillInterGapZero)
    fillInter <- TRUE

  registerDoMC(cores=nCpu)

  if(fillInter & meanWithInter)
    stop("not meaningful")

  if (file.info(as.character(file1))$size == 0) {
    set1DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
  } else {
    set1DT <- read.delim(file1, header = FALSE, stringsAsFactors = FALSE)
    colnames(set1DT) <- c("chromo", "start", "end")
  }
  if (file.info(as.character(file2))$size == 0) {
    set2DT <- data.frame(chromo=character(0), start = numeric(0), end =numeric(0))
  } else {
    set2DT <- read.delim(file2, header=FALSE, stringsAsFactors = FALSE)
    colnames(set2DT) <- c("chromo", "start", "end")
  }
  if(fillInterGapZero){
    set1DT_nofilled <- set1DT
    set2DT_nofilled <- set2DT
  }

  if(fillInter) {
    set1DT <- fill_part(set1DT, chrSize)
    set2DT <- fill_part(set2DT, chrSize)
  }

  if(fillInterGapZero){
    if(nrow(set1DT) > 0)
      set1DT$partType <- unlist(sapply(1:nrow(set1DT), function(x)
        ifelse(any(set1DT_nofilled$start == set1DT$start[x] & set1DT_nofilled$end == set1DT$end[x]), "domain", "gap")))
    if(nrow(set2DT) > 0)
      set2DT$partType <- unlist(sapply(1:nrow(set2DT), function(x)
        ifelse(any(set2DT_nofilled$start == set2DT$start[x] & set2DT_nofilled$end == set2DT$end[x]), "domain", "gap")))
  }

  if(nrow(set1DT) > 0)
    rownames(set1DT) <- paste0("P", 1:nrow(set1DT))
  if(nrow(set2DT) > 0)
    rownames(set2DT) <- paste0("Q", 1:nrow(set2DT))

  if(nrow(set1DT) == 0 & nrow(set2DT) > 0){
    MoC_score <- 0
  } else if(nrow(set1DT) > 0 & nrow(set2DT) == 0){
    MoC_score <- 0
  } else if( (nrow(set1DT) == nrow(set2DT))  & ( all(set1DT$start == set2DT$start) & all(set1DT$end == set2DT$end) ) ) {
    MoC_score <- 1
  } else if( (nrow(set1DT) == nrow(set2DT)) & (nrow(set1DT) == 1) ) {
    MoC_score <- 1
  } else {
    all_fragmentDT <- foreach(i = 1:nrow(set1DT), .combine='rbind') %dopar% {
      ref_start <- set1DT$start[i]
      ref_end <- set1DT$end[i]
      # if exact match
      all_matches <- which((set2DT$start == ref_start & set2DT$end == ref_end)  |
                           # nested
                           (set2DT$start >= ref_start & set2DT$end <= ref_end) |
                           # overlap left
                           (set2DT$start <= ref_start & set2DT$end >= ref_start) |
                           # overlap right
                           (set2DT$start <= ref_end & set2DT$end >= ref_end))
      all_matches_2 <- which(set2DT$end >= ref_start & set2DT$start <= ref_end)
      stopifnot(all(all_matches == all_matches_2))

      fragmentDT <- foreach(i_match = all_matches, .combine='rbind') %dopar% {
        # "fragment" size
        tmp_range <- c(set2DT$start[i_match] :set2DT$end[i_match])
        tmp_range <- tmp_range[tmp_range >= ref_start & tmp_range <= ref_end]
        frag_overlap <- tmp_range[length(tmp_range)] - tmp_range[1] + 1
        c(rownames(set1DT)[i], rownames(set2DT)[i_match], frag_overlap)
      }
      fragmentDT
    }

    # there is no intersect -> 0
    if(is.null(all_fragmentDT)) {
      MoC_score <- 0
    } else {
      # through matrix otherwise drop if nrow=1
      all_fragmentDT <- as.data.frame(matrix(all_fragmentDT, ncol=3), stringsAsFactors=F)
      rownames(all_fragmentDT) <- NULL
      colnames(all_fragmentDT) <- c("set1", "set2", "intersect_size")
      all_fragmentDT$intersect_size <- as.numeric(as.character(all_fragmentDT$intersect_size))

      MoC_score <- foreach(i = 1:nrow(all_fragmentDT), .combine = 'sum') %dopar% {
        d1 <- all_fragmentDT$set1[i]
        d2 <- all_fragmentDT$set2[i]
        # get the size of the domain from P
        Pi <- set1DT[d1, "end"]  - set1DT[d1, "start"] + 1
        # get the size of the domain from Q
        Qi <- set2DT[d2, "end"]  - set2DT[d2, "start"] + 1
        Fij <- all_fragmentDT$intersect_size[i]
        moc <- ((Fij*Fij)/(Pi * Qi))
        if(fillInterGapZero) {
          if(set1DT[all_fragmentDT$set1[i], "partType" ] != set2DT[all_fragmentDT$set2[i], "partType" ])
            moc <- 0
        }
        moc
      }
      # if wanted, correct for the number of clusters
      if(correctClust) {
        # penalize number of clusters / number tot of bins
        penaltyTerm <- nrow(all_fragmentDT) / (chrSize/binSize)
        MoC_score <- MoC_score + (1 - penaltyTerm)
      }
    }

    if(! noMinusOne)
      MoC_score <- MoC_score - 1

    MoC_score <- MoC_score/(sqrt(nrow(set1DT) * nrow(set2DT)) - 1)
  }
  if(meanWithInter) {
    bd_only_set1DT <- bd_only(set1DT, chrSize = chrSize)
    file1_BD_only <- sub("_final_domains.txt", "_final_domains_BD_only.txt", file1)
    write.table(bd_only_set1DT, file = file1_BD_only, col.names = F, row.names = F, quote=F, sep="\t")

    bd_only_set2DT <- bd_only(set2DT, chrSize = chrSize)
    file2_BD_only <- sub("_final_domains.txt", "_final_domains_BD_only.txt", file2)
    write.table(bd_only_set2DT, file = file2_BD_only, col.names = F, row.names = F, quote=F, sep="\t")

    MoC_score_inter <- get_MoC(file1_BD_only, file2_BD_only,chrSize, fillInter=fillInter, meanWithInter = FALSE )
    MoC_score <- (MoC_score + MoC_score_inter)/2
  }
  MoC_score
}
```

Calculate measure of concordance for each chromosome:

```
moc_results <- list()

for (chr in names(chr_sizes)) {
  file1 <- paste0(path1,"CiFi_", chr, "_KR.bed3")
  file2 <- paste0(path1, "HiC_", chr, "_KR.bed3")

  if (!file.exists(file1) || !file.exists(file2)) {
    message(sprintf("Skipping %s — missing file(s)", chr))
    moc_results[[chr]] <- NA
    next
  }

  # Try to compute MoC
  result <- tryCatch({
    score <- get_MoC(
      file1 = file1,
      file2 = file2,
      chrSize = chr_sizes[chr],
      nCpu = 4,
      fillInter = FALSE,
      meanWithInter = FALSE,
      correctClust = TRUE,
      binSize = 50000,
      noMinusOne = FALSE
    )
    score
  }, error = function(e) {
    message(sprintf("Error in %s: %s", chr, e$message))
    NA
  })

  # Store result and print progress
  moc_results[[chr]] <- result
  message(sprintf("Finished %s — MoC: %.4f", chr, result))
}
```

### 3. Calculate jaccard index to measure concordance


bedtools jaccard info available [here](https://bedtools.readthedocs.io/en/latest/content/tools/jaccard.html)

Overall jaccard index:

```
#!/bin/bash

# Define chromosome names
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" \
"chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" \
"chr21" "chr22" "chrX" "chrY" "chrM")

# Loop over each chromosome
for chr in "${chromosomes[@]}"; do
    file1="../PATH/TO/CiFi_topdom/CiFi_KR_${chr}_KR.bed3"
    file2="../PATH/TO/HiC_topdom/HiC_KR_${chr}_KR.bed3"

    # Run bedtools jaccard and extract the Jaccard score (skip the header, get second line, 5th column)
    jaccard=$(bedtools jaccard -a "$file1" -b "$file2" | tail -n 1 | cut -f 3)
    
    echo "${chr}: ${jaccard}"
done
```

Jaccard index requiring 25% reciprocal overlap between TADs:

```
#!/bin/bash

# Define chromosome names
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" \
"chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" \
"chr21" "chr22" "chrX" "chrY" "chrM")

# Loop over each chromosome
for chr in "${chromosomes[@]}"; do
    file1="../PATH/TO/CiFi_topdom/CiFi_KR_${chr}_KR.bed3"
    file2="../PATH/TO/HiC_topdom/HiC_KR_${chr}_KR.bed3"

    # Run bedtools jaccard and extract the Jaccard score (skip the header, get second line, 5th column)
    jaccard=$(bedtools jaccard -f .25 -r -a "$file1" -b "$file2" | tail -n 1 | cut -f 3)
    
    echo "${chr}: ${jaccard}"
done
```

### 4. Compare TAD gaps between CiFi and Hi-C

cat CiFi_merged_mapq1_50000_KR_chr{1..22}_KR.bed > CiFi_merged_all_chr1_22_KR.bed

