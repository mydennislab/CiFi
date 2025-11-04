#!/bin/bash

# Author: Colin Shew
# Last updated: 18 Feb 2020

# Call domains on .hic matrix (Knight-Ruiz normalization) for all chromosomes at give resolution
# Usage:
# bash hic2topdom.sh sample1 /path/to/sample1.hic sample2 /path/to/sample2.hic [resolution in bp]

#module load R/3.6.1 # whichever version of R has TopDom installed

samp1=$1
path1=$2
samp2=$3
path2=$4
r=$5

echo "running TopDom at resolution" $r

for chr in {1..22} X Y; do
        echo "----------------------------------- PROCESSING CHR"$chr "-----------------------------------"

        # dump matrix
        java -jar /share/dennislab/programs/juicer/CPU/common/juicer_tools.1.9.9_jcuda.0.8.jar dump -d observed KR $path1 $chr $chr BP $r $samp1"_chr"$chr"_KR.txt"
        java -jar /share/dennislab/programs/juicer/CPU/common/juicer_tools.1.9.9_jcuda.0.8.jar dump -d observed KR $path2 $chr $chr BP $r $samp2"_chr"$chr"_KR.txt"
        echo "done dumping matrices; chr"$chr

        # generate bins
        nbins1=$(wc -l $samp1"_chr"$chr"_KR.txt" | awk '{print $1}')
        nbins2=$(wc -l $samp2"_chr"$chr"_KR.txt" | awk '{print $1}')
        if [ $nbins1 -ne $nbins2 ]; then
                echo "*** WARNING: unequal bin number in chr"$chr "; N bins:" $nbins1 $nbins2
                #continue
                if [ $nbins1 -gt $nbins2 ]; then # more bins in samp1 --> crop to samp2 dimensions
                        echo ".....cropping" $samp1"_chr"$chr"_KR.txt"
                        head -$nbins2 $samp1"_chr"$chr"_KR.txt" > tmp # trim rows
                        awk -v b=1 -v e=$nbins2 'BEGIN{FS=OFS="\t"} {for (i=b;i<=e;i++) printf "%s%s", $i, (i<e ? OFS : ORS)}' tmp > $samp1"_chr"$chr"_KR.txt" # trim columns
                        nbins1=$nbins2
                elif [ $nbins1 -lt $nbins2 ]; then # more bins in samp2 --> crop to samp1 dimensions
                        echo ".....cropping" $samp2"_chr"$chr"_KR.txt"
                        head -$nbins1 $samp2"_chr"$chr"_KR.txt" > tmp # trim rows
                        awk -v b=1 -v e=$nbins1 'BEGIN{FS=OFS="\t"} {for (i=b;i<=e;i++) printf "%s%s", $i, (i<e ? OFS : ORS)}' tmp > $samp2"_chr"$chr"_KR.txt" # trim columns
                fi
        fi
        for i in $(seq 1 $nbins1); do echo "chr"$chr$'\t'$((($i-1)*$r))$'\t'$(($i*$r)); done > bins.tmp

        # format into (N x N+3) and remove NaN
        paste bins.tmp $samp1"_chr"$chr"_KR.txt" > tmp
        mv tmp $samp1"_chr"$chr"_KR.txt"
        sed -i 's/NaN/0/g' $samp1"_chr"$chr"_KR.txt"

        paste bins.tmp $samp2"_chr"$chr"_KR.txt" > tmp
        mv tmp $samp2"_chr"$chr"_KR.txt"
        sed -i 's/NaN/0/g' $samp2"_chr"$chr"_KR.txt"
        echo "done formatting matrices; chr"$chr

        # run TopDom (uses default w=5)
        Rscript topdom.R $samp1"_chr"$chr"_KR.txt"
        Rscript topdom.R $samp2"_chr"$chr"_KR.txt"
        echo "done running TopDom; chr"$chr

        # create bedpe and bed outputs
        awk -v OFS='\t' '{print "#chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score", "strand1", "strand2", "color"}' $samp1"_chr"$chr"_KR.bed" | head -1 > $samp1"_chr"$chr"_KR.bedpe" # generate header
        awk -v OFS='\t' '{if ($4=="domain") print $1,$2,$3,$1,$2,$3,".",".",".",".","255,255,0"}' $samp1"_chr"$chr"_KR.bed" >> $samp1"_chr"$chr"_KR.bedpe"
        tail -n +2 $samp1"_chr"$chr"_KR.bedpe" | awk -v OFS='\t' '{print $1,$2,$3}' > $samp1"_chr"$chr"_KR.bed3"

        awk -v OFS='\t' '{print "#chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score", "strand1", "strand2", "color"}' $samp2"_chr"$chr"_KR.bed" | head -1 > $samp2"_chr"$chr"_KR.bedpe" # generate header
        awk -v OFS='\t' '{if ($4=="domain") print $1,$2,$3,$1,$2,$3,".",".",".",".","255,255,0"}' $samp2"_chr"$chr"_KR.bed" >> $samp2"_chr"$chr"_KR.bedpe"
        tail -n +2 $samp2"_chr"$chr"_KR.bedpe" | awk -v OFS='\t' '{print $1,$2,$3}' > $samp2"_chr"$chr"_KR.bed3"
        echo "done creating BED outputs; chr"$chr
done

# concatenate all chr
cat $samp1"_chr"*"_KR.bedpe" > $samp1"_KR.bedpe"
cat $samp2"_chr"*"_KR.bedpe" > $samp2"_KR.bedpe"
cat $samp1"_chr"*"_KR.bed3" | sort -k1,1 -k2,2n > $samp1"_KR.bed3"
cat $samp2"_chr"*"_KR.bed3" | sort -k1,1 -k2,2n > $samp2"_KR.bed3"
