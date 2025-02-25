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
