sample=$1
peaks=$2
bins=$3
folder=$4
bedtools subtract -A -a ${folder}/${sample}.bed -b ${peaks}.bed > ${folder}/${sample}_nopeaks.bed
bedtools intersect -a ${bins} -b ${folder}/${sample}_nopeaks.bed -c > ${folder}/${sample}_nopeaks_counts.txt

if [[ ! -f "${folder}/${peaks}_bins.bed" ]]; then
	bedtools intersect -a ${bins} -b ${peaks}.bed -wao | awk '{print $1"\t"$2"\t"$3-5"\t"$7}' - | bedtools merge -i - -c 4 -o sum |  awk '{print $1"\t"$2"\t"$3+5"\t"$4}' > ${folder}/${peaks}_bins.bed
fi

