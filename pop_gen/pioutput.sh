#!/usr/bin/env bash
#SBATCH -J exon_fst_pi_vcftools
#SBATCH --ntasks-per-node=1
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00 # hours
#SBATCH --mem 1G
#SBATCH -o /scratch/med7xdv/popgen/err/exon_vcf.%A_%a.out # Standard output
#SBATCH -e /scratch/med7xdv/popgen/err/exon_vcf.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# start
echo "Start" $(date)

# Exon Functions
exonFunc () {

# Load Modules
module load vcftools
module load tabix

# Working folder is core folder where this pipeline is being run.
wd=/scratch/med7xdv/popgen

# Input file
IN_GZVCF=${wd}/data/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.ann.vcf.gz

# Master Bed file
bed=${wd}/data/miss10.daphnia.pulex.merged.RMoutHiCGM.NsandDepthandChrEnd.final.merge50.bed

# Meta Exon
exon=($( ls ${wd}/exons/* ))

# Move to working directory
cd ${wd}/

# Go through each exon
while read -r i; do

#i="Scaffold_7757_HRSCAF_8726       3472398 3472582 exon    Daphnia04572-RA"

# Start
echo ${i}
date

# Extract metadata
chrom=$( echo ${i} | cut -f1 -d " " )
start=$( echo ${i} | cut -f2 -d " " )
stop=$( echo ${i} | cut -f3 -d " " )
gene=$( echo ${i} | cut -f5 -d " " )

# Calculate length
len=$(echo $i | awk '{$6 = $3-$2 } 1' | cut -f6 -d " ")

# VCF functions
analy="--window-pi ${len} --keep ${wd}/samples/d8samples.txt"

# Output name
out_namey="all_pi_nomiss"

# Subset VCF and Run VCFTools
tabix -h ${IN_GZVCF} ${chrom}:${start}-${stop} |
vcftools \
--vcf - \
--max-missing 1 \
--chr ${chrom} \
--from-bp ${start} \
--to-bp ${stop} \
--exclude-positions ${bed} \
`echo ${analy}` \
-c | \
awk 'NR>1 {print $0 "\t" "'"$start"'" "\t" "'"$stop"'" "\t" "'"$gene"'"}' \
>> ${wd}/piout/${out_namey}_${SLURM_ARRAY_TASK_ID}

# Finish exon
done < ${exon[${SLURM_ARRAY_TASK_ID}]}
}

# Exports function
export -f exonFunc

# Go through each function
exonFunc ${SLURM_ARRAY_TASK_ID}
echo "VCF completed" $(date)

# Submit script:
# 1580 exons
# sbatch --array=0-1580 pifst.sh
