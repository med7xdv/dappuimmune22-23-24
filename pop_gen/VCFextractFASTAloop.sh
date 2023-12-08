#!/usr/bin/env bash
#SBATCH -J exon_fst_pi_vcftools
#SBATCH --ntasks-per-node=1
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00 # hours
#SBATCH --mem 1G
#SBATCH -o /scratch/med7xdv/popgen/err/vcf.%A_%a.out # Standard output
#SBATCH -e /scratch/med7xdv/popgen/err/vcf.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load bcftools
module load samtools

wd="/project/berglandlab/madi/all_samples/vcf"
reference_genome="${wd}/dpugenome.fa"
vcf_file="/project/berglandlab/connor/new_vcf2/combined.filtsnps10bpindels_snps_filthighmiss.qual.pass.vcf.gz"
samp="/scratch/med7xdv/popgen/dbunksamples.txt"
genes="${wd}/list_of_samples/genes.txt"


while IFS=$'\t' read -r name scaffold start stop ; do
	printf "%b\n" "${name}"
	printf "%b\n" "${scaffold}"
	printf "%b\n" "${start}"
	printf "%b\n" "${stop}"


	# read sample names from text file and loop through each sample
	while IFS= read -r sampname; do
	
  		# region of gene
  		reg="$scaffold:$start-$stop"
  
		# set output directory
		out="/scratch/med7xdv/popgen/mk/dbunk_fastas/${name}"
		mkdir -p "$out"
		
		# create reference
		samtools faidx "$reference_genome" "$scaffold:$start-$stop" > ${out}/"${name}.fasta"

  		# create sample fasta with variants
  		cat "${out}/${name}.fasta" | \
    		bcftools consensus $vcf_file \
    		--sample $sampname > \
    		"${out}/${sampname}.${name}.consensus.fa"

  		# rename fasta header
   		rename="$sampname.$scaffold.$start.$stop"
  		sed -i "s/>$reg/>$rename/g" "${out}/${sampname}.${name}.consensus.fa"
	
  		
	done < $samp
	
	
	# make the outgroup file
	outgroup="2018_Pulicaria_Pond21_22"

	cat "${out}/${name}.fasta" | \
  		bcftools consensus $vcf_file \
  		--sample $outgroup > \
  		"${out}/outgroup.fa"

  	sed -i "s/>$reg/>$outgroup/g" "${out}/outgroup.fa"
  	
  	#combine all the fasta files
  	cd $out
  	cat *.fa > \
  		combined_${name}.fasta
  	
	
done < $genes




