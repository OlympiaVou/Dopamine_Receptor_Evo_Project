

# The ID file must contain:
#	1) the protein IDs in the first column,
#	2) the Species name in the second column in a specific format (e.g. Homo_sapiens), the same as in the gff

# The gff path must have each gff file in a folder with the same name the species appear in the ID file:


#id_file=$1
#gff_path=$2
#fasta_path=$3
#species=$4
#prot_of_interest=$5
#window=$6

# Read the variables from the configuration file
source ./microsynteny_preprocessing_input.conf
#source ./keep_longest_from_blast.conf

# Check if the gff directory exists
if [[ ! -d "${gff_path}${species}/" ]]; then               # If the dirrectory of the fasta does not exist
        echo "Directory ${gff_path}${species}/ not found." # Print Error
        exit 1
fi

# Check if the fasta directory exists
if [[ ! -d "${fasta_path}" ]]; then               # If the dirrectory of the fasta does not exist
        echo "Directory ${fasta_path} not found." # Print Error
        exit 1
fi

# Check if the ID file exists
if [[ ! -f "${id_file}" ]]; then               # If the dirrectory of the fasta does not exist
        echo "File ${id_file} not found." # Print Error
        exit 1
fi

echo -e "STARTING THE PREPROCESS\n"
echo -e "===========================================\n"


# Path to the IDs file
echo -e "ID FILE: ${id_file}\n"

# Directory of the gff files where the gff for the specific organism is in a separate directory named with the species name
echo -e "GFF PATH: ${gff_path}\n"

# Directory of the fasta file where the proteome fasta for the specific organism is in a separate directory named with the species name
echo -e "FASTA PATH: ${fasta_path}\n"

# Species name for analysis in the following format (eg. Homo_sapiens)
echo -e "SPECIES: ${species}\n"

# Gene window size for the analysis
echo -e "GENE WINDOW: ${window}\n"

echo -e "===========================================\n"

# Break down the Species name
first=$(echo "${species}" | cut -d'_' -f1)
second=$(echo "${species}" | cut -d'_' -f2)
sp="${first:0:1}${second:0:1}"


# Get the GeneID from the proteinID
gene_IDs=()
for protein in $(grep "${species}" ${id_file} | awk {'print $1'}); do
	if [ ${db} == "NCBI" ]; then
		matches=$(grep "${protein}" ${gff_path}/${species}/* | awk -F';' {'print $3'} | grep -oP "(?<=GeneID:)\d+" | uniq )
	elif [ ${db} == "ENSEMBL" ]; then
		rna=$(grep "${protein}" ${gff_path}/${species}/* | awk -F';' {'print $2'} | awk -F':' {'print $2'})
		matches=$(grep "${rna}" ${gff_path}/${species}/* | awk -F';' {'print $2'} | grep -oP "(?<=Parent=gene:)\w+" | uniq)
		#echo "${matches}"
	fi
	if [[ -n "${matches}" ]]; then
		gene_IDs+=("${matches}")
	fi
done

echo -e "SAVED GENE IDs\n"

# Save the .gff file for the GeneID given by the user
dup=0
for gene in "${gene_IDs[@]}"; do
	((dup+=1))

	echo -e "RUNNING FOR GENE ID: ${gene}\n"

	# If the file is downloaded from NCBI do the following
	if [ ${db} == "NCBI" ]; then
		chromosome=$(grep "GeneID:${gene}" ${gff_path}/${species}/* | head -1 | awk '{print $1}') # Save the ID of the chromosome where the gene of interest is found
		output=$(grep -P "\tgene\t" ${gff_path}/${species}/*| # Save only the rows with the genes
			grep "gene_biotype=protein_coding" | # Save only the rows with protein coding genes
			grep -A ${window} -B ${window} "GeneID:${gene}"| # Save the rows consisting genes before and after the gene of interest in the window selected
			grep -oP "(?<=GeneID:)\d+" | # Save the GeneIDs for these genes
			while read ID; do # For each GeneID do the following
				grep -P "^${chromosome}" ${gff_path}/${species}/* | # Search only the chromosome where the gene of interest is found
				grep -P "GeneID:${ID}\W" # Save all the lines containing the GeneID
			done
			)

		# Count lines before and after the gene of interest
		upstream_count=$(echo "${output}" | grep -P "\tgene\t" | grep -B ${window} "GeneID:${gene}" | grep -v "GeneID:${gene}" | wc -l)
		downstream_count=$(echo "${output}" | grep -P "\tgene\t" | grep -A ${window} "GeneID:${gene}" | grep -v "GeneID:${gene}" | wc -l)
		echo -e "SAVED GENE GAPS UPSTREAM & DOWNSTREAM OF GENE OF INTEREST\n"

		# Save output in a file
		echo "${output}" > ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_for_AGAT.gff
		echo -e "SAVED GENOMIC REGIONS GENE IDs\n"
	
	# If the file is downloaded from ENSEMBL do the following
	elif [ ${db} == "ENSEMBL" ]; then
		chromosome=$(grep "${gene}" ${gff_path}/${species}/* | head -1 | awk '{print $1}') # Save the ID of the chromosome where the gene of interest is found
		output1=$(grep -P "\tgene\t" ${gff_path}/${species}/*| # Save only the rows with the genes
			grep "biotype=protein_coding"| # Save only the rows with protein coding genes
			grep -A ${window} -B ${window} "${gene}" | # Save the rows consisting genes before and after the gene of interest in the window selected
			grep -oP "(?<=ID\=gene:)\w+" | # Save the gene IDs for these genes
			while read ID; do # For each gene ID do the following
				grep -P "^${chromosome}" ${gff_path}/${species}/*| # Search only the chromosome where the gene of interest is found
				grep -P "gene:${ID}\W" # Save all the lines containing the gene ID (gene, mRNA)
			done
			)
		
		rna=$(echo "${output1}" | grep -P "\tmRNA\t" | awk -F';' {'print $1'} | awk -F':' {'print $2'}) # Save all the lines containing info about the mRNAs for the genes saved before
		output2=$(for i in $(echo "${rna}");do  grep -P "Parent=transcript:${i}" ${gff_path}/${species}/*; done) # Save all the entries (exon, CDS) for the mRNAs saved in the previous step
		output="${output1}\n${output2}" # Concatenate the two outputs to have all the info in one variable

		# Count lines before and after the gene of interest
		upstream_count=$(echo "${output}" | grep -P "\tgene\t" | grep -B ${window} "gene:${gene}" | grep -v "gene:${gene}" | wc -l)
		downstream_count=$(echo "${output}" | grep -P "\tgene\t" | grep -A ${window} "gene:${gene}" | grep -v "gene:${gene}" | wc -l)
		echo -e "SAVED GENE GAPS UPSTREAM & DOWNSTREAM OF GENE OF INTEREST\n"

		# Save output in a file		
		echo "${output}" > ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_not_ready_for_AGAT.gff
		sed -E 's/(ID|Parent|geneID)=(gene|transcript|CDS):/\1=/g' ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_not_ready_for_AGAT.gff | sort -k4,4 -n > ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_for_AGAT.gff
		echo -e "SAVED GENOMIC REGIONS GENE IDs\n"

	fi
	
	# Calculate missing genes
	missing_upstream=$(( window - upstream_count ))
	missing_downstream=$(( window - downstream_count ))

	# Save output in a file
	#echo "${output}" > ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_for_AGAT.gff
	echo -e "=========================\nAGAT: START"

	# Search for the longest transcript in the .gff file
	agat_sp_keep_longest_isoform.pl -gff ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_for_AGAT.gff -o ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_long.gff
	echo -e "AGAT: DONE\n========================="

	# Save their protein IDs
	if [ ${db} == "NCBI" ]; then
		grep -oP "(?<=protein_id\=)\w+_\d+\.\d+" ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_long.gff | awk '!seen[$0]++' > ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_long.txt
	elif [ ${db} == "ENSEMBL" ]; then
		grep -P "\tCDS\t" ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_long.gff | grep -oP "(?<=ID\=)\w+\d+" | awk '!seen[$0]++' > ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_long.txt
	fi

	# Download the sequence of each protein ID in a multi-fasta
	#python3 ~/scripts/select_seqs_from_multifasta.py ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_long.txt ${fasta_path}/${species}_ncbi_proteome.fa -o ${species}_paralog${dup}_${prot_of_interest}_${window}_genes_window.fa
	#echo -e "DOWNLOAD SEQUENCES: DONE\n========================="

	# Create a file with the ProteinIDs, the GeneIDs, and the SpeciesNames
   	all_gene_IDs=$(grep -P "\tgene\t" ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_long.gff | grep -oP "(?<=GeneID:)\d+" | awk '!seen[$0]++' )

	# Save only the gene insertions of the longest in a gff file
	grep -P "\tgene\t" ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_long.gff > ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_only_long.gff

	# Concatenate the geneIDs and proteinIDs in a tsv file
	paste ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_long.txt <(printf "%s\n" "${all_gene_IDs[@]}") > ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_IDs.tsv

	# Extract the protein IDs
	protein_ids=$(awk '{print $1}' ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_IDs.tsv)

	# Extract the last two digits before the first dot in the first column of tab-separated file
	chr=$(awk -v sp="${sp}" -F'\t' '{ match($1, /([0-9]{2})\./, a); if (a[1]) print sp a[1] }' ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_only_long.gff)

	# Extract the start of the gene
	start=$(awk -F'\t' '{ print $4 }' ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_only_long.gff)
	end=$(awk -F'\t' '{ print $5 }' ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_only_long.gff)

	# Extract the direction of the gene
	direction=$(awk -F'\t' '{ print $7 }' ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_only_long.gff)

	# Save the data in a gff file in a format compatible for microsynteny analysis
	paste <(echo "${chr}") <(echo "${protein_ids}") <(echo "${start}") <(echo "${end}") <(echo "${direction}") > ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_nogaps.gff
	echo -e "GFF for Microsynteny Analysis: DONE\n========================="

	# Add missing upstream genes
	if (( missing_upstream > 0 )); then
		for ((i = 1; i <= missing_upstream; i++)); do
			echo -e "$(echo "${chr}" | head -n 1)\tEND\t0\t0\t." >> ${species}_${prot_of_interest}${dup}_${window}_window.gff
		done
	fi

	# Add actual genes
	cat ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes_nogaps.gff >> ${species}_${prot_of_interest}${dup}_${window}_window.gff

	# Add missing downstream genes
	if (( missing_downstream > 0 )); then
    	for ((i = 1; i <= missing_downstream; i++)); do
        	echo -e "$(echo "${chr}" | head -n 1)\tEND\t0\t0\t." >> ${species}_${prot_of_interest}${dup}_${window}_window.gff
    	done
	fi

	# Remove the extra files
	rm ${species}_paralog${dup}_${prot_of_interest}_${window}_window_pc_genes*
done
