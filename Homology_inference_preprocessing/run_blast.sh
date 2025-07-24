#!/bin/bash

# Check if the number of arguments is correct
if [ $# -ne 4 ]; then
    echo "Usage: $0 blast_type /path/query /path/selected_species number_of_hits"
    exit 1
fi

blast_type=$1
query_path=$2
db_path=$3
number_of_hits=$4

blast_types=("blastn" "blastx" "blastp" "tblastx" "tblastn") # Initialize the valid types of blast


# Check if the blast type is valid
valid_blast_type=false	# Initialize the valid_blast_type to FALSE
for type in "${blast_types[@]}"; do		# Check all the valid blast types
    if [ "$type" == "$blast_type" ]; then	# If the users blast_type is equal to a valid type from blast_types
        valid_blast_type=true			# Change the valid_blast_type variable to TRUE
        break
    fi
done

if ! $valid_blast_type; then
    echo "Invalid blast type"
    exit 1
fi


# Check if the query file exists
if [ ! -f "$query_path" ]; then
    echo "File $query_path not found."
    exit 1
fi

# Check if the database directory exists
if [ ! -d "$db_path" ]; then
    echo "Directory $db_path not found."
    exit 1
fi


echo "Starting BLAST process with type $blast_type"

mkdir -p blast_results/out  blast_results/tabular blast_results/fasta blast_results/first_hits # Create a directory & subdirectories to save the ouput of blast

# Iterate over files in the directory
for file in $(ls "$db_path"*.* | cut -d. -f1 | sort -u); do
	filename=$(basename "$file")  # Save the name of the database without its path
	echo "START"
	echo "Running BLASTP for $filename"
	$blast_type -query $query_path -db $file -out blast_results/out/$filename.out -max_target_seqs $number_of_hits  # Run blast &  save the results  in .out format
	$blast_type -query $query_path -db $file -out blast_results/tabular/$filename.tsv -max_target_seqs $number_of_hits -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" # Run blast &  save the results  in .tsv format
	awk '{print $2}' blast_results/tabular/$filename.tsv | sort | uniq > blast_results/hit_ids.txt
	blastdbcmd -db $file -entry_batch blast_results/hit_ids.txt -out blast_results/fasta/$filename.fasta
	#$blast_type -query $query_path -db $file -out blast_results/tabular/$filename.tsv -max_target_seqs $number_of_hits -outfmt 6
	if [ -s blast_results/tabular/$filename.tsv ]; then  # Check if results file is not empty
 		#awk 'NR==1 {print ">" $2 "\n" $13}' blast_results/tabular/$filename.tsv > blast_results/first_hits/first_hit_$filename.fasta  # Extract the first hit
		awk -v file_name="$filename" 'NR==1 {print ">" $2 "_" file_name "\n" $13}' blast_results/tabular/$filename.tsv > blast_results/first_hits/first_hit_$filename.fasta  # Extract the first hit
    		echo "First hit saved to blast_results/first_hits/first_hit_$filename.fasta"
	else
    		echo "No hits found"
	fi
	#$blast_type -query $query_path -db $file -out blast_results/xml/$filename.xml -outfmt 5 # Run blast &  save the results in .xml format
	echo "DONE"
done


