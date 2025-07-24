# This script selects the first n # of hits from a blast analysis & prints them as stdout in fasta format

# Check if the number of arguments is correct
if [ $# -ne 3 ]; then
    echo "Usage: $0 file_path_out file_path_fasta ID_file  > new_fasta.fasta"
    exit 1
fi

out_path=$1 # Path to the .out file of BLAST
fasta_path=$2 # Path to the multi-fasta file of BLAST
IDs=$3 # Path to the file with the IDs of seqs you want to keep

num_of_seqs=$(grep '>' $fasta_path | wc -l) # Find the # of the sequences in the fasta file


# Initialize the line the sequences' IDs start
start_ind=$(grep -n "Sequences" $out_path| awk '{print $1}'| awk -F ':' '{print $1}') # Find the line that starts with "Sequences" and save its index
start_ind=$((start_ind + 2)) # Add 2 to the index to have the index of he first hit's ID



current=1 # Initialize the line to start saving the IDs


for id in $(cat $IDs); do
	echo ${id}
	index_line=$(grep -n "${id}" $fasta_path | awk '{print $1}'| awk -F ':' '{print $1}') # Save the index of the line of the fasta with the saved ID
        index_line=$(( index_line + 1)) # Add 1 to the index in order to save the sequence
        awk -v start=$index_line 'NR >= start && /^>/ {exit} NR >= start {print}' "$fasta_path" # Print the lines from the index_line until the line that starts with '>'
done


#while [ $current -le $number_of_hits ]; do # While the current # is lower or equal to number_of_hits do the following
#	total_ID=$(awk -v start=$start_ind 'NR >= start && /^>/ {exit} NR >= start {print}' $out_path | head -$current | tail -1) # Save the line with info about the current hit
#	ID="${total_ID%% *}" # Save only the ID
#	grep ">${ID}" $fasta_path # Print the line of the fasta with the saved ID
#	index_line=$(grep -n ">${ID}" $fasta_path | awk '{print $1}'| awk -F ':' '{print $1}') # Save the index of the line of the fasta with the saved ID
#	index_line=$(( index_line + 1)) # Add 1 to the index in order to save the sequence
#	awk -v start=$index_line 'NR >= start && /^>/ {exit} NR >= start {print}' "$fasta_path" # Print the lines from the index_line until the line that starts with '>'
#	current=$(( current + 1)) # Add 1 to the current variable to continue
#done
