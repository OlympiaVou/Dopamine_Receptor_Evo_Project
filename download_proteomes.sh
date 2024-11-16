#!/bin/bash


# Check if the number of arguments provided is exactly one
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 path/to/Uniprot_ID_file"
    exit 1
fi

UP_file="$1"

# Check if the file specified by UP_file exists
if [ ! -f "$UP_file" ]; then
    echo "File not found: $UP_file"
    exit 1
fi

while read -r word; do
    # Use wget to download the fasta file for each UniProt ID
    wget -O "${word}.fasta.gz" "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28proteome%3A${word}%29%29"
done < "$UP_file"

# Unzip all downloaded gzip files
gunzip UP*.gz
