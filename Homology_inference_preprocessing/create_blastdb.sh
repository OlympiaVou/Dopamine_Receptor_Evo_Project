#!/bin/bash

# Check if the argument is provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 directory dbtype"
    exit 1
fi

# Check if the directory exists
if [ ! -d "$1" ]; then
    echo "Directory $1 not found."
    exit 1
fi


path=$1
type=$2 # This argument can be prot or nucl

# Iterate over files in the directory
for file in "$path"/*.fa; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        species="${filename%.fa}"
        makeblastdb -in "$file" -dbtype $type -out "${species}_db" -parse_seqids
    fi
done


