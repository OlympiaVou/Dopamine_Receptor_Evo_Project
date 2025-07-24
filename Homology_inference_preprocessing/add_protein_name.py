'''
This script adds the name of the protein to the ID of each sequence of a fasta/multi-fasta file.
INPUT:
1. A string with the name of the protein
2. A fasta/multifasta file
3. The name of the output file

OUTPUT: A new fasta file with the name of the protein.
'''

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("protein_name", type=str)   # Save the protein name
parser.add_argument("fasta", type=str)          # Save the input fasta file
parser.add_argument("output_file", type=str)    # Save the output fasta file
args = parser.parse_args()

input = args.fasta				# Save the input file name as input
out = args.output_file 				# Save the output file name as out
protein = args.protein_name 			# Save the protein name as protein

from Bio import SeqIO # Library to process fasta files

records = list(SeqIO.parse(input, 'fasta')) # Read and save all the fasta entries of the input file

with open(out, 'w') as output:              						# Open the output file
	for record in records:								# For each species in the list of species
		output.write(f'>{record.description}_{protein}\n{record.seq}\n')  	# Write the entry in the output file
