
# Import libraries
from Bio import Entrez
from Bio import SeqIO
import argparse

# Parse the input file with the IDs
parser = argparse.ArgumentParser()
parser.add_argument("IDs", help="Give the path to the input file")
parser.add_argument("input_type",choices=["prot","nucl"], help="Specify the input type: 'prot' OR 'nucl'")
parser.add_argument("output_type",choices=["prot", "nucl"], help="Specify the output type: 'prot' OR 'nucl'")
parser.add_argument("-o","--out", help="Give the path to the output file")
parser.add_argument("email", help="Give your email")
args = parser.parse_args()
IDs_file = args.IDs
in_type = args.input_type
out_type = args.output_type
email = args.email

if args.out:
	out_file = args.out
	print (f"The sequences will be saved to ", args.out)
else:
	out_file = "./IDs_sequences.fasta"
	print (f"The sequences will be saved to ", out_file)



# NCBI requires an email address to monitor tool usage
Entrez.email = email


IDs = [] # Initialize an empty list to save the IDs

with open(IDs_file, 'r') as file:
	for line in file:
		line = line.rstrip('\n')
		if line and line[0] != "E":
			IDs.append(line)

print(f"Number of sequences to DOWNLOAD:", len(IDs))

if out_type == "nucl":
    db = "nucleotide"
elif out_type == "prot":
    db = "protein"

# Fetch the sequences from NCBI
#handle = Entrez.efetch(db="protein", id=",".join(IDs), rettype="fasta", retmode="text")
handle = Entrez.efetch(db=db, id=",".join(IDs), rettype="fasta", retmode="text")


# Save to a file
with open(out_file, "w") as fasta_file:
    fasta_file.write(handle.read())
handle.close()


