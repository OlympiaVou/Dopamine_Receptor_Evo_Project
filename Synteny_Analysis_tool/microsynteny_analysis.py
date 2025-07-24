'''
Microsynteny analysis: All-by-all Pairwise Comparisons
Calculation of syntenic scores and mapping of the synteny in a graph based on the orthogroups computed by OrthoFinder

    Parameters:
        gffs (directory): Containing all the gff files needed for the analysis created by the microsynteny preprocessing
        species order (.txt file): Containing one species name per row (this indicates the order the species will appear in the plots)
	species names (.csv file): Containing one scientific species name per row, and the name it will appear in the plot separated by a ; (eg. Homo_sapiens;Human)
	orthogroups file (.tsv file): Containing the Orthogroup IDs and the protein IDs in each Orthogroup, per species
    Returns:
        Heatmap with Synteny Scores
	Graph of the Syntenic Map (optional)
'''

#syntenic_map_plot = 0 	# 0: NO Syntenic Map, 1: YES Syntenic Map
#window = 50		# Set the gene window using for the analysis
#highlight = window + 1	# Set the highlighted gene in the Syntenic Map
#a = 0.5			# Set the coefficient of the two scores in the combined score (combined = a*pres_abs + (1-a)*lcm_fr)

# Import the libraries that will be needed
import argparse
import pandas as pd
import os
from microsynteny_analysis_functions import *
import matplotlib.patches as mpatches
import seaborn as sns
import matplotlib.pyplot as plt
import configparser

# Load the config file
config = configparser.ConfigParser()
#config.read('~/master_thesis/scripts/Synteny_Analysis_tool/microsynteny_analysis_config.ini')
config_path = '/home/olympia/master_thesis/scripts/Synteny_Analysis_tool/microsynteny_analysis_config.ini'

if not os.path.exists(config_path):
    print("File not found:", config_path)
else:
    config.read(config_path)
    print("Sections found:", config.sections())
# Read values from the [Settings] section
syntenic_map_plot = config.getint('Settings', 'syntenic_map_plot')
window = config.getint('Settings', 'window')
highlight = config.getint('Settings', 'highlight')
a = config.getfloat('Settings', 'a')


# Parse the input files
parser = argparse.ArgumentParser()
parser.add_argument("gffs", type=str)
parser.add_argument("species_order", type=str)
parser.add_argument("species_names", type=str)
parser.add_argument("orthogroups_file", type=str)

args = parser.parse_args()
gff = args.gffs
orthogroups = args.orthogroups_file
sp_names_file = args.species_names
org_order_file = args.species_order


'''
---------- STEP 1 ----------
Read and save the species order, and the Orthogroup IDs in a dictionary
'''
# Read the species order
org_order_file_read = open(org_order_file, "r")
org_order = org_order_file_read.read().splitlines()

# Read the species names
sp_names_file_read = open(sp_names_file, "r")
sp_names = sp_names_file_read.read().splitlines()
sp_names_dict = dict(line.split(';') for line in sp_names)

# Save the protein IDs and the matching Orthogroup IDs in a dictionary
id_to_og = {}
with open(orthogroups, newline='') as f:
	for line in f:
		parts = line.strip().split('\t')
		og_id = parts[0]
		for prot_id in parts[1:]:
			prots = prot_id.strip().split(', ')
			for prot in prots:
				id_to_og[prot] = og_id



'''
---------- STEP 2 ----------
Read the gff files and save them in a DataFrame
'''
# Initialize the DataFrame
data = pd.DataFrame()
sp_direction = pd.DataFrame()

# For each gff file separately do the following
for org in org_order:
	# Save the common species name
	file = gff + org
	species_name = os.path.basename(file).split('.')[0]
	target_species_name = species_name.split('_')[0]
	sp_name = sp_names_dict[target_species_name] 
	paralog = species_name.split('_')[2]
	
	# Initialize the lists for the proteinID and the orientation of the gene
	ids = []
	direction = []

	# Read the file
	with open(file) as f:
		for line in f:
			parts = line.strip().split()
			if len(parts) >= 5:
				ids.append(parts[1])  # Get the second column (protein ID)
				direction.append(parts[4]) # Get the fifth column (gene orientation)

	#match = [id_to_og.get(prot_id,'NoID') for prot_id in ids] 	# Match the protein IDs with the Orthogroup IDs
	match = ["END" if prot_id == "END" else id_to_og.get(prot_id if not prot_id.startswith("ENS") else f"{prot_id}.1","NoID") for prot_id in ids] # Match the protein IDs with the Orthogroup IDs
	g_dir = [ 1 if i=="+" else 0 for i in direction] # Create a list with the orientation of all the genes in binary format

	data[f"{sp_name}_{paralog}"] = match # Add the Orthogroup IDs to the total data DataFrame
	sp_direction[f"{sp_name}_{paralog}"] = g_dir # Add the direction of all the genes to the total direction DataFrame


'''
---------- STEP 3 ----------
Compute the Syntenic Scores
'''
lcs_scores = pd.DataFrame( index=data.columns, columns=data.columns) # Initialize a dataframe to save the lcs scores
pres_abs_scores = pd.DataFrame( index=data.columns, columns=data.columns) # Initialize a dataframe to save the jaccard scores
combined_scores = pd.DataFrame( index=data.columns, columns=data.columns) # Initialize a dataframe to save the combined scores

# Calculate all-by-all pairwise synteny scores
for i in data.columns:
	for j in data.columns:
		scores = lcs_pa(list(data[i]), list(data[j])) # Compute the scores using the functions imported
		lcs = scores[0]
		pres_abs_score = scores[1]

		lcs_scores.loc[i,j] = float(lcs) # Add the lcm score to the lcm matrix
		pres_abs_scores.loc[i,j] = float(pres_abs_score) # Add the jaccard score to the jaccard matrix

		combined_scores.loc[i,j] = a*float(pres_abs_score) + (1-a)*float(lcs) # Compute and add the combined score to the combined matrix

# Transform all the matrices to be matrices of floats
lcs_scores = lcs_scores.astype(float)
pres_abs_scores = pres_abs_scores.astype(float)
combined_scores = combined_scores.astype(float)



'''
---------- STEP 4 ----------
Plot the Syntenic Scores
'''
# Create the heatmap
plt.figure(figsize=(15, 17))
#sns.heatmap(combined_scores, cmap='rocket', annot=True, fmt='.2f', square=True)
sns.heatmap(combined_scores, cmap='rocket', annot=False, fmt='.2f', square=True)
plt.xticks(rotation=45, ha='right')
plt.title('Syntenic Score Heatmap')
plt.subplots_adjust(top=0.92, bottom=0.2)
plt.show()
plt.close()


'''
---------- STEP 5 ----------
Plot the Syntenic Map
'''
if syntenic_map_plot: # Check if the user want the syntenic map

	# Replace OGs that appear in only one species with 'NoID'
	all_values = data.values.ravel()
	og_counts = pd.Series(all_values).value_counts()
	singleton_ogs = og_counts[og_counts == 1].index

	# Replace singleton OGs in the DataFrame with 'NoID'
	for col in data.select_dtypes(include='object'):
		data[col] = data[col].map(lambda x: 'NoID' if x in singleton_ogs else x)

	# Get list of the columns
	target_cols = data.columns[:]

	# Unique RefIDs and color mapping
	all_ids = pd.unique(data[target_cols].values.ravel())
	ref_ids = [rid for rid in all_ids if rid != 'NoID' and rid != 'END']


	palette = sns.color_palette("hls", len(ref_ids))
	id_to_color = {rid: color for rid, color in zip(ref_ids, palette)}

	# Plot setup
	fig, ax = plt.subplots(figsize=(16, 9))

	# Plot the syntenic map using the imported functions
	synteny_blocks(data, sp_direction, target_cols, id_to_color, ax)
	synteny_connecting_lines(data, target_cols, id_to_color, ax)

	# Highlight the gene of interest
	if highlight:
		ax.axvline(x=highlight, color='red', linestyle=':', linewidth=1.5)

	legend_handles = [mpatches.Patch(color=color, label=ref_id) for ref_id, color in id_to_color.items()]
	legend_handles.append(mpatches.Patch(color='lightgray', label='NoID'))

	# Plot visualization
	ax.set_yticks(range(len(target_cols)))
	ax.set_xticks([i for i in range(1, len(data) + 2) if i%2==1])
	ax.set_yticklabels(target_cols)
	ax.set_xlim(0, len(data)+1)
	ax.set_ylabel("Species")
	ax.set_xlabel("Relative Position of Genes")
	ax.set_title("Synteny Map (10-gene window)")
	#ax.legend(handles=legend_handles, title="Reference IDs", loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=10, frameon=False)
	#fig.subplots_adjust(bottom=0.25)  # Adjust this if legend is still clipped
	plt.tight_layout()
	plt.show()
	plt.close()


