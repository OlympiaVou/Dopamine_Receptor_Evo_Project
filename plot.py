''' This script takes a txt file with species names and saves their IDs.
INPUT:  1) .txt file - with species names
        2) the NCBI file speclist.txt
OUTPUT : A .txt file with the selected species IDs
'''


###################### STEP 0 ######################
# Import libraries & Save the files names using argparse

import matplotlib.pyplot as plt
import numpy as np
import argparse
parser = argparse.ArgumentParser()

# Save the scores file name
parser.add_argument("file_name", type=str, help="Display a .txt file that contains what you want to plot")
# Save the index file name
#parser.add_argument('--i',"index_file", type=str, help="Display a .txt file the indeces of the points you want to highlight")
args = parser.parse_args()

score_file = args.file_name
#index_file = args.index_file



###################### STEP 1 ######################
# Save the scores in a list

with open(score_file, 'r') as file:			# Open the scores file in reading mode
	scores = []					# Initialize the list with the scores
	for line in file:				# For each line consisting a score
		line=line.strip()
		line = line.replace(" ", "")
		scores = line.split(',')[1:-1]
		scores = [float(i) for i in scores]
		#scores.append(float(line.strip()))	# Add the score to the list
		#scores.append(line)      # Add the score to the list

#print(scores)



###################### STEP 2 ######################
# Save the indeces to highlight in a list
'''
with open(index_file, 'r') as file:                     # Open the scores file in reading mode
        indeces = []                                     # Initialize the list with the scores
        for line in file:                               # For each line consisting a score
                indeces.append(int(line.strip()))      # Add the score to the list

print(indeces)
'''


###################### STEP 3 ######################
# Plot the scores

#scores_np = np.array(scores)
scores_ind = [i for i in range(len(scores))]
#scores_ind = np.array([i for i in range(len(scores))])

#print(scores_np)
#print(scores_ind)

plt.plot(scores_ind, scores)
plt.title("Alignment scores")
plt.xlabel("Residues")
plt.ylabel("Alignment score")
#for i in indeces:
#	plt.plot(i, scores[i], 'o', color='red')
plt.show()
