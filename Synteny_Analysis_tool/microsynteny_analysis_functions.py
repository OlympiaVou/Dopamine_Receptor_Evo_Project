
def synteny_connecting_lines(data, target_cols, id_to_col, ax):
	# For each RefID (excluding 'NoID'), collect coordinates and connect them
	refid_coords = {}

	# Build reverse lookup of where each RefID appears
	for i, col in enumerate(target_cols):
		for j, ref_id in enumerate(data[col]):
			if ref_id != "NoID":
				refid_coords.setdefault(ref_id, []).append((j, i))  # (x=index, y=species row)

	# Plot lines connecting the same RefID across species
	for ref_id, coords in refid_coords.items():
		if ref_id not in id_to_col:
			continue  # Skip anything unmapped

		if len(coords) >= 2:
			coords = sorted(coords, key=lambda x: x[1])  # Sort by species row (top to bottom)
			for k in range(len(coords) - 1):
				x1, y1 = coords[k]
				x2, y2 = coords[k + 1]
				line_color = id_to_col[ref_id]
				ax.plot([x1+1, x2+1], [y1, y2], color=line_color, linewidth=1.5, alpha=0.9, zorder=0)


def synteny_blocks(data, rel_direction, target_cols, id_to_color, ax):
	for i, col in enumerate(target_cols):
		y_pos = i  # Position of each species in the y-axis
		# Add a black baseline for each row considerring the end of the scaffold/chromosome
		non_end_indices = [i for i, val in enumerate(data[col]) if val != 'END']
		if non_end_indices:
			xmin = min(non_end_indices) + 0.5
			xmax = max(non_end_indices) + 1.5
			ax.hlines(y=y_pos, xmin=xmin, xmax=xmax, color='black', linewidth=1)
		
		#ax.hlines(y=y_pos, xmin=0.5, xmax=len(data[data[col] != 'END'])+0.5, color='black', linewidth=1)
		
		sp_dir = rel_direction[col]
		for j, ref_id in enumerate(data[col]):
			if sp_dir[j]==0:
				shape = ">"
			else:
				shape = "<"
			if ref_id != 'NoID' and ref_id != 'END':
				#ax.plot(j, y_pos, 's', color="blue", markersize=8, mec='black', mew=0.4)
				ax.plot(j+1, y_pos, marker=shape, color=id_to_color[ref_id], markersize=8, mec='black', mew=0.4)
			elif ref_id == 'NoID':
				ax.plot(j+1, y_pos, marker=shape, color="lightgray", markersize=8, mec='black', mew=0.4)
			elif ref_id == 'END':
				continue
	#ax.legend(handles=legend_handles, loc="upper right", bbox_to_anchor=(1.05, 1))


def lcs(sp1, sp2):
	'''
    Compute a modified Longest Common Subsequence (LCS) score between two sequences,
    considering both forward and reverse matches, excluding "NoID" elements.

    Parameters:
        sp1 (list): First sequence of genes (using their Orthogroup IDs).
        sp2 (list): Second sequence of genes (using their Orthogroup IDs).

    Returns:
        float: Combined LCS score, giving full weight to forward matches and half to reverse matches.
    '''
	if sp1 == sp2: # If the two sequences of genes are completely the same 
		forward_len_lcs = len(sp1) # The Longest Common Subsequence is equal to the length of the sequence
		reverse_len_lcs = 0 # Initialize the reverse Longest Common Subsequence as 0 to avoid mistakes
	else:
		# Initialize two matrices to store the lengths of the Longest Common Subsequence (LCS)
		forward_lcs_matrix = [[0 for x in range(len(sp1)+1)] for x in range(len(sp2)+1)]
		reverse_lcs_matrix = [[0 for x in range(len(sp1)+1)] for x in range(len(sp2)+1)]

		# Iterate through each element in sp1 and sp2
		for i in range(len(sp1)):
			for j in range(len(sp2)):
				# ---------- Forward LCS Computation ----------
				if sp1[i] == sp2[j] and sp1[i] != "NoID": # If the elements match and are  not "NoID", they are considered part of the LCS
					if i == 0 or j == 0:
						forward_lcs_matrix[i+1][j+1] = 1 # Initialize: if at the first row or column, LCS is 1
					else:
						forward_lcs_matrix[i+1][j+1] = forward_lcs_matrix[i][j] + 1 # Extend the previous LCS by 1
				else:
					# Carry forward the max value from the adjacent computed entries
					forward_lcs_matrix[i+1][j+1] = max(forward_lcs_matrix[i][j+1], forward_lcs_matrix[i+1][j])
				
				# ---------- Reverse LCS Computation ----------
				if sp1[i] == sp2[-j] and sp1[i] != "NoID": # If the elements match and are  not "NoID", they are considered part of the LCS
					if i == 0 or j == 0:
						reverse_lcs_matrix[i+1][j+1] = 0 # Initialize: if at the first row or column, LCS is 0, to avoid short matches
					else:
						reverse_lcs_matrix[i+1][j+1] = reverse_lcs_matrix[i][j] + 1 # Extend the previous LCS by 1
				else:
					# Carry forward the max value from the adjacent computed entries
					reverse_lcs_matrix[i+1][j+1] = max(reverse_lcs_matrix[i][j+1], reverse_lcs_matrix[i+1][j])
		
		# Get the final LCS lengths from both matrices
		forward_len_lcs = forward_lcs_matrix[-1][-1]
		reverse_len_lcs = reverse_lcs_matrix[-1][-1]

	# Combine the forward and reverse scores
	if reverse_len_lcs <= 1: # If the length of the reverse Longest Common Subsequence is 1 it needs to be ignored
		len_lcs = forward_len_lcs
	else: # If the length of the reverse Longest Common Subsequence is >1 it needs to be included but with a 0.75 coefficient
		len_lcs = forward_len_lcs  + reverse_len_lcs*0.75

	return len_lcs



def lcs_pa(s1,s2):
	# Calculate the lcm score & normalize it
	lcs_score = lcs(s1,s2)
	norm_lcs = lcs_score/max(len(s1),len(s2))

	# Calculate the jaccard score
	if s1 == s2:
		intersection = s1
		total = s1
	else:
		intersection = set(s1) & set(s2) # Find the union
		total = set(s1) | set(s2)
	pres_abs_score = len(intersection)/max(len(s1),len(s2)) # Compute the jaccard score normalizing by the maximum length
	#pres_abs_score = len(intersection)/len(total) # Compute the jaccard score normalizing by the maximum length

	return (norm_lcs, pres_abs_score)



