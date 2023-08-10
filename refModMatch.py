

import pysam 
import sys
import numpy as np 
import gzip
import subprocess
import shlex
import tqdm
from array import array


bam = pysam.AlignmentFile(sys.argv[1],'rb')
header = bam.header.to_dict()

output_bam_name = sys.argv[2]
output_bam = pysam.AlignmentFile(output_bam_name,"wb", header = header)

def getMatchedModBases(read):

	reverseComplement = {'C':'G',
						 'A':'T',
						 'G':'C',
						 'T':'A'}

	forward_mods = read.modified_bases_forward
	# get modified bases

	aligned_pairs = np.array(read.get_aligned_pairs(with_seq=True,matches_only=False))
	# get the sequence to find mismatches, and we want non-matches (indels) as well 
	aligned_pairs = aligned_pairs[aligned_pairs[:,0] != None ,:]
	none_idx = aligned_pairs[:,2] == None 
	# index of Nones (indel)

	aligned_pairs[:,2][none_idx] = "n" 
	# replace None with n for string operations 

	seq_len = read.infer_read_length()

	# print(aligned_pairs[:,0])

	if read.is_reverse:
		none_idx = aligned_pairs[:,0] == None
		aligned_pairs[:,0][~none_idx] = np.abs(aligned_pairs[:,0][~none_idx] - seq_len) - 1


	aligned_pairs = aligned_pairs[np.argsort(aligned_pairs[:,0])]

	aligned_pairs_forward = aligned_pairs[:,0].T

	orientation = "+" 
	
	if read.is_reverse:
		orientation = "-"

	modification_val_dict = {}	
	
	for mod in forward_mods:


		mod_forward_positions = np.array(forward_mods[mod])[:,0].T.astype(int)

		values, mod_idx, aligned_pairs_idx = np.intersect1d(mod_forward_positions, 
															aligned_pairs_forward, 
															assume_unique=True, 
															return_indices = True)
		
		correctBase = mod[0]
		if read.is_reverse:
			correctBase = reverseComplement[mod[0]]

		# Aligned pairs idx represents all the values encoded in the MM Tag
		# we want the index of aligned_pairs_idx where the correct base is encoded

		correct_idx = np.argwhere(aligned_pairs[:,2][aligned_pairs_idx] == correctBase).T[0]

		modification_val_dict[mod[0] + "+" + mod[-1]] = correct_idx
		# correct_idx will let us directly mask the delta-encoded MM tag
		# and quickly modify the positions

	return modification_val_dict

def convertMM(encoding,correct_idx):

	 
	if len(correct_idx) > 0:

		cumsum_correct = np.cumsum(encoding + 1)[correct_idx]

		new_encoding = np.insert(np.diff(cumsum_correct),0,cumsum_correct[0]) - 1
		return new_encoding
	
	else:
		
		return None
	




def modifyMMAndMLTags(read, mod_idx_dict):

	original_mm = read.get_tag('MM').split(';') [:-1]
	# gets us the actual array we will need to index, MM
	
	original_ml = np.array(read.get_tag('ML'))
	# get us the actual array of probabilities we need to index, ML
		
	new_mm_list = []

	ml_start_counter = 0
	new_ml = []
	
	for mod_str in original_mm:
		
		mod_arr = mod_str.split(',')
		if len(mod_arr) > 1:
			

			correct_idx = mod_idx_dict[mod_arr[0]]
			
			encoding = np.array(mod_arr[1:]).astype(int)

			new_ml.append(original_ml[correct_idx + ml_start_counter])

			ml_start_counter += len(encoding)

			new_encoding = convertMM(encoding, correct_idx)

			if new_encoding is not None:

				new_mm_list.append(','.join([mod_arr[0]] + list(new_encoding.astype(str))))
			else:
				new_mm_list.append(mod_arr[0])
		else:
			new_mm_list.append(mod_arr[0])



	new_mm = ';'.join(new_mm_list) + ';'
	if len(new_ml) == 0:
		new_ml = None
	
	else:
		new_ml = np.concatenate(new_ml)	 
	
	return new_mm, new_ml

for read in tqdm.tqdm(bam):

	if read.has_tag('MM'):

		mod_idx_dict = getMatchedModBases(read)

		new_mm, new_ml = modifyMMAndMLTags(read, mod_idx_dict)
		
		if new_ml is not None:
			read.set_tag('ML' ,array('B', new_ml))

		read.set_tag('MM', new_mm, "Z")
		# print('---')


	output_bam.write(read)

		# tag_idx = getIndexOfTagValsToRemove(read, bad_read_positions)

		# For generating actual MM array we can try fast mm manipulation or if doesen't work 
		# we can do quick sequence manipulation from forward_sequence
