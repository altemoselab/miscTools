

import pysam 
import sys
import numpy as np 
import gzip
import subprocess
import shlex
import tqdm

# remove modifications that are not reference matched
# requires MD tag
# all base modifications must be in reference to the forward strand (0)


bam = pysam.AlignmentFile(sys.argv[1],'rb')
header = bam.header.to_dict()

output_bam_name = sys.argv[2]
output_bam = pysam.AlignmentFile(output_bam_name,"wb", header = header)



for read in tqdm.tqdm(bam.fetch()):

	if read.has_tag('MM'):

		aligned_pairs = np.array(read.get_aligned_pairs(with_seq=True,matches_only=False))

		# find all indel/padding positions
		none_idx = aligned_pairs[:,2]==None

		# replace those with lowercase string for vectorized operation
		aligned_pairs[:,2][none_idx] = "n"
		# sub idx now reperesents all bad positions based on lowercase search, as per pysam spec
		sub_idx = np.char.islower(aligned_pairs[:,2].astype(str)) 
		
		bad_read_positions = aligned_pairs[:,0][sub_idx].astype(int)

		seq_len = read.infer_read_length()

		if read.is_reverse:
			bad_read_positions = np.sort(np.abs(bad_read_positions - seq_len) - 1)

		# bad read positions is based on query alignment 
		# positions so these would need to be adjusted for alignment orientation

		# can use pysam's read_modified_bases forward to get the forward 
		# positions for every modification
		# filter those and recieve the index's that pass
		# to then grab the positions of 

		forward_mods = read.modified_bases_forward # gets the forward coord of the modifications
		# print(read.get_tag('MM'))
		original_mm = read.get_tag('MM')   # gets us the actual array we will need to index, MM
		original_ml = read.get_tag('ML')  # get us the actual array of probabilities we need to index, ML
		
		original_mm = original_mm.split(';')[:-1]

		for mm in original_mm:
			split_mm = mm.split(',')
			base_mod = split_mm[0]
			forward_mod_positions = np.array(forward_mods.get((base_mod[0],0,base_mod[-1])))
			mm_positions = np.array(split_mm[1:]).astype(int)
			if forward_mod_positions.size > 1 and len(mm_positions) > 0:

				forward_mod_positions = forward_mod_positions[:,0]

				forward_mod_to_remove = np.isin(forward_mod_positions, bad_read_positions) # removing True, keeping False
				if np.sum(forward_mod_to_remove) > 0:

					forward_mod_to_remove_diff_encoding = np.diff(np.hstack([[False],forward_mod_to_remove,[False]]).astype(int)) # value of 1 indicates beginning of stretch
																					# value of -1 indicates end of stretch (remove becomes False)
																					# value of 0 means continuation
			
					# this is in "python index" so just use [start:stop], stop is also 
					# the value that you add the new delta too
					idx_remove_start = np.where(forward_mod_to_remove_diff_encoding == 1)[0]
					idx_remove_stop = np.where(forward_mod_to_remove_diff_encoding == -1)[0]

					if idx_remove_stop[-1] == len(forward_mod_to_remove):
						idx_remove_start = idx_remove_start[:-1]
						idx_remove_stop = idx_remove_stop[:-1]


					tmp_mm = mm_positions.copy()

					ranges_idx = range(len(idx_remove_start))
					print(forward_mod_to_remove)
					for i in ranges_idx:

						start, stop = int(idx_remove_start[i]), int(idx_remove_stop[i])
						new_val = tmp_mm[stop] + np.sum(mm_positions[start:stop] + 1)

						tmp_mm[stop] = mm_positions[stop] + np.sum(mm_positions[start:stop] + 1)
					new_mm = tmp_mm[forward_mod_to_remove==False]

					print(mm_positions)
					print(new_mm)
					print('---')




				# to change the MM positiong it is sum(lost values + 1) of positions because of 0 based delta calculation







				








		# break















		# need to map the MM format onto reference / genomic coordinates
		# what would be fastest is to get all lowercase values + None values 
		# and then remove MM coords corresponding to those values

		






		











	
	output_bam.write(read)






	# for read in tqdm(infile.fetch()):
	
	# if read.has_tag('MM'):
		
	# 	for_seq = read.get_forward_sequence().upper()

	# 	T_count = for_seq.count("T")

	# 	mod_read = read.modified_bases_forward

	# 	analog_count = 0

	# 	for mod in mod_read:
	# 		if mod[0] == "T":
	# 			if mod[2] == 'e' or mod[2] == 'b':

	# 				mm_ml = np.vstack(mod_read[mod]).T

	# 				analog_count += np.sum(mm_ml[1] > prob_cutoff)

	# 	print(analog_count / T_count)



	# 	# outfile.write(read)




