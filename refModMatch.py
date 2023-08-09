

import pysam 
import sys
import numpy as np 
import gzip
import subprocess
import shlex
import tqdm
from array import array

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
		none_idx = aligned_pairs[:,2] == None

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

		new_ml = []
		ml_counter = 0
		all_mm = []
		mod_type_tracker = []
		for mm in original_mm:
			
			split_mm = mm.split(',')
			base_mod = split_mm[0]
			forward_mod_positions = np.array(forward_mods.get((base_mod[0],0,base_mod[-1])))
			mm_positions = np.array(split_mm[1:]).astype(int)
			mod_type_tracker.append(base_mod)

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
					for i in ranges_idx:

						start, stop = int(idx_remove_start[i]), int(idx_remove_stop[i])
						new_val = tmp_mm[stop] + np.sum(mm_positions[start:stop] + 1)

						tmp_mm[stop] = mm_positions[stop] + np.sum(mm_positions[start:stop] + 1)
					new_mm = tmp_mm[forward_mod_to_remove==False]
					if len(new_mm) == 0:
						output_bam.write(read)
						continue
					all_mm.append(new_mm)
					

					## new_MM is good and represents the new mms we want to use
					# now we need to index the ML tag using the forward_mod_to_remove 
					# and keep count of how many positions we consume 
					# to then filter for the next modification type

					tmp_ml = np.array(original_ml[ml_counter:ml_counter + len(forward_mod_to_remove)])
					new_ml.append(tmp_ml[forward_mod_to_remove==False])
					
					ml_counter =+ len(forward_mod_to_remove)


		if len(new_ml) < 2:
			output_bam.write(read)
			continue
		new_ml = np.concatenate(new_ml)

		output_ml = list(new_ml.astype(int))
		read.set_tag('ML' ,array('B', output_ml))
		
		output_mm_list = []
		
		for mod_type, mod_mm_positions in zip(mod_type_tracker, all_mm):
			output_mm_list.append(mod_type + "," + ','.join(mod_mm_positions.astype(str)))


		output_mm_string = ';'.join(output_mm_list) + ";"

		read.set_tag('MM', output_mm_string, "Z")


					## new_MM is good and represents the new mms we want to use
					# now we need to index the ML tag using the forward_mod_to_remove 
					# and keep count of how many positions we consume 
					# to then filter for the next modification type





				# to change the MM positiong it is sum(lost values + 1) of positions because of 0 based delta calculation







				








		# break















		# need to map the MM format onto reference / genomic coordinates
		# what would be fastest is to get all lowercase values + None values 
		# and then remove MM coords corresponding to those values

		






		











	
	output_bam.write(read)






