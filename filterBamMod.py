import pysam 
from pybedtools import BedTool
import sys
import numpy as np 
import re
from array import array
from tqdm import tqdm
# filter MM / ML tags based on scaled probability distribution


input_bam = sys.argv[1] # sorted bam w/ index
output_bam = sys.argv[2] # unsorted bam output

# argv[4] is 0 - 255 inclusive, decimal value to be transformed to 255 space

mod_tags = [str(x) for x in sys.argv[3].split(',')]
probability_cutoffs = [int(x) for x in sys.argv[4].split(',')]


infile = pysam.AlignmentFile(input_bam, 'rb', 
								threads = 1, check_sq=False)

outfile = pysam.AlignmentFile(output_bam, 
								mode = 'wb',
								template = infile, 
								threads = 1)

def getDefinedMod(mod_tags,read):	
	
	defined_mods = []
	for mod_tag in mod_tags:	
		if mod_tag == 'EdU':
			base = "T"
			mod  = "e"
		elif mod_tag == "BrdU":
			base = "T"
			mod = "b"
		elif mod_tag == "5mC":
			base = "C"
			mod = "m"
		elif mod_tag == "6mA":
			base = "A"
			mod ="a"
		else:
			print('''Modified Base not supported. \n Supported base-mods: Edu,BrdU,5mC,6mA. \n 
				Proper usage: python3 filterBaMod.py 
				bam output_name modification probability_cutoff( int 0-255 )''' )



		# if read.is_reverse:
		# 	direction = '-'
		# else:
		# 	direction = "+"

		defined_mod = base + "+" + mod
		defined_mods.append(defined_mod)

	return defined_mods



for read in tqdm(infile.fetch(until_eof=True)):
	if read.has_tag('MM'):
		
		mods = read.get_tag("MM", with_value_type=True)
		mods_prob = read.get_tag("ML", with_value_type=True)

		defined_mods = getDefinedMod(mod_tags,read)


		ml = mods_prob[0]
		mm = mods[0]

		ml = np.array(ml)

		region_specific_start_index = 0

		new_mods = []
		new_mods_prob = []


		for mm_specific_mod in mm.split(';')[:-1]:


			mm_specific_mod_arr = mm_specific_mod.split(',')
			
			mm_specification = mm_specific_mod_arr[0]
			
			mm_mod_positions = np.array(mm_specific_mod_arr[1:]).astype(int)


			number_mods = len(mm_mod_positions)


			#(start, stop) 0 - based, start-included, end excluded
			mod_specific_ml_coords = (region_specific_start_index,region_specific_start_index + number_mods)
			if mm_specification in defined_mods:

			# 	# cumulative_sum(MM_Tag + 1 ) represents the number of preceding Ts
			# 	# take the numpy difference there and then you get the number since 

				probability_cutoff = probability_cutoffs[defined_mods.index(mm_specification)]

				idx = ml[mod_specific_ml_coords[0]:mod_specific_ml_coords[1]] >= probability_cutoff

				if np.sum(idx) > 0:

					new_mod_positions_cumsum = np.cumsum((mm_mod_positions + 1))[idx]

					new_mod_positions = np.diff(new_mod_positions_cumsum) - 1 
					new_mod_positions = np.insert(new_mod_positions, 0, new_mod_positions_cumsum[0] - 1)

					new_ml = ml[mod_specific_ml_coords[0]:mod_specific_ml_coords[1]][idx]

					
					# generate new string to be added to new mods and new mods prob 

					new_mods.append(mm_specification + "," + ','.join(new_mod_positions.astype(str)))
					new_mods_prob.append(new_ml.astype(int))

			else:
				
				new_mods.append(mm_specific_mod)
				new_mods_prob.append(ml[mod_specific_ml_coords[0]:mod_specific_ml_coords[1]])

				# now we know where we are in the arr, and we clip down 
				# the ml array to be only our region of interest
			region_specific_start_index += number_mods



		new_mods_string = ';'.join(new_mods) + ';'

		if len(new_mods_prob) == 0:
			read.set_tag("MM",None)
			read.set_tag("ML",None)
		# new_tags = [('MM', new_mods, "Z"), ('ML' ,list(np.concatenate(new_mods_prob)), 'C')]
		else:	
			output_ml = list(np.concatenate(new_mods_prob).astype(int))
			# print(output_ml)
			# exit()
			read.set_tag('MM', new_mods_string, "Z")
			read.set_tag('ML' ,array('B', output_ml))


		outfile.write(read)
