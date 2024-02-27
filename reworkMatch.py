

import pysam 
import sys
import numpy as np 
import gzip
import subprocess
import shlex
import tqdm
from array import array
import argparse
import pyfaidx

def inputArgs():
	'''Parse in arguments. '''
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-b','--bam', 
		type = str, 
		help = 'bam file')
	
	parser.add_argument('-o','--output_bam', 
		type = str, 
		help = 'output bam file')
	
	parser.add_argument('-m','--motif',
		type=str,
		help = ' string with motif, position of base. supports single base, CG, GC currently \n example: A,0;CG,0;GC,1')

	args = parser.parse_args()

	return args.bam, args.output_bam, args.motif


def parseMotif(input_motifs):

	split_input = input_motifs.split(";")

	all_motifs = []
	base_pos_in_motif = []


	for motif in split_input:

		split_motif = motif.split(',')

		all_motifs.append(split_motif[0])
		base_pos_in_motif.append(int(split_motif[1]))
	return all_motifs, base_pos_in_motif

def generateMMTag(mm_labels,mm_tags,ml_vals,read,tag_idxs):

	mm_string = str(read.get_tag("MM"))

	output_MM = []
	output_ML = []


	for mm_single_tag_string in mm_string.split(';')[:-1]:
		single_tag_list = mm_single_tag_string.split(',')
		tag_label = single_tag_list[0]
		og_mm_tag_cumsum = np.cumsum(np.array(single_tag_list[1:]).astype(int) + 1 )


		if 'C+h' in tag_label: # no hemi-methylation currently
			continue
		
		try:
			idx = mm_labels.index(tag_label)
		except:
			continue

		new_label = mm_labels[idx]

		# new_forward_pos = mm_tags[idx]

		new_ml = ml_vals[idx]

		mm_tag_idx = tag_idxs[idx]
		
		idxed_cumsum = og_mm_tag_cumsum[mm_tag_idx]

		if len(idxed_cumsum) == 0:
			continue


		new_mm = [str(idxed_cumsum[0])] + list((np.diff(idxed_cumsum) - 1).astype(str))



		tmp_mm_vals = ','.join(new_mm)
		tmpMM = str(tag_label) + ',' + tmp_mm_vals + ';'

		# print(tmpMM, new_ml)
		output_MM.append(tmpMM)
		output_ML.append(new_ml)
	
	if len(output_ML) == 0:
		return None, None

	MM = ''.join(output_MM)
	ML  = array('B', np.concatenate(output_ML))

	return MM, ML





def parseBam(bam, motifs, base_pos_in_motif, output_bam):


	parseA = False
	parseCG = False
	parseGC = False

	if 'A' in motifs:
		parseA = True
	if 'CG' in motifs:
		parseCG = True
	if 'GC' in motifs:
		parseGC = True

	for read in tqdm.tqdm(bam.fetch()):


		aligned_pairs = np.array(read.get_aligned_pairs(with_seq=True,matches_only=True))
		# get the sequence to find mismatches, and we want non-matches (indels) as well 
		aligned_pairs = aligned_pairs[aligned_pairs[:,0] != None ,:]
		none_idx = aligned_pairs[:,2] == None 
		# index of Nones (indel)
		aligned_pairs[:,2][none_idx] = "n" 

		# seq = np.array(list(read.get_forward_sequence()))
		seq_len = len(read.get_forward_sequence())
		
		mod_bases = read.modified_bases_forward
		# seq = np.array(list(read.get_forward_sequence()))

		reverse_state = read.is_reverse

		# here we will store the modified MM tags + ML tags
		# in the order they are processed and output
		# in order to properly format the MM and ML tag
		mm_labels = []
		mm_tags = []
		ml_vals = []
		tag_idxs = []

		for m in mod_bases:
			
			if 'h' in m: # not supporting hemi-methylation currently
				continue
			
			mod_stack = np.vstack(mod_bases[m]).T
			
			if reverse_state:
				mod_stack[0]= np.abs(mod_stack[0] - seq_len) - 1

			ap_idx_of_all_mods = np.isin(aligned_pairs[:,0], mod_stack[0])


			# then we generate the forward idx, which we will use to parse the sequence
			# to define the 
			
			if m[0] == "A" and parseA:
				base = "A"
				if reverse_state:
					
					base = "T"
				
				forward_idx = aligned_pairs[:,0][ap_idx_of_all_mods][aligned_pairs[:,2][ap_idx_of_all_mods] == base].astype(int)

				if reverse_state:
					forward_idx = np.abs(seq_len - (forward_idx+ 1)) 

			
				new_label = 'A+a'


			if m[0] == "C":
				
				if parseCG or parseGC:

					base = "C"
					neighbor_base = "G"
					
					if reverse_state:
					
						base = "G"
						neighbor_base = "C"
					
					base_match_idx = aligned_pairs[:,2][ap_idx_of_all_mods] == base
					# need to clip 1 base of the edge of base_match_idx before combining					
					if parseCG :
						
						ap_idx_of_all_mods_plus_one = (np.argwhere(ap_idx_of_all_mods).T[0] + 1)[1:-1]
						
						neighbor_idx = aligned_pairs[:,2][ap_idx_of_all_mods_plus_one] == neighbor_base

						CG_idx = (base_match_idx[1:-1]) & (neighbor_idx)

						ap_CG_idx = ap_idx_of_all_mods_plus_one[CG_idx] - 1 

					if parseGC : 
					
						ap_idx_of_all_mods_minus_one = (np.argwhere(ap_idx_of_all_mods).T[0] - 1)[1:-1]
						
						neighbor_idx = aligned_pairs[:,2][ap_idx_of_all_mods_minus_one] == neighbor_base

						GC_idx = (base_match_idx[1:-1]) & (neighbor_idx)

						ap_GC_idx = ap_idx_of_all_mods_minus_one[GC_idx] + 1 


					if parseCG and parseGC:

						forward_idx = aligned_pairs[:,0][np.unique(np.concatenate([ap_CG_idx,ap_GC_idx]))]
					elif parseCG:
						forward_idx = aligned_pairs[:,0][ap_CG_idx]
					elif parseGC:
						forward_idx = aligned_pairs[:,0][ap_GC_idx]
					else:
						print('Modification not supported')
						return None

					forward_idx = forward_idx.astype(int)
					if reverse_state:
						forward_idx = np.abs(seq_len - (forward_idx+ 1)) 
					
					new_label = 'C+m'

			forward_idx.sort()
			
			# forward idx is the index of the base in the forward sequence
			# of the read ( original read as it comes of sequencer )

			tag_idx = np.isin(mod_stack[0],forward_idx) 
			new_ml_tags = mod_stack[1][tag_idx]

			mm_labels.append(new_label)
			mm_tags.append(forward_idx)
			ml_vals.append(new_ml_tags)
			tag_idxs.append(tag_idx)
		
		MM, ML = generateMMTag(mm_labels,mm_tags,ml_vals,read,tag_idxs)

		if MM == None or ML == None:
			# read with no MM or ML tag
			read.set_tag('ML' , None)

			read.set_tag('MM', None)
		else:

			read.set_tag('ML' , ML)

			read.set_tag('MM', MM, "Z")

		# generateMMTag(mm_labels,mm_tags,ml_vals,read,seq)
			#### 
			# forward index will no give us correct position for every base we want to keep in the MM tag
			# we now can use to filter MM tag and grab the correct ML tags 
			# it is in the forward coordinate space 
			# so it is always encoded as "+" in the MM spec 

			output_bam.write(read)






				








def main():

	bam_file, output_bam_name, input_motifs = inputArgs()

	motifs,base_pos_in_motif = parseMotif(input_motifs)
	motifs,base_pos_in_motif = ['A','CG','GC'],[0,0,1]

	# motifs,base_pos_in_motif = ['A'],[0]

	bam = pysam.AlignmentFile(bam_file,'rb')
	header = bam.header.to_dict()
	
	output_bam = pysam.AlignmentFile(output_bam_name,"wb", header = header)
	
	parseBam(bam, motifs, base_pos_in_motif, output_bam)





if __name__=="__main__":
	main()


