

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
		help = 'motif to parse, current options A,(C)G,G(C),(C)C \n Provide multiple contexts as comma seperated list, example: A,CG,GC,CC')

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

def generateMMVals(read,forward_idx, base):

	# because forward idx is adjusted to the forward coordinates we can use it
	seq = np.array(list(read.get_forward_sequence()))

	binary = np.zeros(shape=len(seq),dtype=int)
	binary[seq==base] = 1
	
	cumsum_forward = np.cumsum(binary)[forward_idx]	
	try:
		new_mm = list([(cumsum_forward[0] - 1).astype(str)]) + list((np.diff(cumsum_forward) - 1).astype(str))

		return new_mm

	except:

		return None
	# in the sequence will then convert 

def createTag(mm,ml,tag):
	outMM = tag + ',' +','.join(mm) + ';'
	outML = ml 	
	return outMM, outML

def writeRead(file,read,raw_mm,raw_ml):


	if len(raw_ml) > 0:

		final_mm = ''.join(raw_mm)
		final_ml = array('B', np.concatenate(raw_ml))

	else:

		final_mm = None
		final_ml = None

	read.set_tag('MM', final_mm, "Z")
	read.set_tag('ML' ,final_ml)

	file.write(read)


def parseBam(bam, motifs, output_bam):


	# reverse parsing is wrong

	parseA = False
	parseCG = False
	parseGC = False
	parseCC = False


	if 'A' in motifs:
		parseA = True
	if 'CG' in motifs:
		parseCG = True
	if 'GC' in motifs:
		parseGC = True
	if 'CC' in motifs:
		parseCC = True



	# counter = 0
	for read in tqdm.tqdm(bam):
		# counter += 1
		# if counter > 2500:
		# 	return
		# if not read.is_reverse:
		# 	continue

		aligned_pairs = np.array(read.get_aligned_pairs(with_seq=True,matches_only=False))
		# get the sequence to find mismatches, and we want non-matches (indels) as well 
		none_idx = aligned_pairs[:,0] != None
		aligned_pairs = aligned_pairs[aligned_pairs[:,0] != None ,:]
		none_idx = aligned_pairs[:,2] == None 
		# index of Nones (indel)
		aligned_pairs[:,2][none_idx] = "n" 

		seq = read.get_forward_sequence()
		seq_len = len(seq)


		
		mod_bases = read.modified_bases_forward
		# modified bases in the original forward coordinates

		reverse_state = read.is_reverse


		if reverse_state:
			aligned_pairs[:,0] = (seq_len-1) - aligned_pairs[:,0].astype(int)
			aligned_pairs = aligned_pairs[np.argsort(aligned_pairs[:,0]),:] #### LINE THAT FIXED PROBLEM OF VERY FEW T


		all_new_mm = []
		all_new_ml = []

		for m in mod_bases:
			
			if 'h' in m: # not supporting hemi-methylation currently
				continue
			
			mod_stack = np.vstack(mod_bases[m])

			ap_idx_of_all_mods = np.isin(aligned_pairs[:,0].astype(int), mod_stack[:,0].astype(int),assume_unique=True,kind='table')

			# handling reverse strand is incorrect

			if m[0] == "A" and parseA:
				base = "A"
				new_label = 'A+a'

				if read.is_reverse:
					base = "T"		

					# print(seq[mod_stack[0]])		


				forward_idx = aligned_pairs[:,0][ap_idx_of_all_mods][ aligned_pairs[:,2][ap_idx_of_all_mods] == base ].astype(int)
				
				forward_idx.sort()
				
				tag_idx = np.isin( mod_stack[:,0].astype(int),forward_idx,assume_unique=True,kind='table')

				# everything is correct up to this point 
				# need to change how MM Tag is generated i think
				new_ml_tags = mod_stack[:,1][tag_idx]
				

				new_mm = generateMMVals(read, forward_idx,"A")

				if type(new_mm) != type(None):

					output_MM, output_ML = createTag(new_mm,new_ml_tags,'A+a') 
		
					all_new_mm.append(output_MM)
					all_new_ml.append(output_ML)


			if m[0] == "C" and m[2]=="m":
				
				if parseCG or parseGC:

					base = "C"
					neighbor_base = "G"
					
					if read.is_reverse:
					
						base = "G"
						neighbor_base = "C"
					


					base_match_idx = aligned_pairs[:,2][ap_idx_of_all_mods].T == base
					# all instances of where bases match at the position


					# the way to do this is to actually calculate the forward position -+ 1 from base_match_idx and then check which

					# so take base_match_idx 
					# keep it in the same space as base match idx



					aligned_pairs_scalar_index = np.argwhere(ap_idx_of_all_mods).T[0][base_match_idx].astype(int)

					if len(aligned_pairs_scalar_index) < 1:
						continue

					if aligned_pairs_scalar_index[0] == 0:
						aligned_pairs_scalar_index = aligned_pairs_scalar_index[1:]
					
					try:
						if aligned_pairs_scalar_index[-1] >= len(aligned_pairs[:,0]) - 1:
							aligned_pairs_scalar_index = aligned_pairs_scalar_index[:-1]
					except:
						continue
					# need to clip 1 base of the edge of base_match_idx before combining					
					
					if parseCG :

						ap_idx_of_all_mods_plus_one = aligned_pairs_scalar_index + 1
						neighbor_idx = aligned_pairs[:,2][ap_idx_of_all_mods_plus_one] == neighbor_base

						
						ap_CG_idx = aligned_pairs_scalar_index[neighbor_idx]

					if parseGC : 

						ap_idx_of_all_mods_minus_one = aligned_pairs_scalar_index - 1
						neighbor_idx = aligned_pairs[:,2][ap_idx_of_all_mods_minus_one] == neighbor_base

						ap_GC_idx = aligned_pairs_scalar_index[neighbor_idx]
						


					if parseCG and parseGC:

						forward_idx = aligned_pairs[:,0][np.unique(np.concatenate([ap_CG_idx,ap_GC_idx]))]
					elif parseCG:
						forward_idx = aligned_pairs[:,0][ap_CG_idx]
					elif parseGC:
						forward_idx = aligned_pairs[:,0][ap_GC_idx]
					else:
						print('Modification not supported')
						return None


				if parseCC:
					base = "C"
					neighbor_base = "C"
					
					if read.is_reverse:
					
						base = "G"
						neighbor_base = "G"

					base_match_idx = aligned_pairs[:,2][ap_idx_of_all_mods].T == base
					aligned_pairs_scalar_index = np.argwhere(ap_idx_of_all_mods).T[0][base_match_idx].astype(int)

					if len(aligned_pairs_scalar_index) < 1:
						continue

					if aligned_pairs_scalar_index[0] == 0:
						aligned_pairs_scalar_index = aligned_pairs_scalar_index[1:]
					
					try:
						if aligned_pairs_scalar_index[-1] >= len(aligned_pairs[:,0]) - 1:
							aligned_pairs_scalar_index = aligned_pairs_scalar_index[:-1]
					except:
						continue
					# need to clip 1 base of the edge of base_match_idx before combining					
					

					ap_idx_of_all_mods_plus_one = aligned_pairs_scalar_index + 1
					neighbor_idx = aligned_pairs[:,2][ap_idx_of_all_mods_plus_one] == neighbor_base

					
					ap_CC_idx = aligned_pairs_scalar_index[neighbor_idx]

					#if :
					#	forward_idx = aligned_pairs[:,0][ap_CC_idx]
					try:
						forward_idx = aligned_pairs[:,0][np.unique(np.concatenate([forward_idx,ap_CC_idx]))]
					except:
						forward_idx = aligned_pairs[:,0][ap_CC_idx]
				forward_idx = np.unique(forward_idx.astype(int))
				forward_idx.sort()
				# forward idx is in the correct orientation

				new_label = 'C+m'

				tag_idx = np.isin(mod_stack[:,0].astype(int),forward_idx,assume_unique=True,kind='table')
				
				
				new_ml_tags = mod_stack[:,1][tag_idx]
				

				new_mm = generateMMVals(read, forward_idx,"C")
				
				if type(new_mm) != type(None):

					output_MM, output_ML = createTag(new_mm,new_ml_tags,'C+m') 

					all_new_mm.append(output_MM)
					all_new_ml.append(output_ML)

		writeRead(output_bam,read,all_new_mm,all_new_ml)




def main():

	bam_file, output_bam_name, input_motifs = inputArgs()

	motifs = input_motifs.split(',')

	bam = pysam.AlignmentFile(bam_file,'rb')
	header = bam.header.to_dict()
	
	output_bam = pysam.AlignmentFile(output_bam_name,"wb", header = header)
	
	parseBam(bam, motifs, output_bam)

if __name__=="__main__":
	main()


