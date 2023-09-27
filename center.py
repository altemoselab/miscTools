# Centering tool for single molecule reads overlapping bed file 

# Input modBam, bed with genomic coordinates
# By Default it searches for match to reference, can be turned off with -n (speed up)
# Outputs: Centered molecule bed12 (-), centered to beginning of interval
# 		option: genome_tabular center output Default and standard)
				# chr		start	stop 	strand	centering_pos	A_acount mA_count 

#output will be printed to stdout
#output is not adjusted to orientat

import pysam 
import sys
import numpy as np 
import gzip
import subprocess
import shlex
import argparse
import tqdm
from pybedtools import BedTool
import regex as re
import pyfaidx
from Bio.Seq import Seq


def inputArgs():
	'''Parse in arguments. '''
	
	parser = argparse.ArgumentParser()

	parser.add_argument('--bam', 
		type = str, 
		default = '-',
		help = 'bam file, with index')

	parser.add_argument('--bed',
		type=str,
		help = 'output bam name')

	parser.add_argument('-b','--bases',
		type=str,
		default='A',
		help = 'modification base, options = [ A,CG ]. \n Can be multiple but must be seperated by comma: A,CG [NOT SUPPORTED YET] ')

	parser.add_argument('-m','--min_ml_score',
		type=int,
		default=125,
		help='min-ml-score for modification')
	
	parser.add_argument('-g','--genome',
		type=str,
		help='reference genome, with index, not gzipped')

	parser.add_argument('-d','--distance',
		type=int,
		default=10000,
		help='distance away from center to look at mods, sequences are dervied from read alignment')


	parser.add_argument('--no_match', action='store_true')
	args = parser.parse_args()

	return args.bam, args.bed, args.bases, args.min_ml_score, args.no_match, args.distance, args.genome



def outputCenter(bam,bed,ref,bases,min_ml_score,distance):



	if ',' in bases:
		split_bases = bases.split(',')
		base='list'
	else:
		if bases == "A":
			mod_code = ('A', 0, 'a')
			base = bases
			revbase='T'
			bytebase = bytes(base, encoding='utf-8')
			byterevbase = bytes(revbase, encoding='utf-8')
		elif bases == 'CG':
			mod_code = ('C', 0, 'm')
			base = 'CG'
			bytebaseC = bytes('C', encoding='utf-8')
			bytebaseG = bytes('G', encoding='utf-8')
		else:
			sys.stderr.write('invalid base {}, options are A, CpG, or A,CpG'.format(bases))
			exit()
	


	for interval in tqdm.tqdm(bed,mininterval=0.5):
		if len(interval) >= 6:
			if "+" == str(interval[5]) or "-" == str(interval[5]):
				strand = interval[-1]
			elif "+" == str(interval[3]) or "-" == str(interval[3]):
				strand = interval[3]
			else:
				strand = "."
		elif len(interval) == 4:

			if "+" == str(interval[3]) or "-" == str(interval[3]):
				strand = interval[3]
			else:
				strand = "."
		else:
			strand = "."

		ref_start = int(interval[1]) - distance
		ref_end = int(interval[1]) + distance
		read_length = ref_end - ref_start


		seq  = ref[str(interval[0])][int(interval[1]) - distance:int(interval[1]) + distance]


		base_idx = [m.start(0) for m in re.finditer(base, seq)]
		if base == "A":
			rev_comp_base_idx = [m.start(0) for m in re.finditer(revbase, seq)]
		elif base == "CG":
			rev_comp_base_idx = [m+1 for m in base_idx] # count Gs
		
		sorted_pos = np.sort(np.array(base_idx + rev_comp_base_idx)) 

		base_count = np.zeros(shape=len(sorted_pos))
		modification_count = np.zeros(shape=len(sorted_pos))
		# the reverse complement also
		 # has potential mod sites so we must search for them on the rev_seq sequence 

		# base_idx will be the positions in which we can have a modification that we are tracking
		# we will have a counter based on an index linked to base_idx



		for read in bam.fetch(str(interval[0]),ref_start,ref_end):

			ap = np.array(read.get_aligned_pairs(matches_only=True,with_seq=False)).T.astype(int)
			# if read.is_reverse:
			# 	for_seq = np.frombuffer(bytes(str(Seq(read.get_forward_sequence()).reverse_complement()), "utf-8"), dtype="S1")
			# else:
			for_seq = np.frombuffer(bytes(read.get_forward_sequence(), "utf-8"), dtype="S1")

			# reverse read can contribute at T 
			# forwrad read can contribute at A

			if base == 'A':
				# if read.is_reverse:
				# 	for_seq_match_idx = np.argwhere(for_seq==bytebase).T[0]
				# else:
				for_seq_match_idx = np.argwhere(for_seq==bytebase).T[0]
			
			elif base == 'CG':
			# if base is CG , C is what we count for sense strand
			# G is option for antisense strand
					C_match = np.argwhere(for_seq==bytebaseC).T[0]
					if C_match[-1] == len(for_seq)-1:
						C_match = C_match[:-1]
					CG_idx = np.argwhere(for_seq[C_match+1]==bytebaseG).T[0]
					G_match = C_match[CG_idx] + 1
					C_match = C_match[CG_idx]
					if read.is_reverse:
						for_seq_match_idx = G_match
					else:
						for_seq_match_idx = C_match
					# print(C_match,G_match)

			# NEED TO ADJUST COORDS FOR REVERSE MAPPING READS BEFORE ACCESSING ALIGNED PAIRS

			correct_coord = lambda x: np.abs(x - read.infer_query_length()) + 1

			match_ref_coord = ap[1][np.isin(ap[0],correct_coord(for_seq_match_idx))]
			match_ref_coord = match_ref_coord[(match_ref_coord >= ref_start) & (match_ref_coord < ref_end)] - ref_start
			# print(match_ref_coord)

			#### UP TOTAL BASE COUNT ###
			# print(match_ref_coord)
			base_count[np.isin(sorted_pos,match_ref_coord)] += 1
			###

			# need map positions to reference that will count 
			# fastest way to do this is to identify all positions in read that are A/T or CG
			# then use aligned pairs to pull out ref coords that match 
			# adjust by ( ref_coord >= ref_start and ref_coord < ref_stop- ref_start ) - ref_start
			# and this gives you relative postision that pass and then you can populate array 

			# we have the reference positions in seq 
			# we have the alignment information in ap
			# we can now use these two to sieve modifications that much both


			### modified bases
			mod_bases_forward = np.array(read.modified_bases_forward[mod_code]).T
 		
			accepted_mods = mod_bases_forward[0][mod_bases_forward[1] >= min_ml_score]
			# if read.is_reverse:
			# print(for_seq[accepted_mods])
				# 	# coords are 0 based, length is inherently 1-based, so we must add 1
			# 	accepted_mods  = np.abs(accepted_mods - read.infer_query_length()) + 1
				

			mod_match_ref_coord = ap[1][np.isin(ap[0],correct_coord(accepted_mods[np.isin(accepted_mods,for_seq_match_idx)]))]
			mod_match_ref_coord = mod_match_ref_coord[(mod_match_ref_coord >= ref_start) & (mod_match_ref_coord < ref_end)] - ref_start

			### UP MODIFICATION COUNT
			modification_count[np.isin(sorted_pos,mod_match_ref_coord)] += 1
			###

		# print(modification_count)
		# print(base_count)
		sorted_pos_idx = range(len(sorted_pos))
		sorted_pos = (sorted_pos + ref_start) - distance
		for i in sorted_pos_idx:
			out_interval = [interval[0],
							str(sorted_pos[i]),
							str(sorted_pos[i] + 1),
							strand,
							str(interval[1]),
							str(int(base_count[i])),
							str(int(modification_count[i]))]
			print('\t'.join(out_interval))







def main():
	
	bam_file, bed_file, bases, min_ml_score,no_match,distance,genome = inputArgs()

	bed = BedTool(bed_file)
	bam = pysam.AlignmentFile(bam_file,'rb',check_sq=False) 
	ref = pyfaidx.Fasta(genome, sequence_always_upper=True, as_raw=True)


	outputCenter(bam,bed,ref,bases,min_ml_score,distance)




if __name__=="__main__":
	main()



		# if len(interval) >= 6:
		# 	if "+" == str(interval[5]) or "-" == str(interval[5]):
		# 		strand = interval[-1]
		# 	elif "+" == str(interval[3]) or "-" == str(interval[3]):
		# 		strand = interval[3]
		# 	else:
		# 		strand = "."
		# elif len(interval) == 4:

		# 	if "+" == str(interval[3]) or "-" == str(interval[3]):
		# 		strand = interval[3]
		# 	else:
		# 		strand = "."
		# else:
		# 	strand = "."
