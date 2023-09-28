
import pysam 
import sys
import numpy as np 
import argparse
from tqdm import tqdm 
from pybedtools import BedTool


def inputArgs():
	'''Parse in arguments. '''
	
	parser = argparse.ArgumentParser()

	parser.add_argument('--bam', 
		type = str, 
		default = '-',
		help = 'bam file, with index')

	parser.add_argument('--bed',
		type=str,
		help = 'bed intervals to calculate refMatched bedGraph, make sure elements do not overlap')

	parser.add_argument('-m','--min_ml_score',
		type=int,
		default=125,
		help='min-ml-score for modification')

	args = parser.parse_args()

	return args.bam, args.bed, args.min_ml_score



def outputCenter(bam,bed,min_ml_score):

	output_dict = {}#chrom and then in each chrom it is 
	for interval in tqdm(bed,mininterval=0.5):

		output_dict[str(interval[0])] = {}

		# for pos in range(int(interval[1]),int(interval[2])+1):
		# 	output_dict[str(interval[0])][pos] = [0,0] #overlap,modification 

		for read in bam.fetch(str(interval[0]),int(interval[1]),int(interval[2])):

			ap = np.array(read.get_aligned_pairs(matches_only=True,with_seq=True)).T.astype(str)
			# if read.is_reverse:
			# 	for_seq = np.frombuffer(bytes(str(Seq(read.get_forward_sequence()).reverse_complement()), "utf-8"), dtype="S1")
			# else:
			
			if read.is_reverse:
				ap_ref_match = ap[:,ap[2]=="T"]
			else:
				ap_ref_match = ap[:,ap[2]=="A"]

			ref_align_matches = ap_ref_match[1].astype(int)
			forward_align_matches = ap_ref_match[0].astype(int)
			
			within_window_ref_align_matches = ref_align_matches[(ref_align_matches >= int(interval[1])) & (int(interval[2]) > ref_align_matches)]
			within_window_ref_align_matches_list = list(within_window_ref_align_matches)


			# only forward positions within genomic window
			forward_align_matches = forward_align_matches[(ref_align_matches >= int(interval[1])) & (int(interval[2]) > ref_align_matches)] 


			for pos in within_window_ref_align_matches_list:
				if pos not in output_dict[str(interval[0])]:
					output_dict[str(interval[0])][pos] = [0,0]
				output_dict[str(interval[0])][pos][0] += 1

			mods = np.array(read.modified_bases_forward[('A',0,'a')]).T.astype(int)
			valid_mods = mods[0][mods[1] >= min_ml_score]

			# Not here is where we would have to convert
			if read.is_reverse:
				valid_mods = np.abs(valid_mods - read.infer_query_length()) + 1

			mod_counter_ref_matches = list(within_window_ref_align_matches[np.isin(forward_align_matches,valid_mods)])
			# print(mod_counter_ref_matches) 
			for pos in mod_counter_ref_matches:
				# print(str(interval[0]),pos)
				output_dict[str(interval[0])][pos][1] += 1



	for c in output_dict:
		for pos in output_dict[c]:
			mod_count = output_dict[c][pos][1]
			base_count = output_dict[c][pos][0]
			if base_count == 0:
				print(c,"\t",pos,"\t",pos+1,"\t",0,"\t",0,"\t",0)
			else:
				print(c,"\t",pos,"\t",pos+1,"\t",output_dict[c][pos][1]/output_dict[c][pos][0],"\t",output_dict[c][pos][1],"\t",output_dict[c][pos][0])









def main():
	
	bam_file, bed_file, min_ml_score = inputArgs()

	bed = BedTool(bed_file)
	bam = pysam.AlignmentFile(bam_file,'rb',check_sq=False) 



	outputCenter(bam,bed,min_ml_score)




if __name__=="__main__":
	main()




