import sys
import argparse
import numpy as np
import gzip
from tqdm import tqdm


def main():
	args = parse()
	process(args)

def parse():
	parser = argparse.ArgumentParser(description='FT --all extract output to bed files (for DiMeLo, by DD)')

# Data and Run Directories
	parser.add_argument('--tsv', dest='tsv',type=str,help='tsv file from ft extract --all')
	parser.add_argument('--clean',dest='clean', action='store_true',help='if set, trims overlapping bed12 ends for compatibility with bigBed conversion')

	args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
	
	return args


def str2Arr(a):
	a = np.array(a.split(',')[:-1])
	if len(a) < 1:
		return None
	else:
		return a.astype(np.uint32)

def processLine(line,idxs):
	
	split_line = line.strip().split()
	# -1 mC_qual -2 ref_mC_pos -3 mA_qual -4 ref_mA_pos 

	fields = [split_line[i] for i in idxs]

	fields[-4] = str2Arr(fields[-4])
	fields[-3] = str2Arr(fields[-3])
	fields[-2] = str2Arr(fields[-2])
	fields[-1] = str2Arr(fields[-1])

	return fields 


def generateBed12(fields,block_start_index, block_qual_index,cutoff = 225,label_add = None ):

	entry = [fields[0], fields[1], fields[2], fields[3] + label_add, '1000', fields[4],fields[1],fields[2],'128,0,128']

	starts = fields[block_start_index]
	quals = fields[block_qual_index]

	mods = fields[block_start_index]
	
	mod_quals = fields[block_qual_index]

	if type(starts) == type(None):
		blockCounts = 2
		blockStarts = '0,' + str((int(fields[2]) - int(fields[1])) - 1)
		blockSizes = '1,1'
		blockCounts = str(blockCounts)
	elif np.sum(mod_quals >= cutoff) == 0:
		blockCounts = 2
		blockStarts = '0,' + str((int(fields[2]) - int(fields[1])) - 1)
		blockSizes = '1,1'
		blockCounts = str(blockCounts)
	else:	
		idx = mod_quals >= cutoff
		blockCounts = 2 + int(np.sum(idx))
		blockStarts = '0,' + ','.join(list(((mods[idx] - int(fields[1]))).astype(str))) + ',' + str((int(fields[2]) - int(fields[1])) - 1)
		blockSizes = '1,' + ','.join(['1' for i in range(blockCounts - 1) ])
		blockCounts = str(blockCounts)

	entry.append(blockCounts)
	entry.append(blockSizes)
	entry.append(blockStarts)

	entry = '\t'.join(entry)

	return entry


	# chrom, start, stop,name, score, strand, start, stop, itemrgb, block coutn, block size , block start

def process(args):
	
	idxs = np.array([0,1,2,3,5,28,29,31,32]).astype(int)
	



	if 'gz' in args.tsv:
		with gzip.open(args.tsv,'rt') as handle:
			for line in tqdm(handle,miniters=100):
				if "#" in line:
					continue
				
				fields = processLine(line,idxs)

				# first generate mA bed
				# then generate mC bed
				# create mA bed
				# -1 mC_qual -2 ref_mC_pos -3 mA_qual -4 ref_mA_pos 

				# process m6A 

				m6a_entry = generateBed12(fields, block_start_index = -4 ,block_qual_index = -3, cutoff=225,label_add='_m6A')

				# process mCpG 

				mcpg_entry = generateBed12(fields, block_start_index = -2 ,block_qual_index = -1, cutoff=200,label_add='_mCpG')

				print(m6a_entry)
				print(mcpg_entry)

	# else:
	# 	with open(args.tsv,'r') as handle:
	# 		for line in handle:
	# 			if "#" in line:
	# 				continue
				
	# 			# do

if __name__=="__main__":
	main()
