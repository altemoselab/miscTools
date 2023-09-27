
import sys
import numpy as np 
import itertools

fasta = sys.argv[1] 



# def allPermutations(alphabet,n):

# 	for i in n:






alphabet = "ACGT"

with open(fasta,'r') as file:
	for line in file:
		if ">" in line:
			read_name = line.strip()
		else:
			size = 0 
			seq = np.array(list(line.strip()))
			N_idx = seq == "N"
			# x = itertools.product("A","C","G","T",repeat = 8)
			for c in itertools.product(alphabet,repeat = np.sum(N_idx)):
				
				seq[N_idx] = np.array(c)
				print(read_name+"_"+''.join(c))
				print(''.join(list(seq)))