import sys


meme_header="MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"

def writeOut(current_meme,file_name):
	with open(file_name+"_meme.txt",'w') as out_handle:
		out_handle.write(current_meme)
	return

with open(sys.argv[1],'r') as handle:
	
	current_meme=meme_header
	first=True
	for line in handle:
		
		if line[:5] == "MOTIF":

			if first:
				file_name = line.strip().split('.')[3]
				continue
			else:
				writeOut(current_meme, file_name)
				file_name = line.strip().split('.')[3]

		current_meme =+ line

