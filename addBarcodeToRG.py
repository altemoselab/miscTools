#guppy output and guppy barcoder add barcode to bam file



import sys
import pysam 
from tqdm import tqdm


inbam = pysam.AlignmentFile(sys.argv[1],'rb')
barcode_file = sys.argv[2]
barcodes_to_keep = sys.argv[3].split(',') # comma seperated list of barcodes including leading 0s
outbam = sys.argv[1][:-3]+'barcoded.bam'

#RG_sub_dict = {}

dict_header = inbam.header.to_dict()
dict_header['RG']= [] # list of dictionaries for every barcode there is ID and SM
for barcode in barcodes_to_keep:
	new_dict = {'ID': barcode, 'SM': 'barcode'+barcode }
	dict_header['RG'].append(new_dict)

outbam = pysam.AlignmentFile(outbam, 'wb', header = dict_header)


readname_to_bc_dict = {}

with open(barcode_file,'r') as handle:
	handle.readline()
	for line in handle:

		split_line = line.split("\t")
		read_id, barcode = split_line[0],split_line[1]

		barcode = barcode.replace('barcode','')
		if barcode in barcodes_to_keep:
			readname_to_bc_dict[read_id] = barcode


	
for read in tqdm(inbam):
#	print(str(read.query_name))
	barcode_to_add = readname_to_bc_dict.get(str(read.query_name))
	if barcode_to_add:
		read.set_tag(tag="RG",value=barcode_to_add,value_type="Z")

		outbam.write(read)






















