#guppy output and guppy barcoder add barcode to bam file



import sys
import pysam 



inbam = pysam.AlignmentFile(sys.argv[1],'rb')
barcode_file = sys.argv[2]
barcodes_to_keep = sys.argv[3].split(',') # comma seperated list of barcodes including leading 0s


readname_to_bc_dict = {}


with open(barcode_file,'r') as handle:
	for line in handle:

		split_line = line.split("\t")
		read_id, barcode = split_line[0],split_line[1]

		barcode.replace('barcode','')
	
		if barcode in barcodes_to_keep:
			readname_to_bc_dict[read_id] = barcode


for read in inbam:

	barcode_to_add = readname_to_bc_dict[str(read.query_name)]

	read.set_tag(tag="RG",value=barcode_to_add,value_type="Z")

	print(read)






















