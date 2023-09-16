# split bed12-mols into input window size 



import numpy as np 
from pybedtools import BedTool
import sys

bed = BedTool(sys.argv[1])
window = int(sys.argv[2])


for interval in bed:

	interval_start, interval_stop = int(interval.start), int(interval.stop)

	feature_starts = np.array(interval[-1].split(',')[1:-1]).astype(int)

	feature_starts_relative_to_genome = feature_starts + interval_start

	windowed_interval_start = interval_start 
	


	while windowed_interval_start < interval_stop:

		windowed_interval_stop = windowed_interval_start + window
		new_interval_feature_starts = feature_starts_relative_to_genome[
																		(feature_starts_relative_to_genome >= windowed_interval_start) & 
																		(feature_starts_relative_to_genome < windowed_interval_stop)
																		] -  windowed_interval_start

		num_new_feature = len(new_interval_feature_starts)
		size_str_arr = ','.join(['1'] * num_new_feature)
		new_interval_feature_starts_str_arr= ','.join(list(new_interval_feature_starts.astype(str)))





		new_interval = [str(interval[0]),
						str(windowed_interval_start),
						str(windowed_interval_stop),
						str(interval[3]),
						str(interval[4]),
						str(interval[5]),
						str(windowed_interval_start),
						str(windowed_interval_stop),
						str(interval[8]),
						str(num_new_feature),
						str(size_str_arr),
						str(new_interval_feature_starts_str_arr)] 

		print('\t'.join(new_interval))



		windowed_interval_start = windowed_interval_stop


