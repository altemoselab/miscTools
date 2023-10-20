


pod5_file=$1 #pod5 file or files 
subset_size=$2 # integer, recommended between 20-50  


pod5 view --ids --no-header $pod5_file > 

split -n l/$subset_size subset_

for file in subset_* ; do 
	mv $file ${file}.txt
	mkdir $file
	mv ${file}.txt ${file}/
done


for dir in subset_* ; do 
	
	pod5 filter -i ${dir}/${dir}.txt -o ${dir}/${dir}.pod5 &

done 

exit