hicConvertFormat -m $1 -o file.ginteractions --inputFormat cool \
--outputFormat ginteractions

awk -F "\t" '{print 0, $1, $2, 0, 0, $4, $5, 1, $7}' file.ginteractions.tsv | \
sort -k2,2d -k6,6d > file.ginteractions.tsv.short.sorted

alias juicer='java -Xms512m -Xmx2048m -jar /oak/stanford/groups/altemose/dubocd/juicer/scripts/juicer_tools.2.20.00.jar'

juicer pre -r 2000,4000,8000,16000,32000,64000,128000 file.ginteractions.tsv.short.sorted file.ginteractions.tsv.short.sorted.hic $2
