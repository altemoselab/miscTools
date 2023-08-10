# miscTools
Miscellaneous scripts for working with DiMeLo-seq + other sequencing data types.

### addBarcodeToRG.py 
This is a script that will take the output of the guppy barcoder and add the barcode to the read group in a bam file.
Usage: 
```addBarcodeToRG.py <bam> <guppy_barcoder_summary_tsv> 
        <comma_seperated_list_of_barcodes_to_keep, including leading 0s (example 01,02,03,08,11,12) > ```

### refModMatch.py 
This is a script which will only keep bases in the MM tag if they are a match to the reference genome.
Acceps stdin or bam file name as input.
Usage:
```refModMatch.py < input.bam > < output.bam > ```
or (with stdin)
``` refModMatch.py - <output.bam> ```
This should operate on about 1.5k reads per sec/thread. You can increase throughput by subsampling bams with samtools
and providing each a thread.
Here's an example:
``` samtools view -b input.bam chr6 | refModMatch.py - output.chr6.bam ```


