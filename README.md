###bamStats
Tool to generate bam statistics and pdf/html reports. To be used on the UMC hpc cluster.

####Usage
`perl bamStats.pl [options] -bam <bamfile1.bam> -bam <bamfile2.bam>`

#####Required:
- bam

#####Optional (some have default settings):
######Whole genome sequencing statistics
- wgs
- coverage_cap 250
    
######RNA sequencing statistics
- rna
- ref_flat refFlat.gz
- strand [NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND]
    
######Capture sequencing statistics (exome)
- capture
- targets targets.bed
- baits baits.bed
    
######Other
- output_dir = ./bamStats
- genome genome.fasta
- queue veryshort
- queue_threads 1
- queue_mem 8
- picard_path /path/to/picard/


