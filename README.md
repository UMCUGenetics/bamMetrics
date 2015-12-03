# bamMetrics
Tool to generate bam statistics and pdf/html reports. To be used on the UMC hpc cluster.

## Usage
`perl bamMetrics.pl [options] -bam <bamfile1.bam> -bam <bamfile2.bam>`

### Required:
- bam

### Optional (some have default settings):

#### Whole genome sequencing statistics
- wgs
- coverage_cap 250

#### RNA sequencing statistics
- rna
- ref_flat refFlat.gz
- strand [NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND]
- ribosomal_intervals rRNA.intervallist

#### Capture sequencing statistics (exome)
- capture
- targets targets.bed
- baits baits.bed

#### Other
- output_dir ./bamMetrics
- run_name bamMetrics
- genome genome.fasta
- queue all.q
- queue_time    2:0:0
- queue_threads 1
- queue_mem 8
- picard_path /path/to/picard/

## Dependencies
- Perl
- R
- sge
- Picard >= 1.119

### Perl Modules
- strict
- Getopt::Long
- Pod::Usage
- POSIX
- Cwd
- File::Basename

### R packages
- ggplot2
- knitr
- markdown
- reshape
- xtable
- tools
- brew
