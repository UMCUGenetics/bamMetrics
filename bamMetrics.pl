#!usr/bin/perl
### Robert Ernst
### bamMetrics
### Tool for generating html and/or pdf reports with bam statistics

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(tmpnam);
use Cwd qw(cwd abs_path);
use File::Basename qw( dirname );

### Input options ###

# Declare options
my @bams;
my $single_end = "";

my $wgs = "";
my $coverage_cap = 250;

my $rna = "";
my $ref_flat = "/hpc/cog_bioinf/data/annelies/RNA_Seq/hg19.refFlat.gz";
my $strand = "SECOND_READ_TRANSCRIPTION_STRAND";
my $rib_interval = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.rRNA.intervallist";

my $capture = "";
my $targets = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/sorted_Homo_sapiens.GRCh37.74_nopseudo_noRNA_CDS_picard.bed";
my $baits = "/hpc/cog_bioinf/ENRICH/PICARD/sorted_SS_exome_v5_S04380110_Covered_picard.bed";

#my $pdf = "";
#my $html = "";
my $output_dir = cwd()."/bamMetrics";
my $run_name = "bamMetrics";

# Picard and Cluster settings
my $queue = "veryshort";
my $queue_threads = 1;
my $queue_mem = 8;
my $picard_path = "/hpc/cog_bioinf/common_scripts/picard-tools-1.119";
my $genome = "/hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta";

# Development settings
my $debug = "";

# Parse options
GetOptions ("bam=s" => \@bams,
	    "single_end" => \$single_end,
	    "wgs" => \$wgs,
	    "coverage_cap=i" => \$coverage_cap,
	    "rna" => \$rna,
	    "ref_flat=s" => \$ref_flat,
	    "strand=s" => \$strand,
	    "ribosomal_intervals=s" => \$rib_interval,
	    "capture" => \$capture,
	    "targets=s" => \$targets,
	    "baits=s" => \$baits,
#	    "pdf" => \$pdf,
#	    "html" => \$html,
	    "output_dir=s" => \$output_dir,
	    "run_name=s", => \$run_name,
	    "queue=s" => \$queue,
	    "queue_threads=i" =>  \$queue_threads,
	    "queue_mem=i" => \$queue_mem,
	    "picard_path=s" => \$picard_path,
	    "genome=s" => \$genome,
	    "debug" => \$debug
	    ) or pod2usage(1);

# Check user input
if( !@bams ) { pod2usage(1) };
	### check file existence

### Create output dirs
if(! -e $output_dir){
    mkdir($output_dir) or die "Could not create directory: $output_dir";
}
my $tmp_dir = $output_dir."/tmp";
if(! -e $tmp_dir){
    mkdir($tmp_dir) or die "Could not create directory: $tmp_dir";
}

### Run picard tools ###
my @picardJobs;
my @wgsmetrics;
my @hsmetrics;
my @rnametrics;
my @bam_names;
my $javaMem = $queue_threads * $queue_mem;
my $picard = " java -Xmx".$javaMem."G -jar ".$picard_path;

foreach my $bam (@bams) {
    #Parse bam file name
    $bam = abs_path($bam);
    my $bam_name = (split("/",$bam))[-1];
    $bam_name =~ s/.bam//;
    push(@bam_names, $bam_name);
    print "\n$bam_name \t $bam \n";

    #Create picard output folder per bam file
    my $bam_dir = $output_dir."/".$bam_name;
    if(! -e $bam_dir){
	mkdir($bam_dir) or die "Could not create directory: $bam_dir";
    }

    # Multiple metrics
    my $output = $bam_dir."/".$bam_name."_MultipleMetrics.txt";
    if( $single_end ){
	if(! (-e $output.".alignment_summary_metrics" && -e $output.".base_distribution_by_cycle_metrics" && -e $output.".quality_by_cycle_metrics" && -e $output.".quality_distribution_metrics") ) {
	    my $command = $picard."/CollectMultipleMetrics.jar R=".$genome." ASSUME_SORTED=TRUE INPUT=".$bam." OUTPUT=".$output." PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=QualityScoreDistribution PROGRAM=QualityScoreDistribution";
	    my $jobID = bashAndSubmit(
		command => $command,
		jobName => "MultipleMetrics_".$bam_name."_".get_job_id(),
		tmpDir => $tmp_dir,
		outputDir => $bam_dir,
		queue => $queue,
		queueThreads => $queue_threads,
		);
	    push(@picardJobs, $jobID);
	}
    } else { #paired
	if(! (-e $output.".alignment_summary_metrics" && -e $output.".base_distribution_by_cycle_metrics" && -e $output.".insert_size_metrics" && -e $output.".quality_by_cycle_metrics" && -e $output.".quality_distribution_metrics") ) {
	    my $command = $picard."/CollectMultipleMetrics.jar R=".$genome." ASSUME_SORTED=TRUE INPUT=".$bam." OUTPUT=".$output." PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=QualityScoreDistribution";
	    my $jobID = bashAndSubmit(
		command => $command,
		jobName => "MultipleMetrics_".$bam_name."_".get_job_id(),
		tmpDir => $tmp_dir,
		outputDir => $bam_dir,
		queue => $queue,
		queueThreads => $queue_threads,
		);
	    push(@picardJobs, $jobID);
	}
    }
    # Library Complexity
    $output = $bam_dir."/".$bam_name."_LibComplexity.txt";
    if(! -e $output) {
	my $command = $picard."/EstimateLibraryComplexity.jar INPUT=".$bam." OUTPUT=".$output;
	my $jobID = bashAndSubmit(
	    command => $command,
	    jobName => "LibComplexity_".$bam_name."_".get_job_id(),
	    tmpDir => $tmp_dir,
	    outputDir => $bam_dir,
	    queue => $queue,
	    queueThreads => $queue_threads,
	    );
	push(@picardJobs, $jobID);
    }
    # WGS
    if($wgs){
	my $output = $bam_dir."/".$bam_name."_WGSMetrics.txt";
	push(@wgsmetrics, $output);
	if(! -e $output) {
	    my $command = $picard."/CollectWgsMetrics.jar R=".$genome." INPUT=".$bam." OUTPUT=".$output." MINIMUM_MAPPING_QUALITY=1 COVERAGE_CAP=".$coverage_cap;
	    my $jobID = bashAndSubmit(
		command => $command,
		jobName => "WGSMetrics_".$bam_name."_".get_job_id(),
		tmpDir => $tmp_dir,
		outputDir => $bam_dir,
		queue => $queue,
		queueThreads => $queue_threads,
		);
	    push(@picardJobs, $jobID);
	}
	## Add sambamba stats?
    }

    # RNA
    if($rna){
	my $output = $bam_dir."/".$bam_name."_RNAMetrics.txt";
	push(@rnametrics, $output);
	if(! -e $output) {
	    my $command = $picard."/CollectRnaSeqMetrics.jar R=".$genome." REF_FLAT=".$ref_flat." ASSUME_SORTED=TRUE INPUT=".$bam." OUTPUT=".$output." STRAND_SPECIFICITY=".$strand." RIBOSOMAL_INTERVALS=".$rib_interval;
	    my $jobID = bashAndSubmit(
		command => $command,
		jobName => "RNAMetrics_".$bam_name."_".get_job_id(),
		tmpDir => $tmp_dir,
		outputDir => $bam_dir,
		queue => $queue,
		queueThreads => $queue_threads,
		);
	    push(@picardJobs, $jobID);
	}
    }

    # CAPTURE
    if($capture){
	my $output = $bam_dir."/".$bam_name."_HSMetrics.txt";
	push(@hsmetrics, $output);
	if(! -e $output) {
	    my $command = $picard."/CalculateHsMetrics.jar R=".$genome." INPUT=".$bam." OUTPUT=".$output." BAIT_INTERVALS=".$baits." TARGET_INTERVALS=".$targets." METRIC_ACCUMULATION_LEVEL=SAMPLE";
	    my $jobID = bashAndSubmit(
		command => $command,
		jobName => "HSMetrics_".$bam_name."_".get_job_id(),
		tmpDir => $tmp_dir,
		outputDir => $bam_dir,
		queue => $queue,
		queueThreads => $queue_threads,
		);
	    push(@picardJobs, $jobID);
	}
    }
}

### Parse HSMetrics or WGSMetrics
my $root_dir = dirname(abs_path($0));

if( @wgsmetrics ) {
    print "\n Parse WGSMetrics \n";
    my $command = "perl $root_dir/parsePicardOutput.pl -output_dir ".$output_dir." ";
    foreach my $wgsmetric (@wgsmetrics) { $command .= "-wgsmetrics $wgsmetric " }
    my $jobID = bashAndSubmit(
	command => $command,
	jobName => "parse_wgsmetrics_".get_job_id(),
	tmpDir => $tmp_dir,
	outputDir => $output_dir,
	queue => $queue,
	queueThreads => $queue_threads,
	holdJobs => join(",",@picardJobs),
	);
    push(@picardJobs, $jobID);
}
if( @rnametrics ) {
    print "\n Parse RNAMetrics \n";
    my $command = "perl $root_dir/parsePicardOutput.pl -output_dir ".$output_dir." ";
    foreach my $rnametric (@rnametrics) { $command .= "-rnametrics $rnametric " }
    my $jobID = bashAndSubmit(
	command => $command,
	jobName => "parse_rnametrics_".get_job_id(),
	tmpDir => $tmp_dir,
	outputDir => $output_dir,
	queue => $queue,
	queueThreads => $queue_threads,
	holdJobs => join(",",@picardJobs),
	);
    push(@picardJobs, $jobID);
}
if( @hsmetrics ) {
    print "\n Parse HSMetrics \n";
    my $command = "perl $root_dir/parsePicardOutput.pl -output_dir ".$output_dir." ";
    foreach my $hsmetric (@hsmetrics) { $command .= "-hsmetrics $hsmetric " }
    my $jobID = bashAndSubmit(
	command => $command,
	jobName => "parse_hsmetrics_".get_job_id(),
	tmpDir => $tmp_dir,
	outputDir => $output_dir,
	queue => $queue,
	queueThreads => $queue_threads,
	holdJobs => join(",",@picardJobs),
	);
    push(@picardJobs, $jobID);
}

### Run Rplots ###
my $command = "Rscript $root_dir/bamMetrics.R -output_dir $output_dir -root_dir $root_dir -run_name $run_name -samples ".join(" -samples ", @bam_names);
my $jobID = bashAndSubmit(
    command => $command,
    jobName => "bamMetrics_report_".$run_name,
    tmpDir => $tmp_dir,
    outputDir => $output_dir,
    queue => $queue,
    queueThreads => $queue_threads,
    holdJobs => join(",",@picardJobs),
);
push(@picardJobs, $jobID);

### Clean tmp
if(! $debug){
    my $command = "rm -r $tmp_dir";
    my $jobID = bashAndSubmit(
	command => $command,
	jobName => "bamMetrics_clean_".get_job_id(),
	tmpDir => $tmp_dir,
	outputDir => $output_dir,
	queue => $queue,
	queueThreads => $queue_threads,
	holdJobs => join(",",@picardJobs),
	log_output => "/dev/null",
    );
    push(@picardJobs, $jobID);
}

### Functions ###
sub bashAndSubmit {
    my %args = (
	jobName => "bamMetrics",
	holdJobs => "",
	@_);

    my $jobID = $args{jobName};
    my $bashFile = $args{tmpDir}."/".$jobID.".sh";
    my $log_output = $args{tmpDir};
    
    if($args{log_output}) {
	$log_output = $args{log_output};
    }
    
    open BASH, ">$bashFile" or die "cannot open file $bashFile\n";
    print BASH "#!/bin/bash\n\n";
    print BASH "cd $args{outputDir}\n";
    print BASH "$args{command}\n";
    close BASH;
    
    if( $args{holdJobs} ){
	system "qsub -q $args{queue} -pe threaded $args{queueThreads} -o $log_output -e $log_output -N $jobID -hold_jid $args{holdJobs} $bashFile";
    } else {
	system "qsub -q $args{queue} -pe threaded $args{queueThreads} -o $log_output -e $log_output -N $jobID $bashFile";
    }
    return $jobID;
}

sub get_job_id {
    my $id = tmpnam();
    $id =~ s/\/tmp\/file//;
    return $id;
}

__END__

=head1 SYNOPSIS

$ perl bamMetrics.pl [options] -bam <bamfile1.bam> -bam <bamfile2.bam>

    Required:
     -bam

=head1 OPTIONS

    Whole genome sequencing statistics
     -wgs
     -coverage_cap <250>

    RNA sequencing statistics
     -rna
     -ref_flat </hpc/cog_bioinf/data/annelies/RNA_Seq/hg19.refFlat.gz>
     -strand [NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND]
     -ribosomal_intervals </hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.rRNA.intervallist>

    Capture sequencing statistics (exome)
     -capture
     -targets </hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/sorted_Homo_sapiens.GRCh37.74_nopseudo_noRNA_CDS_picard.bed>
     -baits </hpc/cog_bioinf/ENRICH/PICARD/sorted_SS_exome_v5_S04380110_Covered_picard.bed>

    Other:
     -single_end (default is paired end)
     -output_dir <./bamMetrics>
     -run_name <bamMetrics>
     -genome </hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta>
     -queue <veryshort>
     -queue_threads 1
     -queue_mem 8
     -picard_path </hpc/cog_bioinf/common_scripts/picard-tools-1.119>

=cut

