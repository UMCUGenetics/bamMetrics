#!usr/bin/perl
### Robert Ernst
### bamStats
### Tool for generating html and/or pdf reports with bam statistics

use strict;
use warnings;
use Getopt::Long;

### Input options ###

# Declare options
my @bams;

my $wgs = "";
my $rna = "";
my $exome = "";
my $targets;
my $baits;

my $pdf = "";
my $html = "";

# Parse options
GetOptions ("bam=s" => \@bams,
	    "wgs" => \$wgs,
	    "rna" => \$rna,
	    "exome" => \$exome,
	    "targets=s" => \$targets,
	    "baits=s" => \$baits,
	    "pdf" => \$pdf,
	    "html" => \$html
	    ) or die usage();

# Check user input
if( !@bams || !($wgs || $rna || ($exome && $targets && $baits)) || !($pdf || $html)) { usage() };

### Run picard tools ###
foreach my $bam (@bams) {
    my $bamName = (split("/",$bam))[-1];
    $bamName =~ s/.bam//;
    print $bam ."\t".$bamName."\n";
}

### Functions ###
sub usage{
    warn "Usage: perl bamStats.pl -bam <bamfile> -bam <bamfile2> [-wgs | -rna | (-exome -targets <target_file.txt>  -baits <baits_file.txt>)] -html -pdf";
    exit;
}