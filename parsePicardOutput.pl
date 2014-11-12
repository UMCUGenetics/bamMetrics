#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

### Input options ###

# Declare options
my @hsmetrics;
my @wgsmetrics;
my @rnametrics;
my $output_dir = "";

# Parse options
GetOptions ("hsmetrics=s" => \@hsmetrics,
	    "wgsmetrics=s" => \@wgsmetrics,
	    "rnametrics=s" => \@rnametrics,
	    "output_dir=s" => \$output_dir
	    ) or pod2usage(1);

if(! (@hsmetrics || @wgsmetrics || @rnametrics) || ! $output_dir) { pod2usage(1) };

### Parse hsmetrics files
if( @hsmetrics ){
    my $filename = $output_dir."/HSMetrics_summary.txt";
    open(SUMFILE, ">", $filename) || die ("Can't create $filename");
    my $printed_header = 0;
    
    foreach my $file (@hsmetrics) {
	print "\t Parsing: ". $file . "\n";
	open(FILE, $file) || die ("Can't open $file");
	my $bait_intervals;
	my $target_intervals;
	
	#Processing headerlines -> beginning with # or empty lines.
	while(<FILE> =~ /(^\s*#)(.*)/ || <FILE> eq "") {
	    my $header_line = $2;
	    if ($header_line =~ m/BAIT_INTERVALS=(\S*).TARGET_INTERVALS=(\S*)/){ #grep bait and target interval file settings.
		$bait_intervals = $1;
		$target_intervals = $2;
	    }
	}
	
	#Processing table
	my $table_header = <FILE>;
	unless ($printed_header) { #print table header once.
	    print SUMFILE "sample \t baitIntervals \t targetIntervals \t". $table_header;
	    $printed_header = 1;
	}
	
	my $line = <FILE>; #grep statistics, here we assume one row with statistics.
	my ($sample) = ($file =~ m/.+\/(.+)_HSMetrics.txt/);
	print SUMFILE $sample ."\t". $bait_intervals ."\t". $target_intervals ."\t". $line;
    }
    close SUMFILE;
}

### Parse wgsmetrics files
if( @wgsmetrics ){
    my $filename = $output_dir."/WGSMetrics_summary.txt";
    open(SUMFILE, ">", $filename) || die ("Can't create $filename");
    my $printed_header = 0;
    
    print "Parsing WGSMetric files \n";
    foreach my $file (@wgsmetrics) {
	print "\t Parsing: ". $file . "\n";
	open(FILE, $file) || die ("Can't open $file");
    
	#Processing headerlines -> beginning with # or empty lines.
	while(<FILE> =~ /(^\s*#)(.*)/ || <FILE> eq "") {}
    
	my $table_header = <FILE>;
	unless ($printed_header) { #print table header once.
	    print SUMFILE "sample \t". $table_header;
	    $printed_header = 1;
	}

	my $line = <FILE>; #grep statistics, here we assume one row with statistics.
	my ($sample) = ($file =~ m/.+\/(.+)_WGSMetrics.txt/);
	print SUMFILE $sample ."\t". $line;
    }
    close SUMFILE
}

### Parse rnametrics files
if( @rnametrics ){
    my $filename = $output_dir."/RNAMetrics_summary.txt";
    open(SUMFILE, ">", $filename) || die ("Can't create $filename");
    my $printed_header = 0;
    
    print "Parsing RNAMetric files \n";
    foreach my $file (@rnametrics) {
	print "\t Parsing: ". $file . "\n";
	open(FILE, $file) || die ("Can't open $file");
    
	#Processing headerlines -> beginning with # or empty lines.
	while(<FILE> =~ /(^\s*#)(.*)/ || <FILE> eq "") {}
    
	my $table_header = <FILE>;
	unless ($printed_header) { #print table header once.
	    print SUMFILE "sample \t". $table_header;
	    $printed_header = 1;
	}

	my $line = <FILE>; #grep statistics, here we assume one row with statistics.
	my ($sample) = ($file =~ m/.+\/(.+)_RNAMetrics.txt/);
	print "\n\n $sample \n\n";
	print SUMFILE $sample ."\t". $line;
    }
    close SUMFILE
}

__END__

=head1 SYNOPSIS

$ perl parsePicardOutput.pl -output_dir /path/to/output -hsmetrics <hsmetrics.txt> -wgsmetrics <wgsmetrics.txt>

=cut