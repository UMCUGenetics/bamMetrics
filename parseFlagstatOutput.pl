#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

### Input options ###

# Declare options
my @flagstat_files;
my $output_dir = "";

# Parse options
GetOptions ("flagstat=s" => \@flagstat_files,
	    "output_dir=s" => \$output_dir
	) or pod2usage(1);

if(! @flagstat_files || ! $output_dir) { pod2usage(1) };

### Parse flagstat files
my $filename = $output_dir."/flagstat_summary.txt";
open(SUMFILE, ">", $filename) || die ("Can't create $filename");

my @samples;
my @read_counts;
my @mapped_percentages;

print "Parsing flagstat files \n";
foreach my $file (@flagstat_files) {
    print "\t Parsing: ". $file . "\n";
    open(FILE, $file) || die ("Can't open $file");
    
    my ($sample) = ($file =~ m/.+\/(.+).flagstat/);
    push(@samples, $sample);
    while (my $line = <FILE>){
	if ($. == 1 ) {
	    my $total_reads = (split(' ',$line))[0];
	    push(@read_counts, $total_reads);
	}
	if ($. == 5 ) { 
	    my ($mapped_percentage) = $line =~ /(\d{2}\.\d{2}\%)/;
	    push(@mapped_percentages, $mapped_percentage);
	}
    }
    close FILE
}

print SUMFILE "sample\t". join("\t",@samples) ."\n";
print SUMFILE "Total number of reads\t". join("\t",@read_counts) ."\n";
print SUMFILE "Percentage reads mapped\t". join("\t",@mapped_percentages) ."\n";

close SUMFILE

__END__

=head1 SYNOPSIS

$ perl parsePicardOutput.pl -output_dir /path/to/output -flagstat <my_file.flagstat>

=cut