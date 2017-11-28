#!/usr/bin/env perl

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
my @diff_chr_percentages;
my @diff_chr_mapq_percentages;

print "Parsing flagstat files \n";
foreach my $file (@flagstat_files) {
    print "\t Parsing: ". $file . "\n";
    open(FILE, $file) || die ("Can't open $file");

    my ($sample) = ($file =~ m/.+\/(.+).flagstat/);
    push(@samples, $sample);
    my $total_reads;
    my $mapped_reads;
    
    while (my $line = <FILE>){
	if ($. == 1) {
	    $total_reads = (split(' ',$line))[0];
	    push(@read_counts, $total_reads);
	}
	if ($. == 5) {
	    $mapped_reads = (split(' ',$line))[0];
	    my $mapped_percentage = ($mapped_reads / $total_reads) * 100;
	    push(@mapped_percentages, sprintf('%.2f%%', $mapped_percentage));
	}
	if ($. == 12){
	    my $diff_chr = (split(' ',$line))[0];
	    my $diff_chr_percentage = ($diff_chr / $total_reads) * 100;
	    push(@diff_chr_percentages, sprintf('%.2f%%', $diff_chr_percentage));
	    
	}
	if ($. == 13){
	    my $diff_chr_mapq = (split(' ',$line))[0];
	    my $diff_chr_mapq_percentage = ($diff_chr_mapq / $total_reads) * 100;
	    push(@diff_chr_mapq_percentages, sprintf('%.2f%%', $diff_chr_mapq_percentage));
	}
    }
    close FILE
}

print SUMFILE "sample\t". join("\t",@samples) ."\n";
print SUMFILE "Total number of reads\t". join("\t",@read_counts) ."\n";
print SUMFILE "Percentage reads mapped\t". join("\t",@mapped_percentages) ."\n";
print SUMFILE "Percentage reads with mate mapped to different chr\t". join("\t",@diff_chr_percentages) ."\n";
print SUMFILE "Percentage reads with mate mapped to different chr (mapQ>=5)\t". join("\t",@diff_chr_mapq_percentages) ."\n";


close SUMFILE

__END__

=head1 SYNOPSIS

$ perl parseFlagstatOutput.pl -output_dir /path/to/output -flagstat <my_file.flagstat>

=cut
