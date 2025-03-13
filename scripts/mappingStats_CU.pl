#!/usr/bin/perl -w

use strict;

if (scalar(@ARGV) != 2)
{
	print "Returns mapping stats assuming existence of basename.sam[.gz] and basename.bed.\n";
	print "Results are appended to outfile, so can specify same outfile for multiple simultaneous samples.\n";
	print "Usage: mappingStats.pl basename outfile\n";
	exit (1);
}


my ($basename, $outfile) = @ARGV;

my $gz="";
$gz=".gz" if (-e "$basename.sam.gz");

die "Could not find $basename.sam[.gz]\n" if ! -e "$basename.sam$gz";
die "Could not find $basename.bed\n" if ! -e "$basename.bed"; 

# total reads
print "Counting total reads...\n";
my $total = trim(`samtools view -S $basename.sam$gz | wc -l`);

# mapped reads
print "Counting mapped reads...\n";
my $mapped = trim(`samtools view -S -F 4 $basename.sam$gz | wc -l`);
my $mappedPercent = sprintf "%.3f", $mapped/$total;
my $mappedPercentPretty = 100*$mappedPercent . "%";

# unique mapped reads
print "Counting unique mapped reads...\n";
my $mappedUniq = trim(`cat $basename.bed | wc -l`);
my $mappedUniqPercent = sprintf "%.3f", $mappedUniq/$total;
my $mappedUniqPercentPretty = 100*$mappedUniqPercent . "%";

# junction reads ($6 is CIGAR string)
# my $junction = trim(`awk '{if (match(\$6, /N/)){print \$0}}' $basename.sam$gz | wc -l`);

# unique junction reads ($11 is block count):
print "Counting unique junction reads...\n";
my $junctionUniq = trim(`awk '{if (match(\$11, /,/)) {print \$0}}' $basename.bed | wc -l`);
my $junctionUniqPercent = sprintf "%.3f", $junctionUniq/$total;
my $junctionUniqPercentPretty = 100*$junctionUniqPercent . "%";

my $out = "$basename\t$total\t$mapped\t$mappedPercentPretty\t$mappedUniq\t$mappedUniqPercentPretty\t$junctionUniq\t$junctionUniqPercentPretty\n";
print $out;
open (OUTFILE, ">>$outfile");
print OUTFILE $out;
close (OUTFILE);
print "Done.\n";

sub trim {
    (my $s = $_[0]) =~ s/^\s+|\s+$//g;
    return $s;        
}
