#!/usr/bin/perl -w
use strict;
use warnings;

use File::Basename;
use Carp;
use Getopt::Long;

use MyConfig;
use Data::Dumper;

my $progDir = "~/codes/splice-site-qtl/";

my $prog = basename ($0);

my $evo = 0;
my $mkUniqName = 0;
my $overrideSS = 0;
my $minMaxInc = 0;
my $excludeLow = 0;
my $excludeHigh = 1;

my $cache = getDefaultCache ($prog);

my $verbose = 0;
my $naString = 'NaN';

GetOptions (
	'evo'=>\$evo,
	'uniq-name'=>\$mkUniqName,
	'override-ss'=>\$overrideSS,
	'min-max-inc:f'=>\$minMaxInc,
	'exclude-low:f'=>\$excludeLow,
	'exclude-high:f'=>\$excludeHigh,
	"na-string:s"=>\$naString,
	'cache:s'=>\$cache,
	'v'=>\$verbose);

if (@ARGV != 3)
{
	print "convert canonical transcripts file to bed format\n";
	print "$prog [options] <canonical.txt> <ss.usage.txt> <out.txt\n";

	print " --evo         : the transcript file has evolutionary columns\n";
	print " --override-ss : override splice sites in the input transcript file (need to use --uniq-name if names are not unique)\n";
	print " --uniq-name   : add suffix to make unique transcript name\n";
	print " --min-max-inc  [float]  : min max (constitutive) inclusion level ($minMaxInc)\n";
	print " --exclude-low  [float]  : exclude psi's lower than the specified value ($excludeLow)\n";
	print " --exclude-high [float]  : exclude psi's higher than the specified value ($excludeHigh)\n";
	print " --na-string    [string] : na string (default:$naString)\n";
	print " -v                : verbose\n";
	exit (1);
}


my ($inFile, $ssFile, $outFile) = @ARGV;


my %geneHash;

if ($overrideSS)
{
	print "matching genes and splice sites ...\n" if $verbose;

	my $ret = system ("mkdir $cache");
	Carp::croak "cannot mkdir $cache\n" if $ret != 0;

	my $ssBedFile = "$cache/ss.bed";
	my $cmd = "cut -f 1-6 $ssFile > $ssBedFile";
	print "CMD=$cmd\n" if $verbose;
	$ret = system ($cmd);
	Carp::croak "CMD=$cmd failed\n" if $ret != 0;

	my $tsBedFile = "$cache/ts.bed";
	$cmd = "awk '{print \$3\"\\t\"\$5-1\"\\t\"\$6\"\\t\"\$1\"\\t\"\$2\"\\t\"\$4}' $inFile > $tsBedFile";
	$cmd = "awk 'BEGIN{i=0} {print \$3\"\\t\"\$5-1\"\\t\"\$6\"\\t\"\$1\".\"i\"\\t\"\$2\"\\t\"\$4; i=i+1}' $inFile > $tsBedFile" if $mkUniqName;

	if ($evo)
	{
		$cmd = "cut -f 3- $inFile | awk '{print \$3\"\\t\"\$5-1\"\\t\"\$6\"\\t\"\$1\"\\t\"\$2\"\\t\"\$4}' > $tsBedFile";
    	$cmd = "cut -f 3- $inFile | awk 'BEGIN{i=0} {print \$3\"\\t\"\$5-1\"\\t\"\$6\"\\t\"\$1\".\"i\"\\t\"\$2\"\\t\"\$4; i=i+1}' > $tsBedFile" if $mkUniqName;
	}

	print "CMD=$cmd\n" if $verbose;
    $ret = system ($cmd);
	Carp::croak "CMD=$cmd failed\n" if $ret != 0;

	my $ss2tsBedFile = "$cache/ss_vs_ts.bed";
	$cmd = "perl $progDir/scripts/tagoverlap.pl -ss -d \"##\" -region $tsBedFile $ssBedFile $ss2tsBedFile --keep-cache";
	print "CMD=$cmd\n" if $verbose;
    $ret = system ($cmd);
	Carp::croak "CMD=$cmd failed\n" if $ret != 0;
	
	my $iter1 = 0;
	my $iter2 = 0;

	my $fin;
	open ($fin, "<$ss2tsBedFile") || Carp::croak "cannot open file $ss2tsBedFile to read\n";
	while (my $line = <$fin>)
	{
		chomp $line;
		my ($chrom, $chromStart, $chromEnd, $name, $score, $strand) = split ("\t", $line);
	   	my ($ssId, $tsId) = split("##", $name);

		my $junctionType = substr ($ssId, -3);
	
		if (($strand eq '+' && $junctionType eq '5ss') ||
			($strand eq '-' && $junctionType eq '3ss'))
		{
			push @{$geneHash{$tsId}{'jStarts'}}, $chromEnd;
			$iter1++;
		}
		else
		{
			push @{$geneHash{$tsId}{'jEnds'}}, $chromEnd;
			$iter2++;
		}
	}
	close ($fin);

	print "$iter1 junction starts, $iter2 junction ends matched to genes\n" if $verbose;
}



my ($fin, $fout);


#my $tsBedFile = "$cache/ts.bed";

print "loading splice site usage ...\n" if $verbose;

my (%junctionStartHash, %junctionEndHash);

my $sampleNum = 0;
my $iter1 = 0;
my $iter2 = 0;
open ($fin, "<$ssFile") || Carp::croak "cannot open file $ssFile to read\n";
while (my $line = <$fin>)
{
	chomp $line;
	next if $line =~/^\#/;
	next if $line =~/^\s*$/;

	my ($chrom, $chromStart, $chromEnd, $name, $score, $strand, $type, $isoforms, @cols) = split ("\t", $line);

	if ($sampleNum == 0)
	{
		$sampleNum = @cols;
	}
	else
	{
		Carp::croak "inconsistent number of samples: $line\n" if $sampleNum != @cols;
	}

	my $junctionType = substr ($name, -3); 

	#Carp::croak "strand = $strand, jtype=$junctionType\n";

	if (($strand eq '+' && $junctionType eq '5ss') ||
		($strand eq '-' && $junctionType eq '3ss'))
	{
		$junctionStartHash{$chrom}{$strand}{$chromEnd} = \@cols;
		$iter1++;
	}
	else
	{
		$junctionEndHash{$chrom}{$strand}{$chromEnd} = \@cols;
		$iter2++;
	}
}


close ($fin);

print "$iter1 junction starts, $iter2 junction ends loaded\n" if $verbose;

print "add splice site usage ...\n" if $verbose;

open ($fin, "<$inFile") || Carp::croak "cannot open file $inFile to read\n";
open ($fout, ">$outFile") || Carp::croak "cannot open file $outFile to write\n";

my $iter = 0;
while (my $line = <$fin>)
{
	chomp $line;
	next if $line=~/^\s*$/;

	my @cols = split ("\t", $line);

	my ($db, $chromRef) = ("", "");
	if ($evo)
	{
		$db = shift @cols; 
		$chromRef = shift @cols;
	}

	my ($name0, $paralog, $chrom, $strand, $tsStart, $tsEnd, $jStartsStr, $jEndsStr) = @cols;

	my $name = $name0;
	$name .= "." . $iter if $mkUniqName;
	$iter++;
	
	my (@jStarts, @jEnds);

	if ($overrideSS)
	{
		if (exists $geneHash{$name} && exists $geneHash{$name}{'jStarts'} && exists $geneHash{$name}{'jEnds'})
		{
			@jStarts = sort {$a <=> $b} @{$geneHash{$name}{'jStarts'}};
			@jEnds = sort {$a <=> $b} @{$geneHash{$name}{'jEnds'}};
		}
		else
		{
			warn "Either 5' or 3' splice sites not detected for $name, skip\n";
			next;
		}
	}
	else
	{
		chop $jStartsStr if $jStartsStr=~/\,$/;
		@jStarts = split (",", $jStartsStr);

		chop $jEndsStr if $jEndsStr=~/\,$/;
		@jEnds = split (",", $jEndsStr);
	}


	my ($maxJStart, $jStartUsage) = retrieveSpliceSiteUsage (\%junctionStartHash, $chrom, $strand, \@jStarts);
	my ($maxJEnd, $jEndUsage) = retrieveSpliceSiteUsage (\%junctionEndHash, $chrom, $strand, \@jEnds);

	if ($evo)
	{
		print $fout join ("\t", $db, $chromRef), "\t";
	}

	print $fout join ("\t", $name0, $paralog, $chrom, $strand, $tsStart, $tsEnd);
	print $fout "\t", join ("\t", join(",", @jStarts) . ",", join(",", @jEnds) . ",");

	for (my $i = 0; $i < @$jStartUsage; $i++)
	{
		if ($maxJStart >= $minMaxInc || $maxJEnd >= $minMaxInc)
		{
			print $fout "\t", join("\t", $jStartUsage->[$i], $jEndUsage->[$i]);
		}
		else
		{
			my ($jStartStr, $jEndStr) = ("", "");
			for (my $i = 0; $i < @jStarts; $i++)
			{
				$jStartStr .= $naString . ",";
			}
			for (my $i = 0; $i < @jEnds; $i++)
			{
				$jEndStr .= $naString . ",";
			}
			print $fout "\t", join("\t", $jStartStr, $jEndStr);
		}
	}
	print $fout "\n";
}

close ($fin);
close ($fout);

#system ("rm -rf $cache");

sub retrieveSpliceSiteUsage
{
	my ($jHash, $chrom, $strand, $jPos) = @_;
	my @usage;

	my $maxInc = 0;
	foreach my $p (@$jPos)
	{
		for (my $j = 0; $j < $sampleNum; $j++)
		{
			my $psi = exists $jHash->{$chrom}{$strand}{$p} ? $jHash->{$chrom}{$strand}{$p}[$j] : $naString;
			$psi = $naString if $psi eq 'NA' || $psi eq 'NaN';
			
			if ($psi ne $naString)
			{
				$maxInc = $psi if $psi >= $maxInc;
				$psi = $naString if ($psi < $excludeLow || $psi > $excludeHigh);
				#Carp::croak "psi=$psi, excludeLow=$excludeLow, excludeHigh=$excludeHigh\n";
			}

			#$psi = sprintf("%.3f", $psi) if $psi ne $naString;

			$usage[$j] .= $psi . ",";
		}
	}
	return ($maxInc, \@usage);
}

