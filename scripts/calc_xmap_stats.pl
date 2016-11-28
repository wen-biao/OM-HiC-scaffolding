# $Id: calc_xmap_stats.pl 3200 2014-09-05 14:03:58Z psheth $

#!/usr/bin/perl -w

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);

# This adds "${CURRENT_SCRIPT_PATH}/perl5/" direcory at run time to the @INC array
# This script sould sit one level above the additional Perl modules directory.
BEGIN {
	my $script_path = abs_path(dirname($0));
	my $module_path = abs_path($script_path . "/perl5");
	my $lib = glob("~/scripts/HybridScaffold/scripts/perl5");
	unshift @INC, $module_path;
	unshift @INC, $lib;
	my $lib2; my $lib3;
	if ($] >= 5.010000 && $] <= 5.011000) {
		$module_path = $module_path."/5.10.1"; 
		$lib = $lib."/5.10.1";
		$lib2 = $lib."/x86_64-linux-thread-multi";
		$lib3 = $module_path."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.014000 && $] <= 5.015000) {
		$module_path = $module_path."/5.14.4"; 
		$lib = $lib."/5.14.4";
		$lib2 = $lib."/x86_64-linux-thread-multi";
		$lib3 = $module_path."/x86_64-linux-thread-multi";}
	else {
		print "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X\n";
		die "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X\n"; 
		exit; }
	unshift @INC, $module_path;
	unshift @INC, $lib;
	unshift @INC, $lib2;
	unshift @INC, $lib3;
	#print "$0 library paths:\n\t"; print join("\n\t", @INC); print "\n";
}

use BNG::Utility;

#my $xmapDir = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results/v1_ngs_to_hybrid";
my $xmapDir = $ARGV[0];
#my $xmapFile = "v1_ngs_to_hybrid.xmap";
my $xmapFile = $ARGV[1];

#my $refDir = "/mnt/bionf_tmp/apang/triobppAdj2/NA12878/ngs_alignment_pac_bio_mtSinai/rerun_coassembly_050814/v2_output_T10/output_results";
#my $refFile = "step2.hybrid.cmap";
my $refDir = $ARGV[0];
my $refFile = $ARGV[2];

#my $outDir = "$xmapDir";
#my $outFile = "ngs_coverage_040514.txt";
my $outDir = $ARGV[0];
my $outFile = $ARGV[3];

my $cmapLengthRef = getCmapLength($refDir, $refFile);
$cmapLengthRef = getAlignment($xmapDir, $xmapFile, $cmapLengthRef);
$cmapLengthRef = calcCoverage($cmapLengthRef);
printCoverage($outDir, $outFile, $cmapLengthRef);

			#$cmapLength{$content[0]} = {size => $content[1], overlapRef => [], sumCoverage => 0, nrCoverage => 0, nrPercentCoverage => 0};
sub printCoverage	{
	my ($dir, $file, $cmapLengthRef) = @_;
	open(OUT, ">$dir/$file") or dieLog ("ERROR: printCoverage: cannot write to $dir/$file: $!\n");
my $count = 0;	
#print "printCoverage: writing to $file\n";
	print OUT "refContig\trefContig_size\tquery_cov_bp\tquery_unique_cov_bp\tquery_unique_cov_percent\n";
	print "refContig\trefContig_size\tquery_cov_bp\tquery_unique_cov_bp\tquery_unique_cov_percent\n";
	foreach my $id (keys %$cmapLengthRef)	{
		my $cRef = $cmapLengthRef->{$id};
		my $sNrPercentCoverage = sprintf("%.2f", $cRef->{nrPercentCoverage});
		print OUT join("\t", $id, $cRef->{size}, $cRef->{sumCoverage}, $cRef->{nrCoverage}, $sNrPercentCoverage)."\n";
		print join("\t", $id, $cRef->{size}, $cRef->{sumCoverage}, $cRef->{nrCoverage}, $sNrPercentCoverage)."\n";
$count += 1;
	} # foreach id
#print "\twritten $count records\n";
	close OUT;
} # printCoverage

sub calcCoverage	{
	my ($cmapLengthRef) = @_;
	foreach my $id (keys %$cmapLengthRef)	{
		next if (scalar(@{$cmapLengthRef->{$id}{overlapRef}}) == 0);
		$cmapLengthRef->{$id}{overlapRef} = sortByCoord($cmapLengthRef->{$id}{overlapRef}, "start");

		my $oStart = $cmapLengthRef->{$id}{overlapRef}[0]{start};	
		my $oEnd = $cmapLengthRef->{$id}{overlapRef}[0]{end};
		my $oSizeSum = 0;
		my $redundantSizeSum = $oEnd - $oStart;
		for (my $i = 1; $i < scalar(@{$cmapLengthRef->{$id}{overlapRef}}); $i += 1)	{
			my $iRef = $cmapLengthRef->{$id}{overlapRef}[$i];
			# update redundantSizeSum
			$redundantSizeSum += ($iRef->{end} - $iRef->{start});
			
			if ($oStart <= $iRef->{start} && $iRef->{start} <= $oEnd)	{
				# overlap, iterate
				$oEnd = $iRef->{end};
			} else	{
				# non overlap, update
				$oSizeSum += ($oEnd - $oStart);
				# reset
				$oStart = $iRef->{start};
				$oEnd = $iRef->{end};
			} # if overlap
		} # for i
		# update the last interval
		$oSizeSum += ($oEnd - $oStart);

		$cmapLengthRef->{$id}{sumCoverage} = $redundantSizeSum;
		$cmapLengthRef->{$id}{nrCoverage} = $oSizeSum;
		$cmapLengthRef->{$id}{nrPercentCoverage} = ($oSizeSum / $cmapLengthRef->{$id}{size} * 100);

	} # foreach id
	return $cmapLengthRef;
} # calcCoverage

sub sortByCoord	{
	my ($hRef, $coord) = @_;
	@$hRef = sort	{
		$a->{$coord}	<=>	$b->{$coord}
	} @$hRef;
	return $hRef;
} # sortByCoord

sub getAlignment	{
	my ($dir, $file, $cmapLengthRef) = @_;
	open(IN, "$dir/$file") or dieLog ("ERROR: getAlignment: cannot read $dir/$file: $!\n");
my $count = 0;	#print "getAlignment: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		next if ($line =~ /^#/);
		$line =~ s/^\s+//;	$line =~ s/\s+/\t/g;
		my @content = split(/\t/, $line);
		my ($qId, $refId, $qStart, $qEnd, $start, $end) = @content[1..6];
		dieLog ("ERROR: getAlignment: cannot get reference contig '$refId' information\n") if (! exists $cmapLengthRef->{$refId});
		push(@{$cmapLengthRef->{$refId}{overlapRef}}, {start => $start, end => $end});
$count += 1;
	} # while line
#print "\tread in $count records\n";
	close IN;
	return $cmapLengthRef;
} # getAlignment

sub getCmapLength	{
	my ($dir, $file) = @_;
	my %cmapLength = ();
	open(IN, "$file") or dieLog ("ERROR: getCmapLength: cannot read in $file: $!\n");
my $count = 0;	#print "getCmapLength: reading in $file\n";
	while (my $line = <IN>)	{
		chomp $line;	$line =~ s/\r//g;
		next if ($line =~ /^#/);
		my @content = split(/\t/, $line);
		if ($content[4] == 0)	{
			# labelChannel is 0
			$cmapLength{$content[0]} = {size => $content[1], overlapRef => [], sumCoverage => 0, nrCoverage => 0, nrPercentCoverage => 0};
$count += 1;
		} # if content
	} # while line
#print "\tread in $count records\n";
	close IN;
	return \%cmapLength;
} # getCmapLength
