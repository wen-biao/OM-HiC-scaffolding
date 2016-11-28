# $Id: findNGS.pl 3200 2014-09-05 14:03:58Z psheth $

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

# this script identifies which NGS contigs (id < 100000) went into the hybrid scaffold

my $homeDir = $ARGV[0];
my $inFile = "step1.merge.pairs.txt";

my $outFile = $ARGV[1];

open(IN, "$homeDir/$inFile") or dieLog ("ERROR: Cannot open $homeDir/$inFile: $!\n");
my $count = 0;	print "reading in $inFile\n";
my $skip = <IN>;
my %nrNgs = ();
while (my $line = <IN>)	{
	chomp $line;
	$line =~ s/\r//g;
	$line =~ s/"//g;
	$line =~ s/\s+/\t/g;

	my @content = split(/\t/, $line);
	my ($id1, $id2, $hybridId, $stage) = @content[1..3];
$count += 1;
	
	if ($id1 < 100000)	{
		$nrNgs{$id1} = 1;
	} # if id1

	if ($id2 < 100000)	{
		$nrNgs{$id2} = 1
	} # if id2
} # while line	
print "\tread in $count records\n";
close IN;

open(OUT, ">$homeDir/$outFile") or dieLog ("ERROR: Cannot write to $homeDir/$outFile: $!\n");
print "writing to $outFile\n";	my $outCount = 0;
#print OUT "ngsContig"."\n";
foreach my $id (sort numeric keys %nrNgs)	{
	print OUT $id."\n";	
$outCount += 1;
} # foreach 
print "\twritten $outCount records\n";
close OUT;

sub numeric	{$a <=> $b}
