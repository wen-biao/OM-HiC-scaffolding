#!/usr/bin/perl

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
my $cmap_file = $ARGV[0];

  my ($cmap_data, $numContig, $contigLength_ref) = readCMap($cmap_file);
  my @contigLength=@$contigLength_ref;

  my ($min, $max, $mean, $median, $n50value, $total) = getContigStat(@contigLength);
  print "N contigs = $numContig\n";
  if ($numContig == 0) {
	  print "Min contig length = N/A Mb\n";
	  print "Median contig length = N/A Mb\n";
	  print "Mean contig length = N/A Mb\n";
	  print "Contig N50 = N/A Mb\n";
	  print "Max contig length = N/A Mb\n";
	  print "Total contig length = N/A Mb\n";  	
  } else {
	  print "Min contig length = " . $min/1000000 . " Mb\n";
	  print "Median contig length = " . $median/1000000 . " Mb\n";
	  print "Mean contig length = " . $mean/1000000 . " Mb\n";
	  print "Contig N50 = " . $n50value/1000000 . " Mb\n";
	  print "Max contig length = " . $max/1000000 . " Mb\n";
	  print "Total contig length = " . $total/1000000 . " Mb\n";
  }

=pod

=head1 NAME

calc_cmpa_stats.pl - Calculate stats about CMAPs

=head1 SYNOPSIS

calc_cmpa_stats.pl <CMAP_File>

Examples:

    calc_cmap_stat.pl Ler0_GM0327.cmap

=head1 DESCRIPTION

This script is used to caculate statitical data about a CMAP.

[1] "N contigs = 173"
[1] "Min contig length = 0.1436865 Mb"
[1] "Median contig length = 0.481299 Mb"
[1] "Mean contig length = 0.624446769364162 Mb"
[1] "Contig N50 = 0.75731195 Mb"
[1] "Max contig length = 2.683027 Mb"
[1] "Total contig length = 108.0292911 Mb"


=head1 ARGUMENTS

calc_cmap_stat.pl takes the following arguments:

=over

=item <CMAP_File>

A CMAP File with version 0.1

=back

=cut
