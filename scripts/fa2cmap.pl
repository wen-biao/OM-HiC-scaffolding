# $Id: fa2cmap.pl 3281 2014-09-28 22:57:49Z psheth $

#!/usr/bin/perl -w
########################################################################
# File: fa2cmap.pl                                                     #
# Date: 08/27/2014                                                     #
# Purpose: Transform fasta file to BioNano cmap file format            #
#                                                                      #
# Author: Xiang Zhou, Computational Biologist                          #
# Email : xzhou@bionanogenomics.com                                    #
# Affiliation: Research Department, BioNano Genomics Inc.              #
#                                                                      #
# Usage:                                                               #
#   fa2cmap.pl <-h> <-v> <-i fasta_file> <-n|-s enzyme_name|enzyme_seq>#
#              <-m min_labels_filter> <-M min_size_filter(Kb)>         #
#     -h    : This help message                                        #
#     -v    : Verbose output (Default: OFF)                            #
#     -i    : Input fasta file                                         #
#     -n    : Using enzyme by name                                     #
#     -s    : Using enzyme by sequence                                 #
#     -m    : Filter: Minimum labels  (Default: 0)                     #
#     -M    : Filter: Minimum size (Kb)  (Default: 0)                  #
#                                                                      #
#   NOTE: CMAP index is 1-based                                        #
########################################################################

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

use POSIX;
use Getopt::Std;
use File::Slurp qw(edit_file_lines);
use BNG::Utility;
#$Getopt::Std::STANDARD_HELP_VERSION = 1;

sub Init;
sub Usage;
sub Find_enzymes;
sub Print_cmap_header;
sub Generate_cmap;
sub Uniq;

my %enzyme = (
	"BspQI" => "GCTCTTC",
	"BbvCI" => "CCTCAGC",
	"BsmI"  => "GAATGC",
	"BsrDI" => "GCAATG",
	"bseCI" => "ATCGAT"
	
	# You can add more enzymes here ...
	
);
my $current_enzyme;
my ($min_labels, $min_length) = (0, 0);
my ($FASTA, $CMAP, $KEY);
my ($filename, $filename_key);

my %opts;
my @loci = ();
my ($count, $count_valid) = (0, 0);
my @density;
my $density;
my ($total_length_original, $total_length_filtered, $total_nicks, $num_cmaps) = (0) x 4;

my ($command, $seq);
my ($fastaHeader, $fastaHeader_previous);
my ($A, $C, $G, $T, $N, $global_N, $global_GC, $global_ACGT, $global_ACGTN);
my ($N_percentage, $GC_percentage, $global_N_percentage, $global_GC_percentage);

Init();

open($FASTA, $opts{i}) || dieLog ("ERROR: Can't open $opts{i}: $!\n");
$filename = $opts{i};
$filename =~ s/(\S+)\.\w+$/$1_$current_enzyme/;
$filename = $filename."_".$min_length."Kb_".$min_labels."labels.cmap";
$filename_key = $opts{i};
$filename_key =~ s/(\S+)\.\w+$/$1_$current_enzyme/;
$filename_key = $filename_key."_".$min_length."Kb_".$min_labels."labels_key.txt";

open($KEY, ">$filename_key") || dieLog ("ERROR: Can't open $filename_key: $!\n");
print $KEY("# CMAP = $filename\n");
print $KEY("# filter: Minimum Labels = $min_labels\n");
print $KEY("# filter: Minimum Size (Kb) = $min_length\n");
print $KEY("CompntId\tCompntName\tCompntLength\n");

# Verbose output?
if(defined($opts{v})){
	print("Input file:\t$opts{i}\n");
	print("Output file:\t$filename\n");
	print("ID keyfile:\t$filename_key\n");
	print("Enzyme used:\t$current_enzyme($enzyme{$current_enzyme})\n\n");

	print("Contigs\t             Length\t         Frequency\t           GC%\t  N%\n");
	print("-----------------------------------------------------------------------\n");
}

Print_cmap_header($filename);

while(my $line = <$FASTA>){
	chomp $line;
	$line =~ s/\r//g;
	
	if($line =~ /^>/){
		$fastaHeader_previous = $fastaHeader;
		$fastaHeader = substr($line, 1);
		
		if($count != 0){
			$A = ($seq =~ tr/A/A/);
			$C = ($seq =~ tr/C/C/);
			$G = ($seq =~ tr/G/G/);
			$T = ($seq =~ tr/T/T/);
			$N = ($seq =~ tr/N/N/);
			$global_N += $N;
			$global_GC += ($C+$G);
			$global_ACGT += ($A+$C+$G+$T);
			$global_ACGTN += ($A+$C+$G+$T+$N);
			
			$total_length_original += length($seq);
			
			@loci = Find_enzymes($seq, $enzyme{$current_enzyme});
			
			# Filter by $min_labels and $min_length
			if(scalar(@loci) >= $min_labels && length($seq) >= $min_length * 1000){
				Generate_cmap($filename, $count, \@loci, length($seq));
				print $KEY("$count\t$fastaHeader_previous\t", length($seq), "\n");
				
				$density[$count] = sprintf("%7.3f", @loci/length($seq) * 100000);
				my $length_tmp = sprintf("%11.3f", length($seq)/1000000);
				
				if($A+$C+$G+$T == 0){
					$GC_percentage = sprintf("%5.2f", 0);
				}
				else{
					$GC_percentage = sprintf("%5.2f", ($C+$G)/($A+$C+$G+$T)*100);
				}
				if($N == 0){
					$N_percentage = sprintf("%6.2f", 0);
				}
				else{
					$N_percentage = sprintf("%6.2f", $N/length($seq)*100);
				}
				
				# Verbose output?
				if(defined($opts{v})){
					print("Contig[$count]:\t$length_tmp (MB)\t$density[$count] nick(s)/100KB\t$GC_percentage%\t$N_percentage%\n");
				}
				
				$total_length_filtered += length($seq);
				$total_nicks += @loci;
				
				$count_valid++;
			}
		}
		
		$seq = "";
		$count++;
	}
	else{
		$seq .= uc($line);
	}
}
$A = ($seq =~ tr/A/A/);
$C = ($seq =~ tr/C/C/);
$G = ($seq =~ tr/G/G/);
$T = ($seq =~ tr/T/T/);
$N = ($seq =~ tr/N/N/);
$global_N += $N;
$global_GC += ($C+$G);
$global_ACGT += ($A+$C+$G+$T);
$global_ACGTN += ($A+$C+$G+$T+$N);

$total_length_original += length($seq);

@loci = Find_enzymes($seq, $enzyme{$current_enzyme});

# Filter by $min_labels and $min_length
if(scalar(@loci) >= $min_labels && length($seq) >= $min_length * 1000){
	Generate_cmap($filename, $count, \@loci, length($seq));
	print $KEY("$count\t$fastaHeader\t", length($seq), "\n");
	
	$density[$count] = sprintf("%7.3f", @loci/length($seq) * 100000);
	my $length_tmp = sprintf("%11.3f", length($seq)/1000000);
	if($A+$C+$G+$T == 0){
		$GC_percentage = sprintf("%5.2f", 0);
	}
	else{
		$GC_percentage = sprintf("%5.2f", ($C+$G)/($A+$C+$G+$T)*100);
	}
	if($N == 0){
		$N_percentage = sprintf("%6.2f", 0);
	}
	else{
		$N_percentage = sprintf("%6.2f", $N/length($seq)*100);
	}
	# Verbose output?
	if(defined($opts{v})){
		print("Contig[$count]:\t$length_tmp (MB)\t$density[$count] nick(s)/100KB\t$GC_percentage%\t$N_percentage%\n");
	}
	
	$total_length_filtered += length($seq);
	$total_nicks += @loci;
	
	$count_valid++;
}
close($FASTA);
close($KEY);

#$command = "perl -pi -e 's|N/A|$num_cmaps| if " . '$.' . " == 5' $filename";
#print("$command\n");
#system("$command");

edit_file_lines {s |N/A|$num_cmaps| } $filename;

# Verbose output?
if(defined($opts{v})){
	print("\n================================Summary================================\n");
	print("Total contigs processed: $count\n");
	print("Total contigs generated: $count_valid\n");
	print("Total length of the original contigs: $total_length_original\n");
	print("Total length of the filtered contigs: $total_length_filtered\n");
	$density = sprintf("%6.3f", $total_nicks/$total_length_filtered * 100000);
	$global_GC_percentage = sprintf("%5.3f", $global_GC/$global_ACGT*100);
	$global_N_percentage = sprintf("%5.3f", $global_N/$global_ACGTN*100);
	print("Global nick frequency:\t$density nick(s) /100KB\n");
	print("Global GC percentage:\t$global_GC_percentage%\n");
	print("Global N percentage:\t$global_N_percentage%\n");
	print("=======================================================================\n");
}

######################################################################
#                           Subroutines                              #
######################################################################
sub Init{
	my $opt_string = 'hvi:n:s:m:M:';
	if(!getopts("$opt_string", \%opts)){
		print("ERROR: Invalid parameter(s)! Try -h for more information.\n");
		Usage();
	}
	Usage() if $opts{h};
	
	if(defined($opts{n})){
		if(uc(substr($opts{n}, 0, 4)) eq "BSPQ"){
			$current_enzyme = "BspQI";
		} elsif(uc(substr($opts{n}, 0, 4)) eq "BBVC"){
			$current_enzyme = "BbvCI";
		} elsif(uc(substr($opts{n}, 0, 3)) eq "BSM"){
			$current_enzyme = "BsmI";
		} elsif(uc(substr($opts{n}, 0, 4)) eq "BSRD"){
			$current_enzyme = "BsrDI";
		} elsif(uc(substr($opts{n}, 0, 4)) eq "BSEC"){
			$current_enzyme = "bseCI";
		} else	{
			dieLog ("ERROR: Unknown enzyme specified.\n");
		}
		
		# You may add more enzymes ...
	}
	elsif(defined($opts{s})){

		if(uc($opts{s}) eq $enzyme{BspQI}){
			$current_enzyme = "BspQI";
		} elsif(uc($opts{s}) eq $enzyme{BbvCI}){
			$current_enzyme = "BbvCI";
		} elsif(uc($opts{s}) eq $enzyme{BsmI}){
			$current_enzyme = "BsmI";
		} elsif(uc($opts{s}) eq $enzyme{BsrDI}){
			$current_enzyme = "BsrDI";
		} elsif(uc($opts{s}) eq $enzyme{bseCI}){
			$current_enzyme = "bseCI";
		} else	{
			dieLog ("ERROR: Unknown enzyme specified.\n");
		} 
	}
	else{
		print("ERROR: Invalid parameter(s)! Try -h for more information.\n");
		Usage();
	}
	if(!defined($opts{i})){
		print("ERROR: Invalid parameter(s)! Try -h for more information.\n");
		Usage();
	}
	if(defined($opts{m})){
		$min_labels = $opts{m};
		if($min_labels < 0){
			print("ERROR: Invalid input value(s)! Try -h for more information.\n");
			Usage();
		}
	}
	if(defined($opts{M})){
		$min_length = $opts{M};
		if($min_length < 0){
			print("ERROR: Invalid input value(s)! Try -h for more information.\n");
			Usage();
		}
	}
}

sub Usage{
	print << "EOF";

Usage: $0 <-h> <-v> <-i fasta_file> <-n|-s enzyme_name|enzyme_seq>
       <-m min_labels_filter> <-M min_size_filter(Kb)>
  -h    : This help message
  -v    : Verbose output [Default: OFF]
  -i    : Input fasta file
  -n    : Using enzyme by name
  -s    : Using enzyme by sequence
  -m    : Filter: Minimum labels  (Default: 0)
  -M    : Filter: Minimum size (Kb)  (Default: 0)

NOTE: CMAP index is 1-based
EOF
	exit;
}

sub Find_enzymes{
	my ($seq, $enzyme) = @_;
	my @result;
	
	# Find the enzymes in the forward strand, staring from the first nucleotide!!!
	my $current_loc = index($seq, $enzyme, 0);
	while ($current_loc != -1){
		push(@result, $current_loc+1);
		$current_loc = index($seq, $enzyme, $current_loc + 1);
	}
	
	my $enzyme_rc = reverse($enzyme);
	$enzyme_rc =~ tr/ACGTUN/TGCAAN/;
	
	# Find the rc(enzymes) in the forward strand, staring from the first nucleotide!!!
	$current_loc = index($seq, $enzyme_rc, 0);
	while ($current_loc != -1){
		push(@result, $current_loc+1);
		$current_loc = index($seq, $enzyme_rc, $current_loc + 1);
	}
	
	return Uniq(@result);
}

sub Print_cmap_header{
	my ($filename) = @_;
	my $OUT;
	
	open($OUT, ">$filename") || dieLog ("ERROR: Can't open $filename: $!\n");
	
	my $str = << "EOF";
# CMAP File Version:	0.1
# Label Channels:	1
# Nickase Recognition Site 1:	$enzyme{$current_enzyme}
# Enzyme1:	Nt.$current_enzyme
# Number of Consensus Nanomaps:	N/A
#h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence
#f int	float	int	int	int	float	float	int	int
EOF
	print $OUT($str);
	close($OUT);
}

sub Generate_cmap{
	my ($filename, $ID, $loci_ref, $length) = @_;
	
	my $OUT;
	my $i;
	my $length_float = sprintf("%.1f", $length);
	my @sorted_loci = sort {$a <=> $b} @$loci_ref;
	
	open($OUT, ">>$filename") || dieLog ("ERROR: Can't open $filename: $!\n");
	
	for($i = 0; $i < @sorted_loci; $i++){
		my $loci_float = sprintf("%.1f", $sorted_loci[$i]);
		print $OUT("$ID\t$length_float\t", scalar(@sorted_loci), "\t", $i+1, "\t1\t$loci_float\t1.0\t1\t1\n");
	}
	if(scalar(@sorted_loci) != 0){
		print $OUT("$ID\t$length_float\t", scalar(@sorted_loci), "\t", $i+1, "\t0\t$length_float\t0.0\t1\t0\n");
		$num_cmaps++;
	}
	close($OUT);
}

sub Uniq{
	my %seen;
	grep {!$seen{$_}++} @_;
}

__END__

./fa2cmap.pl -i test.fa -n bspq -v -m 0 -M 0
./fa2cmap.pl -i test.fa -s GCTCTTC -v -m 0 -M 0

