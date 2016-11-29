# $Id: hybridScaffold.pl 3516 2015-01-17 22:33:36Z psheth $

#!/usr/bin/perl -w
#########################################################################################################################################
# File: hybridScaffold.pl                                                                  												#
# Date: 07/24/2014                                                                         												#
# Purpose: Merge NGS scaffolds with BioNano CMAP                                           												#
#                                                                                          												#
#                                                                                          												#
# Usage:                                                                                   												#
#  hybridScaffold.pl <-h> <-n ngs_file> <-b bng_cmap_file> <-c config_file> <-o output_folder> <-B> <-N> <-f> <-m molecules_bnx> <-v>	#
#     -h    : This help message                                                            												#
#     -n    : Input NGS FASTA or CMAP file [required]                                      												#
#     -b    : Input BioNano CMAP [required]                                                												#
#     -c    : Merge configuration file [required]                                          												#
#	  -o    : Output folder [required]                                                     												# 
#     -B    : Use unfiltered BioNano CMAP                                               												#
#     -N    : Use unfiltered NGS CMAP                                                   												#
#	  -f    : Force output and overwrite any existing files                                												#
#     -m    : Input BioNano molecules BNX for molecules to hybrid scaffold alignment [optional]    										#
#	  -v    : Print pipeline version information																					    #
#																			               												#
#########################################################################################################################################

use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);

# This adds "${CURRENT_SCRIPT_PATH}/scripts/perl5/" direcory at run time to the @INC array
# This script sould sit two levels above the additional Perl modules directory
BEGIN {
	my $script_path = abs_path(dirname($0));
	my $module_path = abs_path($script_path . "/scripts/perl5");
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

use IO::Handle;
use PerlIO::Util;
use Getopt::Std;
use Config::Simple;
use File::Path qw(mkpath rmtree);
use File::Slurp;
use File::Copy qw(copy move);
use File::Copy::Recursive qw(fcopy);
use Scalar::Util qw(looks_like_number);
use Data::Dumper;
use XML::Simple;
use DateTime;
use DateTime::Format::Human::Duration;
use File::Basename;
use threads;
use IO::Select;
use IPC::Open3;
use File::Spec;
use File::Find;
sub Init;
sub Usage;
sub CHECK_Config;
sub split_xmap;
sub Version;
sub find_refaligner;
sub log10;
use List::MoreUtils qw(uniq);
use BNG::Utility;

my %opts;
my @cmd = ();	
my $cmdRef = \@cmd;
my ($outResults, $errResults) = ("", "");
my $plDir=abs_path(dirname($0));

#read comand line args
Init();

if(!defined($opts{n}) || !-e $opts{n}){
	 #print("Invalid input NGS file: $opts{n}! \n");
	 print("ERROR: Invalid input NGS file: $opts{n}! \n");
	 Usage(); }
elsif(!defined($opts{b}) || !-e $opts{b}) {
	 #print("Invalid input BNG file: $opts{b}! \n");
	 print("ERROR: Invalid input BNG file: $opts{b}! \n");
	 Usage(); }
elsif(!defined($opts{c}) || !-e $opts{c}){
	 #print("Invalid input configuration file: $opts{c}! \n");
	 print("ERROR: Invalid input configuration file: $opts{c}! \n");
	 Usage(); }
if(!defined($opts{o})){
	#print "No output folder defined!";
	print "ERROR: No output folder defined!";
	Usage(); }
if ($opts{f}) {
	if (-d $opts{o}) {
		rmtree($opts{o}); } }
elsif (-d $opts{o}) {
	dieLog( "ERROR: Output folder: $opts{o} already exists! Use -f option to overwrite.\n"); }
	
eval { mkpath($opts{o}) };
if ($@) {
	print("ERROR: Couldn't create $opts{o}: $@"); 
	print("ERROR: Output folder invalid or does not exist! \n");
	Usage(); }
	
my (@BNGpath) = split('/', $opts{b});
my $n = scalar(@BNGpath) - 1;
my $BNGfile = $BNGpath[$n];
$BNGfile =~ s/\./_/g ; 
my (@NGSpath) = split('/', $opts{n}); 
$n = scalar(@NGSpath) - 1;
my $NGSfile = $NGSpath[$n];
$NGSfile =~ s/\./_/g ;
my $prefix = $BNGfile."_".$NGSfile."_HYBRID_SCAFFOLD";

#my $log_file = $opts{o}."/".$prefix."_log.txt";
my @s = split ('/', $opts{o});
my $log_file_2 = $opts{o}."/".$s[-1]."_log.txt";

#my $log_file_2 = $BNGpath[-1] . "_HYBRID_SCAFFOLD_log.txt";
$prefix = abs_path("$prefix");
for (*STDOUT, *STDERR) {
  #$_->autoflush; $_->push_layer(tee => "$log_file");
  $_->autoflush; $_->push_layer(tee => "$log_file_2");
}
	
#print header
print "**************************\n";
print "*****BioNano Genomics*****\n";
print "******BNG-NGS Merge*******\n";
print "**************************\n\n";

Version();
print "\n";

my $hostname = `hostname`;
print "Running on host: $hostname\n";

my $dtStart = DateTime->now;
print "Start time: "; print join ' ', $dtStart->ymd, $dtStart->hms; print "\n\n";

print qx/ps -o args $$/;

#make sure all input files exist 
print "\n";
if (-e $opts{n} && -e $opts{b} && -e $opts{c}) {
	
	print "NGS file: $opts{n} \n" ;
	print "BNG file: $opts{b} \n";
	print "Configuration file: $opts{c} \n"; 
	print "Output folder: $opts{o}\n"; 
	if (defined($opts{m})) {
		print "Molecules file: $opts{m} \n";}
	}
else {
	print("One or more input files do not exist. \n");
	Usage(); }

#make sure input BNG CMAP file contains at least one contig
my (undef, $numContig, undef) = readCMap($opts{b});
if (!($numContig > 0)) {
	#print "ERROR: Input BNG file $opts{b} does not contain any contigs!\n\n";
	dieLog( "ERROR: Input BNG file $opts{b} does not contain any contigs!\n\n");
	exit;
}

#load config file
my $XML = new XML::Simple(KeyAttr=>[]);
my $config = $XML->XMLin($opts{c});

# parse config XML file
my %config = Parse($config);

#my $cfg = new Config::Simple($opts{c});
#%config = $cfg->vars();
#tie %config, "Config::Simple", $opts{c};

#check config file
CHECK_Config(\%config);

copy "$opts{c}", "$opts{o}" ;

#run fasta to cmap conversion (if needed)
my $enzyme = $config{'fasta2cmap.enzyme'};
my $minLabels = $config{'fasta2cmap.minLabels'};
my $minLength = $config{'fasta2cmap.minLength'};

my @file = split(/\./, $opts{n});
my $filename = $opts{n};
$filename =~ s/(\S+)\.\w+$/$1_$enzyme/;
#my $ngs_cmap = $filename.".cmap";
my $ngs_cmap = $filename."_".$minLength."Kb_".$minLabels."labels.cmap";

if ($file[$#file] =~ m/cmap/i) {
	print "\nInput NGS file seems to be already in CMAP format!\nSkipping FASTA to CMAP conversion.\n";
	$ngs_cmap = $opts{n} ; }
elsif (! -e $ngs_cmap) {
	my $dtStartStage = DateTime->now;
	print "\nBeginning FASTA to CMAP conversion...\n";
	print "Using Enzyme: $enzyme\n";
	print "Minimum Length: $minLength Kb\nMinimum Labels: $minLabels\n";
	#system("$^X fa2cmap.pl -v -i $opts{n} -n $enzyme > $opts{o}/fa2cmap.log");
	
	chdir $plDir or dieLog( "ERROR: Cannot change directory to $plDir: $!\n");	
	@$cmdRef = ($^X, "scripts/fa2cmap.pl", "-v", "-i", $opts{n}, "-n", $enzyme, "-m", $minLabels, "-M", $minLength);
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand ($cmdRef);
	eval { mkpath("$opts{o}/fa2cmap") };
	if ($@) {
	  print "Couldn't create $opts{o}/fa2cmap: $@"; }
	open (OUT, ">$opts{o}/fa2cmap/fa2cmap.log"); print OUT $outResults."\n"; close OUT;
	open (ERR, ">$opts{o}/fa2cmap/fa2cmap.errlog"); print ERR $errResults."\n"; close ERR;
	
	#copy "$ngs_cmap", "$opts{o}" ;
	my $dtEndStage = DateTime->now;
	my $span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);	
	
	# error check
	errCheck($outResults, $errResults, "ERROR:", "FASTA to CMAP conversion complete in $span.", "ERROR: FASTA to CMAP conversion cannot be completed.");
	
	$dtStartStage = DateTime->now;
	print "\nBeginning FASTA header conversion...\n";
	chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
	
	my $filename_key = $opts{n};
	$filename_key =~ s/(\S+)\.\w+$/$1_$enzyme/;
	$filename_key = $filename_key."_".$minLength."Kb_".$minLabels."labels_key.txt";	
	my $filename = $opts{n};$filename =~ s/(\S+)\.\w+$/$1_$enzyme/;
	$filename = $filename."_".$minLength."Kb_".$minLabels."labels.cmap";
	
	@$cmdRef = ($^X, "scripts/fa_key_convert.pl", $opts{n}, $filename_key);
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand ($cmdRef);
	eval { mkpath("$opts{o}/fa2cmap") };
	if ($@) {
	  print "Couldn't create $opts{o}/fa2cmap: $@"; }
	open (OUT, ">$opts{o}/fa2cmap/faHeader_to_cmapId.log"); print OUT $outResults."\n"; close OUT;
	open (ERR, ">$opts{o}/fa2cmap/faHeader_to_cmapId.errlog"); print ERR $errResults."\n"; close ERR;	
	
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);	
	
	# error check
	errCheck($outResults, $errResults, "ERROR:", "FASTA header conversion complete in $span.", "ERROR: FASTA header conversion cannot be completed.");
	my $faPath = "$opts{o}/fa2cmap";
	eval { mkpath($faPath) };
	copy "$filename_key", "$faPath" ; 
	#@s = split(/\./,$opts{n}); 
	#$filename_key = $opts{n};
	$filename =~ s/(\S+)\.\w+$/$1/;
	my $temp = $filename."_CmapIdHeaders.fa"; print "\nNew FASTA with CMAP Ids as headers: $temp\n";
	copy "$temp", "$faPath" ;
	copy "$opts{n}", "$faPath" ;
	}
else {
	print "\nIn-silico digested $ngs_cmap already exists!\nSkipping FASTA to CMAP conversion.\n"; }
$ngs_cmap = abs_path($ngs_cmap);

#make sure input NGS CMAP file contains at least one contig
(undef, $numContig, undef) = readCMap($ngs_cmap);
if (!($numContig > 0)) {
	#print "ERROR: Input NGS file $ngs_cmap does not contain any contigs!\n\n";
	dieLog( "ERROR: Input NGS file $ngs_cmap does not contain any contigs!\n\n");
	exit;
}

#perform initial alignment 
my $dtStartStage = DateTime->now;
print "\nBeginning initial NGS CMAP to BioNano CMAP alignment...\n";
my $outputDir = $opts{o}."/align1";
eval { mkpath($outputDir) };
if ($@) {
  print "Couldn't create $outputDir: $@"; }

my $refaligner = $config{'global.refaligner'};
if ($refaligner =~ /~/) {
	 $refaligner = glob($refaligner); }
else {
	if (-e $refaligner) {
		$refaligner = abs_path($refaligner); } }
#make sure RefAligner exists
if (! -e $refaligner) {
	warn "WARNING: RefAligner binary does not exist at $refaligner as defined in config file. Trying to find it...\n";
	my $script_path = abs_path(dirname($0));
	my @s = split('/',$script_path);
	my $val = pop(@s); $val = pop(@s);
	my $home_path = join('/',@s);
	$refaligner = $home_path."/tools/RefAligner"; 
	if (! -e $refaligner) {	
		$refaligner = glob("~/tools/RefAligner");
		if (! -e $refaligner) {		
			#dieLog ("ERROR: RefAligner binary does not exist at $refaligner\n"); }
			dieLog ("ERROR: RefAligner binary cannot be found!\n"); }
		else {
			print "RefAligner binary found at $refaligner\n"; } }
	else {
			print "RefAligner binary found at $refaligner\n"; } }	  
else {
			print "RefAligner binary found at $refaligner\n"; } 

my $T=$config{'align1.T'};
my $endoutlier=$config{'align1.endoutlier'};
my $outlier=$config{'align1.outlier'};
my $biaswt=$config{'align1.biaswt'};
my $sd=$config{'align1.sd'};
my $res=$config{'align1.res'};
my $sf=$config{'align1.sf'};


#$cmd = "cd '$outputDir'; $refaligner -f -ref $ngs_cmap -i $opts{b} -o align1 -endoutlier $endoutlier -outlier $outlier -extend $extend -FN $FN -FP $FP -sf $sf -sd $sd -sr $sr -res $res -resSD $resSD -mres $mres -A $A -biaswt $biaswt -M $M -Mfast $Mfast -maxmem $maxmem -maxthreads $maxthreads -deltaX $deltaX -deltaY $deltaY -xmapchim $xmapchim -RepeatMask $RepeatMask -RepeatRec $RepeatRec -T $T -stdout -stderr";
chdir $outputDir or dieLog( "ERROR: Cannot change directory to $outputDir: $!\n");
@$cmdRef = ($refaligner, "-f", "-ref", $ngs_cmap, "-i", $opts{b}, "-o", "align1", "-T", $T, "-endoutlier", $endoutlier, "-outlier", $outlier, "-biaswt", $biaswt, "-sd", $sd, "-res", $res, "-sf", $sf,  "-stdout", "-stderr", "-hashgen", split(/\s+/, "5 3 2.4 1.5 0.05 5.0 1 1 1"), "-hash", "-hashdelta", 50);
#########################################################@$cmdRef = ($refaligner, "-f", "-ref", $ngs_cmap, "-i", $opts{b}, "-o", "align0", "-endoutlier", $endoutlier, "-outlier", $outlier, "-extend", $extend, "-FN", $FN, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-sr", $sr, "-res", $res, "-resSD", $resSD, "-mres", $mres, "-A", $A, "-biaswt", $biaswt, "-M", $M, "-Mfast", $Mfast, "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-deltaX", $deltaX, "-deltaY", $deltaY, "-xmapchim", $xmapchim, "-xmapUnique", $xmapUnique,"-RepeatMask", split(/\s+/,$RepeatMask), "-RepeatRec", split(/\s+/, $RepeatRec), "-T", $T, "-stdout", "-stderr");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
my $dtEndStage = DateTime->now;
my $span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheckRefAligner("align1.stdout", "END of output", "Initial alignment complete in $span.", "ERROR: Initial alignment cannot be completed.");

my $align1_xmap = $outputDir."/align1.xmap";
my $align1_r_cmap = $outputDir."/align1_r.cmap";
my $align1_q_cmap = $outputDir."/align1_q.cmap";


#check to make sure that intial alignment XMAP contains alignments
my @int_xmap = read_file($align1_xmap);
my $headLines=0;
my $alignLines=0;
foreach (@int_xmap) {
	if ($_ =~ /^#/ ) {
		$headLines++; }
	else {
		$alignLines++; } }
if ($alignLines < 1) {
	dieLog ("\nERROR: No intial alignments found between $opts{n} and $opts{b}\n"); }
else {
	print "\n$alignLines alignments found between $opts{n} and $opts{b}\n"; }

#run AssignAlignType.R script
$dtStartStage = DateTime->now;
print "\nBeginning AssignAlignType...\n";
$outputDir = $opts{o}."/assignAlignType";
eval { mkpath($outputDir) };
if ($@) {
  print "Couldn't create $outputDir: $@"; }
my $assignAlignType_xmap = $outputDir."/assignAlignType.xmap";
my $assignAlignType_r_cmap = $outputDir."/assignAlignType_r.cmap";
my $assignAlignType_q_cmap = $outputDir."/assignAlignType_q.cmap";
my $T_cutoff=$config{'assignAlignType.T_cutoff'};
$T_cutoff = -log10($T_cutoff);
my $max_overhang=$config{'assignAlignType.max_overhang'};

chdir $plDir or dieLog( "ERROR: Cannot change directory to $plDir: $!\n");
#@$cmdRef = ("Rscript", "./scripts/AssignAlignType.R", $align1_xmap, $align1_r_cmap, $align1_q_cmap, $assignAlignType_xmap, $assignAlignType_r_cmap, $assignAlignType_q_cmap, $T_cutoff, $max_overhang, $ngs_cmap, $bppAdjust;
#*****@$cmdRef = ("$^X", "./scripts/AssignAlignType.pl", $align1_xmap, $align1_r_cmap, $align1_q_cmap, $assignAlignType_xmap, $assignAlignType_r_cmap, $assignAlignType_q_cmap, $T_cutoff, $max_overhang, $ngs_cmap, $bppAdjust);
@$cmdRef = ("$^X", "./scripts/AssignAlignType.pl", $align1_xmap, $align1_r_cmap, $align1_q_cmap, $assignAlignType_xmap, $assignAlignType_r_cmap, $assignAlignType_q_cmap, $T_cutoff, $max_overhang, $ngs_cmap, $opts{b});

print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
open(OUT, ">$outputDir/assignAlignType.log") or dieLog ("ERROR: Cannot write to $outputDir/assignAlignType.log: $!\n");	print OUT "$outResults";	close OUT;
open(ERR, ">$outputDir/assignAlignType.errlog") or dieLog ("ERROR: Cannot write to $outputDir/assignAlignType.errlog: $!\n");	print ERR "$errResults";	close ERR;
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheck($outResults, $errResults, "ERROR:", "AssignAlignType complete in $span.", "ERROR: AssignAlignType cannot be completed.");
#$cmd = "cd $plDir; Rscript ./scripts/AssignAlignType.R $align1_xmap $align1_r_cmap $align1_q_cmap $assignAlignType_xmap $assignAlignType_r_cmap $assignAlignType_q_cmap $T_cutoff $max_overhang $ngs_cmap $opts{b} &> $outputDir/assignAlignType.log";
#print "Running command: $cmd\n";
#system($cmd) == 0 
#	or die "system $cmd failed: $?";

#print number of BioNano and NGS contigs that have been flagged as conflicting
my @bn; my @ngs;
my @sticky = read_file($assignAlignType_xmap);
foreach (@sticky) {
	my @s = split("\t", $_);
	push (@bn, $s[1]);
	push (@ngs, $s[2]);
}
@bn = uniq @bn;
@ngs = uniq @ngs;
print scalar(@bn)." BNG contigs have been flagged as conflicting and were filtered out\n";
print scalar(@ngs)." NGS contigs have been flagged as conflicting and were filtered out\n";

#check AssignAlign outputs and Merge inputs contain at least one contig
my $assignAlign_r = abs_path("$opts{o}/assignAlignType/assignAlignType_r.cmap");
(undef, $numContig, undef) = readCMap($assignAlign_r);
if (!($numContig > 0)) {;
	print "WARNING: AssignAlign output file $assignAlign_r does not contain any contigs! All NGS contigs have flagged as conflicting.\n";
	warn "WARNING: AssignAlign output file $assignAlign_r does not contain any contigs! All NGS contigs have flagged as conflicting.\n";
}
my $assignAlign_q = abs_path("$opts{o}/assignAlignType/assignAlignType_q.cmap");
(undef, $numContig, undef) = readCMap($assignAlign_q);
if (!($numContig > 0)) {
	print "WARNING: AssignAlign output file $assignAlign_q does not contain any contigs! All BNG contigs have flagged as conflicting.\n";
	warn "WARNING: AssignAlign output file $assignAlign_q does not contain any contigs! All BNG contigs have flagged as conflicting.\n";
}

#run MergeNGS_BN.R
$dtStartStage = DateTime->now;
print "\nBeginning MergeNGS_BN...\n";
$outputDir = $opts{o}."/mergeNGS_BN";
eval { mkpath($outputDir) };
if ($@) {
  print "Couldn't create $outputDir: $@"; }

my $merge_Tvalue=$config{'mergeNGS_BN.merge_Tvalue'};
my $id_shift=$config{'mergeNGS_BN.id_shift'};
my $max_merge_rounds=$config{'mergeNGS_BN.max_merge_rounds'};
my $maxthreads=$config{'global.maxthreads'};
$endoutlier=$config{'mergeNGS_BN.endoutlier'};
$outlier=$config{'mergeNGS_BN.outlier'};
$biaswt=$config{'mergeNGS_BN.biaswt'};
$sd=$config{'mergeNGS_BN.sd'};
$res=$config{'mergeNGS_BN.res'};
$sf=$config{'mergeNGS_BN.sf'};
my $RepeatMask=$config{'mergeNGS_BN.RepeatMask'};
my $RepeatRec=$config{'mergeNGS_BN.RepeatRec'};
############################################################my $pairmerge=$config{'mergeNGS_BN.pairmerge'};
#*****
my $readparameters="$opts{o}/align1/align1.errbin";

#if using unfiltered BioNano CMAP and/or unfiltered NGS CMAP
if ($opts{B} && !$opts{N}) {
	#make sure AssignAlign outputs and Merge inputs contain at least one contig
	(undef, $numContig, undef) = readCMap($assignAlign_r);
	if (!($numContig > 0)) {
		#print "ERROR: AssignAlign output $assignAlign_r does not contain any contigs! All NGS contigs have flagged as conflicting. You can try using unfiltered NGS option.\n\n";
		dieLog ("ERROR: AssignAlign output $assignAlign_r does not contain any contigs! All NGS contigs have flagged as conflicting. You can try using unfiltered NGS option.\n\n");
		exit;
	}

	chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
#	@$cmdRef = ("Rscript", "./scripts/MergeNGS_BN.R", "$opts{o}/mergeNGS_BN", $refaligner, $merge_Tvalue, "./", "mergeNGS_BN.log", "../assignAlignType/assignAlignType_r.cmap", $opts{b}, $id_shift, $max_merge_rounds, $maxthreads, $endoutlier, $outlier, $biaswt, $sd, $res, $sf, split(/\s+/, $RepeatMask), split(/\s+/, $RepeatRec));
	@$cmdRef = ("$^X", "./scripts/MergeNGS_BN.pl", "$opts{o}/mergeNGS_BN", $refaligner, $merge_Tvalue, "./", "mergeNGS_BN_script.log", "../assignAlignType/assignAlignType_r.cmap", $opts{b}, $id_shift, $max_merge_rounds, $maxthreads, $endoutlier, $outlier, $biaswt, $sd, $res, $sf, split(/\s+/, $RepeatMask), split(/\s+/, $RepeatRec)); #######################split(/\s+/, $pairmerge), $readparameters);

	#$cmd = "cd $plDir; Rscript ./scripts/MergeNGS_BN.R $opts{o}/mergeNGS_BN $refaligner $merge_Tvalue ./ mergeNGS_BN.log ../assignAlignType/assignAlignType_r.cmap $opts{b} $id_shift $max_merge_rounds $maxthreads $endoutlier $outlier $biaswt $sd $res $sf $RepeatMask $RepeatRec";
	print "*Using unfiltered BioNano CMAP and filtered NGS CMAP*\n"; }
elsif (!$opts{B} && !$opts{N}) {
	#make sure AssignAlign outputs and Merge inputs contain at least one contig
	(undef, $numContig, undef) = readCMap($assignAlign_r);
	if (!($numContig > 0)) {
		#print "ERROR: AssignAlign output $assignAlign_r does not contain any contigs! All NGS contigs have flagged as conflicting. You can try using unfiltered NGS option.\n\n";
		dieLog( "ERROR: AssignAlign output $assignAlign_r does not contain any contigs! All NGS contigs have flagged as conflicting. You can try using unfiltered NGS option.\n\n");
		exit;
	}
	(undef, $numContig, undef) = readCMap($assignAlign_q);
	if (!($numContig > 0)) {
		#print "ERROR: AssignAlign output $assignAlign_q does not contain any contigs! All BNG contigs have flagged as conflicting. You can try using unfiltered BNG option.\n\n";
		dieLog ("ERROR: AssignAlign output $assignAlign_q does not contain any contigs! All BNG contigs have flagged as conflicting. You can try using unfiltered BNG option.\n\n");
		exit;
	}
	
	chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
    # @$cmdRef = ("Rscript", "./scripts/MergeNGS_BN.R", "$opts{o}/mergeNGS_BN", $refaligner, $merge_Tvalue, "./", "mergeNGS_BN.log", "../assignAlignType/assignAlignType_r.cmap", "../assignAlignType/assignAlignType_q.cmap", $id_shift, $max_merge_rounds, $maxthreads, $endoutlier, $outlier, $biaswt, $sd, $res, $sf, split(/\s+/, $RepeatMask), split(/\s+/, $RepeatRec));
	@$cmdRef = ("$^X", "./scripts/MergeNGS_BN.pl", "$opts{o}/mergeNGS_BN", $refaligner, $merge_Tvalue, "./", "mergeNGS_BN_script.log", "../assignAlignType/assignAlignType_r.cmap", "../assignAlignType/assignAlignType_q.cmap", $id_shift, $max_merge_rounds, $maxthreads, $endoutlier, $outlier, $biaswt, $sd, $res, $sf, split(/\s+/, $RepeatMask), split(/\s+/, $RepeatRec)); ###################split(/\s+/, $pairmerge), $readparameters);

	#$cmd = "cd $plDir; Rscript ./scripts/MergeNGS_BN.R $opts{o}/mergeNGS_BN $refaligner $merge_Tvalue ./ mergeNGS_BN.log ../assignAlignType/assignAlignType_r.cmap ../assignAlignType/assignAlignType_q.cmap $id_shift $max_merge_rounds $maxthreads $endoutlier $outlier $biaswt $sd $res $sf $RepeatMask $RepeatRec";
	print "*Using filtered BioNano CMAP and filtered NGS CMAP*\n"; }
elsif (!$opts{B} && $opts{N}) {
	#make sure AssignAlign outputs and Merge inputs contain at least one contig
	(undef, $numContig, undef) = readCMap($assignAlign_q);
	if (!($numContig > 0)) {
		#print "ERROR: AssignAlign output .$assignAlign_q does not contain any contigs! All BNG contigs have flagged as conflicting. You can try using unfiltered BNG option.\n\n";
		dieLog( "ERROR: AssignAlign output .$assignAlign_q does not contain any contigs! All BNG contigs have flagged as conflicting. You can try using unfiltered BNG option.\n\n");
		exit;
	}
	
	chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
	# @$cmdRef = ("Rscript", "./scripts/MergeNGS_BN.R", "$opts{o}/mergeNGS_BN", $refaligner, $merge_Tvalue, "./", "mergeNGS_BN.log", $ngs_cmap, "../assignAlignType/assignAlignType_q.cmap", $id_shift, $max_merge_rounds, $maxthreads, $endoutlier, $outlier, $biaswt, $sd, $res, $sf, split(/\s+/, $RepeatMask), split(/\s+/, $RepeatRec));
	@$cmdRef = ("$^X", "./scripts/MergeNGS_BN.pl", "$opts{o}/mergeNGS_BN", $refaligner, $merge_Tvalue, "./", "mergeNGS_BN_script.log", $ngs_cmap, "../assignAlignType/assignAlignType_q.cmap", $id_shift, $max_merge_rounds, $maxthreads, $endoutlier, $outlier, $biaswt, $sd, $res, $sf, split(/\s+/, $RepeatMask), split(/\s+/, $RepeatRec)); #############################split(/\s+/, $pairmerge), $readparameters);

	#$cmd = "cd $plDir; Rscript ./scripts/MergeNGS_BN.R $opts{o}/mergeNGS_BN $refaligner $merge_Tvalue ./ mergeNGS_BN.log $opts{n} ../assignAlignType/assignAlignType_q.cmap $id_shift $max_merge_rounds $maxthreads $endoutlier $outlier $biaswt $sd $res $sf $RepeatMask $RepeatRec"; 
	print "*Using filtered BioNano CMAP and unfiltered NGS CMAP*\n";}
elsif ($opts{B} && $opts{N}) {
	chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
	# @$cmdRef = ("Rscript", "./scripts/MergeNGS_BN.R", "$opts{o}/mergeNGS_BN", $refaligner, $merge_Tvalue, "./", "mergeNGS_BN.log", $ngs_cmap, $opts{b}, $id_shift, $max_merge_rounds, $maxthreads, $endoutlier, $outlier, $biaswt, $sd, $res, $sf, split(/\s+/, $RepeatMask), split(/\s+/, $RepeatRec));
	@$cmdRef = ("$^X", "./scripts/MergeNGS_BN.pl", "$opts{o}/mergeNGS_BN", $refaligner, $merge_Tvalue, "./", "mergeNGS_BN_script.log", $ngs_cmap, $opts{b}, $id_shift, $max_merge_rounds, $maxthreads, $endoutlier, $outlier, $biaswt, $sd, $res, $sf, split(/\s+/, $RepeatMask), split(/\s+/, $RepeatRec)); #####################################split(/\s+/, $pairmerge), $readparameters);

	#$cmd = "cd $plDir; Rscript ./scripts/MergeNGS_BN.R $opts{o}/mergeNGS_BN $refaligner $merge_Tvalue ./ mergeNGS_BN.log $opts{n} $opts{b} $id_shift $max_merge_rounds $maxthreads $endoutlier $outlier $biaswt $sd $res $sf $RepeatMask $RepeatRec"; 
	print "*Using unfiltered BioNano CMAP and unfiltered NGS CMAP*\n";}

print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
open(OUT, ">$opts{o}/mergeNGS_BN/MergeNGS_BN.log") or dieLog ("ERROR: Cannot write to $opts{o}/mergeNGS_BN/MergeNGS_BN.log: $!\n");	print OUT "$outResults";	close OUT;
open(ERR, ">$opts{o}/mergeNGS_BN/MergeNGS_BN.errlog") or dieLog ("ERROR: Cannot write to $opts{o}/mergeNGS_BN/MergeNGS_BN.errlog: $!\n");	print ERR "$errResults";	close ERR;
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheck($outResults, $errResults, "ERROR:", "MergeNGS_BN complete in $span.", "ERROR: MergeNGS_BN cannot be completed");
#system($cmd) == 0 
#	or die "system $cmd failed: $?";

# error check


#Align original NGS contigs to Hybrid map 

#find NGS contigs that went into hybrid assembly and make new NGS cmap
$dtStartStage = DateTime->now;
$outputDir = $opts{o}."/mergeNGS_BN";
print "\nBeginning extraction of NGS contigs in hybrid assembly...\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
@$cmdRef = ($^X, "scripts/findNGS.pl", $outputDir, "ngsContigs.txt");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
open(OUT, ">$outputDir/findNGS.log") or dieLog ("ERROR: Cannot write to $outputDir/findNGS.log: $!\n");	print OUT $outResults;	close OUT;
open(ERR, ">$outputDir/findNGS.errlog") or dieLog ("ERROR: Cannot write to $outputDir/findNGS.errlog: $!\n");	print ERR $errResults;	close ERR;
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheck($outResults, $errResults, "ERROR:", "Extraction of NGS contigs complete in $span.", "ERROR: Extraction of NGS contigs cannot be completed.");

#my $results = capture("cd $plDir; $^X scripts/findNGS.pl $outputDir ngsContigs.txt");
#open (OUT, ">$outputDir/findNGS.log"); print OUT $results; close OUT;
my $ngsContigsFile = "$outputDir/ngsContigs.txt";

$dtStartStage = DateTime->now;
print "\nBeginning generation of filtered NGS cmap...\n";
my @ngs_cmap = read_file($ngs_cmap);
my @ngs_contigs_orig = read_file($ngsContigsFile);
my @ngs_contigs;
foreach my $line (@ngs_contigs_orig) {
	chomp($line);
	push (@ngs_contigs, $line); }
my $final_outputDir = $opts{o}."/align_final";
eval { mkpath($final_outputDir) };
if ($@) {
  print "Couldn't create $final_outputDir: $@"; }
open (OUT, ">$final_outputDir/filtered_NGS.cmap");
copy "$outputDir/step2.hybrid.cmap", "$final_outputDir/step2.hybrid.cmap" ;

my %h = map {$_ => 1 } @ngs_contigs;
foreach (@ngs_cmap) {
	my $l = $_;
	if ($_ =~ /^#/ ) {
		print OUT $_; }
	else {
		my ($CMapId, $ContigLength,	$NumSites, $SiteID, $LabelChannel, $Position, $StdDev, $Coverage, $Occurrence) = split("\t",$_); 	
		if (exists $h{$CMapId}) {
			print OUT $l; } } }
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
print "Generation of filtered NGS cmap complete in $span.\n";	

#align just the NGS contigs in hybrid back to hybrid assembly
$dtStartStage = DateTime->now;
print "\nBeginning alignment of filtered NGS cmap to Hybrid CMAP...\n";

(@BNGpath) = split('/', $opts{b});
$n = scalar(@BNGpath) - 1;
$BNGfile = $BNGpath[$n];
$BNGfile =~ s/\./_/g ; 
(@NGSpath) = split('/', $opts{n}); 
$n = scalar(@NGSpath) - 1;
$NGSfile = $NGSpath[$n];
$NGSfile =~ s/\./_/g ;
my $prefixOrig = $BNGfile."_".$NGSfile."_HYBRID_SCAFFOLD";
$prefix = $BNGfile."_".$NGSfile."_NGScontigs_HYBRID_SCAFFOLD";

my $filteredNGS_file = "$final_outputDir/filtered_NGS.cmap";
my $hybrid_cmap = "$final_outputDir/step2.hybrid.cmap";

#$refaligner = $config{'global.refaligner'};
$endoutlier=$config{'align_final.endoutlier'};
$outlier=$config{'align_final.outlier'};
my $extend=$config{'align_final.extend'};
my $FN=$config{'align_final.FN'};
my $FP=$config{'align_final.FP'};
$sf=$config{'align_final.sf'};
$sd=$config{'align_final.sd'};
my $sr=$config{'align_final.sr'};
$res=$config{'align_final.res'};
my $resSD=$config{'align_final.resSD'};
my $mres=$config{'align_final.mres'};
my $A=$config{'align_final.A'};
$biaswt=$config{'align_final.biaswt'};
my $M=$config{'align_final.M'};
my $Mfast=$config{'align_final.Mfast'};
my $maxmem=$config{'global.maxmem'};
$maxthreads=$config{'global.maxthreads'};
my $deltaX=$config{'align_final.deltaX'};
my $deltaY=$config{'align_final.deltaY'};
$RepeatMask=$config{'align_final.RepeatMask'}; 
$RepeatRec=$config{'align_final.RepeatRec'}; 
$T=$config{'align_final.T'};
my $BestRef=$config{'align_final.BestRef'};
my $nosplit=$config{'align_final.nosplit'};

chdir $final_outputDir or dieLog ("ERROR: Cannot change to $final_outputDir: $!\n");
@$cmdRef = ($refaligner, "-f", "-ref", $hybrid_cmap, "-i", $filteredNGS_file, "-o", $prefix, "-endoutlier", $endoutlier, "-outlier", $outlier, "-extend", $extend, "-FN", $FN, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-sr", $sr, "-res", $res, "-resSD", $resSD, "-mres", $mres, "-A", $A, "-biaswt", $biaswt, "-M", $M, "-Mfast", $Mfast, "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-deltaX", $deltaX, "-deltaY", $deltaY, "-RepeatMask", split(/\s+/, $RepeatMask), "-RepeatRec", split(/\s+/, $RepeatRec), "-T", $T, "-BestRef", $BestRef, "-nosplit", $nosplit, "-XmapStatWrite", "$prefix.stats", "-stdout", "-stderr", "-hashgen", split(/\s+/, "5 3 2.4 1.5 0.05 5.0 1 1 1"), "-hash", "-hashdelta", 50);
#####################################@$cmdRef = ($refaligner, "-f", "-ref", $hybrid_cmap, "-i", $filteredNGS_file, "-o", $prefix, "-endoutlier", $endoutlier, "-outlier", $outlier, "-extend", $extend, "-FN", $FN, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-sr", $sr, "-res", $res, "-resSD", $resSD, "-mres", $mres, "-A", $A, "-biaswt", $biaswt, "-M", $M, "-Mfast", $Mfast, "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-deltaX", $deltaX, "-deltaY", $deltaY, "-RepeatMask", split(/\s+/, $RepeatMask), "-RepeatRec", split(/\s+/, $RepeatRec), "-T", $T, "-BestRef", $BestRef, "-nosplit", $nosplit, "-XmapStatWrite", "$prefix.stats", "-stdout", "-stderr");
print "Running command: ", (join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
errCheckRefAligner("$prefix.stdout", "END of output", "Alignment of filtered NGS CMAP to Hybrid CMAP complete in $span.", "ERROR: alignment of filtered NGS CMAP to Hybrid CMAP cannot be completed.");

#$cmd = "cd $final_outputDir; $refaligner -f -ref $hybrid_cmap -i $filteredNGS_file -o $prefix -endoutlier $endoutlier -outlier $outlier -extend $extend -FN $FN -FP $FP -sf $sf -sd $sd -sr $sr -res $res -resSD $resSD -mres $mres -A $A -biaswt $biaswt -M $M -Mfast $Mfast -maxmem $maxmem -maxthreads $maxthreads -deltaX $deltaX -deltaY $deltaY -RepeatMask $RepeatMask -RepeatRec $RepeatRec -T $T -BestRef $BestRef -nosplit $nosplit -XmapStatWrite $prefix.stats -stdout -stderr";
#print "Running command: $cmd\n";
#system($cmd) == 0 
#	or die "system $cmd failed: $?";

#Align original BNG contigs to Hybrid map 

#find BNG contigs that went into hybrid assembly and make new BNG cmap
$dtStartStage = DateTime->now;
$outputDir = $opts{o}."/mergeNGS_BN";
print "\nBeginning extraction of BNG contigs in hybrid assembly...\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
@$cmdRef = ($^X, "scripts/findBNG.pl", $outputDir, "bngContigs.txt");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
open(OUT, ">$outputDir/findBNG.log") or dieLog ("ERROR: Cannot write to $outputDir/findBNG.log: $!\n");	print OUT $outResults;	close OUT;
open(ERR, ">$outputDir/findBNG.errlog") or dieLog ("ERROR: Cannot write to $outputDir/findBNG.errlog: $!\n");	print ERR $errResults;	close ERR;
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheck($outResults, $errResults, "ERROR:", "Extraction of BNG contigs complete in $span.", "ERROR: Extraction of BNG contigs cannot be completed.");

my $bngContigsFile = "$outputDir/bngContigs.txt";

$dtStartStage = DateTime->now;
print "\nBeginning generation of filtered BNG cmap...\n";
my @bng_cmap = read_file($opts{b});
my @bng_contigs_orig = read_file($bngContigsFile);
my @bng_contigs;
foreach my $line (@bng_contigs_orig) {
	chomp($line);
	push (@bng_contigs, $line); }
$final_outputDir = $opts{o}."/align_final";
eval { mkpath($final_outputDir) };
if ($@) {
  print "Couldn't create $final_outputDir: $@"; }
open (OUT, ">$final_outputDir/filtered_BNG.cmap");
copy "$outputDir/step2.hybrid.cmap", "$final_outputDir/step2.hybrid.cmap" ;

%h = map {$_ => 1 } @bng_contigs;
foreach (@bng_cmap) {
	my $l = $_;
	if ($_ =~ /^#/ ) {
		print OUT $_; }
	else {
		my ($CMapId, $ContigLength,	$NumSites, $SiteID, $LabelChannel, $Position, $StdDev, $Coverage, $Occurrence) = split("\t",$_); 	
		if (exists $h{$CMapId}) {
			print OUT $l; } } }
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
print "Generation of filtered BNG cmap complete in $span.\n";	

#align just the BNG contigs in hybrid back to hybrid assembly
$dtStartStage = DateTime->now;
print "\nBeginning alignment of filtered BNG cmap to Hybrid CMAP...\n";

(@BNGpath) = split('/', $opts{b});
$n = scalar(@BNGpath) - 1;
$BNGfile = $BNGpath[$n];
$BNGfile =~ s/\./_/g ; 
(@NGSpath) = split('/', $opts{n}); 
$n = scalar(@NGSpath) - 1;
$NGSfile = $NGSpath[$n];
$NGSfile =~ s/\./_/g ;
my $prefixBNG = $BNGfile."_".$NGSfile."_BNGcontigs_HYBRID_SCAFFOLD";

my $filteredBNG_file = "$final_outputDir/filtered_BNG.cmap";
$hybrid_cmap = "$final_outputDir/step2.hybrid.cmap";

#$refaligner = $config{'global.refaligner'};
$endoutlier=$config{'align_final.endoutlier'};
$outlier=$config{'align_final.outlier'};
$extend=$config{'align_final.extend'};
$FN=$config{'align_final.FN'};
$FP=$config{'align_final.FP'};
$sf=$config{'align_final.sf'};
$sd=$config{'align_final.sd'};
$sr=$config{'align_final.sr'};
$res=$config{'align_final.res'};
$resSD=$config{'align_final.resSD'};
$mres=$config{'align_final.mres'};
$A=$config{'align_final.A'};
$biaswt=$config{'align_final.biaswt'};
$M=$config{'align_final.M'};
$Mfast=$config{'align_final.Mfast'};
$maxmem=$config{'global.maxmem'};
$maxthreads=$config{'global.maxthreads'};
$deltaX=$config{'align_final.deltaX'};
$deltaY=$config{'align_final.deltaY'};
$RepeatMask=$config{'align_final.RepeatMask'}; 
$RepeatRec=$config{'align_final.RepeatRec'}; 
$T=$config{'align_final.T'};
$BestRef=$config{'align_final.BestRef'};
$nosplit=$config{'align_final.nosplit'};

chdir $final_outputDir or dieLog ("ERROR: Cannot change to $final_outputDir: $!\n");
@$cmdRef = ($refaligner, "-f", "-ref", $hybrid_cmap, "-i", $filteredBNG_file, "-o", $prefixBNG, "-endoutlier", $endoutlier, "-outlier", $outlier, "-extend", $extend, "-FN", $FN, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-sr", $sr, "-res", $res, "-resSD", $resSD, "-mres", $mres, "-A", $A, "-biaswt", $biaswt, "-M", $M, "-Mfast", $Mfast, "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-deltaX", $deltaX, "-deltaY", $deltaY, "-RepeatMask", split(/\s+/, $RepeatMask), "-RepeatRec", split(/\s+/, $RepeatRec), "-T", $T, "-BestRef", $BestRef, "-nosplit", $nosplit, "-XmapStatWrite", "$prefix.stats", "-stdout", "-stderr", "-hashgen", split(/\s+/, "5 3 2.4 1.5 0.05 5.0 1 1 1"), "-hash", "-hashdelta", 50);
####################################@$cmdRef = ($refaligner, "-f", "-ref", $hybrid_cmap, "-i", $filteredBNG_file, "-o", $prefixBNG, "-endoutlier", $endoutlier, "-outlier", $outlier, "-extend", $extend, "-FN", $FN, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-sr", $sr, "-res", $res, "-resSD", $resSD, "-mres", $mres, "-A", $A, "-biaswt", $biaswt, "-M", $M, "-Mfast", $Mfast, "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-deltaX", $deltaX, "-deltaY", $deltaY, "-RepeatMask", split(/\s+/, $RepeatMask), "-RepeatRec", split(/\s+/, $RepeatRec), "-T", $T, "-BestRef", $BestRef, "-nosplit", $nosplit, "-XmapStatWrite", "$prefix.stats", "-stdout", "-stderr");
print "Running command: ", (join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
errCheckRefAligner("$prefixBNG.stdout", "END of output", "Alignment of filtered BNG CMAP to Hybrid CMAP complete in $span.", "ERROR: alignment of filtered BNG CMAP to Hybrid CMAP cannot be completed.");









#merge hybrid cmap with step1.NGS.naive.cmap
$dtStartStage = DateTime->now;
print "\nMerging Hybrid CMAP with naive NGS CMAP...\n";
my $NGS_naive_cmap = "$outputDir/step1.NGS.naive.cmap";
#my $BNG_naive_cmap = "$outputDir/step1.BN.naive.cmap";
my $hybrid_naiveNGS = "HYBRID_SCAFFOLD_naiveNGS_merged.cmap";
#my $hybrid_naiveBNG = "HYBRID_SCAFFOLD_naiveBNG_merged.cmap";
chdir $final_outputDir or dieLog ("ERROR: Cannot change to $final_outputDir: $!\n");
@$cmdRef = ($refaligner, "-f", "-i", $hybrid_cmap, "-i", $NGS_naive_cmap, "-o", "HYBRID_SCAFFOLD_naiveNGS_merged", "-merge", "-stdout", "-stderr");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheckRefAligner("HYBRID_SCAFFOLD_naiveNGS_merged.stdout", "END of output", "Merging Hybrid CMAP with naive NGS CMAP complete in $span.", "ERROR: Merging Hybrid CMAP with naive NGS CMAP cannot be completed.");

#merge hybrid cmap with step1.BN.naive.cmap
$dtStartStage = DateTime->now;
print "\nMerging Hybrid CMAP with naive BNG CMAP...\n";
#my $NGS_naive_cmap = "$outputDir/step1.NGS.naive.cmap";
my $BNG_naive_cmap = "$outputDir/step1.BN.naive.cmap";
#my $hybrid_naiveNGS = "HYBRID_SCAFFOLD_naiveNGS_merged.cmap";
my $hybrid_naiveBNG = "HYBRID_SCAFFOLD_naiveBNG_merged.cmap";
chdir $final_outputDir or dieLog ("ERROR: Cannot change to $final_outputDir: $!\n");
@$cmdRef = ($refaligner, "-f", "-i", $hybrid_cmap, "-i", $BNG_naive_cmap, "-o", "HYBRID_SCAFFOLD_naiveBNG_merged", "-merge", "-stdout", "-stderr");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
# error check
errCheckRefAligner("HYBRID_SCAFFOLD_naiveBNG_merged.stdout", "END of output", "Merging Hybrid CMAP with naive BNG CMAP complete in $span.", "ERROR: Merging Hybrid CMAP with naive BNG CMAP cannot be completed.");





#$cmd = "cd $final_outputDir; $refaligner -f -i $hybrid_cmap -i $NGS_naive_cmap -o HYBRID_SCAFFOLD_naiveNGS_merged -merge -stdout -stderr";
#print "Running command: $cmd\n";
#system($cmd) == 0 
#	or die "system $cmd failed: $?";
$hybrid_naiveNGS = abs_path("$final_outputDir/$hybrid_naiveNGS");
$hybrid_naiveBNG = abs_path("$final_outputDir/$hybrid_naiveBNG");
#print "Merging Hybrid CMAP with naive NGS CMAP complete.\n";

##align just the NGS contigs in hybrid back to hybrid + naiveNGS assembly
#$dtStartStage = DateTime->now;
#print "\nBeginning alignment of filtered NGS cmap to Hybrid + naive NGS cmap...\n";

#my $prefix2 = "HYBRID_SCAFFOLD_naiveNGS_merged";
#chdir $final_outputDir or dieLog ("ERROR: Cannot change to $final_outputDir: $!\n");	
#@$cmdRef = ($refaligner, "-f", "-ref", $hybrid_naiveNGS, "-i", $filteredNGS_file, "-o", $prefix2, "-endoutlier", $endoutlier, "-outlier", $outlier, "-extend", $extend, "-FN", $FN, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-sr", $sr, "-res", $res, "-resSD", $resSD, "-mres", $mres, "-A", $A, "-biaswt", $biaswt, "-M", $M, "-Mfast", $Mfast, "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-deltaX", $deltaX, "-deltaY", $deltaY, "-RepeatMask", split(/\s+/, $RepeatMask), "-RepeatRec", split(/\s+/, $RepeatRec), "-T", $T, "-BestRef", $BestRef, "-nosplit", $nosplit, "-XmapStatWrite", "$prefix2.stats", "-stdout", "-stderr");
#print "Running command: ".(join(" ", @$cmdRef))."\n";
#($outResults, $errResults) = runCommand($cmdRef);
#$dtEndStage = DateTime->now;
#$span = DateTime::Format::Human::Duration->new();
#$span = $span->format_duration_between($dtEndStage, $dtStartStage);
## error check
#errCheckRefAligner("$prefix2.stdout", "END of output", "Alignment of filtered NGS CMAP to Hybrid + naive NGS CMAP complete to $span.", "ERROR: Alignment of filtered NGS CMAP to Hybrid + naive NGS CMAP cannot be completed");
##$cmd = "cd $final_outputDir; $refaligner -f -ref $hybrid_naiveNGS -i $filteredNGS_file -o $prefix2 -endoutlier $endoutlier -outlier $outlier -extend $extend -FN $FN -FP $FP -sf $sf -sd $sd -sr $sr -res $res -resSD $resSD -mres $mres -A $A -biaswt $biaswt -M $M -Mfast $Mfast -maxmem $maxmem -maxthreads $maxthreads -deltaX $deltaX -deltaY $deltaY -RepeatMask $RepeatMask -RepeatRec $RepeatRec -T $T -BestRef $BestRef -nosplit $nosplit -XmapStatWrite $prefix2.stats -stdout -stderr";
##print "Running command: $cmd\n";
##system($cmd) == 0 
##	or die "system $cmd failed: $?";

#align molecules BNX to Hybrid CMAP if BNX file provided
my $prefix3;
my $prefix4;
my $mqrDir;
if(defined($opts{m})) {
	$dtStartStage = DateTime->now;
	print "\nBeginning alignment of molecules BNX to Hybrid CMAP...\n";
	$endoutlier=$config{'mqr.endoutlier'};
	$outlier=$config{'mqr.outlier'};
	$FP=$config{'mqr.FP'};
	$sf=$config{'mqr.sf'};
	$sd=$config{'mqr.sd'};
	$sr=$config{'mqr.sr'};
	$A=$config{'mqr.A'};
	$biaswt=$config{'mqr.biaswt'};
	$M=$config{'mqr.M'};
	$Mfast=$config{'mqr.Mfast'};
	$maxmem=$config{'global.maxmem'};
	$maxthreads=$config{'global.maxthreads'};
	$T=$config{'mqr.T'};
	$nosplit=$config{'mqr.nosplit'};
	$BestRef=$config{'mqr.BestRef'};
	my $resbias=$config{'mqr.resbias'};
	my $minlen=$config{'mqr.minlen'};
	my $S=$config{'mqr.S'};
	my $randomize=$config{'mqr.randomize'};
	my $subset=$config{'mqr.subset'};
	my $hash = $config{'mqr.hash'};
	#$prefix3 = $prefixOrig."_MoleculeQualityReport";
	my @s = split(m/\//, $prefixOrig);
	$prefixOrig = pop(@s);
	$prefix3 = $prefixOrig."_MoleculeQualityReport";
	$mqrDir = $opts{o}."/mqr";
	
	@$cmdRef=();
	
	eval { mkpath($mqrDir) };
	if ($@) {
		print "Couldn't create $mqrDir: $@"; }
	if ((defined($randomize)) && (defined($subset)) && ($randomize ne '0') && ($subset ne '0') && (looks_like_number($randomize)) ) {
		chdir $mqrDir or dieLog ("ERROR: Cannot change to $mqrDir: $!\n");
		@$cmdRef = ($refaligner, "-f", "-ref", $hybrid_cmap, "-i", $opts{m}, "-o", $prefix3, "-nosplit", $nosplit, "-BestRef", $BestRef, "-biaswt", $biaswt, "-Mfast", $Mfast, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-A", $A, "-outlier", $outlier, "-endoutlier", $endoutlier, "-S", $S, "-sr", $sr, "-resbias", split(/\s+/, $resbias), "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-M", $M, "-minlen", $minlen, "-T", $T, "-XmapStatWrite", "$prefix3.stats", "-randomize", $randomize, "-subset", split(/\s+/, $subset), split(/\s+/, $hash), "-stdout", "-stderr");
		#$cmd = "cd $mqrDir; $refaligner -f -ref $hybrid_cmap -i $opts{m} -o $prefix3 -nosplit $nosplit -BestRef $BestRef -biaswt $biaswt -Mfast $Mfast -FP $FP -sf $sf -sd $sd -A $A -outlier $outlier -endoutlier $endoutlier -S $S -sr $sr -resbias $resbias -maxmem $maxmem -maxthreads $maxthreads -M $M -minlen $minlen -T $T -XmapStatWrite $prefix3.stats -randomize $randomize -subset $subset -stdout -stderr"; }
	} else {
		chdir $mqrDir or dieLog ("ERROR: Cannot change to $mqrDir: $!\n");	
		@$cmdRef = ($refaligner, "-f", "-ref", $hybrid_cmap, "-i", $opts{m}, "-o", $prefix3, "-nosplit", $nosplit, "-BestRef", $BestRef, "-biaswt", $biaswt, "-Mfast", $Mfast, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-A", $A, "-outlier", $outlier, "-endoutlier", $endoutlier, "-S", $S, "-sr", $sr, "-resbias", split(/\s+/, $resbias), "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-M", $M, "-minlen", $minlen, "-T", $T, "-XmapStatWrite", "$prefix3.stats", split(/\s+/, $hash), "-stdout", "-stderr");
		#$cmd = "cd $mqrDir; $refaligner -f -ref $hybrid_cmap -i $opts{m} -o $prefix3 -nosplit $nosplit -BestRef $BestRef -biaswt $biaswt -Mfast $Mfast -FP $FP -sf $sf -sd $sd -A $A -outlier $outlier -endoutlier $endoutlier -S $S -sr $sr -resbias $resbias -maxmem $maxmem -maxthreads $maxthreads -M $M -minlen $minlen -T $T -XmapStatWrite $prefix3.stats -stdout -stderr"; 
	}
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	chdir $mqrDir or dieLog ("ERROR: Cannot change to $mqrDir: $!\n");
	($outResults, $errResults) = runCommand($cmdRef);
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	# error check
	errCheckRefAligner("$prefix3.stdout", "END of output", "Alignment of molecules BNX to Hybrid CMAP complete in $span.", "ERROR: Alignment of molecules BNX to Hybrid CMAP cannot be completed.");
	#print "Running command: $cmd\n";
	#system($cmd) == 0 
	#	or die "system $cmd failed: $?"; 
		
	
	#$dtStartStage = DateTime->now;
	#print "\nBeginning alignment of molecules BNX to Hybrid + naive NGS CMAP...\n";
	#$prefix4 = $prefix2."_MoleculeQualityReport";
	#if ((defined($randomize)) && (defined($subset)) && ($randomize ne '0') && ($subset ne '0') && (looks_like_number($randomize)) ) {
		#chdir $mqrDir or dieLog ("ERROR: Cannot change to $mqrDir: $!\n");	
		#@$cmdRef = ($refaligner, "-f", "-ref", $hybrid_naiveNGS, "-i", $opts{m}, "-o", $prefix4, "-nosplit", $nosplit, "-BestRef", $BestRef, "-biaswt", $biaswt, "-Mfast", $Mfast, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-A", $A, "-outlier", $outlier, "-endoutlier", $endoutlier, "-S", $S, "-sr", $sr, "-resbias", split(/\s+/, $resbias), "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-M", $M, "-minlen", $minlen, "-T", $T, "-XmapStatWrite", "$prefix3.stats", "-randomize", $randomize, "-subset", split(/\s+/, $subset), "-stdout", "-stderr");
		##$cmd = "cd $mqrDir; $refaligner -f -ref $hybrid_naiveNGS -i $opts{m} -o $prefix4 -nosplit $nosplit -BestRef $BestRef -biaswt $biaswt -Mfast $Mfast -FP $FP -sf $sf -sd $sd -A $A -outlier $outlier -endoutlier $endoutlier -S $S -sr $sr -resbias $resbias -maxmem $maxmem -maxthreads $maxthreads -M $M -minlen $minlen -T $T -XmapStatWrite $prefix3.stats -randomize $randomize -subset $subset -stdout -stderr"; }
	#} else {
		#chdir $mqrDir or dieLog ("ERROR: Cannot change to $mqrDir: $!\n");	
		#@$cmdRef = ($refaligner, "-f", "-ref", $hybrid_naiveNGS, "-i", $opts{m}, "-o", $prefix4, "-nosplit", $nosplit, "-BestRef", $BestRef, "-biaswt", $biaswt, "-Mfast", $Mfast, "-FP", $FP, "-sf", $sf, "-sd", $sd, "-A", $A, "-outlier", $outlier, "-endoutlier", $endoutlier, "-S", $S, "-sr", $sr, "-resbias", split(/\s+/, $resbias), "-maxmem", $maxmem, "-maxthreads", $maxthreads, "-M", $M, "-minlen", $minlen, "-T", $T, "-XmapStatWrite", "$prefix3.stats", "-stdout", "-stderr");
		##$cmd = "cd $mqrDir; $refaligner -f -ref $hybrid_naiveNGS -i $opts{m} -o $prefix4 -nosplit $nosplit -BestRef $BestRef -biaswt $biaswt -Mfast $Mfast -FP $FP -sf $sf -sd $sd -A $A -outlier $outlier -endoutlier $endoutlier -S $S -sr $sr -resbias $resbias -maxmem $maxmem -maxthreads $maxthreads -M $M -minlen $minlen -T $T -XmapStatWrite $prefix3.stats -stdout -stderr"; 
	#}
	#print "Running command: ".(join(" ", @$cmdRef))."\n";
	#($outResults, $errResults) = runCommand($cmdRef);
	#$dtEndStage = DateTime->now;
	#$span = DateTime::Format::Human::Duration->new();
	#$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	## error check
	#errCheckRefAligner("$prefix4.stdout", "END of output", "Alignment of molecules BNX to Hybrid + naive NGS CMAP complete in $span.", "ERROR: Alignment of molecules BNX to Hybrid + naive NGS CMAP cannot be completed.");
	##print "Running command: $cmd\n";
	##system($cmd) == 0 
	##	or die "system $cmd failed: $?"; 
	
	$dtStartStage = DateTime->now;
	print "\nBeginning spliting of XMAPS...\n";
		
	#split MoleculeQualityReport.xmap into individual contig xmaps	
	my $xmap = abs_path($mqrDir."/".$prefix3.".xmap");
	my $rcmap = abs_path($mqrDir."/".$prefix3."_r.cmap");
	my $qcmap = abs_path($mqrDir."/".$prefix3."_q.cmap");
	my $contig_prefix = $prefix3."_contig";
	my $molPath = $mqrDir."/".$prefixOrig."_molecules_to_maps";
	
	#my $xmap2 = abs_path($mqrDir."/".$prefix4.".xmap");
	#my $rcmap2 = abs_path($mqrDir."/".$prefix4."_r.cmap");
	#my $qcmap2 = abs_path($mqrDir."/".$prefix4."_q.cmap");
	#my $contig_prefix2 = $prefix4."_contig";
	#my $molPath2 = $mqrDir."/".$prefix2."_molecules_to_maps";
	
	my $firstThread = threads->create(\&split_xmap,$xmap, $rcmap, $qcmap, $contig_prefix, $molPath);	
	#my $secondThread = threads->create(\&split_xmap,$xmap2, $rcmap2, $qcmap2, $contig_prefix2, $molPath2);
	
	#$_->join() foreach ( $firstThread, $secondThread );
	$_->join() foreach ( $firstThread );
	
	$dtEndStage = DateTime->now;
	$span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	print "\nXMAP splitting complete in $span.\n";
	
}

#Copy final output files to top level of output directory
my $file = $BNGfile."_".$NGSfile."_HYBRID_SCAFFOLD.cmap";
copy "$final_outputDir/step2.hybrid.cmap", "$opts{o}/$file" ; 
copy "$final_outputDir/$prefix.xmap", "$opts{o}" ;
copy "$final_outputDir/$prefix"."_q.cmap", "$opts{o}" ;
copy "$final_outputDir/$prefix"."_r.cmap", "$opts{o}" ;
copy "$final_outputDir/$prefixBNG.xmap", "$opts{o}" ;
copy "$final_outputDir/$prefixBNG"."_q.cmap", "$opts{o}" ;
copy "$final_outputDir/$prefixBNG"."_r.cmap", "$opts{o}" ;
#copy "$hybrid_naiveNGS", "$opts{o}" ;
#copy "$final_outputDir/$prefix2.xmap", "$opts{o}" ;
#copy "$final_outputDir/$prefix2"."_q.cmap", "$opts{o}" ;
#copy "$final_outputDir/$prefix2"."_r.cmap", "$opts{o}" ;
if(defined($opts{m}) && defined($mqrDir) && defined($prefix3) && defined($prefix4)) {
	copy "$mqrDir/$prefix3.xmap", "$opts{o}" ;
	copy "$mqrDir/$prefix3"."_q.cmap", "$opts{o}" ;
	copy "$mqrDir/$prefix3"."_r.cmap", "$opts{o}" ;
	
	copy "$mqrDir/$prefix4.xmap", "$opts{o}" ;
	copy "$mqrDir/$prefix4"."_q.cmap", "$opts{o}" ;
	copy "$mqrDir/$prefix4"."_r.cmap", "$opts{o}" ;
}
my $faPath = "$opts{o}/fa2cmap";
eval { mkpath($faPath) };
copy "$ngs_cmap", "$faPath" ;
copy "$opts{b}", "$faPath" ;

#Calculate some stats about CMAPs

$dtStartStage = DateTime->now;
print "\nCalculating statistics...\n\n";

print "Original BioNano CMAP statistics:\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
# @$cmdRef = ("Rscript", "./scripts/calc_cmap_stats.R", $opts{b});
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $opts{b});
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# errorCheck
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for the original BioNano CMAP file");
# print statistics
print "$outResults\n";

#$cmd = "cd $plDir; Rscript ./scripts/calc_cmap_stats.R $opts{b}"; 
#system($cmd) == 0 
#	or die "system $cmd failed: $?";

print "Original NGS CMAP statistics:\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
# @$cmdRef = ("Rscript", "./scripts/calc_cmap_stats.R", $ngs_cmap);
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $ngs_cmap);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for original NGS CMAP file");
# print statistics
print "$outResults\n";
#$cmd = "cd $plDir; Rscript ./scripts/calc_cmap_stats.R $ngs_cmap"; 
#system($cmd) == 0 
#	or die "system $cmd failed: $?";

print "Filtered BioNano CMAP statistics:\n";
if (!$opts{B}) {	
	$assignAlignType_q_cmap = abs_path($assignAlignType_q_cmap);
	chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
	# @$cmdRef = ("Rscript", "./scripts/calc_cmap_stats.R", $assignAlignType_q_cmap);
	@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $assignAlignType_q_cmap);
	print "Running command: ".(join(" ", @$cmdRef))."\n";
	($outResults, $errResults) = runCommand($cmdRef);
	# error check
	errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for Filtered BioNano CMAP file");
	# print statistics
	print "$outResults\n";
	#$cmd = "cd $plDir; Rscript ./scripts/calc_cmap_stats.R $assignAlignType_q_cmap"; 
	#system($cmd) == 0 
	#	or die "system $cmd failed: $?";
}
else {
	print "Using unfiltered original BioNano CMAP\n\n";
}

print "Filtered NGS CMAP statistics:\n";
if (!$opts{N}) {
	$assignAlignType_r_cmap = abs_path($assignAlignType_r_cmap);
	chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
	# @$cmdRef = ("Rscript", "./scripts/calc_cmap_stats.R", $assignAlignType_r_cmap);
	@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $assignAlignType_r_cmap);
	print "Running command: ".(join(" ", @$cmdRef))."\n";	
	($outResults, $errResults) = runCommand($cmdRef);
	# error check
	errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for Filtered NGS CMAP file");
	# print statistics
	print "$outResults\n";
	#$cmd = "cd $plDir; Rscript ./scripts/calc_cmap_stats.R $assignAlignType_r_cmap"; 
	#system($cmd) == 0 
	#	or die "system $cmd failed: $?";
}
else {
	print "Using unfiltered original NGS CMAP\n\n";
}

print "Step1 naive (leftovers) BioNano CMAP statistics:\n";
$BNG_naive_cmap = abs_path($BNG_naive_cmap);
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
# @$cmdRef = ("Rscript", "./scripts/calc_cmap_stats.R", $BNG_naive_cmap);
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $BNG_naive_cmap);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for Step1 naive (leftovers) BioNano CMAP file");
# print statistics
print "$outResults\n";
#$cmd = "cd $plDir; Rscript ./scripts/calc_cmap_stats.R $BNG_naive_cmap"; 
#system($cmd) == 0 
#	or die "system $cmd failed: $?";

print "Step1 naive (leftovers) NGS CMAP statistics:\n";
$NGS_naive_cmap = abs_path($NGS_naive_cmap);
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
# @$cmdRef = ("Rscript", "./scripts/calc_cmap_stats.R", $BNG_naive_cmap);
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $NGS_naive_cmap);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for Step1 naive (leftovers) NGS CMAP file");
# print statistics
print "$outResults\n";
#$cmd = "cd $plDir; Rscript ./scripts/calc_cmap_stats.R $NGS_naive_cmap"; 
#system($cmd) == 0 
#	or die "system $cmd failed: $?";

print "BNG contigs in hybrid CMAP statistics:\n";
$filteredBNG_file = abs_path($filteredBNG_file);
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
# @$cmdRef = ("Rscript", "./scripts/calc_cmap_stats.R", $filteredNGS_file);
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $filteredBNG_file);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for CMAP of BNG contigs used in hybrid scaffolding");
# print statistics
print "$outResults\n";
#$cmd = "cd $plDir; Rscript ./scripts/calc_cmap_stats.R $filteredNGS_file"; 
#system($cmd) == 0 
#	or die "system $cmd failed: $?";

print "NGS contigs in hybrid CMAP statistics:\n";
$filteredNGS_file = abs_path($filteredNGS_file);
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
# @$cmdRef = ("Rscript", "./scripts/calc_cmap_stats.R", $filteredNGS_file);
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $filteredNGS_file);
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for CMAP of NGS contigs used in hybrid scaffolding");
# print statistics
print "$outResults\n";
#$cmd = "cd $plDir; Rscript ./scripts/calc_cmap_stats.R $filteredNGS_file"; 
#system($cmd) == 0 
#	or die "system $cmd failed: $?";

print "Hybrid CMAP statistics:\n";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
# @$cmdRef = ("Rscript", "./scripts/calc_cmap_stats.R", "$opts{o}/$file");
@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", "$opts{o}/$file");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for Hybrid CMAP");
# print statistics
print "$outResults\n";
#$cmd = "cd $plDir; Rscript ./scripts/calc_cmap_stats.R $opts{o}/$file"; 
#print "Running command: $cmd\n\n";
#system($cmd) == 0 
#	or die "system $cmd failed: $?";
	
#print "\nHybrid + naive NGS CMAP statistics:\n";
#chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
## @$cmdRef = ("Rscript", "./scripts/calc_cmap_stats.R", $hybrid_naiveNGS);
#@$cmdRef = ("$^X", "./scripts/calc_cmap_stats.pl", $hybrid_naiveNGS);
#print "Running command: ".(join(" ", @$cmdRef))."\n";
#($outResults, $errResults) = runCommand($cmdRef);
## error check
#errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for Hybrid + naive NGS CMAP");
## print statistics
#print "$outResults\n";

#$cmd = "cd $plDir; Rscript ./scripts/calc_cmap_stats.R $hybrid_naiveNGS";
#print "Running command: $cmd\n\n";
#system($cmd) == 0 
#	or die "system $cmd failed: $?";

#calculate XMAP stats and output to .stats file

$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
print "\nCalculating CMAP statistics complete in $span.\n";

$dtStartStage = DateTime->now;
print "\nCalculating XMAP statistics...\n";

$file = "$opts{o}/$file";
chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
@$cmdRef = ($^X, "scripts/calc_xmap_stats.pl", $opts{o}, "$prefix.xmap", $file, "$prefix.xmap.stats");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
open(OUT, ">$final_outputDir/calc_xmap_stats.log") or dieLog ("ERROR: Cannot write to $final_outputDir/calc_xmap_stats.log: $!\n");	print OUT $outResults; close OUT;
open(ERR, ">$final_outputDir/calc_xmap_stats.errlog") or dieLog ("ERROR: Cannot write to $final_outputDir/calc_xmap_stats.errlog: $!\n");	print ERR $errResults;	close ERR;
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for xmap $prefix.xmap");
# print statistics
print "$outResults\n";
#$cmd = "cd $plDir; $^X scripts/calc_xmap_stats.pl $opts{o} $prefix.xmap $file $prefix.xmap.stats &> $final_outputDir/calc_xmap_stats.log";
#print "\n\nRunning command: $cmd\n\n";
#$results = capture($cmd);

chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
@$cmdRef = ($^X, "scripts/calc_xmap_stats.pl", $opts{o}, "$prefixBNG.xmap", $file, "$prefixBNG.xmap.stats");
print "Running command: ".(join(" ", @$cmdRef))."\n";
($outResults, $errResults) = runCommand($cmdRef);
open(OUT, ">$final_outputDir/calc_xmap_stats.log") or dieLog ("ERROR: Cannot write to $final_outputDir/calc_xmap_stats.log: $!\n");	print OUT $outResults; close OUT;
open(ERR, ">$final_outputDir/calc_xmap_stats.errlog") or dieLog ("ERROR: Cannot write to $final_outputDir/calc_xmap_stats.errlog: $!\n");	print ERR $errResults;	close ERR;
# error check
errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for xmap $prefixBNG.xmap");
# print statistics
print "$outResults\n";
#$cmd = "cd $plDir; $^X scripts/calc_xmap_stats.pl $opts{o} $prefix.xmap $file $prefix.xmap.stats &> $final_outputDir/calc_xmap_stats.log";
#print "\n\nRunning command: $cmd\n\n";
#$results = capture($cmd);

#chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");	
#@$cmdRef = ($^X, "scripts/calc_xmap_stats.pl", $opts{o}, "$prefix2.xmap", $hybrid_naiveNGS, "$prefix2.xmap.stats");
#print "Running command: ".(join(" ", @$cmdRef))."\n";
#($outResults, $errResults) = runCommand($cmdRef);
#open(OUT, ">>$final_outputDir/calc_xmap_stats.log") or dieLog ("ERROR: Cannot write to $final_outputDir/calc_xmap_stats.log: $!\n");	print OUT $outResults; close OUT;
#open(ERR, ">>$final_outputDir/calc_xmap_stats.errlog") or dieLog ("ERROR: Cannot write to $final_outputDir/calc_xmap_stats.errlog: $!\n");	print ERR $errResults;	close ERR;
## error check
#errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for xmap $prefix2.xmap");
## print statistics
#print "$outResults\n";
##$cmd = "cd $plDir; $^X scripts/calc_xmap_stats.pl $opts{o} $prefix2.xmap $hybrid_naiveNGS $prefix2.xmap.stats &>> $final_outputDir/calc_xmap_stats.log";
##print "\n\nRunning command: $cmd\n\n";
##$results = capture($cmd);

#if(defined($opts{m})) {
	
	#chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
	#@$cmdRef = ($^X, "scripts/calc_xmap_stats.pl", $opts{o}, "$prefix3.xmap", $file, "$prefix3.xmap.stats");
	#print "Running command: ".(join(" ", @$cmdRef))."\n";
	#($outResults, $errResults) = runCommand($cmdRef);
	#open(OUT, ">$mqrDir/calc_xmap_stats.log") or dieLog ("ERROR: Cannot write to $mqrDir/calc_xmap_stats.log: $!\n");	print OUT $outResults;	close OUT;
	#open(ERR, ">$mqrDir/calc_xmap_stats.errlog") or dieLog ("ERROR: Cannot write to $mqrDir/calc_xmap_stats.errlog: $!\n");	print ERR $errResults;	close ERR;
	## error check
	#errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for xmap $prefix3.xmap");
	## print statistics
	#print "$outResults\n";


	##$cmd = "cd $plDir; $^X scripts/calc_xmap_stats.pl $opts{o} $prefix3.xmap $file $prefix3.xmap.stats &> $mqrDir/calc_xmap_stats.log";
	##print "\n\nRunning command: $cmd\n\n";
	##$results = capture($cmd);

	#chdir $plDir or dieLog ("ERROR: Cannot change directory to $plDir: $!\n");
	#@$cmdRef = ($^X, "scripts/calc_xmap_stats.pl", $opts{o}, "$prefix4.xmap", $hybrid_naiveNGS, "$prefix4.xmap.stats");
	#print "Running command: ".(join(" ", @$cmdRef))."\n";
	#($outResults, $errResults) = runCommand($cmdRef);
	#open(OUT, ">>$mqrDir/calc_xmap_stats.log") or dieLog ("ERROR: Cannot write to $mqrDir/calc_xmap_stats.log: $!\n");	print OUT $outResults;	close OUT;
	#open(ERR, ">>$mqrDir/calc_xmap_stats.errlog") or dieLog ("ERROR: Cannot write to $mqrDir/calc_xmap_stats.errlog: $!\n");	print ERR $errResults;	close ERR;
	## error check
	#errCheck($outResults, $errResults, "ERROR:", "", "ERROR: cannot calculate statistics for xmap $prefix4.xmap");
	## print statistics
	#print "$outResults\n";	

	##$cmd = "cd $plDir; $^X scripts/calc_xmap_stats.pl $opts{o} $prefix4.xmap $hybrid_naiveNGS $prefix4.xmap.stats &>> $mqrDir/calc_xmap_stats.log";
	##print "\n\nRunning command: $cmd\n\n";
	##$results = capture($cmd);
#}

$dtEndStage = DateTime->now;
$span = DateTime::Format::Human::Duration->new();
$span = $span->format_duration_between($dtEndStage, $dtStartStage);
print "XMAP statistics calculation complete in $span.\n";

print "\n";
my $dtEnd = DateTime->now;
print "End time: "; print join ' ', $dtEnd->ymd, $dtEnd->hms; print "\n\n";

$span = DateTime::Format::Human::Duration->new();
print 'Total elapsed time: ', $span->format_duration_between($dtEnd, $dtStart); print "\n";

print "\n\nMerging of $opts{b} with $opts{n} is complete.\n\n";


print "END of output\n";



#print "Modules used:\n";
#print join("\n", map { s|/|::|g; s|\.pm$||; $_ } keys %INC);
#print "\n\n";
#print  "Files used:\n";
#print join("\n", map { s|/|/|g; $_ } keys %INC);
#print "\n\n";
#print "Absolute path:\n";
#chdir($plDir);
#my @files;

#find(\&d, "$plDir");

#my $allModules = join("#", map { s|/|/|g; $_ } keys %INC);
#my @modules = split("#", $allModules);

#foreach my $file (@files) {
	#foreach my $module (@modules) {
		#if ($file =~ /$module/) {
			#print "$file\n"; 
			#fcopy("$file", "perl_modules/$file"); } } }

#sub d {
#-f and -r and push  @files, $File::Find::name;
#}











######################################################################
#                           Subroutines                              #
######################################################################

# this subroutine is used to check RefAligner processes to see if it finishes successfully
# it assumes that -stderr and -stdout was used to run RefAligner
sub errCheckRefAligner	{
	my ($file, $completeString, $completeMsg, $dieMsg) = @_;
	open(IN, "$file") or dieLog ("ERROR: Cannot open $file: $!\n");
	my $presence = 0;
	while (my $line = <IN>)	{
		if ($line =~ /$completeString/)	{
			# if the line contains the string that indicates successful completion
			$presence = 1;
		} # if line
	} # while line
	if ($presence == 0)	{
		dieLog ("ERROR: $dieMsg\n");
	} else	{
		print "$completeMsg\n";
	} # if presence
	close IN;
} # errCheckRefAligner

# this subrountine is used to call non-RefAligner processes to see if there is error in the child
# it assumes that the child process prints out "ERROR:" tag to all messages
sub errCheck	{
	my ($outResults, $errResults, $errorString, $completeMsg, $dieMsg) = @_;
	if ($outResults =~ /$errorString/)	{
		dieLog ("ERROR: $dieMsg\n"); } 
	elsif ($errResults =~ /$errorString/)	{
		dieLog ("ERROR: $dieMsg\n"); } 
	else {
		print "$completeMsg\n"; } 
}

sub runCommand  {
        my ($argsRef) = @_;
#print "runCommand: before open3\n";
        my $pid = open3(my $CMD_IN, my $out, my $err, @$argsRef);
#print "runCommand: after open3\n";

        close($CMD_IN);

        my $outResults = "";
        my $errResults = "";
        my $sel = new IO::Select;
        $sel->add($out, $err);
        while(my @fhs = $sel->can_read) {
                foreach my $fh (@fhs) {
                        my $line = <$fh>;
                        unless(defined $line) {
                                $sel->remove($fh);
                                next;
                        } # unless line
                        if($fh == $out) {
                                $outResults .= "$line";
                                #print "$line";
                        }elsif($fh == $err) {
                                $errResults .= "$line";
                                #print "$line";
                        }else{
                                dieLog ("ERROR: This should never execute!");
                        } # if fh
                } # foreach fh
        } # while
#print "runCommand: before waitpid\n";
        my $ret=waitpid ($pid, 0); # reap the exit code
#print "runCommand: after waitpid\n";
        return ($outResults, $errResults);
} # runCommand

sub Init{
	my $opt_string = 'hn:b:c:o:BNfm:v';
	if(!getopts("$opt_string", \%opts)){
		print("Invalid parameter(s)! Try -h for more information.\n");
		Usage();
	}
	Usage() if $opts{h};
	if ($opts{v}) {
		Version();
		exit; }
	$opts{b} = abs_path($opts{b});
	$opts{n} = abs_path($opts{n}); 
	$opts{c} = abs_path($opts{c});
	if(defined($opts{o})) {
		$opts{o} = abs_path($opts{o}); }
	if(defined($opts{m})) {
		if(!-e $opts{n}) {
			print "\nInput molecules BNX file $opts{m} does not exist! Skipping molecules to Hybrid CMAP alignment.\n";
			undef $opts{m}; }
		else {
			$opts{m} = abs_path($opts{m}); } }
}

sub Version{
	my $dtStartStage = DateTime->now;

	my @revs;
	my @pipelineFiles = process_files($plDir);
	foreach (@pipelineFiles) {
		if (!-d $_) {
			#if ($_ !~ /svn/ && $_ !~ /perl5/) {
				open (FILE,$_);
				my @file = <FILE>; 
				close FILE;
				my $start = '$Id:'; my $end = '$';
				#print "$_\n";
				foreach (@file) {
					if ($_ =~ m/\$Id/) { 
						chomp($_);
						my ($wanted) = $_ =~ /\$Id:(.*)\$/;
						$wanted =~ s/^\s+|\s+$//g;
						push @revs, $wanted;
						#print $wanted."\n"; 
						last; 
						} } } } 
						#}
	
	my $revLine; my $revNum=0;
	@revs = uniq @revs;
	foreach (@revs) {
		my @split = split(/\s+/);
		#print "Version: $split[2]  Line: $_\n";
		if ($split[1] >= $revNum) {
			$revNum = $split[1];
			$revLine = $_; } }
	
	print "$revLine\n";  

	my $dtEndStage = DateTime->now;
	my $span = DateTime::Format::Human::Duration->new();
	$span = $span->format_duration_between($dtEndStage, $dtStartStage);
	#print "Total time: $span\n";
}

sub process_files {
	my $path = shift;

    opendir (DIR, $path)
        or dieLog ("ERROR: Unable to open $path: $!");

    # We are just chaining the grep and map
    # This is the same as: LIST = map(EXP, grep(EXP, readdir()))
    my @files =
        # Third: Prepend the full path
        map { $path . '/' . $_ }
        # Second: take out '.' and '..'
        grep { !/^\.{1,2}$/ }
        # First: get all files
        readdir (DIR);

    closedir (DIR);

    for (@files) {
        if (-d $_) {
            # Add all of the new files from this directory
            # (and its subdirectories, and so on... if any)
            if (abs_path($_) !~ /5.10.1/ && abs_path($_) !~ /5.14.4/ && abs_path($_) !~ /svn/ && abs_path($_) =~ /scripts/) {
				push @files, process_files ($_); } }

        else {
            # Do whatever you want here =) .. if anything.
        }
    }
    # NOTE: we're returning the list of files
    return @files;
}

sub Usage{
	print << "EOF";
	
Usage: perl hybridScaffold.pl <-h> <-n fasta_file> <-b bng_cmap_file> <-c config_file> <-o output_folder> <-B> <-N> <-f> <-m molecules_bnx> <-v>
      -h    : This help message         
      -n    : Input NGS FASTA or CMAP file [required]
      -b    : Input BioNano CMAP  [required]
      -c    : Merge configuration file [required]
      -o    : Output folder [required]
      -B    : Use unfiltered BioNano CMAP 
      -N    : Use unfiltered NGS CMAP 
      -f    : Force output and overwrite any existing files
      -m    : Input BioNano molecules BNX for molecules to hybrid scaffold alignment [optional]
      -v    : Print pipeline version information
      
EOF
	exit;
} # Usage

sub Parse{
	my ($config) = @_;
	my $tmp_key;
	my $tmp_val_1;
	my $tmp_val_2="";

	my %hash;

	foreach my $node (keys %$config){
		#print $node, "\n";  # align_final

		foreach my $item (keys %{$config->{$node}}){
			#print $item, "\n";  # flag

			if(ref($config->{$node}{$item}) eq 'ARRAY'){
				foreach my $attr (@{$config->{$node}{$item}}){
					foreach my $value (keys %$attr){
						if($value eq 'attr'){
							$tmp_key = $node . "." . $attr->{$value};
						}
						if($value eq 'val0'){
							$tmp_val_1 = $attr->{$value};
						}
						if($value eq 'val1'){
							$tmp_val_2 = $attr->{$value};
						}
					}

					if(length($tmp_val_2)){
						#print $tmp_key, " = ", $tmp_val_1, " ", $tmp_val_2, "\n";
						$hash{$tmp_key} = $tmp_val_1 . " " . $tmp_val_2;
						$tmp_val_2 = "";
					}
					else{
						#print $tmp_key, " = ", $tmp_val_1, "\n";
						$hash{$tmp_key} = $tmp_val_1;
					}
				}
			}
			elsif(ref($config->{$node}{$item}) eq 'HASH'){
				foreach my $attr (keys %{$config->{$node}{$item}}){
					if($attr eq 'attr'){
						$tmp_key = $node . "." . $config->{$node}{$item}{$attr};
					}
					if($attr eq 'val0'){
						$tmp_val_1 = $config->{$node}{$item}{$attr};
					}
					if($attr eq 'val1'){
						$tmp_val_2 = $config->{$node}{$item}{$attr};
					}
				}

				if(length($tmp_val_2)){
					#print $tmp_key, " = ", $tmp_val_1, " ", $tmp_val_2, "\n";
					$hash{$tmp_key} = $tmp_val_1 . " " . $tmp_val_2;
					$tmp_val_2 = "";
				}
				else{
					#print $tmp_key, " = ", $tmp_val_1, "\n";
					$hash{$tmp_key} = $tmp_val_1;
				}
			}
			else{
				dieLog ("ERROR: in parsing XML configuration file!\n");
			}
		}
	}
	return (%hash);
} # Parse

sub CHECK_Config{
	my ($config) = @_;
	foreach my $item (keys %$config){
		if($item =~ 'RepeatMask' || $item =~ 'RepeatRec' || $item =~ 'subset' || $item =~ 'resbias'){
			my @data = split(" ", $config->{$item});
			# Test if the option exists and the value is defined
			if(@data != 2 || !defined($data[0]) || !length($data[0])){
				dieLog ("ERROR: Configuration file error, please make sure that $item has two values.\n");
			}
			if(@data != 2 || !defined($data[1]) || !length($data[1])){
				dieLog ("ERROR: Configuration file error, please make sure that $item has two values.\n");
			}
		}
		else{
			# Test if the option exists and the value is defined
			if(!defined($config->{$item}) || !length($config->{$item})){
				dieLog ("ERROR: Configuration file error, please make sure that $item has a value.\n");
			}
		}
		
		# Test if the required numerical values are numerical
		my @non_num_params = qw(global.refaligner fasta2cmap.enzyme);
		if(!$item ~~ @non_num_params){
			if($item =~ 'RepeatMask' || $item =~ 'RepeatRec' || $item =~ 'subset' || $item =~ 'resbias'){
				my @data = split(" ", $config->{$item});
				if(!looks_like_number($data[0]) || !looks_like_number($data[1])){
					dieLog ("ERROR: Configuration file error, please make sure that $item is numerical.\n");
				}
			}
			else{
				if(!looks_like_number($config->{$item})){
					dieLog ("ERROR: Configuration file error, please make sure that $item is numerical.\n");
				}
			}
		}
	}
}


sub split_xmap {
	my ($xmap, $rcmap, $qcmap, $contig_prefix, $molPath) = @_;
	 	
	my @xmap_in = read_file($xmap);	
	my @rcmap_in = read_file($rcmap);
	my @qcmap_in = read_file($qcmap);
	
	my %rcmap;
	my $rcmap_header="";
	my %qcmap;
	my $qcmap_header="";
	
	my @xmapheader;
	
	print "\nSplitting $xmap...\n";
	
	for (my $i=0; $i <= 4; $i++) {
		push (@xmapheader, $xmap_in[$i]); }
		
	eval { mkpath($molPath) };
	if ($@) {
		dieLog ("ERROR: Couldn't create $molPath: $@"); }
		
	
	foreach my $rline (@rcmap_in) {
	if ($rline =~ /^#/ ) {
		$rcmap_header = $rcmap_header . $rline; }
	else {
		my ($CMapId,$ContigLength,$NumSites,$SiteID,$LabelChannel,$Position,$StdDev,$Coverage,$Occurrence,$GmeanSNR,$lnSNRsd) = split("\t",$rline);
		if ($CMapId) {
			push @{ $rcmap{$CMapId} }, $rline; }}}
			
	foreach my $qline (@qcmap_in) {
	if ($qline =~ /^#/ ) {
		$qcmap_header = $qcmap_header . $qline;  }
	else {
		my ($CMapId,$ContigLength,$NumSites,$SiteID,$LabelChannel,$PositionStdDev,$Coverage,$Occurrence) = split("\t",$qline);
		if ($CMapId) {
			push @{ $qcmap{$CMapId} }, $qline; }}}
	
	
	my $prevContigID = 0;
	foreach my $xline (@xmap_in) {
		if ($xline !~ /^#/ ) {			
			my ($XmapEntryID, $QryContigID, $RefContigID, $QryStartPos, $QryEndPos, $RefStartPos, $RefEndPos, $Orientation, $Confidence, $HitEnum) = split("\t",$xline);  
			if ($RefContigID != $prevContigID) {
				
				my $xmapName = $contig_prefix.$RefContigID.".xmap";
				my $rcmapName = $contig_prefix.$RefContigID."_r.cmap";
				my $qcmapName = $contig_prefix.$RefContigID."_q.cmap";
				
				open (XMAP, ">>$molPath/$xmapName");
				open (RCMAP, ">>$molPath/$rcmapName");
				open (QCMAP, ">>$molPath/$qcmapName");
				
				foreach my $s (@xmapheader) {
					$s =~ s/^\s+//;
					print XMAP $s; }
				print XMAP "# Reference Maps From:	".$rcmapName."\n";
				print XMAP "# Query Maps From:	".$qcmapName."\n";
				print XMAP "#h XmapEntryID	QryContigID	RefContigID	QryStartPos	QryEndPos	RefStartPos	RefEndPos	Orientation	Confidence	HitEnum\n";
				print XMAP "#f int        	int        	int        	float      	float    	float      	float    	string     	float     	string\n"; 
				print XMAP $xline;
				

				
				# RCMAP output
				print RCMAP $rcmap_header; 
				if($rcmap{$RefContigID}){
					print RCMAP join("", @{ $rcmap{$RefContigID} }); }
				
				#QCMAP output
				print QCMAP $qcmap_header;
				if($qcmap{$QryContigID}){
					print QCMAP join("", @{ $qcmap{$QryContigID} }); }
				
				$prevContigID = $RefContigID;
					
				close(QCMAP);
				close(RCMAP);
				close(XMAP);
				}				
				
			elsif ($RefContigID == $prevContigID) {
				
				my $qcmapName = $contig_prefix.$RefContigID."_q.cmap";
				my $xmapName = $contig_prefix.$RefContigID.".xmap";
				
				open (XMAP, ">>$molPath/$xmapName");
				open (QCMAP, ">>$molPath/$qcmapName");				

				print XMAP $xline; 
				
				if($qcmap{$QryContigID}){
					print QCMAP join("", @{ $qcmap{$QryContigID} }); }

				close(QCMAP);
				close(XMAP);				
				
				} } }
	
	print "Splitting $xmap complete.\n"; 	
}

sub find_refaligner {
	if (! -d $_) {
		if ($_ =~ /RefAligner/i) {
			$refaligner = $File::Find::name; }}
}

sub log10 {
	my $n = shift;
    return log($n)/log(10);
}
