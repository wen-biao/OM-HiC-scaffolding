# $Id: MergeNGS_BN.pl 3489 2015-01-05 21:31:51Z apang $

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

print "\nInfo: Running the command: $0 @ARGV\n";

use File::Path qw(make_path);
use File::Copy;
use BNG::Utility;
use BNG::refAlignerRun;

#declare and read in command line args
my $outputDir=$ARGV[0];
my $refaligner=$ARGV[1];
my $merge_Tvalue=$ARGV[2];
my $scratchDir=$ARGV[3];
my $logFile=$ARGV[4];
my $ngs_cmap_fn=$ARGV[5];
my $bng_cmap_fn=$ARGV[6];
my $id_shift=$ARGV[7];
my $max_merge_rounds=$ARGV[8];
my $maxthreads=$ARGV[9];
my $endoutlier=$ARGV[10];
my $outlier=$ARGV[11];
my $biaswt=$ARGV[12];
my $sd=$ARGV[13];
my $res=$ARGV[14];
my $sf=$ARGV[15];
my $RepeatMask = $ARGV[16] . " " . $ARGV[17];
my $RepeatRec = $ARGV[18] . " " . $ARGV[19];
#####################################my $pairmerge = $ARGV[20] . " " . $ARGV[21];
#######################################
#my $readparameters = $ARGV[22];
my $readparameters="$outputDir/../align1/align1.errbin";
print "readpara\t$readparameters\t20\t$ARGV[20]\t$ARGV[22]\n";

#my $readparameters = $ARGV[20];
# C:/Users/ZhanyangZhu/workspace/BioNano/dev/Export_HybridScaffold/data/arabidopsis/mergeNGS_BN 
# "C:/Program Files/BioNano Genomics/SetupRefAligner/WindowsRefAligner" 
# 1e-10 ./ mergeNGS_BN_pl.log ../assignAlignType/assignAlignType_r.cmap 
# ../assignAlignType/assignAlignType_q.cmap 100000 40 64 0 1e-4 0 0.1 2.9 0.2 2 0.01 0.7 0.6
#print "\nRepeatMask=$RepeatMask;\nRepeatRec=$RepeatRec\noutlier=$outlier\nbiaswt=$biaswt\n"; 
##################################my $doDebug = 1; 
my $doDebug = 0; 

#####################
# Step 0. Initiation.
#####################
chdir $outputDir;	# change current directory to $outputDir
print "cwd=$outputDir\n";
make_path($scratchDir); # make sure $scratchDir exisit
# 0.1. Read in NGS CMap
open(LOG, ">$logFile");
my $logFH = *LOG;

my ($ngs_cmap, $nc_ngs_cmap, undef) = readCMap($ngs_cmap_fn);
logMessage($logFH, "Read NGS contig '" .$ngs_cmap_fn . "' completed with ", $nc_ngs_cmap, " cmaps.");  

# 0.2. Read in BioNano Assembled Cmap    
my ($bng_cmap, $nc_bng_cmap, undef) = readCMap($bng_cmap_fn);
logMessage($logFH, "Read BioNano contig '" .$bng_cmap_fn . "' completed with ", $nc_bng_cmap, " cmaps."); 

# 0.3. Shift BioNano ContigId by offset to prevent ContigID collision
shiftCMapIds($bng_cmap, $id_shift);
my $bionano_idshifted_cmap_fn = $scratchDir . "bionano.idshift.cmap";
writeCMapFile($bng_cmap, $bionano_idshifted_cmap_fn, 0);
logMessage($logFH, "Output BioNano contig with ID shift of '" . $id_shift . "' completed and output to '" . $bionano_idshifted_cmap_fn, "'.");
logMessage($logFH, "Step 0. Initiation Completed");
# merge round ids - max 130 (5*26) rounds; starts from 1	
my @mrg_rounds_ids = (" ", "A".."Z", "AA".."AZ", "BA".."BZ", "CA".."CZ", "DA".."DZ");
if ($max_merge_rounds > 130 ) {
	$max_merge_rounds = 130; # make sure max rounds no more than 130
} # if max_merge_rounds

#####################
# Step 1.
# we first run merges between BN/NGS, to minimize within group merge,
# we only merge between merged and leftover from round2
##################### 
my $step1_hybrid_cmap_name=""; 
my $all_step1_mrg_pairs = {};

logMessage($logFH, "Step 1. Pair Merge between NGS\/BioNano started.");

my $rounds = 1;
my $do_mrg = 1;

# keep track whether there is any left over files (both NGS and BN), if this is true, that means the merge has exhausted all input records. if that is the case, the last successful merge is actually this.mrg.round, not prev.mrg.round
my $no_left_over = 0;	
my $this_mrg_prefix = "";
my $prev_mrg_prefix = "";
my $this_mrg_round_id = "";
my $prev_mrg_round_id = "";
my $pairmerge_single_round_2files_run = 
	  new BNG::refAlignerRun({binPath 	=> $refaligner,
						 # i			=> [$inputFile1,$inputFile2],
						 # o			=> $mrg_output_prefix,
						 T			=> $merge_Tvalue,
						 pairmerge	=> 1, 					##########################################################pair merge
						 preferNGS	=> 1,
						 endoutlier	=> $endoutlier,
						 outlier	=> $outlier,
						 biaswt		=> $biaswt,
						 maxthreads => $maxthreads,
						 sd 		=> $sd,                 #noise.params
						 res 		=> $res,
						 sf			=> $sf,
						 f			=> 1,                   # overwrite
						 stdout		=> 1,
						 stderr		=> 1,
						 first		=> -1,
						 RepeatMask => "$RepeatMask",
						 RepeatRec 	=> "$RepeatRec",	  
						 ############################################ pairmerge  => "$pairmerge",
						 readparameters => "$readparameters",
                          NoBpp		=> 1        					 
	  } );
#my ($outResults, $errResults, $job_status);
while ($do_mrg) { 
	$this_mrg_round_id = $mrg_rounds_ids[$rounds];
	if ($rounds > 1) { 
		$prev_mrg_round_id = $mrg_rounds_ids[$rounds-1];
	} # if rounds
	$this_mrg_prefix = "Mrg" . $this_mrg_round_id;
	$prev_mrg_prefix = "Mrg" . $prev_mrg_round_id;

	# for each round, we take BioNano/NGS merged contig and merge against leftover NGS contig
	# for first round, we take all BioNano contig and merge against all NGS contig
	my  $inputFile1 = $scratchDir . $prev_mrg_prefix . "_mrged.cmap";
	my  $inputFile2 = $scratchDir . $prev_mrg_prefix . "_leftover.cmap";       
	if ($rounds == 1) { 
		$inputFile1 = $bionano_idshifted_cmap_fn;
		$inputFile2 = $ngs_cmap_fn;	    	
	} # if rounds
	logMessage($logFH, "\tRound " . $rounds . " " . $this_mrg_prefix . " started between '" . $inputFile1 . "' and '" . $inputFile2 . "'");

	my $mrg_output_prefix = $scratchDir . '/' . $this_mrg_prefix;
	$pairmerge_single_round_2files_run->setParams({i => [$inputFile1,$inputFile2], o => $mrg_output_prefix} );	
	my $cmd = $pairmerge_single_round_2files_run->getCMD();	
	logMessage($logFH,$cmd);
	my ($outResults, $errResults, $job_status) = $pairmerge_single_round_2files_run->runCMD();
	exit $job_status if ($job_status !=0);
	my (@cleanup_file_list) = find_cmap_list_w_prefix($mrg_output_prefix);
	push(@cleanup_file_list, $mrg_output_prefix . '.stdout');
	push(@cleanup_file_list, $mrg_output_prefix . '.align');

	my ($numMrgP, $this_mrg_pairs)= parsingMrgStdOut($mrg_output_prefix . '.stdout');
	my $MgrResult = {};
	$MgrResult->{mergeCnt} =$numMrgP;
	$MgrResult->{mergePairs} = $this_mrg_pairs;
	#       read in the two cmaps to only save the number of cmap ids:
	my ($cmap_m_1, $num1, undef) = readCMap($inputFile1);
	my ($cmap_m_2, $num2, undef) = readCMap($inputFile2);
	if ($numMrgP == 0) { 
		logMessage($logFH, "No merge found for " . $mrg_output_prefix);
		# define $result when no merge was found:
		#        # there was no merge at all.
		# $MgrResult->{mergeCnt} =0;
		$MgrResult->{leftoverInput1} = $num1;
		$MgrResult->{leftoverInput2} = $num2;
		# $MgrResult->{mergePairs} = {};          
	} else {	
		my $leftover_contigIds_1 = uniqueIDs(getCMapIds($cmap_m_1), $this_mrg_pairs->{ContigID1});
		my $l1 = @$leftover_contigIds_1;
		my $leftover_contigIds_2 = uniqueIDs(getCMapIds($cmap_m_2), $this_mrg_pairs->{ContigID2});
		my $l2 = @$leftover_contigIds_2;
		logMessage($logFH, "\tAfter merge, LeftOver input1 include " . $l1 . " contigs");
		logMessage($logFH, "\tAfter merge, LeftOver input2 include " . $l2 . " contigs");    

		$MgrResult->{leftoverInput1} = $l1;
		$MgrResult->{leftoverInput2} = $l2;

		# Now we need to merge cmap files into final products.
		# first we do merged

		my $this_output_mrged_list_prefix = $mrg_output_prefix . "_mrged";
		my $this_output_mrged_cmap_file = $this_output_mrged_list_prefix . ".cmap";
		my $this_output_mrged_list = getArrayFileNames($mrg_output_prefix . "_contig",  $this_mrg_pairs->{ResultContigID}, ".cmap");

		my $run_merge = 
			new BNG::refAlignerRun({binPath 	=> $refaligner,
							 i			=> \@$this_output_mrged_list,
							 o			=> $this_output_mrged_list_prefix,
							 merge		=> 1,
							 f			=> 1,                   			# overwrite
							 stdout		=> 1,
							 stderr		=> 1	          					 
			} );	
		logMessage($logFH, $run_merge->getCMD());
		my ($outResults, $errResults, $job_status) = $run_merge->runCMD();
		exit $job_status if ($job_status !=0);
		push(@cleanup_file_list, $this_output_mrged_list_prefix . ".idmap");
		push(@cleanup_file_list, $this_output_mrged_list_prefix . ".stdout");

		my $this_output_leftover_list_prefix = $mrg_output_prefix . "_leftover";
		my $this_output_leftover_cmap_file = $this_output_leftover_list_prefix . ".cmap";
		my $this_output_leftover_txt_file = $this_output_leftover_list_prefix . ".txt";
		my @comb_leftoverIds = (@$leftover_contigIds_1, @$leftover_contigIds_2);
		# check if any ids in @comb_leftoverIds:
		my $NcIds=@comb_leftoverIds;
		if ($NcIds !=0) {
			my $this_output_leftover_list = getArrayFileNames($mrg_output_prefix . "_contig", \@comb_leftoverIds, ".cmap");
			open TMP_FH, ">$this_output_leftover_txt_file" or dieLog ("ERROR: Unable to writet to file " + $this_output_leftover_txt_file + " - " + $!);
			foreach my $fn (@$this_output_leftover_list) { 
				print TMP_FH "$fn\n";
				push(@cleanup_file_list, $fn);
			} # foreach fn
			close (TMP_FH);
			$run_merge = 
			new BNG::refAlignerRun({binPath 	=> $refaligner,
							 if			=> $this_output_leftover_txt_file,
							 o			=> $this_output_leftover_list_prefix,
								 merge		=> 1,
							 f			=> 1,                   			# overwrite
							 stdout		=> 1,
							 stderr		=> 1	          					 
			} );	
			logMessage($logFH, $run_merge->getCMD());
			my ($outResults, $errResults, $job_status) = $run_merge->runCMD();
			exit $job_status if ($job_status !=0);	  	        
			# we cleanup now
			push(@cleanup_file_list, $this_output_leftover_list_prefix, ".idmap");
			push(@cleanup_file_list, $this_output_leftover_list_prefix, ".stdout");
		} else {
			#			# there are merges and no more leftover remains, then stop; this round is the last successful round
			$do_mrg = 0;
			$no_left_over = 1;		  	
		} # if NcIds
	} # if numMrgP
	# decide if to clean the temp files:
	if ($doDebug == 0) {
		foreach my $tmpfn (@cleanup_file_list) { 
			unlink($tmpfn) or warn "Could not unlink $tmpfn: $!";;
		}  # foreach tmpfn
		my $ntf = @cleanup_file_list; 
		logMessage($logFH, "\tClean up completed for " . $ntf . " files."); 	        
	} # if doDebug
	# update the merged pairs information
	$all_step1_mrg_pairs->{$rounds} = $MgrResult;
	# decide if to continue to next round:  
	if ($numMrgP == 0) {  
		# no merges are possible, then stop; previous round is the last successful round
		$do_mrg = 0;
	} # if numMrgP
	# check if maximum allowed merge rounds has been reached, then stop; (depending on the merge result above, either this round or previous round is the last successful round)
	if ($rounds >= $max_merge_rounds) {
		$do_mrg = 0; 
	} # if rounds
	logMessage($logFH, "Finished round " . $rounds . " merge of " . $this_mrg_prefix . " with " . $numMrgP . " pair merges");    	

	if ($do_mrg) {
		$rounds = $rounds+1;
	}  # if do_mrg
}  # while do_mrg
my $total_1Rounds = $rounds;
if ($rounds == 1)	{
	logMessage($logFH, "1.1 No Step1 merge was possible. Stop the merge program.");
	logMessage($logFH, "ERROR: 1.1. There are no merge possible between NGS and BioNano map in step1.");
	dieLog ("ERROR: 1.1. There are no merge possible between NGS and BioNano map in step1.");
} else {
	logMessage($logFH, "1.1. now we combine things from last successful merge into step 1. ");
} # if prev_mrg_prefix
my $step1_mrg = "";
if ($no_left_over)	{
	$step1_mrg = $this_mrg_prefix;	# we capture last successful merge
} else { 
	$step1_mrg = $prev_mrg_prefix  # we capture last successful merge
} # if no_left_over
my $step1_output = $scratchDir . $step1_mrg . "_combined";
my $to_be_merged_list = [$scratchDir . $step1_mrg . "_mrged.cmap"];
if (! $no_left_over) {
	$to_be_merged_list= [$to_be_merged_list->[0], $scratchDir .  $prev_mrg_prefix . "_leftover.cmap"]  # .1 and .2 if needed
}  # if no_left_over
logMessage($logFH, "1.1. Capture " . $step1_mrg . " as last successful merge and now merge");
my $run_merge11 = new BNG::refAlignerRun({binPath 	=> $refaligner,
				 i			=> $to_be_merged_list,
				 o			=> $step1_output,
				 merge		=> 1,
				 f			=> 1,                   			# overwrite
				 stdout		=> 1	          					 
} );
logMessage($logFH, $run_merge11->getCMD());
my ($outResults, $errResults, $job_status) = $run_merge11->runCMD();
exit $job_status if ($job_status !=0);
my @allMrgIds = ();
for (my $i=1; $i<=$total_1Rounds; $i++) { 
	push(@allMrgIds, @{$all_step1_mrg_pairs->{$i}->{mergePairs}->{ContigID1}});
	push(@allMrgIds, @{$all_step1_mrg_pairs->{$i}->{mergePairs}->{ContigID2}});
} # for i
my $naive_BN_contigs = uniqueIDs(getCMapIds($bng_cmap), \@allMrgIds); # these are the contigs never participated in merge
my $naive_NGS_contigs = uniqueIDs(getCMapIds($ngs_cmap), \@allMrgIds); # ditto
logMessage($logFH, "Naive BioNano Contig count =  " .  (scalar @$naive_BN_contigs));
logMessage($logFH, "Naive NGS Contig count =  " . (scalar @$naive_NGS_contigs));

# hybrid.contigs =! naive
my ($tmp_cmap, undef, $tmp_cmap_length) = readCMap($step1_output . ".cmap"); 
my @all_naive_ids = @$naive_NGS_contigs;
push(@all_naive_ids, @$naive_BN_contigs);
my $hybrid_cmap = getSubsetCMap($tmp_cmap, \@all_naive_ids, 'exclude');
my $hybrid_contigs = getCMapIds($hybrid_cmap);
logMessage($logFH, "Step1 Hybrid Contig count =  " . (scalar @$hybrid_contigs));

{ # write out cmap files: 
	$step1_hybrid_cmap_name = $scratchDir . "step1.hybrid.cmap";
	writeCMapFile($hybrid_cmap, $step1_hybrid_cmap_name);
	logMessage($logFH, "Step1 Hybrid Cmap exported as '" . $step1_hybrid_cmap_name . "'");

	my $naive_BN_cmap = getSubsetCMap($tmp_cmap, $naive_BN_contigs, 'include');
	my $step1_BN_cmap_name = $scratchDir . "step1.BN.naive.cmap";
	writeCMapFile($naive_BN_cmap, $step1_BN_cmap_name);
	logMessage($logFH, "Step1 BioNano Cmap exported as '" . $step1_BN_cmap_name . "'");

	my $naive_NGS_cmap = getSubsetCMap($tmp_cmap, $naive_NGS_contigs, 'include'); 
	my $step1_NGS_cmap_name = $scratchDir . "step1.NGS.naive.cmap";
	writeCMapFile($naive_NGS_cmap, $step1_NGS_cmap_name);
	logMessage($logFH, "Step1 NGS Cmap exported as '" . $step1_NGS_cmap_name . "'");
	# run some stats:
	my ($min, $max, $mean, $median, $n50value, $total) = getContigStat(@$tmp_cmap_length);
	logMessage($logFH, "Step1 merged Cmap N=" . $tmp_cmap->{nContigs} . "; N50=" . $n50value/1000000 . " Mb; total=" . $total/1000000 . " Mb");
	# keep record of which pairs of fragments were merged
	writeAllMrgPairs($all_step1_mrg_pairs, $scratchDir . "step1.merge.pairs.txt", $total_1Rounds, \@mrg_rounds_ids); 		    	
}        		
            			
# Step2. Now we perform merge between merged cmaps.
logMessage($logFH, "Step2. Now we perform merge between merged hybrid cmaps.");
my $all_step2_mrg_pairs = {};
$do_mrg = 1;
$rounds = 1;
$this_mrg_prefix = "";
$prev_mrg_prefix = "";
my $pairmerge_single_round_1file_run = $pairmerge_single_round_2files_run; # re-use the parameter sets 
$pairmerge_single_round_1file_run->rmParam("first"); # but remove the first option
# $pairmerge_single_round_1file_run->rmParam("i");  # i/o will be re-set every round:
# $pairmerge_single_round_1file_run->rmParam("o");
my $prev_mrg_round = "";
while ($do_mrg) {
	my $this_mrg_round = $mrg_rounds_ids[$rounds];  # 1 = A
	#my $prev_mrg_round = "";
	if ($rounds > 1) { 
		$prev_mrg_round = $mrg_rounds_ids[$rounds-1];  # 1=@ (does not exist), 2=A
	} # if rounds

	$this_mrg_prefix = "Mrg2" . $this_mrg_round;
	$prev_mrg_prefix = "Mrg2" . $prev_mrg_round;
	my 	$inputFile1 = $scratchDir . $prev_mrg_prefix . "_combined.cmap";
	if ($rounds == 1) {
	    $inputFile1 = $step1_hybrid_cmap_name
	}  # if rounds
	if (! -e $inputFile1) { 
		logMessage($logFH, "ERROR: Unable to open hybrid cmap file '" . $inputFile1 . "'");
		dieLog ("ERROR: Unable to open hybrid cmap file '" . $inputFile1 . "'");
	}  # if inputFile1

	logMessage($logFH, "\tRound " . $rounds . " " . $this_mrg_prefix . " started between '" . $inputFile1);
	my $mrg_output_prefix = $outputDir . '/' . $this_mrg_prefix; 
	$pairmerge_single_round_1file_run->setParams({i => $inputFile1, o => $mrg_output_prefix});
	my $cmd = $pairmerge_single_round_1file_run->getCMD();	
	logMessage($logFH,$cmd);
	my ($outResults, $errResults, $job_status) = $pairmerge_single_round_1file_run->runCMD();
	exit $job_status if ($job_status !=0);
	my @result_cmap_list = find_cmap_list_w_prefix($mrg_output_prefix);
	my @cleanup_file_list = @result_cmap_list; 
	push(@cleanup_file_list, $mrg_output_prefix . '.stdout');
	push(@cleanup_file_list, $mrg_output_prefix . '.align');  
#	 then we find out which contig has been merged by parsing out output to find out which pairs got merged.
	my ($numMrgP, $this_mrg_pairs)= parsingMrgStdOut($mrg_output_prefix . '.stdout');
	my $MgrResult = {};
	$MgrResult->{mergeCnt} =$numMrgP;
	$MgrResult->{mergePairs} = $this_mrg_pairs;
	if ($numMrgP == 0) {
		# there was no merge at all.
		logMessage($logFH, "No merge found for " . $mrg_output_prefix);  
	} else {			       
		# Now we need to merge cmap files into final products.
		# In this case we just lump everything into second round, not necessary to classify into merged and leftover
		# since we are not trying to prevent within group merge now.   
		my $this_output_mrged_list_prefix = $mrg_output_prefix . "_combined";
		my $this_output_mrged_txt_file = $this_output_mrged_list_prefix . ".txt";
		my $this_output_mrged_cmap_file = $this_output_mrged_list_prefix . ".cmap";
		#        write(result.cmap.list, file=this.output.mrged.txt.file, ncolumn=1)
		open TMP_FH, ">$this_output_mrged_txt_file" or dieLog ("ERROR: Unable to writet to file " + $this_output_mrged_txt_file + " - " + $!);
		foreach my $fn (@result_cmap_list) { 
			print TMP_FH "$fn\n";
		} # foreach fn
		close (TMP_FH);
		my $run_merge = 
			new BNG::refAlignerRun({binPath 	=> $refaligner,
							 if			=> $this_output_mrged_txt_file,
							 o			=> $this_output_mrged_list_prefix,
								 merge		=> 1,
							 f			=> 1,                   			# overwrite
							 stdout		=> 1,
							 stderr		=> 1	          					 
			} );		    
		logMessage($logFH, $run_merge->getCMD());
		my ($outResults, $errResults, $job_status) = $run_merge->runCMD();
		exit $job_status if ($job_status !=0);
		push(@cleanup_file_list, $this_output_mrged_list_prefix . ".idmap");
		push(@cleanup_file_list, $this_output_mrged_list_prefix . ".stdout");		    
		push(@cleanup_file_list, $this_output_mrged_txt_file);
		$all_step2_mrg_pairs->{$rounds} = $MgrResult;
	} #  $numMrgP != 0
	# decide if to clean the temp files:
	if ($doDebug == 0) {
		foreach my $tmpfn (@cleanup_file_list) { 
			unlink($tmpfn) or warn "Could not unlink $tmpfn: $!";;
		}  # foreach tmpfn
		my $ntf = @cleanup_file_list; 
		logMessage($logFH, "\tClean up completed for " . $ntf . " files."); 
	} # if doDebug
	# decide if to continue to next round:  
	if ($numMrgP == 0 || $rounds >= $max_merge_rounds) { 
		$do_mrg = 0; 
	} # if numMrgP
	logMessage($logFH, "Finished round " . $rounds . " merge of " . $this_mrg_prefix . " with " . $numMrgP . " pair merges");    	

	if ($do_mrg) {
		$rounds = $rounds+1;
	} # if do_mrg
}; #  while loop
my $total_2Rounds = $rounds;
my $step2_hybrid_cmap_name="";	
if($prev_mrg_round eq "")	{
	# there is no merge possible in step2, that is Mrg2A failed
	logMessage($logFH, "2.2 No Step2 merge was possible. We take step1 output as the final output.");
	# that is we copy step1.hybrid.cmap ($step1_hybrid_cmap_name) to step2.hyrbid.cmap ($step2_hybrid_cmap_name2)
	$step2_hybrid_cmap_name=$scratchDir . "step2.hybrid.cmap";
	copy($step1_hybrid_cmap_name, $step2_hybrid_cmap_name) or dieLog ("ERROR: Copy failed: from $step1_hybrid_cmap_name to $step2_hybrid_cmap_name: $!");            
} else	{
	# some merge was possible
	logMessage($logFH, "2.2. We combine things from last successful merge into step 2. ");
	my $step2_mrg = $prev_mrg_prefix;  # we capture last successful merge
	my $step2_output = $scratchDir . "step2.hybrid";
	$to_be_merged_list = $scratchDir . $prev_mrg_prefix . "_combined.cmap";

	my $run_merge21 = new BNG::refAlignerRun({binPath 	=> $refaligner,
							 i			=> $to_be_merged_list,
							 o			=> $step2_output,
								 merge		=> 1,
							 f			=> 1,                   			# overwrite
							 stdout		=> 1	          					 
		  } );
	logMessage($logFH, $run_merge21->getCMD());
	my ($outResults, $errResults, $job_status) = $run_merge21->runCMD();
	exit $job_status if ($job_status !=0);
	logMessage($logFH, "Step 2 completed.");
	# keep record of which pairs of fragments were merged
	my $a2p = $scratchDir . "step2.merge.pairs.txt";
	writeAllMrgPairs($all_step2_mrg_pairs, $a2p, $total_2Rounds, \@mrg_rounds_ids); 	     
	if (! -e $a2p) { 
		logMessage($logFH, "ERROR: Unable to write file '". $a2p);
		dieLog ("ERROR: Unable to write file '". $a2p);
	} # if a2p
	# run some stats:
	my ($tmp_cmap, undef, $tmp_cmap_length) = readCMap($step2_output . ".cmap");
	my ($min, $max, $mean, $median, $n50value, $total) = getContigStat(@$tmp_cmap_length);
	logMessage($logFH, "Step2 merged Cmap N=" . $tmp_cmap->{nContigs} . "; N50=" . $n50value/1000000 . " Mb; total=" . $total/1000000 . " Mb"); 
} # if prev.mrg.prefix
    
# At the end of step2, the super contigs include step2.hybrids + [NGS.naive + BN.naive]

###############
# Step 3. We align NGS to supercontigs. 
# not implemented yet.
close($logFH);
print "$0 finished successfully.";
