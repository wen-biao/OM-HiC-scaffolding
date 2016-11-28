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

use BNG::Utility;

	my $xmap_in 		 = $ARGV[0];
	my $ngs_cmap_in 	 = $ARGV[1];
	my $bng_cmap_in		 = $ARGV[2];
	my $sticky_xmap_fn 	 = $ARGV[3];
	my $filtered_ngs_fn 	 = $ARGV[4];
	my $filtered_bng_fn 	 = $ARGV[5];
	my $T_cutoff 		 = $ARGV[6];
	my $max_overhang 	 = $ARGV[7];
	my $orig_ngs_cmap_in = $ARGV[8];
	my $orig_bng_cmap_in = $ARGV[9];
	
#	print "orig_bng_cmap_in=$orig_bng_cmap_in\n";
#	$xmap_in = 'align1.xmap';
	my $xmap = readXMap($xmap_in);
    my ($ngs_cmap, undef, undef) = readCMap($ngs_cmap_in);
    my ($bn_cmap, undef, undef) = readCMap($bng_cmap_in);
    my ($orig_ngs_cmap, undef, undef) = readCMap($orig_ngs_cmap_in);
    my ($orig_bn_cmap, undef, undef) = readCMap($orig_bng_cmap_in);
# 	print "$orig_bn_cmap->{FileName}\n";
 	my @xmap_data_name = @{ $xmap->{dataName} }; 
		
## we shift 20bp to find out only labels extending outside alignment	
	my @refLeftLabCnt;
	my @refRightLabCnt;
	my @qryLeftLabCnt;
	my @qryRightLabCnt;	
	for (my $i=0; $i < $xmap->{totalHits}; $i++) {		
		my ($numL, $numR) = countOverangLables( $ngs_cmap, $xmap->{hits}->{RefContigID}->[$i], 
		                                        ($xmap->{hits}->{RefStartPos}->[$i]) - 20, 
		                                        ($xmap->{hits}->{RefEndPos}->[$i]) + 20       );
		push(@refLeftLabCnt, $numL);
		push(@refRightLabCnt, $numR);
#        print "ref: $numL, $numR; ";
		if (($xmap->{hits}->{Orientation}->[$i]) eq "+") { 
           ($numL, $numR) = countOverangLables( $bn_cmap, $xmap->{hits}->{QryContigID}->[$i], 
                                                ($xmap->{hits}->{QryStartPos}->[$i]) - 20, 
                                                ($xmap->{hits}->{QryEndPos}->[$i]) + 20      );			
		} else { 
           ($numR, $numL) = countOverangLables( $bn_cmap, $xmap->{hits}->{QryContigID}->[$i], 
                                                ($xmap->{hits}->{QryEndPos}->[$i]) - 20, 
                                                ($xmap->{hits}->{QryStartPos}->[$i]) + 20      );			
		}
		push(@qryLeftLabCnt, $numL);
		push(@qryRightLabCnt, $numR);
	}
	
# create sticky xmap:
    my $stickyXMap={};
	# - no header - $stickyXMap->{headers}
    $stickyXMap->{dataName} = $xmap->{dataName};
    $stickyXMap->{dataType} = $xmap->{dataType};
    $stickyXMap->{MAPFileVersion} = $xmap->{MAPFileVersion}; 
	my $numc = @xmap_data_name;
	for (my $i=0; $i<$numc; $i++) { 
		$stickyXMap->{hits}->{$xmap_data_name[$i]} = [];
	}  	
	my $tg = 0;
	my $stickyNGSContigs={};
	my $stickyBNContigs={};
	for (my $i=0; $i < $xmap->{totalHits}; $i++) {
	  if ((($xmap->{hits}->{Confidence}->[$i]) >= $T_cutoff) && (
	       ($refLeftLabCnt[$i] > $max_overhang &&
	        $qryLeftLabCnt[$i] > $max_overhang) ||
	       ($refRightLabCnt[$i] > $max_overhang &&
   	        $qryRightLabCnt[$i] > $max_overhang) ) ) { 
	        $tg++;
#	        print "($xmap->{hits}->{Confidence}->[$i]), $refLeftLabCnt[$i], $refRightLabCnt[$i],$qryLeftLabCnt[$i], $qryRightLabCnt[$i], $i\n";
		    for (my $m=0; $m<$numc; $m++) { 
			  my $a = $stickyXMap->{hits}->{$xmap_data_name[$m]};
        	  my $b = $xmap->{hits}->{$xmap_data_name[$m]};
        	  my $c = $b->[$i];
        	  push(@$a, $c);
        	  $stickyXMap->{hits}->{$xmap_data_name[$m]}=$a;
		    }
		    $stickyNGSContigs->{ $xmap->{hits}->{RefContigID}->[$i] } = 1;
		    $stickyBNContigs->{ $xmap->{hits}->{QryContigID}->[$i] } = 1;
        } 		
	}
#	print "tg=$tg\n";
	$stickyXMap->{totalHits} = $tg;
	writeXMapFile($sticky_xmap_fn, $stickyXMap, "noheader" );
	
	my $filtered_ngs_cmap = getSubsetCMap($orig_ngs_cmap, $stickyNGSContigs, 'exclude');
	writeCMapFile($filtered_ngs_cmap, $filtered_ngs_fn, 1);
	my $filtered_bionano_cmap = getSubsetCMap($orig_bn_cmap, $stickyBNContigs, 'exclude');
	writeCMapFile($filtered_bionano_cmap, $filtered_bng_fn, 1);
	print "Info: $0 finished successfully.\n\n";
	exit 0;
	
	
