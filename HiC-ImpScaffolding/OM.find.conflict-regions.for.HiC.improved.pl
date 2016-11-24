#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: OM.find.conflict.breaks.pl
#
#        USAGE: ./OM.find.conflict.breaks.pl  
#
#  DESCRIPTION: only for finding position of conflicting OM and Seq in initial alignment
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wen-Biao Jiao (), 
# ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
#      VERSION: 1.0
#      CREATED: 08/01/2016 02:48:20 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Getopt::Long;

##modified from pacbio.ass.err.check.by.bionano.pl

my ($seqKeyFile,$align1Dir,$outdir,$help);
GetOptions (
  'seqKey|k=s'    => \$seqKeyFile,
  'align1Dir|d=s' => \$align1Dir,
  'outDir|o=s'    => \$outdir,
  'help|h+'      => \$help,
);

my $usage="
Usage: $0 seqKeyFile align1_directory out_dir
Options:
  -seqKey|-k the in silico CMAP key file of sequence assembly
  -align1Dir|-d the directory path of initial alignment output 
  -outDir|-o output directory
  -help|-h print this help
";

die $usage if (defined $help);
die $usage if ((!$seqKeyFile) or (!$align1Dir));
 
if (! -e $seqKeyFile) {
  print "$seqKeyFile can't be found\n";
  die $usage; 	
}

print "##seqCmapKey\t$seqKeyFile\n";
print "##align1Dir\t$align1Dir\n";
print "##outdir\t$outdir\n";
print "\n";

system("mkdir -p $outdir");
##step0: hash the cmap and get conflicting seq/cmap
my $seqCmap="$align1Dir/align1/align1_r.cmap";  ###Don't use the original CMAP from $NGSfa, some close two nick sites were merged into one site after running RefAligner  
my $seqCmaphash=hashCmap($seqCmap);
my $seqCmapTrans=hashNGScmapKey($seqKeyFile);

my $seqAlnType="$align1Dir/assignAlignType/assignAlignType_r.cmap";
my $seqConflicts="$outdir/seq.assembly.conflicts.contigs";
my $seqConf=getConflictCtgs($seqAlnType,$seqCmaphash,$seqConflicts);

my $OMcmap="$align1Dir/align1/align1_q.cmap";
my $OMcmaphash=hashCmap($OMcmap);

my $OMalnType="$align1Dir/assignAlignType/assignAlignType_q.cmap";
my $OMconflicts="$outdir/OM.assembly.conflicts.contigs";
my $OMconf=getConflictCtgs($OMalnType,$OMcmaphash,$OMconflicts);

###step1:get initial breaks, and check those "2-breaks" alignment
my $xmap="$align1Dir/align1/align1.xmap";
my $flt=xmapFilter($xmap);
getBreaks($seqCmaphash,$OMcmaphash,$seqConf,$OMconf,$xmap,$seqCmapTrans,$outdir,$flt);

###step2:get breakpoints, and check them whether they are supported "double-contig" and HiRise 
my $seqBreaksM="$outdir/seq.assembly.conflicts.breaks";
my $seqBreaks2="$outdir/seq.assembly.conflicts.breaks2";  
mergeBreaks($seqBreaksM,$seqCmap,$seqBreaks2,$seqCmaphash,0);
my $OMbreaksM="$outdir/OM.assembly.conflicts.breaks";
my $OMbreaks2="$outdir/OM.assembly.conflicts.breaks2";
mergeBreaks($OMbreaksM,$OMcmap,$OMbreaks2,$OMcmaphash,1);


my $seqBreaks3="$outdir/seq.assembly.conflicts.breaks3";
my $OMbreaks3="$outdir/OM.assembly.conflicts.breaks3";
add2supportInf($OMbreaks2,$seqBreaks2,$seqBreaks3,$OMbreaks3);

my $seqBreaks3f="$outdir/seq.assembly.conflicts.breaks3.flt";
my $OMbreaks3f="$outdir/OM.assembly.conflicts.breaks3.flt";
breaksFlt($seqBreaks3,$OMbreaks3,$seqBreaks3f,$OMbreaks3f);

sub getConflictCtgs {
 my ($assCmap,$cmap,$out)=@_;
 my %cmap=%{$cmap};
 open IN,$assCmap;
 my %confCmap;
 open OUT,">$out";
 while (<IN>) {
   chomp;
  	if (/exclude/) {
  	  my ($exc)=$_=~m/exclude: (.*)/;
  	  my @exc=(split /\,/,$exc);
  	  foreach my $id (@exc) {
         $confCmap{$id}=$cmap{$id};
         print OUT "$id\t$cmap{$id}{'leng'}\t$cmap{$id}{'ns'}\n";
  	  }	
  	}	
  }
  close OUT;close IN;
  return \%confCmap;
}


sub getBreaks {
  my ($seqCmap,$OMcmap,$seqConf,$OMconf,$xmap,$seqCmapTrans,$oudir,$flt)=@_;
  
  my %seqCmap=%{$seqCmap};
  my %OMcmap=%{$OMcmap};  
  my %seqConf=%{$seqConf};  
  my %OMconf=%{$OMconf};
  
  my %seqCmapTrans=%{$seqCmapTrans};  
  my $maxHang=2;
  my %flt=%{$flt};    
  
  my $seqBreaks="$outdir/seq.assembly.conflicts.breaks";
  my $pSeqBreaks="$outdir/potential.seq.assembly.conflicts.breaks";
  my $seqBreaksM="$outdir/seq.assembly.conflicts.more.breaks";      
  open OUT1,">$seqBreaks";  
  print OUT1 "#cmapID\tAssID\tLength\tNumLabels\tAlnEndLabel\tAlnEndPos\tNextUnAlnedLabel\tNextUnAlnedPos\tAlignedBNG\t1/2\n";
  open OUT2,">$pSeqBreaks";
  print OUT2 "#cmapID\tAssID\tLength\tNumLabels\tAlnEndLabel\tAlnEndPos\tNextUnAlnedLabel\tNextUnAlnedPos\tAlignedBNG\t1/2\n";  
  open OUT3,">$seqBreaksM";
  print OUT3 "#cmapID\tAssID\tLength\tNumLabels\tAlnEndLabel\tAlnEndPos\tNextUnAlnedLabel\tNextUnAlnedPos\tAlignedBNG\t1/2\n";
    
  my $OMbreaks="$outdir/OM.assembly.conflicts.breaks";
  my $pOMbreaks="$outdir/potential.OM.assembly.conflicts.breaks";
  my $OMbreaksM="$outdir/OM.assembly.conflicts.more.breaks";
  open OUT4,">$OMbreaks";  
  print OUT4 "#cmapID\tAssID\tLength\tNumLabels\tAlnEndLabel\tAlnEndPos\tNextUnAlnedLabel\tNextUnAlnedPos\tAlignedBNG\t1/2\n";
  open OUT5,">$pOMbreaks";
  print OUT5 "#cmapID\tAssID\tLength\tNumLabels\tAlnEndLabel\tAlnEndPos\tNextUnAlnedLabel\tNextUnAlnedPos\tAlignedBNG\t1/2\n";  
  open OUT6,">$OMbreaksM";
  print OUT6 "#cmapID\tAssID\tLength\tNumLabels\tAlnEndLabel\tAlnEndPos\tNextUnAlnedLabel\tNextUnAlnedPos\tAlignedBNG\t1/2\n";
      
  open IN,"grep -v '#' $xmap|sort -k3,3n -k6,6n -k7,7nr|";
  open FLT,">$outdir/xmap.flt";
  my ($ins,$wrong)=(0,0);
  while (<IN>) {
  	chomp;
  	my @t=split /\t/;
  	my @aln;
  	my $line=$_;
  	my ($xmapID,$queID,$refID,$queStartPos,$queEndPos,$refStartPos,$refEndPos,$orientation,$confidence,$hitEnum,$queLen,$refLen,$labelChannel,$aln)=split /\t/;
  	while ($aln=~m/\((\d+)\,(\d+)\)/g) {
      my ($rsite,$qsite)=($1,$2);	
      push @aln,[($rsite,$qsite)];  
    } 
    if ($flt{$t[2].$t[5].$t[6]}) {
        print "inside alignment Seq:$seqCmapTrans{$t[2]}\ninside Alignment: $_\n";
        $ins++;
        if ($seqConf{$t[2]} and $OMconf{$t[1]}) {
          print "##wrong conflict(not unique alignment)\t$line\n";	
          $wrong++;
        }
      	next;
    } 
    print FLT "$line\n";      	        
  	my $seqLeftHang=$aln[0][0]-1;
  	$seqLeftHang=$aln[1][0]-1 if  ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;	
  	my $seqRightHang=$seqCmap{$t[2]}{"ns"}-$aln[$#aln][0]; 
    $seqRightHang=$seqCmap{$t[2]}{"ns"}-$aln[$#aln-1][0]  if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) ); 
    	 	    	  
  	my $OMleftHang=$aln[0][1]-1;
  	$OMleftHang=$aln[1][1]-1 if  ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;  	  
  	my $OMrightHang=$OMcmap{$t[1]}{"ns"}-$aln[$#aln][1];
    $OMrightHang=$OMcmap{$t[1]}{"ns"}-$aln[$#aln-1][1] if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) );  	  	    	    	
  	if ($orientation eq "-") {
  	  $OMleftHang=$OMcmap{$t[1]}{"ns"}-$aln[0][1];
  	  $OMrightHang=$aln[$#aln][1]-1;
  	  
  	  $OMleftHang=$OMcmap{$t[1]}{"ns"}-$aln[1][1] if ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;
  	  $OMrightHang=$aln[$#aln-1][1]-1 if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) );	
  	}    	
  	#next if (($aln[0][0]<=$maxHang) and ($seqConf{$refID}{'ns'}-$aln[$#aln][0]<=$maxHang) ) ; ##good Seq-inside alignment;
  	#next if (($OMrightHang>$maxHang) and ($OMleftHang>$maxHang) and ($#aln+1<3*$maxHang)); ###filter those partial and not unique alignment
  	next if ( ($seqRightHang<=$maxHang) and ($seqLeftHang<=$maxHang));  	  	  	
  	next if ( ($OMleftHang<=$maxHang) and ($OMrightHang<=$maxHang) ) ;
  	next if ( ($seqRightHang<=$maxHang) and ($OMleftHang<=$maxHang) );
  	next if ( ($OMrightHang<=$maxHang) and ($seqLeftHang<=$maxHang) );
  	
  	  	
    if ( ($seqRightHang>$maxHang) and ($OMrightHang>$maxHang) and ( ($seqLeftHang<=$maxHang) or ($OMleftHang<=$maxHang)  )  )  {  ###seq-OM righ part can't be aligned  
      my $br=$aln[$#aln][0];
  	  $br=$aln[$#aln-1][0] if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) );
  	  my $nbr=$br+1;
  	  my $tbr=$aln[$#aln][1];
  	  $tbr=$aln[$#aln-1][1] if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) );
  	  my $ntbr=$tbr+1;
  	  $ntbr=$tbr-1 if ($orientation eq "-");  	  	      	
  	  if ($seqConf{$t[2]} and ($OMconf{$t[1]})) {
  	    print OUT1 "$t[2]\t$seqCmapTrans{$t[2]}\t$seqCmap{$t[2]}{'leng'}\t$seqCmap{$t[2]}{'ns'}\t$br\t$seqCmap{$t[2]}{$br}\t$nbr\t$seqCmap{$refID}{$nbr}\t$t[1]\t$tbr\t$ntbr\t1\n";  	  	      	      	
  	  }else {
  	    print OUT2 "$t[2]\t$seqCmapTrans{$t[2]}\t$seqCmap{$t[2]}{'leng'}\t$seqCmap{$t[2]}{'ns'}\t$br\t$seqCmap{$t[2]}{$br}\t$nbr\t$seqCmap{$refID}{$nbr}\t$t[1]\t$tbr\t$ntbr\t1\n";  	  	      	      	
  	  }
  	  print OUT3 "$t[2]\t$seqCmapTrans{$t[2]}\t$seqCmap{$t[2]}{'leng'}\t$seqCmap{$t[2]}{'ns'}\t$br\t$seqCmap{$t[2]}{$br}\t$nbr\t$seqCmap{$refID}{$nbr}\t$t[1]\t$tbr\t$ntbr\t1\n";
  	  
  	  $br=$aln[$#aln][1];
  	  $br=$aln[$#aln-1][1] if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) );
  	  $nbr=$br+1;
  	  $nbr=$br-1 if ($orientation eq "-");
  	  $tbr=$aln[$#aln][0];
  	  $tbr=$aln[$#aln-1][0] if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) );
  	  $ntbr=$tbr+1;  	    	  	
  	  if ($seqConf{$t[2]} and $OMconf{$t[1]}) {
  	    print OUT4 "$t[1]\t$t[1]\t$OMcmap{$t[1]}{'leng'}\t$OMcmap{$t[1]}{'ns'}\t$br\t$OMcmap{$t[1]}{$br}\t$nbr\t$OMcmap{$t[1]}{$nbr}\t$t[2]\t$tbr\t$ntbr\t1\n";	
  	  }else {
  	    print OUT5 "$t[1]\t$t[1]\t$OMcmap{$t[1]}{'leng'}\t$OMcmap{$t[1]}{'ns'}\t$br\t$OMcmap{$t[1]}{$br}\t$nbr\t$OMcmap{$t[1]}{$nbr}\t$t[2]\t$tbr\t$ntbr\t1\n";	
  	  }  	  	   	  	  	  	    
  	  print OUT6 "$t[1]\t$t[1]\t$OMcmap{$t[1]}{'leng'}\t$OMcmap{$t[1]}{'ns'}\t$br\t$OMcmap{$t[1]}{$br}\t$nbr\t$OMcmap{$t[1]}{$nbr}\t$t[2]\t$tbr\t$ntbr\t1\n";
  	    	  	    	  	  
    }elsif ( ($seqLeftHang>$maxHang) and ($OMleftHang>$maxHang) and ( ($seqRightHang<=$maxHang) or ($OMrightHang<=$maxHang)  ) ) {
      my $br=$aln[0][0];
  	  $br=$aln[1][0] if  ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;
  	  my $nbr=$br-1;  	  	   	  	    
  	  my $tbr=$aln[0][1];  	  	    
  	  $tbr=$aln[1][1] if  ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;
  	  my $ntbr=$tbr-1;
  	  $ntbr=$tbr+1 if ($orientation eq "-");  	  	     	  	      	  	  
  	  if ($seqConf{$t[2]} and $OMconf{$t[1]} ) {
  	    print OUT1 "$t[2]\t$seqCmapTrans{$t[2]}\t$seqCmap{$t[2]}{'leng'}\t$seqCmap{$t[2]}{'ns'}\t$br\t$seqCmap{$t[2]}{$br}\t$nbr\t$seqCmap{$refID}{$nbr}\t$t[1]\t$tbr\t$ntbr\t1\n";  	
  	  }else {
  	    print OUT2 "$t[2]\t$seqCmapTrans{$t[2]}\t$seqCmap{$t[2]}{'leng'}\t$seqCmap{$t[2]}{'ns'}\t$br\t$seqCmap{$t[2]}{$br}\t$nbr\t$seqCmap{$refID}{$nbr}\t$t[1]\t$tbr\t$ntbr\t1\n";	
  	  }  	  	  	  	  	  	    
  	  print OUT3 "$t[2]\t$seqCmapTrans{$t[2]}\t$seqCmap{$t[2]}{'leng'}\t$seqCmap{$t[2]}{'ns'}\t$br\t$seqCmap{$t[2]}{$br}\t$nbr\t$seqCmap{$refID}{$nbr}\t$t[1]\t$tbr\t$ntbr\t1\n";
  	  	  
  	  $br=$aln[0][1];
  	  $br=$aln[1][1] if  ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;
  	  $nbr=$br-1;
  	  $nbr=$br+1 if ($orientation eq "+");
  	  $tbr=$aln[0][0];
  	  $tbr=$aln[1][0] if  ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;
  	  $ntbr=$tbr-1;
  	  if ($seqConf{$t[2]} and $OMconf{$t[1]}) {
  	    print OUT4 "$t[1]\t$t[1]\t$OMcmap{$t[1]}{'leng'}\t$OMcmap{$t[1]}{'ns'}\t$br\t$OMcmap{$t[1]}{$br}\t$nbr\t$OMcmap{$t[1]}{$nbr}\t$t[2]\t$tbr\t$ntbr\t1\n";	
  	  }else {
  	    print OUT5 "$t[1]\t$t[1]\t$OMcmap{$t[1]}{'leng'}\t$OMcmap{$t[1]}{'ns'}\t$br\t$OMcmap{$t[1]}{$br}\t$nbr\t$OMcmap{$t[1]}{$nbr}\t$t[2]\t$tbr\t$ntbr\t1\n";	
  	  }  	  	   	  	  	  	    
  	  print OUT6 "$t[1]\t$t[1]\t$OMcmap{$t[1]}{'leng'}\t$OMcmap{$t[1]}{'ns'}\t$br\t$OMcmap{$t[1]}{$br}\t$nbr\t$OMcmap{$t[1]}{$nbr}\t$t[2]\t$tbr\t$ntbr\t1\n";
  	  	  
    }elsif ( ($seqRightHang>$maxHang) and ($OMrightHang>$maxHang) and ($seqLeftHang>$maxHang) and ($OMleftHang>$maxHang)) {
  	  my $br1=$aln[0][0];
  	  $br1=$aln[1][0] if  ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;
  	  my $nbr1=$br1-1;  	  	  	  	  	  	
  	  my $tbr1=$aln[0][1];
  	  $tbr1=$aln[1][1]  if  ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;
  	  my $ntbr1=$tbr1-1;  
  	  $ntbr1=$tbr1+1 if ($orientation eq "-");
  	  	  	  	  	  	  	
  	  my $br2=$aln[$#aln][0];  	  	  	
  	  $br2=$aln[$#aln-1][0] if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) );  	  	  	
  	  my $nbr2=$br2+1;  	  	  	  	  	  
  	  my $tbr2=$aln[$#aln][1];
  	  $tbr2=$aln[$#aln-1][1] if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) );
  	  my $ntbr2=$tbr2+1;
  	  $ntbr2=$tbr2-1 if ($orientation eq "-");
  	    	  	  	
  	  if ($seqConf{$t[2]} and $OMconf{$t[1]} ) {
  	    print OUT1 "$t[2]\t$seqCmapTrans{$t[2]}\t$seqCmap{$t[2]}{'leng'}\t$seqCmap{$t[2]}{'ns'}\t$br1\t$seqCmap{$t[2]}{$br1}\t$nbr1\t$seqCmap{$refID}{$nbr1}\t$t[1]\t$tbr1\t$ntbr1\t2\n";
  	    print OUT1 "$t[2]\t$seqCmapTrans{$t[2]}\t$seqCmap{$t[2]}{'leng'}\t$seqCmap{$t[2]}{'ns'}\t$br2\t$seqCmap{$t[2]}{$br2}\t$nbr2\t$seqCmap{$refID}{$nbr2}\t$t[1]\t$tbr2\t$ntbr2\t2\n";	
  	  }else {
  	    print OUT2 "$t[2]\t$seqCmapTrans{$t[2]}\t$seqCmap{$t[2]}{'leng'}\t$seqCmap{$t[2]}{'ns'}\t$br1\t$seqCmap{$t[2]}{$br1}\t$nbr1\t$seqCmap{$refID}{$nbr1}\t$t[1]\t$tbr1\t$ntbr1\t2\n";
  	    print OUT2 "$t[2]\t$seqCmapTrans{$t[2]}\t$seqCmap{$t[2]}{'leng'}\t$seqCmap{$t[2]}{'ns'}\t$br2\t$seqCmap{$t[2]}{$br2}\t$nbr2\t$seqCmap{$refID}{$nbr2}\t$t[1]\t$tbr2\t$ntbr2\t2\n";
  	  }  	  	  	  	  	  	  
  	  print OUT3 "$t[2]\t$seqCmapTrans{$t[2]}\t$seqCmap{$t[2]}{'leng'}\t$seqCmap{$t[2]}{'ns'}\t$br1\t$seqCmap{$t[2]}{$br1}\t$nbr1\t$seqCmap{$refID}{$nbr1}\t$t[1]\t$tbr1\t$ntbr1\t2\n";
  	  print OUT3 "$t[2]\t$seqCmapTrans{$t[2]}\t$seqCmap{$t[2]}{'leng'}\t$seqCmap{$t[2]}{'ns'}\t$br2\t$seqCmap{$t[2]}{$br2}\t$nbr2\t$seqCmap{$refID}{$nbr2}\t$t[1]\t$tbr2\t$ntbr2\t2\n";
  	  	 
  	  $br1=$aln[0][1];
  	  $br1=$aln[1][1] if  ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;
  	  $nbr1=$br1-1;  	  	  	  	  	  	
  	  $nbr1=$br1+1 if ($orientation eq "-");  	  
  	  $tbr1=$aln[0][0];
  	  $tbr1=$aln[1][0]  if  ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;  	  	  
  	  $ntbr1=$tbr1-1;    	  
  	  	  	  	  	  	  	  	
  	  $br2=$aln[$#aln][1];  	  	  	
  	  $br2=$aln[$#aln-1][1] if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) );  	  	  	
  	  $nbr2=$br2+1;
  	  $nbr2=$br2-1 if ($orientation eq "-");  	  	  	  	  	  
  	  $tbr2=$aln[$#aln][0];
  	  $tbr2=$aln[$#aln-1][0] if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) );
  	  $ntbr2=$tbr2+1;
  	  if ( $seqConf{$t[2]} and $OMconf{$t[1]} ) {
  	    print OUT4 "$t[1]\t$t[1]\t$OMcmap{$t[1]}{'leng'}\t$OMcmap{$t[1]}{'ns'}\t$br1\t$OMcmap{$t[1]}{$br1}\t$nbr1\t$OMcmap{$t[1]}{$nbr1}\t$t[2]\t$tbr1\t$ntbr1\t2\n";  	  	    
  	    print OUT4 "$t[1]\t$t[1]\t$OMcmap{$t[1]}{'leng'}\t$OMcmap{$t[1]}{'ns'}\t$br2\t$OMcmap{$t[1]}{$br2}\t$nbr2\t$OMcmap{$t[1]}{$nbr2}\t$t[2]\t$tbr2\t$ntbr2\t2\n";	
  	  }else {
  	    print OUT5 "$t[1]\t$t[1]\t$OMcmap{$t[1]}{'leng'}\t$OMcmap{$t[1]}{'ns'}\t$br1\t$OMcmap{$t[1]}{$br1}\t$nbr1\t$OMcmap{$t[1]}{$nbr1}\t$t[2]\t$tbr1\t$ntbr1\t2\n";  	  	    
  	    print OUT5 "$t[1]\t$t[1]\t$OMcmap{$t[1]}{'leng'}\t$OMcmap{$t[1]}{'ns'}\t$br2\t$OMcmap{$t[1]}{$br2}\t$nbr2\t$OMcmap{$t[1]}{$nbr2}\t$t[2]\t$tbr2\t$ntbr2\t2\n";	
  	  }  	  	    	  	  	  	  
  	  print OUT6 "$t[1]\t$t[1]\t$OMcmap{$t[1]}{'leng'}\t$OMcmap{$t[1]}{'ns'}\t$br1\t$OMcmap{$t[1]}{$br1}\t$nbr1\t$OMcmap{$t[1]}{$nbr1}\t$t[2]\t$tbr1\t$ntbr1\t2\n";  	  	    
  	  print OUT6 "$t[1]\t$t[1]\t$OMcmap{$t[1]}{'leng'}\t$OMcmap{$t[1]}{'ns'}\t$br2\t$OMcmap{$t[1]}{$br2}\t$nbr2\t$OMcmap{$t[1]}{$nbr2}\t$t[2]\t$tbr2\t$ntbr2\t2\n";
  	      	  	    	  	     
  	}else {
  	  print "No categories: $line\n";	
  	}  	  	   	  	  	  	 
  }  
  close OUT1;close OUT2;close OUT3;
  close OUT4;close OUT5;close OUT6;
  close FLT; close IN;    
  print "total inside alignment\t$ins\n";
  print "total wrong conflicts\t$wrong\n";
}



sub mergeBreaks {
  my ($break,$cmap,$break2,$hashCmap,$isBNG)=@_;    
  open IN,"grep -v '#' $break|sort -k1,1n -k5,5n|";  
  #cmapID\tAssID\tLength\tNumLabels\tAlnEndLabel\tAlnEndPos\tNextUnAlnedLabel\tNextUnAlnedPos\tAlignedBNG\t1/2\n
  my ($cmapID,$seqID,$n,$prep,$pos,$leng,$target,$nextp,$nextpos,$tbr,$ntbr,$flag);  
  my $br;my $bp;
  open OUT,">$break2";
  my %breaks;my @breaks;
  my %hCmap=%{$hashCmap};
  my $m=0;
  my $f;
  while (<IN>) {
  	chomp;
  	my @t=split /\t/;  	  	
  	if (!$cmapID) {
  	  ($cmapID,$seqID,$leng,$n,$prep,$pos,$nextp,$nextpos,$target,$tbr,$ntbr,$flag)=@t;
  	  $target=$target.",".$prep.",".$pos."\,$tbr\,$ntbr";
  	  $br=($prep+$nextp)/2;
  	  $bp=int( ($pos+$nextpos)/2 ) ;
  	  $m=1;  
  	}elsif ($cmapID ne $t[0]) {  	  
  	  print OUT "$cmapID\t$seqID\t$leng\t$n\t$br\t$bp\t$target\t$flag\t$m\n";  ####to be modified, breakpoint=the mid of ajacent broken labels
  	  push @breaks,[($br,$bp)];
  	  $breaks{$cmapID}=[@breaks];
  	  @breaks=();
  	  ($cmapID,$seqID,$leng,$n,$prep,$pos,$nextp,$nextpos,$target,$tbr,$ntbr,$flag)=@t;
  	  $target=$target.",".$prep.",".$pos."\,$tbr\,$ntbr"; 	  
  	  $br=($prep+$nextp)/2;
  	  $bp=int( ($pos+$nextpos)/2 ) ;
  	  $m=1;  	  
  	}else {
  	  if ($t[4]-$prep<=4) {
  	    $br=($t[4]+$prep)/2;
  	    $bp=int(($t[5]+$pos)/2);  	    
  	    $m++;
  	    $target=$target.";".$t[8].",".$t[4].",".$t[5]."\,$t[9]\,$t[10]";
  	  }else {
  	  	print OUT "$cmapID\t$seqID\t$leng\t$n\t$br\t$bp\t$target\t$flag\t$m\n";
  	  	push @breaks,[($br,$bp)];
  	  	($cmapID,$seqID,$leng,$n,$prep,$pos,$nextp,$nextpos,$target,$tbr,$ntbr,$flag)=@t;
  	  	$target=$target.",".$prep.",".$pos."\,$tbr\,$ntbr";;
  	    $br=($prep+$nextp)/2;
  	    $bp=int(($pos+$nextpos)/2);  	   
  	    $m=1;  	    
  	  }
  	}
  }
  print OUT "$cmapID\t$seqID\t$leng\t$n\t$br\t$bp\t$target\t$flag\t$m\n";
  push @breaks,[($br,$bp)];
  $breaks{$cmapID}=[@breaks];
  close OUT;close IN;
  @breaks=();
  
}

sub xmapFilter {
  my ($xmap)=@_;
  open IN,"grep -v '#' $xmap|sort -k3,3n -k6,6n -k7,7nr|";
  my ($ref,$start,$end);
  my %insideCtg;
  open OUT,">tmp";
  while (<IN>) {
    chomp;
    my @t=split /\t/;
    if ($t[3]>$t[4]) {
      ($t[3],$t[4])=($t[4],$t[3]);
    }
    my $l=join("\t",@t);
    print OUT "$l\n";    
    if (!$ref) {
      ($ref,$start,$end)=@t[2,5,6];      
    }else {
      if ( ($t[2] eq $ref) and ($t[5]>=$start) and ($t[6]<=$end) ) {
      	$insideCtg{$t[2].$t[5].$t[6]}=1;
      	next;
      }else {
      	($ref,$start,$end)=@t[2,5,6];
      } 	         
    }
  } 
  close IN;
  close OUT;
  
  open IN,"grep -v '#' tmp|sort -k2,2n -k4,4n -k5,5nr|";
  my $que;
  while (<IN>) {
    chomp;
    my @t=split /\t/;
    if (!$que) {
      ($que,$start,$end)=@t[1,3,4];      
    }else {
      if ( ($t[1] eq $ref) and ($t[3]>=$start) and ($t[4]<=$end) ) {
      	$insideCtg{$t[2].$t[5].$t[6]}=1;
      	next;
      }else {
      	($ref,$start,$end)=@t[1,3,4];
      } 	         
    }
  } 
  close IN;
  system("rm tmp");
  return \%insideCtg;
}


sub hashCmap {
  my ($in)=@_;
  open IN,$in;
  my %cmap;
  while (<IN>) {
    chomp;
    if (/^\d/) {
      my ($id,$leng,$ns,$site,$pos)=(split /\t/)[0,1,2,3,5];
      $cmap{$id}{"leng"}=$leng;
      $cmap{$id}{"ns"}=$ns;
      $cmap{$id}{$site}=$pos;
    }
  }
  close IN;
  return \%cmap;
}

sub hashNGScmapKey {
  my ($cmapKey)=@_;  
  open IN,$cmapKey;
  my %cmapID;
  while (<IN>) {
    chomp;
    next if (!/^\d/);
    my @t=split /\t/;
    #print "$t[0]\t$t[2]\n";exit;
    $cmapID{$t[0]}=$t[1];
  } 
  close IN;
  return \%cmapID; 	
}

sub breaksFlt {
  my ($seqbreaks,$BNGbreaks,$out1,$out2)=@_;
  open IN,$seqbreaks;
  my ($twoBr,$twoBNG,$goodBNG,$finalBr)=(0,0,0,0,0,0);
  open OUT1,">$out1";
  while (<IN>) {
  	chomp;
  	my @t=split /\t/;
  	if (/#/) {
  		print OUT1 "$_\n" ;
  		next;
  	}
  	if ($t[7]==2) {
  	  $twoBr++;
  	}
  	if ($t[8]==2) {
  	  $twoBNG++;
  	}
  	

  	if ( ($t[8]==2)  ) {
  	  $finalBr++;
  	  print OUT1 "$_\tSeqMiss\n";
  	  if ($t[6]=~m/2seq/) {
  	  	print "$seqbreaks need check:$_\n";  	  	
  	  }
  	}else {
  	  if ($t[6]=~m/2seq/) {
  	  	print "exclude BNG mis-join supported by at least 2 seq\t$_\n";
  	  }else {
  	  	print OUT1 "$_\tUndetermined\n";  	  	
  	  }
  	}   
  }
  close OUT1;close IN;
  print "statistics:----\n";
  print "double breakpoints: $twoBr\n";
  print "two BNG supported:  $twoBNG\n";
  print "Total determined seq mis-joins: $finalBr\n";
  
  open OUT2,">$out2";
  open IN,"$BNGbreaks";
  print "xxxxx\n";
  while (<IN>) {
  	chomp;
  	if (/#/) {
  		print OUT2 "$_\n" ;
  		next;
  	}
  	my @t=split /\t/;
  	if ((!/SeqID/) and (!/2BNG/))  {
  	  print OUT2 "$_\n";
  	}
  }
  close OUT2;close IN;
  
}


sub add2supportInf {
  my ($BNGbreaks,$seqbreaks,$out1,$out2)=@_;
  open IN,$BNGbreaks;
  my %BNGdouble;
  while (<IN>) {
  	chomp;
  	my @t=split /\t/;
  	next if ($t[8]!=2);
    my @tar=split /\;/,$t[6];
    foreach my $target (@tar) {
      my @target=split /\,/,$target;
      $BNGdouble{$target[0]}{$target[1]}="2Seq";
    } 	
  }
  open IN,$seqbreaks;
  my %seqdouble;
  open OUT1,">$out1\n";
  while (<IN>) {
  	chomp;
  	my @t=split /\t/;
  	my @tar=split /\;/,$t[6];
    foreach my $target (@tar) {
      my @target=split /\,/,$target;	
      if ($BNGdouble{$t[0]}{$target[3]}) {
      	$target=$target.",2seq";
      }
      if ($t[8]==2) {
      	$seqdouble{$target[0]}{$target[1]}="2BNG";
      }
    } 
    $t[6]=join("\;",@tar);
    my $l=join("\t",@t);
    print OUT1 "$l\n";   
  }
  
  open OUT2,">$out2";
  open IN,$BNGbreaks; 
  while (<IN>) {
  	chomp;
  	my @t=split /\t/;
  	my @tar=split /\;/,$t[6];
    foreach my $target (@tar) {
      my @target=split /\,/,$target;	
      if ($seqdouble{$t[0]}{$target[3]}) {
      	$target=$target.",2BNG";
      }      
    } 
    $t[6]=join("\;",@tar);
    my $l=join("\t",@t);
    print OUT2 "$l\n";   
  }
}
