#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: pacbio.ass.err.pbcr-falcon.check.by.bionano.pl
#
#        USAGE: ./pacbio.ass.err.pbcr-falcon.check.by.bionano.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wen-Biao Jiao (), 
# ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
#      VERSION: 1.0
#      CREATED: 05/06/2016 10:13:33 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use File::Basename;

my ($pSeqbreaks,$pBNGbreaks,$fSeqbreaks,$fBNGbreaks,$pbcrXmap,$falconXmap,$pbcrKey,$falconKey,$outdir)=@ARGV;


my $pBNGbreaks3_5="$outdir/pbcr.BNG.breaks3.5";
my $pSeqbreaks3_5="$outdir/pbcr.seq.breaks3.5";
add2supportInf($pBNGbreaks,$pSeqbreaks,$pSeqbreaks3_5,$pBNGbreaks3_5);

my $fBNGbreaks3_5="$outdir/falcon.BNG.breaks3.5";
my $fSeqbreaks3_5="$outdir/falcon.seq.breaks3.5";
add2supportInf($fBNGbreaks,$fSeqbreaks,$fSeqbreaks3_5,$fBNGbreaks3_5);

my $pKey=hashNGScmapKey($pbcrKey);
my $fKey=hashNGScmapKey($falconKey);

my $pbcrFlt=xmapFilter($pbcrXmap);
my $pbcrBNGaln=getGoodBNGAln($pbcrXmap,$pbcrFlt);

my $falconFlt=xmapFilter($falconXmap);
my $falconBNGaln=getGoodBNGAln($falconXmap,$falconFlt);


###get PBcR breaks4
my $pBNGbreaks4="$outdir/pbcr.BNG.falcon-BNG-aln.check";
checkBNGbreaks($pBNGbreaks3_5,$falconBNGaln,$fKey,$pBNGbreaks4,"PBcR");
my $pSeqbreaks4="$outdir/pbcr.seq.falcon-BNG-aln.check";
checkSeqConfBNGAln($pSeqbreaks3_5,$falconBNGaln,$fKey,$pSeqbreaks4);


###get FALCON breaks4
my $fBNGbreaks4="$outdir/falcon.BNG.pbcr-BNG-aln.check";
checkBNGbreaks($fBNGbreaks3_5,$pbcrBNGaln,$pKey,$fBNGbreaks4,"FALCON");
my $fSeqbreaks4="$outdir/falcon.seq.pbcr-BNG-aln.check";
checkSeqConfBNGAln($fSeqbreaks3_5,$pbcrBNGaln,$pKey,$fSeqbreaks4);


###get merged.breaks4
my $pMergedBNGbreaks4="$outdir/pbcr.BNG.merged.breaks4";
my $fMergedBNGbreaks4="$outdir/falcon.BNG.merged.breaks4";
mergeBNGbreaks($pBNGbreaks4,$fBNGbreaks4,$pMergedBNGbreaks4,$fMergedBNGbreaks4);

###get breaks5 using breaks4 and HiRise breaks

###get breaks4.flt
my $pSeqbreaks4flt="$outdir/pbcr.seq.breaks4.flt";
my $pBNGbreaks4flt="$outdir/pbcr.BNG.merged.breaks4.flt";
breaksFlt($pSeqbreaks4,$pMergedBNGbreaks4,$pSeqbreaks4flt,$pBNGbreaks4flt);
system("cp $pSeqbreaks4flt $outdir/pbcr.seq.breaks5");
system("cp $pBNGbreaks4flt $outdir/pbcr.BNG.breaks5");

my $fSeqbreaks4flt="$outdir/falcon.seq.breaks4.flt";
my $fBNGbreaks4flt="$outdir/falcon.BNG.merged.breaks4.flt";
breaksFlt($fSeqbreaks4,$fMergedBNGbreaks4,$fSeqbreaks4flt,$fBNGbreaks4flt);
system("cp $fSeqbreaks4flt $outdir/falcon.seq.breaks5");
system("cp $fBNGbreaks4flt $outdir/falcon.BNG.breaks5");


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


sub getGoodBNGAln {
  my ($xmap,$flt)=@_;    
  my %flt=%{$flt};
  my $maxHang=2;
  my %BNGaln;
  open IN,"grep -v '#' $xmap|sort -k2,2n -k4,4n -k5,5n|";
  while (<IN>) {
  	chomp;
  	next if (!/^\d/);
  	my @t=split /\t/;  	
  	if ($flt{$t[2].$t[5].$t[6]}) {
      next;
    }    	    	
  	my ($xmapID,$queID,$refID,$queStartPos,$queEndPos,$refStartPos,$refEndPos,$orientation,$confidence,$hitEnum,$queLen,$refLen,$labelChannel,$aln)=split /\t/;
  	my @aln;
  	while ($aln=~m/\((\d+)\,(\d+)\)/g) {
      my ($rsite,$qsite)=($1,$2);	
      push @aln,[($rsite,$qsite)];  
    }      
  	my $start=$aln[0][1];
  	my $end=$aln[$#aln][1];
  	if ($start>$end) {
  		$start=$aln[$#aln][1];
  		$end=$aln[0][1];
  	}  
  	foreach my $k ($start..$end) {
  		$BNGaln{$queID}{$k}=$refID;
  	} 	  	  	  	  	     
  }
  #exit;
  close OUT;close IN;       
  return \%BNGaln;
}



###check Seq conflicting contig's corresponding BNG contig's good alignments in another assembly 
sub checkSeqConfBNGAln {
  my ($seqbreaks,$goodAln,$labKey,$out)=@_;
  open IN,$seqbreaks;
  my %goodAln=%{$goodAln};
  my %labKey=%{$labKey};
  my $w=2;
  my ($n,$m)=(0,0);
  open OUT1,">$out";
  print OUT1 "#cmapID\tseqID\tLength\tBreaksLabel\tBreakPos\tCMAPAlnInfor\tisDoubleBreaks\tNr-sup-Maps\tPvalue\tGoodBNGAlnInAss2\n";
  while (<IN>) {
    chomp;
    $n++;
    my @t=split /\t/;
    my @ta=split /\;/,$t[6];
    my $check="";
    foreach my $target (@ta) {
      my @target=split /\,/,$target;
      if ($goodAln{$target[0]}{$target[3]}) {
      	my $alnSeq=$goodAln{$target[0]}{$target[3]};
      	my $n=0;
      	foreach my $k ($target[3]-$w..$target[3]+$w) {
      	  if (($goodAln{$target[0]}{$k}) and ($goodAln{$target[0]}{$k} eq $alnSeq))  {
  	  	   $n++;
  	      }  	
      	}
      	if ($n==$w*2+1) {
      	  if (!$check) {
      	  	my $seqID=$labKey{$alnSeq};
      	  	$check="BNG-$target[0]:$seqID";
      	  }else {
      	  	my $seqID=$labKey{$alnSeq};      	  	
      	  	$check=$check.";"."BNG-$target[0]:$seqID";
      	  }
      	}
      }
    }
    if ($check) {
      print OUT1 "$_\t$check\n";	
    }else {
      print OUT1 "$_\t-\n";	
      $m++;
    }
     
  }
  close OUT1;close IN;
  print "$seqbreaks ---- total breaks: $n\t no good Aln: $m\n";
}

###check pbcr/falcon conflicting BNG whether have good alignments with falcon/pbcr in silico maps
sub checkBNGbreaks {
  my ($breaks,$aln,$cmapkey,$out,$f)=@_;
  my %aln=%{$aln};
  open IN,$breaks;
  my %cmapKey=%{$cmapkey};
  my $w=2;
  open OUT,">$out";
  print OUT "#Assembler\tcmapID\tseqID\tLength\tBreaksLabel\tBreakPos\tCMAPAlnInfor\tisDoubleBreaks\tNr-sup-Maps\tPvalue\tGoodAlnInAss2\n";
  my ($n,$m)=(0,0);
  my $k=0;
  while (<IN>) {
  	chomp;
  	next if (/#/);
  	my @t=split /\t/;
  	$n++;
  	my $br=int($t[4]);
  	my $n=0;
  	my $alnRef=$aln{$t[0]}{$br};
  	if (!$alnRef) {
  	  print "$f: no good aln $_\n";	
  	  print OUT "$f\t$_\t-\n";
  	  $m++;
  	  next;
  	}
  	foreach my $k ($br-$w..$br+$w) {
  	  if (($aln{$t[0]}{$k}) and ($aln{$t[0]}{$k} eq $alnRef))  {
  	  	$n++;
  	  }	
  	}
  	if ($n==$w*2+1) {
  		my $seqID=$cmapKey{$alnRef};
  	    print OUT "$f\t$_\tSeqID-$alnRef:$seqID\n";	
  	}else {
  	  my $seqID=$cmapKey{$alnRef};
  	  print OUT "$f\t$_\t-N$n-Seq-$alnRef-$seqID\n";
  	  if ($n>0) {  	  	
  	  	print "$f:N<5:$n\t $_\t$alnRef\t$seqID\n";
  	  	$k++;
  	  } 
  	}
  }
  close OUT;close IN;
  print "$breaks ---- total breaks: $n\t no good Aln: $m\tN<5:$k\n";  
}





###get the final conflicting seq breaks
sub mergeSeqbreaks {
  my ($BNGbreaks,$seqbreaks,$goodAln,$out1,$out2)=@_;
  open IN,$seqbreaks;
  open OUT1,">$out1";
  open OUT2,">$out2";
  while (<IN>) {
  	chomp;
  	my @t=split /\t/;
  	if (($t[8]==2) or ($t[9] eq "HiRise") or ($t[12]==1) )  {
  	  print OUT2 "$_\tD\n";
  	}else {
  	  print OUT2 "$_\tU\n";	
  	}
  }  
  
}

###get the final conflicting BNG breaks
sub mergeBNGbreaks {
  my ($pBNGbreaks,$fBNGbreaks,$out1,$out2)=@_;
  
  my $cmapID;my $br;my $pbreaks; my $pr;
  my %pbcrBreaks;
  my %falconBreaks;
  open IN,"cat $pBNGbreaks $fBNGbreaks|sort -k3,3n -k6,6n|";
  while (<IN>) {
  	chomp;
  	next if (/#/);
  	my @t=split /\t/;  	
  	if (!$cmapID) {
  	  $cmapID=$t[1];
  	  $br=$t[5];
  	  $pr=$t[0];
  	  $pbreaks=join("\;",@t[5..$#t]);
  	  
  	}elsif ($cmapID ne $t[1]) {  	  	  		
  	  $cmapID=$t[1];
  	  $br=$t[5];
  	  $pr=$t[0];
  	  $pbreaks=join("\;",@t[5..$#t]);
  	  	
  	}else {
  	  if ($t[5]-$br<=3) {
  	  	my $breaks=join("\;",@t[5..$#t]);
  	  	if (($t[0] eq "PBcR") and ($pr eq "FALCON"))  {
  	  	  $falconBreaks{$t[1]}{$br}=$breaks;
  	  	  $pbcrBreaks{$t[1]}{$t[5]}=$pbreaks;  	
  	  	}elsif (($t[0] eq "FALCON") and ($pr eq "PBcR"))  {
  	  	  $falconBreaks{$t[1]}{$t[5]}=$pbreaks;
  	  	  $pbcrBreaks{$t[1]}{$br}=$breaks;	
  	  	}  	  	
  	    $br=$t[5];
  	    $pr=$t[0];
  	    $pbreaks=join("\;",@t[5..$#t]);
  	    
  	  }else {  	  	
  	    $br=$t[5];
  	    $pr=$t[0];
  	    $pbreaks=join("\;",@t[5..$#t]);
  	  }
  	}
  }
  close IN;
  
  open IN,$pBNGbreaks;
  open OUT1,">$out1";
  print OUT1 "#Assembler\tcmapID\tseqID\tLength\tBreaksLabel\tBreakPos\tCMAPAlnInfor\tisDoubleBreaks\tNr-sup-Maps\tPvalue\tGoodAlnInAss2\n";
  
  while (<IN>) {
  	chomp;
  	next if (/#/);
  	my @t=split /\t/;
  	if ($pbcrBreaks{$t[1]}{$t[5]}) {
  	  my $breaks=$pbcrBreaks{$t[1]}{$t[5]};
  	  print OUT1 "$_\tFALCON:$breaks\n";
  	}else {
  	  print OUT1 "$_\t-----\n";	
  	}
  }
  close OUT1;close IN;
  
  open IN,$fBNGbreaks;
  open OUT2,">$out2";
  print OUT2 "#Assembler\tcmapID\tseqID\tLength\tBreaksLabel\tBreakPos\tCMAPAlnInfor\tisDoubleBreaks\tNr-sup-Maps\tPvalue\tGoodAlnInAss2\n";
  
  while (<IN>) {
  	chomp;
  	next if (/#/);
  	my @t=split /\t/;
  	if ($falconBreaks{$t[1]}{$t[5]}) {
  	  my $breaks=$falconBreaks{$t[1]}{$t[5]};
  	  print OUT2 "$_\tPBcR:$breaks\n";
  	}else {
  	  print OUT2 "$_\t-----\n";	
  	}
  }
  close OUT2;close IN;
  
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
  	if ($t[10]=~m/BNG/) {
  	  $goodBNG++;
  	}

  	if ( ($t[8]==2) or ($t[10]=~m/BNG/) ) {
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
  print "good BNG-assembly2 alignment: $goodBNG\n";
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








