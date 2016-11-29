#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: OM.breaks.check.ass2.pl
#
#        USAGE: ./OM.breaks.check.ass2.pl  
#
#  DESCRIPTION: check whether conflict CMAP region is fullly aligned by another seq CMAP
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wen-Biao Jiao (), 
# ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
#      VERSION: 1.0
#      CREATED: 08/03/2016 10:48:51 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Getopt::Long;
use File::Basename;


my ($OMbreaks,$seqBreaksBed,$seq2xmap,$seq2KeyFile,$outdir,$help);
GetOptions (
  'seqBreak|s=s'    => \$seqBreaksBed,
  'OMbreaks|c=s' => \$OMbreaks,
  'outDir|o=s'    => \$outdir,
  'ass2xmap|x=s' => \$seq2xmap,
  'ass2key|k=s' => \$seq2KeyFile,
  'help|h+'      => \$help,
);

my $usage="
Usage: $0 seqKeyFile align1_directory out_dir
Options:
  -seqBreak|-s the break file of sequence
  -OMbreaks|-c the break file of optical consensus maps
  -outDir|-o the ouput directory
  -ass2xmap|-x another assembly's alignments with optical consensus maps
  -ass2key|-x another assembly's key file of in silico optical maps
  -help|-h print this help
";

die $usage if (defined $help);
die "can't find the file $seqBreaksBed" if (!$seqBreaksBed);
die "can't find the file $OMbreaks" if (!$OMbreaks);
die "can't find the file $seq2xmap" if (!$seq2xmap); 
die "can't find the file $seq2KeyFile" if (!$seq2KeyFile);


system("mkdir -p $outdir");
my $flt=xmapFilter($seq2xmap);
my $seq2Key=hashNGScmapKey($seq2KeyFile);
my $OM=basename("$OMbreaks");
my $OMbreaksBed="$outdir/$OM.bed";
breaks2Bed($OMbreaks,$OMbreaksBed);

##
my $w=2; ###the minimal aligned nick sites in the upstream/downstream of breaking point to indicate a fully alignment spanning the conflicting region  
my $OMaln=getGoodOMaln($seq2xmap,$flt);

my $seq=basename($seqBreaksBed,".bed");
my $out1="$outdir/$seq.OMGoodAln.bed";
checkSeqConfOMaln($seqBreaksBed,$OMaln,$seq2Key,$out1,$w);
my $out2="$outdir/$OM.OMGoodAln.bed";
checkOMbreaks($OMbreaksBed,$OMaln,$seq2Key,$out2,$w);



###check Seq conflicting contig's corresponding BNG contig's good alignments in another assembly 
sub checkSeqConfOMaln {
  my ($seqbreaks,$goodAln,$labKey,$out,$w)=@_;
  open IN,$seqbreaks;
  my %goodAln=%{$goodAln};
  my %labKey=%{$labKey};  
  my ($n,$m)=(0,0);
  open OUT1,">$out";
  while (<IN>) {
    chomp;
    $n++;
    my @t=split /\t/;
    my @ta=split /\;/,$t[10];
    my $check="";
    foreach my $target (@ta) {
      my @target=split /\,/,$target;
      if ($goodAln{$target[0]}{$target[1]}) {
      	my $alnSeq=$goodAln{$target[0]}{$target[1]};
      	my $n=0;
      	foreach my $k ($target[1]-$w..$target[1]+$w) {
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
      print OUT1 "$_\t###\n";	
      $m++;
    }
     
  }
  close OUT1;close IN;
  print "$seqbreaks ---- total breaks: $n\t no good Aln: $m\n";
}


sub checkOMbreaks {
  my ($breaks,$aln,$cmapkey,$out,$w)=@_;
  my %aln=%{$aln};
  open IN,$breaks;
  my %cmapKey=%{$cmapkey};
  
  open OUT,">$out";
  my $m=0;my $n=0;
  my $k=0;
  my $ttt=$w*2+1;
  while (<IN>) {
  	chomp;
  	next if (/#/);
  	$n++;
  	my @t=split /\t/;
  	my @br=(split /\,/,$t[6]);  	
  	my $check="###";
  	foreach my $br (@br) {
  	  my $t=0;
  	  my $alnRef=$aln{$t[0]}{$br};
  	  if (!$alnRef) {
  	    print "no good aln $_\n";	
  	    #print OUT "$_\t###\n";
  	    $m++;
  	    next;
  	  }	
  	  foreach my $k ($br-$w..$br+$w) {
  	    if (($aln{$t[0]}{$k}) and ($aln{$t[0]}{$k} eq $alnRef))  {
  	  	  $t++;
  	    }	
  	  }
  	  if ($t==$w*2+1) {
  		my $seqID=$cmapKey{$alnRef};
  	    if ($check eq "###") {
  	      $check="CmapID-$alnRef:$seqID";
  	    }else {
  	      $check=$check.";"."CmapID-$alnRef:$seqID";	
  	    }
  	      	      	    	
  	  }else {
  	  	$k++;
  	  	if ($check eq "###") {
  	  	  $check="N<$ttt:$t"; 	
  	  	}else {
  	  	  $check=$check.";"."N<$ttt:$t";	
  	  	}  	  	
  	  }
  	}
  	print OUT "$_\t$check\n";
  }
  close OUT;close IN;
  print "$breaks ---- total breaks: $n\t no good Aln: $m\tN<$ttt:$k\n";  
}

sub getGoodOMaln {
  my ($xmap,$flt)=@_;    
  my %flt=%{$flt};
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
  	###ajust the start and end of aligned region when indel(>=2 nick sites) exists in the start/end of alignment
  	$start=$aln[1][1] if  ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;
  	$end=$aln[$#aln-1][1] if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) ) ;

  	if ($start>$end) {
  		$start=$aln[$#aln][1];
  		$end=$aln[0][1];
  		$start=$aln[$#aln-1][1] if ( ($hitEnum=~m/1M$/) and ($hitEnum!~m/1D1M$/) and ($hitEnum!~m/1I1M$/) ) ;
  		$end=$aln[1][1] if  ( ($hitEnum=~m/^1M/) and ($hitEnum!~m/^1M1D/) and ($hitEnum!~m/^1M1I/)  )  ;  	    
  	}  
  	foreach my $k ($start..$end) {
  		$BNGaln{$queID}{$k}=$refID;
  	} 	  	  	  	  	     
  }
  close OUT;close IN;       
  return \%BNGaln;
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


sub breaks2Bed {
  my ($breaks,$out)=@_;
  open IN,$breaks;
  open OUT,">$out";
  while (<IN>) {
  	chomp;
  	next if (/#/);
  	my @t=split /\t/;
  	if ($t[$#t]==2) {
  	  my ($start,$end)=split /\,/,$t[5];
  	  $start=int($start);
  	  $end=int($end);
  	  my $infor=join("\t",@t[0,2..10]);
  	  print OUT "$t[1]\t$start\t$end\t$infor\n";
  	}else {
  	  my $start=$t[5];
  	  my $end=$t[7];
  	  ($start,$end)=($end,$start) if ($t[4]>$t[6]);  	  	  	  	  	  	
  	  $start=int($start);
  	  $end=int($end);
  	  my $infor=join("\t",@t[0,2..10]);
  	  print OUT "$t[1]\t$start\t$end\t$infor\n";
  	}
  } 
}