#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: OM.merge.low-high.breaks.pl
#
#        USAGE: ./OM.merge.low-high.breaks.pl  
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
#      CREATED: 08/02/2016 04:03:44 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Getopt::Std;

our ($opt_a,$opt_b,$opt_c,$opt_d,$opt_o,$opt_h);

getopts('a:b:c:d:o:h');

my ($seqBrk1,$seqBrk2,$OMbrk1,$OMbrk2,$outdir);

die "Usage: $0 seq_breaks_low_pvalue seq_breaks_high_pvalue OM_breaks_low_pvalue OM_breaks_high_pvalue out_dir\n" if (defined $opt_h);
	
die "Usage: $0 seq_breaks_low_pvalue seq_breaks_high_pvalue OM_breaks_low_pvalue OM_breaks_high_pvalue out_dir\n" if ((!$opt_a) or (! $opt_b) or (! $opt_c) or (! $opt_d));

system("mkdir -p $opt_o");

my $seqBreaks1=$opt_a;
my $seqBreaks2=$opt_b;
$outdir=$opt_o;
my $seqBreaks3="$outdir/seq.assembly.conflicts.breaks3";
my $OMbreaks1=$opt_c;
my $OMbreaks2=$opt_d;
my $OMbreaks3="$outdir/OM.assembly.conflicts.breaks3";

open IN1,$seqBreaks1;
open OUT1,"|sort -k1,1n -k5,5n >$seqBreaks3";
my %seqLow;
while (<IN1>) {
  chomp;
  next if (/#/);
  my @t=split /\t/;
  $seqLow{$t[0]}{$t[4]}{$t[8]}=1;
  print OUT1 "$_\tLP\n";
}
open IN2,$seqBreaks2;
while (<IN2>) {
  chomp;
  next if (/#/);
  my @t=split /\t/;
  if (!$seqLow{$t[0]}{$t[4]}{$t[8]} ) {
  	print OUT1 "$_\tHP\n";
  }
}
close IN1;close IN2;close OUT1;

open IN1,$OMbreaks1;
my %OMLow;
open OUT2,"|sort -k1,1n -k5,5n >$OMbreaks3";
while (<IN1>) {
  chomp;
  next if (/#/);
  my @t=split /\t/;
  $OMLow{$t[0]}{$t[4]}{$t[8]}=1;
  print OUT2 "$_\tLP\n";
}
open IN2,$OMbreaks2;
while (<IN2>) {
  chomp;
  next if (/#/);
  my @t=split /\t/;
  if (!$OMLow{$t[0]}{$t[4]}{$t[8]}) {
  	print OUT2 "$_\tHP\n";
  }
}
close IN1;close IN2;close OUT2;

###step2:get breakpoints, and check them whether they are supported "double-contig" and HiRise 

my $seqBreaks4="$outdir/seq.breaks4";
mergeBreaks($seqBreaks3,$seqBreaks4);

my $OMbreaks4="$outdir/OM.breaks4";
mergeBreaks($OMbreaks3,$OMbreaks4);

sub mergeBreaks {
  my ($break,$break2)=@_;    
  open IN,"grep -v '#' $break|sort -k1,1n -k5,5n|";  
  my ($cmapID,$seqID,$leng,$n,$brkLab,$brkPos,$nextLab,$nextPos,$tar,$tarBrk,$nextTarBrk,$flag,$pval);  
  my ($br,$bp,$nbr,$nbp,$pvalue,$target);
  open OUT,">$break2";
  my $m=0;
  while (<IN>) {
  	chomp;
  	my @t=split /\t/;  	  	
  	if (!$cmapID) {
  	  ($cmapID,$seqID,$leng,$n,$brkLab,$brkPos,$nextLab,$nextPos,$tar,$tarBrk,$nextTarBrk,$flag,$pval)=@t;
  	  $target="$tar\,$tarBrk\,$nextTarBrk";
      $br=$brkLab;
      $bp=$brkPos;
      $nbr=$nextLab;
      $nbp=$nextPos;
  	  $pvalue=$pval;
  	  $m=1;    	  
  	}elsif ($cmapID ne $t[0]) {  	  
  	  print OUT "$cmapID\t$seqID\t$leng\t$n\t$br\t$bp\t$nbr\t$nbp\t$target\t$pvalue\t$m\n";  ####to be modified, breakpoint=the mid of ajacent broken labels
  	  ($cmapID,$seqID,$leng,$n,$brkLab,$brkPos,$nextLab,$nextPos,$tar,$tarBrk,$nextTarBrk,$flag,$pval)=@t;
  	  $target="$tar\,$tarBrk\,$nextTarBrk";
      $br=$brkLab;
      $bp=$brkPos;
  	  $nbr=$nextLab;
      $nbp=$nextPos;
  	  $pvalue=$pval;
  	  $m=1;    	    	  
  	}else {
  	  if ($t[4]-$brkLab<=4) {
  	    $br=$br.",".$t[4];
  	    $bp=$bp.",".$t[5];  	    
  	    $nbr=$nbr.",".$t[6];
  	    $nbp=$nbp.",".$t[7];  	    
  	    $m++;
  	    $target=$target.";"."$t[8]\,$t[9]\,$t[10]";  
  	    $pvalue=$pval.",".$t[12];	    
  	  }else {  	  	
  	    print OUT "$cmapID\t$seqID\t$leng\t$n\t$br\t$bp\t$nbr\t$nbp\t$target\t$pvalue\t$m\n"; 
  	  	($cmapID,$seqID,$leng,$n,$brkLab,$brkPos,$nextLab,$nextPos,$tar,$tarBrk,$nextTarBrk,$flag,$pval)=@t;
  	  	$target="$tar\,$tarBrk\,$nextTarBrk";
        $br=$brkLab;
        $bp=$brkPos;
        $nbr=$nextLab;
        $nbp=$nextPos;
  	    $pvalue=$pval;
  	    $m=1;    	   
  	  }
  	}
  }
  print OUT "$cmapID\t$seqID\t$leng\t$n\t$br\t$bp\t$nbr\t$nbp\t$target\t$pvalue\t$m\n"; 
  close OUT;close IN;
}




