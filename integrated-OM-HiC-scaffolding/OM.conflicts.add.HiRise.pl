#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: OM.misass.add.HR.MP.pl
#
#        USAGE: ./OM.misass.add.HR.MP.pl  
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
#      CREATED: 08/02/2016 05:37:26 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Getopt::Long;
use File::Basename;

my ($seqBreaks,$seqCmap,$HRbreaks,$outdir,$help);
###add HiRise breaks

GetOptions (
  'seqBreak|s=s'    => \$seqBreaks,
  'seqCmap|m=s' => \$seqCmap,
  'outDir|o=s'    => \$outdir,
  'HRbreaks|r=s' => \$HRbreaks,
  'help|h+'      => \$help,
);

my $usage="
Usage: $0 seqKeyFile align1_directory out_dir
Options:
  -seqBreak|-s the break file of sequence
  -seqCmap|-m the sequence in silico cmap in initial alingment with optical cmap
  -outDir|-o output directory
  -HRbreaks|-r hirise breaks bed file
  -help|-h print this help
";

die $usage if (defined $help);
die $usage if ((!$seqBreaks) or (!$seqCmap) or (!$HRbreaks));




system("mkdir -p $outdir");

my $seqCmapHash=hashCmap($seqCmap);
my $name=basename($seqBreaks);
my $seqBreakBed="$outdir/$name.bed";
breaks2Bed($seqBreaks,$seqCmapHash,$seqBreakBed);
my $w=10;
my $win=$w*1000;

system("windowBed -w $win -a $seqBreakBed -b $HRbreaks >$outdir/$name.HR.win$w.bed ");

open IN1,"$outdir/$name.HR.win$w.bed";
my %HRsup;
while (<IN1>) {
  chomp;
  my @t=split /\t/;
  my ($start,$end)=@t[14,15];
  $HRsup{$t[0]}{$t[1]}="$start\,$end";
}
close IN1;

open IN,"$seqBreakBed";

open OUT,">$outdir/$name.HR.bed";
while (<IN>) {
  chomp;
  my @t=split /\t/;
  my $hr="###";
  if ($HRsup{$t[0]}{$t[1]}) {
  	$hr=$HRsup{$t[0]}{$t[1]};
  }
  print OUT "$_\t$hr\n";
}
close OUT;



sub breaks2Bed {
  my ($breaks,$cmapHash,$out)=@_;
  my %cmap=%{$cmapHash};
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
  	}elsif ($t[$#t]==3) {  	  
  	  my ($start,$end)=(split /\,/,$t[5])[0,2];
  	  $start=int($start);
  	  $end=int($end);
  	  my $infor=join("\t",@t[0,2..10]);
  	  print OUT "$t[1]\t$start\t$end\t$infor\n";
  	}else {
  	  my $start=$cmap{$t[0]}{$t[4]};
  	  my $end=$cmap{$t[0]}{$t[6]};
  	  if ($t[4]>$t[6]) {  	  	  	  	
  	  	$start=$cmap{$t[0]}{$t[6]};
  	  	$end=$cmap{$t[0]}{$t[4]};
  	  }
  	  $start=int($start);
  	  $end=int($end);
  	  my $infor=join("\t",@t[0,2..10]);
  	  print OUT "$t[1]\t$start\t$end\t$infor\n";
  	}
  } 
  close IN;
  close OUT;
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