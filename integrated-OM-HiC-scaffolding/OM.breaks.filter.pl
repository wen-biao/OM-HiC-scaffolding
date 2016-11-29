#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: OM.breaks.filter.pl
#
#        USAGE: ./OM.breaks.filter.pl  
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
#      CREATED: 08/03/2016 10:49:21 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use File::Basename;
use Getopt::Long;

my ($OMbreaksBed,$seqBreaksBed,$outdir,$help);


GetOptions (
  'seqBreak|s=s'    => \$seqBreaksBed,
  'OMbreaks|c=s' => \$OMbreaksBed,
  'outDir|o=s'    => \$outdir,
  'help|h+'      => \$help,
);

my $usage="
Usage: $0 seqKeyFile align1_directory out_dir
Options:
  -seqBreak|-s the break file of sequence
  -OMbreaks|-c the break file of optical consensus maps  
  -HRbreaks|-r hirise breaks bed file
   -outDir|-o output directory
  -help|-h print this help
";

die $usage if (defined $help);


open IN,$seqBreaksBed;
my %seqSup;
my ($cmap2,$MP,$HR,$ass2,$r)=(0,0,0,0,0);
my $name1=basename($seqBreaksBed,".bed");
open OUT1,">$outdir/$name1.check.bed";
open OUT2,">$outdir/$name1.flt.bed";
my %seqRemove;
while (<IN>) {
  chomp;
  my @t=split /\t/;
  if ($t[12]==2) {
  	$cmap2++;
  }
  if ($t[13] ne "###") {
  	$HR++;
  }
  if ($t[14] ne "###") {
  	$ass2++;
  }
  if ( ($t[12]==2) or ($t[13] ne "###") or ($t[14] ne "###") or ($t[14] ne "###") ) {
  	print OUT1 "$_\tSeqMis\n";
  	my @OM=split /\;/,$t[10];
  	foreach my $omap (@OM) {
  	  my @tt=split /\,/,$omap;
  	  $seqSup{$tt[0]}{$tt[1]}{$t[3]}=1;
  	}
  	print OUT2 "$_\tSeqMis\n";
  }else {
  	print OUT1 "$_\tUndetermined\n";
  	if ($t[11] ne "HP") {
  	  print OUT2 "$_\tUndetermined\n";  	 
  	}else {
  	   my @tt=split /\,/,$t[10];
  	  $seqRemove{$tt[0]}{$tt[1]}{$t[3]}=1;
  	  $r++;
  	}
  } 
}
print "2CMAP\t$cmap2\n";
print "HiRise\t$HR\n";
print "Assembly2\t$ass2\n";
print "Seq:Undetermined breaks in High Pvalue alignment\t$r\n";

close OUT1;close OUT2;
close IN;

open IN,"$OMbreaksBed";
my $name2=basename($OMbreaksBed,".bed");
open OUT1,">$outdir/$name2.check.bed";
open OUT2,">$outdir/$name2.flt.bed";
$r=0;
while (<IN>) {
  chomp;
  my @t=split /\t/;
  my @br=split /\,/,$t[6];
  my $check=0;
  foreach my $br (@br) {
    if ($seqSup{$t[0]}{$br}) {
      $check++;
    }
  }
  if ($check<$#br+1) {
  	print OUT1 "$_\tUndetermined\n";
  	my @tt=split /\,/,$t[10];
  	if ($seqRemove{$t[0]}{$t[6]}{$tt[0]}) {
  	  $r++;	
  	}else {
  	  print OUT2 "$_\tUndetermined\n";
  	}
  }else {
  	print OUT1 "$_\tSeqMiss\n";
  }
}
close OUT1;close OUT2;close IN;

print "OM:Undetermined breaks in High Pvalue alignment\t$r\n";




