#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: OM.scaffolding.pl
#
#        USAGE: ./OM.scaffolding.pl  
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
#      CREATED: 08/05/2016 10:20:23 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Bio::Seq;
use Bio::SeqIO;
use File::Basename;
use Getopt::Long;
use Cwd;

my ($xmap,$rmap,$qmap,$keyfile,$seq,$outdir,$sID,$nick,$help);

GetOptions (
  'xmap|x=s' => \$xmap,
  'qmap|q=s' => \$qmap,
  'rmap|r=s' => \$rmap,
  'key|k=s'    => \$keyfile,
  'seq|s=s' => \$seq,
  'sID|i=s' => \$sID,
  'nick|n=s' => \$nick,  
  'outDir|o=s'    => \$outdir,
  'help|h+'      => \$help,
);

$nick="GCTCTTC" if (!$nick);

my $usage="
Usage: $0 xmap rmap qmap key sequence scafID nick_site outputdir 
Options:
  -xmap|-x File of cmap alignment between in silico maps of sequences and hybrid-scaffolded optical consensus maps 
  -qmap|-q File of in silico maps generated from assembly sequences 
  -rmap|-r File of hybrid-scaffolded optical consensus maps
  -key|-k In silico CMAP key file of sequence assembly
  -sID|-i Prefix of exported scaffold ID
  -nick|-n Sequence of nick site used in the generation of optical maps  
  -outDir|-o output directory
  -help|-h print this help
";

die $usage if (defined $help);


die "cann't find the file $xmap\n" if (! -e $xmap);
die "cann't find the file $rmap\n" if (! -e $rmap);
die "cann't find the file $qmap\n" if (! -e $qmap);
die "cann't find the file $seq\n" if (! -e $seq);
die "cann't find the file $keyfile\n" if (! -e $keyfile);

if (!$outdir) {
  $outdir=getcwd()."scaffolding";
}
system("mkdir -p $outdir");

my $qmapHash=hashCmap($qmap);
my %qmap=%{$qmapHash};

my $rmapHash=hashCmap($rmap);

my $ctgID=hashAssKey($keyfile);
my %ctgID=%{$ctgID};

my $refMapHash=hashHybridCmap($rmap);
my %refMapHash=%{$refMapHash};

my $name=basename($xmap,".xmap");
my $xmapFlt="$outdir/$name.flt.xmap";

my $insideCtg=xmapFilter($qmapHash,$rmapHash,$xmap,$xmapFlt);
my %insideCtg=%{$insideCtg};
my %inCtg;

my $insideLeng=0;
my $k=0;
foreach (sort {$a<=>$b} keys %insideCtg) {
  print "inside contig: $_\t";
  my $c=$ctgID{$_};
  $inCtg{$c}=1;
  print "$c\tLength:$qmap{$_}{'leng'}\n  xmap:$insideCtg{$_}\n";#exit;
  $insideLeng+=$qmap{$_}{'leng'};
  $k++;
}
print "total number of inside contigs:$k\n";
print "total length of inside contigs:$insideLeng\n";
#exit;

my $assSeq=hashAssSeq($seq);
my %assSeq=%{$assSeq};
my %rawAssSeq=%{$assSeq};

open OUT1,">$outdir/contig.scaffolding.infor";
open OUT3,">$outdir/scaffolds.contig.position.txt";
print OUT1 "hybridScaffoldID(Ref)\thybridScafLeng\tin-silico-ctg1\tctg1ID\tLength\tAlnStart\tAlnEnd\tOrientation\tAlnRefStart\tAlnRefEnd\tCtg1-tail\t";
print OUT1                                        "in-silico-ctg2\tctg2ID\tLength\tAlnStart\tAlnEnd\tOrientation\tAlnRefStart\tAlnRefEnd\tCtg2-tail\tOMdistance\tGap\tBlastGap\n";
system("mkdir -p $outdir/overlap");
open IN,"grep -v '#' $xmapFlt|sort -k3,3n -k6,6n -k7,7nr|";
my %ctgUsed; my %scafseq;
my @prev=();
my ($scafGap,$scafCmapOvp,$scafTailOvp,$scafGapLeng,$scafCmapOvpLeng,$scafTailOvpLeng)=(0,0,0,0,0,0);
while (<IN>) {
  chomp;
  #my ($xmapID,$queID,$refID,$queStartPos,$queEndPos,$refStartPos,$refEndPos,$orientation,$confidence,$hitEnum,$queLen,$refLen,$labelChannel,$aln)=split /\t/;
  my @t=split /\t/;
  if (!@prev) {
  	$scafseq{$t[2]}="";
  	@prev=@t;
  	$ctgUsed{$ctgID{$t[1]}}=1;
  }else {
  	$ctgUsed{$ctgID{$t[1]}}=1;
  	my $ctgID1=$ctgID{$prev[1]};
  	my $ctgID2=$ctgID{$t[1]};
  	if ($t[2] eq $prev[2]) {
  	  if ($t[5]>=$prev[6]) {
  	    ##Hybrid:  --|----------------|-----------\-------------\-----------
  	  	
  	    ##AssCtg1: --|----------------|----....
  	  
  	    ##AssCtg2:                         ....---\--------------\-------  	  
  	    my $hybridGap=int($t[5]-$prev[6]);
  	    my $tailCtg1=int($prev[10]-$prev[4]);
  	    my $tailCtg2=int($t[3]);
  	    $tailCtg1=int($prev[4]-1) if ($prev[7] eq "-");
  	    $tailCtg2=int($t[10]-$t[3]) if ($t[7] eq "-"); 	  
  	    my $gap=int($hybridGap-$tailCtg1-$tailCtg2);
        my $omDist=int($t[5]-$prev[6]);
        print OUT1 "$t[2]\t$prev[11]\t$prev[1]\t$ctgID1\t$prev[10]\t$prev[3]\t$prev[4]\t$prev[7]\t$prev[5]\t$prev[6]\t";
  	  	print OUT1                      "$t[1]\t$ctgID2\t$t[10]\t$t[3]\t$t[4]\t$t[7]\t$t[5]\t$t[6]\t$tailCtg1\t$tailCtg2\t$omDist\t";
  	    if ($hybridGap>=2*($tailCtg1+$tailCtg2)) {  ###
  	  	  print OUT1 "$gap\t$gap\n";
  	  	  print "no overlap: $t[2]\t$ctgID1\t$ctgID2\t$omDist\t$tailCtg1\t$tailCtg2\t$gap\n";
  	  	  $scafGap++;
  	  	  $scafGapLeng+=$gap;
  	      
  	      if (!$scafseq{$prev[2]}) {
  	      	my $scafEnd=length($assSeq{$ctgID1});
  	      	my $rawLeng=length($rawAssSeq{$ctgID1});
  	      	my $sleng=length($assSeq{$ctgID1});
  	      	print OUT3 "$sID-$prev[2]\t$ctgID1\t$rawLeng\t$sleng\t1\t$scafEnd\t";
  	      }else {
  	      	my $scafStart=length($scafseq{$prev[2]})+1;
  	      	my $scafEnd=$scafStart+length($assSeq{$ctgID1})-1;
  	      	my $rawLeng=length($rawAssSeq{$ctgID1});
  	      	my $sleng=length($assSeq{$ctgID1});
  	      	print OUT3 "$sID-$prev[2]\t$ctgID1\t$rawLeng\t$sleng\t$scafStart\t$scafEnd\t";  	      	
  	      }
  	      my %refHash=%{$refMapHash{$t[2]}};
          my $nick1=$refHash{"pos"}{$prev[6]};
          my $nick2=$refHash{"pos"}{$t[5]};
          my $gapN;
          if ($nick1<=$nick2-2) {          	
          	  $gapN=nickScafGaps(\%refHash,$prev[6],$t[5],$tailCtg1,$tailCtg2,$nick);   #($refHash,$start,$end,$tail1,$tail2,$nick)
          	  my $dist1=$refHash{"site"}{$nick1+1}-$prev[6];  	
              my $dist2=$t[5]-$refHash{"site"}{$nick2-1};
          	#print "GapSeq $ctgID1\t$ctgID2\t$prev[6]\t$t[5]\t$nick1\t$nick2\t$dist1\t$dist2: $gapN\n";
          }else {
          	$gapN="N" x $gap;
          }
  	      if ($prev[7] eq "+") {  	      	
  	        $scafseq{$prev[2]}=$scafseq{$prev[2]}.$assSeq{$ctgID1}.$gapN;
  	        print OUT3 "+\n";	
  	      }else {  	      	
  	      	$scafseq{$prev[2]}=$scafseq{$prev[2]}.recomSeq($assSeq{$ctgID1}).$gapN;
  	      	print OUT3 "-\n";
  	      }
  	      
  	        	        	      	  	  	
  	    }else {
 	  	  my ($overlapSeq1,$overlapSeq2);	  	  	  	   
 	  	  if ($prev[7] eq "+") {
  	        $overlapSeq1=substr($rawAssSeq{$ctgID1},length($rawAssSeq{$ctgID1})-$tailCtg1-$tailCtg2,$tailCtg1+$tailCtg2);   	           
  	      }else {
  	      	$overlapSeq1=recomSeq(substr($rawAssSeq{$ctgID1},0,$tailCtg1+$tailCtg2));  	      	
  	      }  	      
  	      if ($t[7] eq "+") {
  	      	$overlapSeq2=substr($rawAssSeq{$ctgID2},0,$tailCtg1+$tailCtg2);	
  	      }else {
  	      	$overlapSeq2=recomSeq(substr($rawAssSeq{$ctgID2},length($rawAssSeq{$ctgID2})-$tailCtg1-$tailCtg2,$tailCtg1+$tailCtg2));
  	      }
  	      
  	      my ($nctgID1,$nctgID2)=($ctgID1,$ctgID2);
  	      $nctgID1=~s/\|quiver//;
  	      $nctgID2=~s/\|quiver//;

  	      my $outfa="$outdir/overlap/$t[2]\_$nctgID1\_$prev[7]_$nctgID2\_$t[7].fa";
  	      open OUT2,">$outfa";
  	      print OUT2 ">$t[2]\_$nctgID1\_$prev[7]\n$overlapSeq1\n";
  	      print OUT2 ">$t[2]\_$nctgID2\_$t[7]\n$overlapSeq2\n";
  	      close OUT2;
  	      system("makeblastdb -in $outfa -dbtype nucl");
  	      my $blastout="$outdir/overlap/$t[2]\_$nctgID1\_$prev[7]_$nctgID2\_$t[7].blastout";
  	      system("blastn -query $outfa -db $outfa -outfmt 6 -out $blastout -evalue 0.00001"); 
  	      
  	      ###refine the overlapping region based on the blast results
  	      open IN2,$blastout or die "no $blastout";
  	      my $ctg2start=0;
  	      while (<IN2>) {
  	        chomp;
  	        my @tt=split /\t/;
  	        next if ($tt[0] eq $tt[1]);
  	        if ( ($tt[2]>90) and (length($overlapSeq1)-$tt[7]<5) and ($tt[8]<5)) {
  	           $ctg2start=$tt[9]; 
  	           print "blast shows tail-overlap:$_\n";
  	           last; 	
  	        }  	      
  	      }
  	      close IN2;

  	      my $l1=length($overlapSeq1);
  	      my $l2=length($overlapSeq2);  	      
  	      if (($ctg2start==0) and ($gap>0))  {
  	      	if (!$scafseq{$prev[2]}) {
  	      	  my $scafEnd=length($assSeq{$ctgID1});
  	      	  my $rawLeng=length($rawAssSeq{$ctgID1});
  	      	  my $sleng=length($assSeq{$ctgID1});
  	      	  print OUT3 "$sID-$prev[2]\t$ctgID1\t$rawLeng\t$sleng\t1\t$scafEnd\t";  	      	  
  	        }else {
  	      	  my $scafStart=length($scafseq{$prev[2]})+1;
  	      	  my $scafEnd=$scafStart+length($assSeq{$ctgID1})-1;
  	      	  my $rawLeng=length($rawAssSeq{$ctgID1});
  	      	  my $sleng=length($assSeq{$ctgID1});
  	      	  print OUT3 "$sID-$prev[2]\t$ctgID1\t$rawLeng\t$sleng\t$scafStart\t$scafEnd\t";  	      
  	        } 
  	      	
  	        my %refHash=%{$refMapHash{$t[2]}};
            my $nick1=$refHash{"pos"}{$prev[6]};
            my $nick2=$refHash{"pos"}{$t[5]};
            my $gapN;
            if ($nick1<=$nick2-2) {
          	  $gapN=nickScafGaps(\%refHash,$prev[6],$t[5],$tailCtg1,$tailCtg2,$nick)   #($refHash,$start,$end,$tail1,$tail2,$nick)
            }else {
          	  $gapN="N" x $gap;
            }
  	        if ($prev[7] eq "+") {
  	          $scafseq{$prev[2]}=$scafseq{$prev[2]}.$assSeq{$ctgID1}.$gapN;
  	          print OUT3 "+\n";	
  	        }else {  	      	
  	      	  $scafseq{$prev[2]}=$scafseq{$prev[2]}.recomSeq($assSeq{$ctgID1}).$gapN;
  	      	  print OUT3 "-\n";
  	        }
  	        print OUT1 "No-tail-overlap:$gap\t$gap\n";
  	        print "blast result shows no tail-overlap: $t[2]\t$ctgID1\t$ctgID2\t$omDist\t$tailCtg1\t$tailCtg2\t$gap\n";
  	        $scafGap++;
  	  	    $scafGapLeng+=$gap;	
  	      }elsif (($ctg2start==0) and ($gap<=0) )  {
  	      	print OUT1 "Tail-overlap:$gap\t0\n";  ###may use NNNNNNNNNNN instead 0
  	      	print "blast result shows no tail-overlap(but gap<0): $t[2]\t$ctgID1\t$ctgID2\t$omDist\t$tailCtg1\t$tailCtg2\t$gap\n";
  	      	$scafTailOvp++;
  	  	    $scafTailOvpLeng+=$ctg2start;
  	  	    
  	  	    if (!$scafseq{$prev[2]}) {
  	      	  my $scafEnd=length($assSeq{$ctgID1});
  	      	  my $rawLeng=length($rawAssSeq{$ctgID1});
  	      	  my $sleng=length($assSeq{$ctgID1});
  	      	  print OUT3 "$sID-$prev[2]\t$ctgID1\t$rawLeng\t$sleng\t1\t$scafEnd\t";  	      	  
  	        }else {
  	      	  my $scafStart=length($scafseq{$prev[2]})+1;
  	      	  my $scafEnd=$scafStart+length($assSeq{$ctgID1})-1;
  	      	  my $rawLeng=length($rawAssSeq{$ctgID1});
  	      	  my $sleng=length($assSeq{$ctgID1});
  	      	  print OUT3 "$sID-$prev[2]\t$ctgID1\t$rawLeng\t$sleng\t$scafStart\t$scafEnd\t";
  	       }   	  	    
  	  	    my $gapN="N" x 10;
  	        if ($prev[7] eq "+") {
  	          $scafseq{$prev[2]}=$scafseq{$prev[2]}.$assSeq{$ctgID1}.$gapN;
  	          print OUT3 "+\n";	
  	        }else {  	      	
  	      	  $scafseq{$prev[2]}=$scafseq{$prev[2]}.recomSeq($assSeq{$ctgID1}).$gapN;
  	      	  print OUT3 "-\n";
  	        }  	        
  	        #$scafGap++;
  	  	    #$scafGapLeng+=$gap;
  	  	    
  	      }elsif (($ctg2start>0) and ($gap>=0)) {
  	      	print OUT1 "Tail-overlap/tandem:$gap\t$gap\n"; ###maybe a tandem repeat here
  	      	print "blast result shows tail-overlap(but gap>0): $t[2]\t$ctgID1\t$ctgID2\t$omDist\t$tailCtg1\t$tailCtg2\t$gap\tsubstr-Ctg2start:$ctg2start\n";
  	      	#$scafTailOvp++;
  	  	    #$scafTailOvpLeng+=$ctg2start;
  	      	#print OUT1 "No-tail-overlap:$gap\t$gap\n";
  	        #print "blast result shows no tail-overlap: $t[2]\t$ctgID1\t$ctgID2\t$omDist\t$tailCtg1\t$tailCtg2\t$gap\n";
  	        
  	      	$ctg2start=0;  	   
  	      	if (!$scafseq{$prev[2]}) {
  	      	  my $scafEnd=length($assSeq{$ctgID1});
  	      	  my $rawLeng=length($rawAssSeq{$ctgID1});
  	      	  my $sleng=length($assSeq{$ctgID1});
  	      	  print OUT3 "$sID-$prev[2]\t$ctgID1\t$rawLeng\t$sleng\t1\t$scafEnd\t";  
  	        }else {
  	      	  my $scafStart=length($scafseq{$prev[2]})+1;
  	      	  my $scafEnd=$scafStart+length($assSeq{$ctgID1})-1;
  	      	  my $rawLeng=length($rawAssSeq{$ctgID1});
  	      	  my $sleng=length($assSeq{$ctgID1});
  	      	  print OUT3 "$sID-$prev[2]\t$ctgID1\t$rawLeng\t$sleng\t$scafStart\t$scafEnd\t";
  	        }    	
  	      	my $gapN="N" x $gap;
  	        if ($prev[7] eq "+") {
  	          $scafseq{$prev[2]}=$scafseq{$prev[2]}.$assSeq{$ctgID1}.$gapN;
  	          print OUT3 "+\n";	
  	        }else {  	      	
  	      	  $scafseq{$prev[2]}=$scafseq{$prev[2]}.recomSeq($assSeq{$ctgID1}).$gapN;
  	      	  print OUT3 "-\n";
  	        }  	        
  	        $scafGap++;
  	  	    $scafGapLeng+=$gap;	
  	      	 
  	      }else { ##$ctg2start>0,gap<0
  	        if (!$scafseq{$prev[2]}) {
  	      	  my $scafEnd=length($assSeq{$ctgID1});
  	      	  my $rawLeng=length($rawAssSeq{$ctgID1});
  	      	  my $sleng=length($assSeq{$ctgID1});
  	      	  print OUT3 "$sID-$prev[2]\t$ctgID1\t$rawLeng\t$sleng\t1\t$scafEnd\t";  
  	        }else {
  	      	  my $scafStart=length($scafseq{$prev[2]})+1;
  	      	  my $scafEnd=$scafStart+length($assSeq{$ctgID1})-1;
  	      	  my $rawLeng=length($rawAssSeq{$ctgID1});
  	      	  my $sleng=length($assSeq{$ctgID1});
  	      	  print OUT3 "$sID-$prev[2]\t$ctgID1\t$rawLeng\t$sleng\t$scafStart\t$scafEnd\t";
  	        } 
  	      	if ($prev[7] eq "+") {
  	          $scafseq{$prev[2]}.=$assSeq{$ctgID1};
  	          print OUT3 "+\n";  	              	           
  	        }else {
  	      	  $scafseq{$prev[2]}.=recomSeq($assSeq{$ctgID1});
  	      	  print OUT3 "-\n";  	      	    	      	
  	        }  	        
  	        print OUT1 "Tail-overlap:$gap\t-$ctg2start\n";
  	        print "blast result shows tail-overlap: $t[2]\t$ctgID1\t$ctgID2\t$omDist\t$tailCtg1\t$tailCtg2\t$gap\tsubstr-Ctg2start:$ctg2start\n";
  	        $scafTailOvp++;
  	  	    $scafTailOvpLeng+=$ctg2start;
  	        
  	      }  	          	      
  	      if ($t[7] eq "+") {
  	        $assSeq{$ctgID2}=substr($assSeq{$ctgID2},$ctg2start,length($assSeq{$ctgID2})-$ctg2start);
  	      }else {
  	      	$assSeq{$ctgID2}=substr($assSeq{$ctgID2},0,length($assSeq{$ctgID2})-$ctg2start);
  	      }  	       	       	        	  	
  	    }  		  		
  	  }else {
  	  ###contig CMAP overlap
  	  ##Hybrid:  --|----------------\---|----------\-----------------------  	  	
  	  ##AssCtg1: --|----------------\---|-  	  
  	  ##AssCtg2:                   -\---|----------\---  	    	  
  ##SuperContigs: -|------------------
  	   	 
  	   	my ($overlapSeq1,$overlapSeq2);    
  	    ###remove the overlaped region
  	    my $hybridGap=int($t[5]-$prev[6]);  	    
  	    my $tailCtg1=int($prev[10]-$prev[4]);
  	    my $tailCtg2=int($t[3]);
  	    $tailCtg1=int($prev[4]-1) if ($prev[7] eq "-");
  	    $tailCtg2=int($t[10]-$t[3]) if ($t[7] eq "-"); 
  	    
  	    if (!$scafseq{$prev[2]}) {
  	      my $scafEnd=length($assSeq{$ctgID1});
  	      my $rawLeng=length($rawAssSeq{$ctgID1});
  	      my $sleng=length($assSeq{$ctgID1});
  	      print OUT3 "$sID-$prev[2]\t$ctgID1\t$rawLeng\t$sleng\t1\t$scafEnd\t";  
  	    }else {
  	      my $scafStart=length($scafseq{$prev[2]})+1;
  	      my $scafEnd=$scafStart+length($assSeq{$ctgID1})-1;
  	      my $rawLeng=length($rawAssSeq{$ctgID1});
  	      my $sleng=length($assSeq{$ctgID1});
  	      print OUT3 "$sID-$prev[2]\t$ctgID1\t$rawLeng\t$sleng\t$scafStart\t$scafEnd\t";  	      
  	    } 
  	    
  	    if ($prev[7] eq "+") {  	     
  	      $scafseq{$prev[2]}.=$assSeq{$ctgID1};
  	      print OUT3 "+\n";	  
  	    }else {
  	      $scafseq{$prev[2]}.=recomSeq($assSeq{$ctgID1});
  	      print OUT3 "-\n";
  	    }  	       	    
  	    
  	    my $overlap=int($prev[6]-$t[5]+1); 
  	    
  	    my %alnSites1;
  	    while ($prev[13]=~m/\((\d+)\,(\d+)\)/g) {
          my ($rsite,$qsite)=($1,$2);	                
          $alnSites1{$rsite}=$qsite;                
        }                
        my %alnSites2;
  	    while ($t[13]=~m/\((\d+)\,(\d+)\)/g) {
          my ($rsite,$qsite)=($1,$2);	                
          $alnSites2{$rsite}=$qsite;                
        }        
        
        my @rsites1=sort {$a<=>$b} keys %alnSites1;
        my @rsites2=sort {$a<=>$b} keys %alnSites2;  
        
        my ($olapEnd1,$olapEnd2);
        if ( ($prev[7] eq "+") and ($t[7] eq "+") )  {
           if ($alnSites1{$rsites2[0]}) {
             $olapEnd1=$alnSites1{$rsites2[0]}-1;
             if ($olapEnd1>0) {
               $olapEnd1=$qmap{$prev[1]}{$olapEnd1};
             }             	
             $overlapSeq1=substr($rawAssSeq{$ctgID1},$olapEnd1,length($rawAssSeq{$ctgID1})-$olapEnd1);
           }else {
           	 $overlapSeq1=substr($rawAssSeq{$ctgID1},length($rawAssSeq{$ctgID1})-$tailCtg2-$tailCtg1-$overlap,$tailCtg2+$tailCtg1+$overlap);
           }
           if ($alnSites2{$rsites1[$#rsites1]}) {
           	 $olapEnd2=$alnSites2{$rsites1[$#rsites1]}+1;
           	 $olapEnd2=$qmap{$t[1]}{$olapEnd2};
           	 $overlapSeq2=substr($rawAssSeq{$ctgID2},0,$olapEnd2); 
           }else {
           	 $overlapSeq2=substr($rawAssSeq{$ctgID2},0,$tailCtg2+$overlap+$tailCtg1);
           }
        }elsif ( ($prev[7] eq "+") and ($t[7] eq "-") )  {
           if ($alnSites1{$rsites2[0]}) {
             $olapEnd1=$alnSites1{$rsites2[0]}-1;
             if ($olapEnd1>0) {
               $olapEnd1=$qmap{$prev[1]}{$olapEnd1};
             }
             $overlapSeq1=substr($rawAssSeq{$ctgID1},$olapEnd1,length($rawAssSeq{$ctgID1})-$olapEnd1);
           }else {
           	 $overlapSeq1=substr($rawAssSeq{$ctgID1},length($rawAssSeq{$ctgID1})-$tailCtg2-$tailCtg1-$overlap,$tailCtg2+$tailCtg1+$overlap);
           }
           if ($alnSites2{$rsites1[$#rsites1]}) {
           	 $olapEnd2=$alnSites2{$rsites1[$#rsites1]}-1;
           	 if ($olapEnd2>0) {
               $olapEnd2=$qmap{$t[1]}{$olapEnd2};
             }           	 
           	 $overlapSeq2=recomSeq(substr($rawAssSeq{$ctgID2},$olapEnd2,length($rawAssSeq{$ctgID2})-$olapEnd2));
           }else {
           	 $overlapSeq2=recomSeq(substr($rawAssSeq{$ctgID2},length($rawAssSeq{$ctgID2})-$tailCtg2-$tailCtg1-$overlap,$tailCtg2+$tailCtg1+$overlap));
           }
        	        	
        }elsif ( ($prev[7] eq "-") and ($t[7] eq "+") ) {
          if ($alnSites1{$rsites2[0]}) {
             $olapEnd1=$alnSites1{$rsites2[0]}+1;             
             $olapEnd1=$qmap{$prev[1]}{$olapEnd1};	
             $overlapSeq1=recomSeq(substr($rawAssSeq{$ctgID1},0,$olapEnd1));
           }else {
           	 $overlapSeq1=recomSeq(substr($rawAssSeq{$ctgID1},0,$tailCtg2+$tailCtg1+$overlap));
           }
           if ($alnSites2{$rsites1[$#rsites1]}) {
           	 $olapEnd2=$alnSites2{$rsites1[$#rsites1]}+1;
           	 $olapEnd2=$qmap{$t[1]}{$olapEnd2};
           	 $overlapSeq2=substr($rawAssSeq{$ctgID2},0,$olapEnd2); 
           }else {
           	 $overlapSeq2=substr($rawAssSeq{$ctgID2},0,$tailCtg2+$overlap+$tailCtg1);
           }
        }else {
          if ($alnSites1{$rsites2[0]}) {
             $olapEnd1=$alnSites1{$rsites2[0]}+1;             
             $olapEnd1=$qmap{$prev[1]}{$olapEnd1};	
             $overlapSeq1=recomSeq(substr($rawAssSeq{$ctgID1},0,$olapEnd1));
           }else {
           	 $overlapSeq1=recomSeq(substr($rawAssSeq{$ctgID1},0,$tailCtg2+$tailCtg1+$overlap));
           }
           if ($alnSites2{$rsites1[$#rsites1]}) {
           	 $olapEnd2=$alnSites2{$rsites1[$#rsites1]}-1;           	 
           	 if ($olapEnd2>0) {
               $olapEnd2=$qmap{$t[1]}{$olapEnd2};
             }           	 
           	 $overlapSeq2=recomSeq(substr($rawAssSeq{$ctgID2},$olapEnd2,length($rawAssSeq{$ctgID2})-$olapEnd2));
           }else {
           	 $overlapSeq2=recomSeq(substr($rawAssSeq{$ctgID2},length($rawAssSeq{$ctgID2})-$tailCtg2-$tailCtg1-$overlap,$tailCtg2+$tailCtg1+$overlap));
           }
        }       
  	    
  	    my $gap=int($hybridGap-$tailCtg1-$tailCtg2);  	
  	    my $omDist=int($t[5]-$prev[6]);    
  	    print OUT1 "$t[2]\t$prev[11]\t$prev[1]\t$ctgID1\t$prev[10]\t$prev[3]\t$prev[4]\t$prev[7]\t$prev[5]\t$prev[6]\t";
  	  	print OUT1                      "$t[1]\t$ctgID2\t$t[10]\t$t[3]\t$t[4]\t$t[7]\t$t[5]\t$t[6]\t$tailCtg1\t$tailCtg2\t$omDist\t";
  	  	print OUT1 "CMAP-overlap:$gap\t";
  	  	
  	  	my ($nctgID1,$nctgID2)=($ctgID1,$ctgID2);
  	    $nctgID1=~s/\|quiver//;
  	    $nctgID2=~s/\|quiver//;
  	  	
  	  	my $ov=abs($gap);
  	  	my $outfa="$outdir/overlap/$t[2]\_$nctgID1\_$prev[7]_$nctgID2\_$t[7]_cmapOverlap_$ov.fa";
  	    open OUT2,">$outfa";
  	    print OUT2 ">$t[2]\_$nctgID1\_$prev[7]\n$overlapSeq1\n";
  	    print OUT2 ">$t[2]\_$nctgID2\_$t[7]\n$overlapSeq2\n";
  	    close OUT2;
  	    system("makeblastdb -in $outfa -dbtype nucl");
  	    my $blastout="$outdir/overlap/$t[2]\_$nctgID1\_$prev[7]_$nctgID2\_$t[7]_cmapOverlap_$ov.blastout";
  	    system("blastn -query $outfa -db $outfa -outfmt 6 -out $blastout -evalue 0.00001");    
  	    
  	    open IN2,$blastout or die "no $blastout";
  	    my $ctg2start=0;
  	    while (<IN2>) {
  	      chomp;
  	      my @tt=split /\t/;
  	      next if ($tt[0] eq $tt[1]);
  	      #print "$tt[0]\t$nctgID1\t$tt[2]\t$tt[7]\t$tt[8]\n";exit;
  	      if ( ($tt[2]>90) and (length($overlapSeq1)-$tt[7]<5) and ($tt[8]<5) ) {
  	        $ctg2start=$tt[9]; 
  	        print "blast result shows CMAP-overlap: $_\n";
  	        last; 	
  	      }  	      
  	    }
  	    close IN2;
  	    print OUT1 "-$ctg2start\n";
  	    my $l1=length($overlapSeq1);
  	    my $l2=length($overlapSeq2);
  	    print "CMAP-overlap: $t[2]\t$ctgID1\t$ctgID2\t$omDist\t$tailCtg1\t$tailCtg2\t$gap\tsubstr-Ctg2start:$ctg2start\n";
  	    $scafCmapOvp++;
  	  	$scafCmapOvpLeng+=$ctg2start;
  	    if ($t[7] eq "+") {
  	      $assSeq{$ctgID2}=substr($assSeq{$ctgID2},$ctg2start,length($assSeq{$ctgID2})-$ctg2start);
  	      
  	    }else {
  	      $assSeq{$ctgID2}=substr($assSeq{$ctgID2},0,length($assSeq{$ctgID2})-$ctg2start);  	      
  	    } 	  	 	      	  	
  	  }	
  	  
  	}else {
  	  my $preCtg=$ctgID{$prev[1]};
  	  if (!$scafseq{$prev[2]}) {
  	    my $scafEnd=length($assSeq{$preCtg});
  	    my $rawLeng=length($rawAssSeq{$preCtg});
  	    my $sleng=length($assSeq{$preCtg});
  	    print OUT3 "$sID-$prev[2]\t$preCtg\t$rawLeng\t$sleng\t1\t$scafEnd\t";  	    
  	  }else {
  	    my $scafStart=length($scafseq{$prev[2]})+1;
  	    my $scafEnd=$scafStart+length($assSeq{$preCtg})-1;
  	    my $rawLeng=length($rawAssSeq{$preCtg});
  	    my $sleng=length($assSeq{$preCtg});
  	    print OUT3 "$sID-$prev[2]\t$preCtg\t$rawLeng\t$sleng\t$scafStart\t$scafEnd\t";  	    
  	  } 
  	  
  	  if ($prev[7] eq "+") {  	  	
  	  	$scafseq{$prev[2]}.=$assSeq{$preCtg};
  	  	print OUT3 "+\n";
  	  }else {
  	  	$scafseq{$prev[2]}.=recomSeq($assSeq{$preCtg});
  	  	print OUT3 "-\n";
  	  }      
  	  $scafseq{$t[2]}="";   	  	
  	}
  	@prev=@t;  	   	  	
  }  
}

my $preCtg=$ctgID{$prev[1]};
if (!$scafseq{$prev[2]}) {
  my $scafEnd=length($assSeq{$preCtg});
  my $rawLeng=length($rawAssSeq{$preCtg});
  my $sleng=length($assSeq{$preCtg});
  print OUT3 "$sID-$prev[2]\t$preCtg\t$rawLeng\t$sleng\t1\t$scafEnd\t";  
}else {
  my $scafStart=length($scafseq{$prev[2]})+1;
  my $scafEnd=$scafStart+length($assSeq{$preCtg})-1;
  my $rawLeng=length($rawAssSeq{$preCtg});
  my $sleng=length($assSeq{$preCtg});
  print OUT3 "$sID-$prev[2]\t$preCtg\t$rawLeng\t$sleng\t$scafStart\t$scafEnd\t";  	 
} 

if ($prev[7] eq "+") {  	  	
  $scafseq{$prev[2]}.=$assSeq{$preCtg};
  print OUT3 "+\n";
}else {
  $scafseq{$prev[2]}.=recomSeq($assSeq{$preCtg});
  print OUT3 "-\n";
}
close OUT3;

open OUT2,">$outdir/hybrid.map.scaffolds.fasta";
my $n=0;my $m=0;
foreach my $id (sort keys %assSeq) {
  #if ($id eq "scaf-185") {
  	#print OUT2 ">$sID-$id\n$assSeq{$id}\n";next;
  #}
  if (!$ctgUsed{$id}) {
  	if ($inCtg{$id}) {
  	  print "inside contig: $id\n";
  	  $m++;
  	  next;
    }
  	print OUT2 ">$sID-$id\n$assSeq{$id}\n";
  }else {
  	$n++;
  }
}
print "scaffolded contigs: $n\n";

foreach my $id (sort keys %scafseq) { 
  print OUT2 ">$sID-$id\n$scafseq{$id}\n";
}

print "scafGapNum\tscafCMAPoverlap\tscafTailOverlap\tscafGapLeng\tscaf-CMAP-overlap-Length\tscaf-Tail-overlap-Length\n";
print "$scafGap\t$scafCmapOvp\t$scafTailOvp\t$scafGapLeng\t$scafCmapOvpLeng\t$scafTailOvpLeng\n";




sub hashAssKey {
  my ($in)=@_;
  open IN,$in;
  my %ctg;
  while (<IN>) {
    chomp;
    next if (!/^\d/);
    my @t=split /\t/;
    $ctg{$t[0]}=$t[1];
  }
  close IN;
  return \%ctg;
}

sub hashAssSeq {
  my ($in)=@_;
  my $seqin=new Bio::SeqIO(-file=>"$in",-format=>"fasta");
  my %ass;
  while (my $seqobj=$seqin->next_seq()) {
  	my $id=$seqobj->id();
    $ass{$id}=$seqobj->seq();
  }
  return \%ass;
}


sub hashCmap {
  my ($in)=@_;
  open IN,$in;
  my %map;
  while (<IN>) {
    chomp;
    if (/^\d/) {
      my ($id,$leng,$ns,$site,$pos)=(split /\t/)[0,1,2,3,5];
      $map{$id}{"leng"}=$leng;
      $map{$id}{"ns"}=$ns;
      $map{$id}{$site}=$pos;
    }
  }
  close IN;
  return \%map;
}

sub hashHybridCmap {
	 my ($in)=@_;
  open IN,$in;
  my %map;
  while (<IN>) {
    chomp;
    if (/^\d/) {
      my ($id,$leng,$ns,$site,$pos)=(split /\t/)[0,1,2,3,5];      
      $map{$id}{"site"}{$site}=$pos;
      $map{$id}{"pos"}{$pos}=$site;
    }
  }
  close IN;
  return \%map;
}


sub recomSeq {
  my ($seq)=@_;
  my $recom=reverse($seq);
  $recom=~tr/ATGCatgc/TACGtacg/;
  return $recom;  
}



sub xmapFilter {
  my ($seqCmap,$OMcmap,$xmap,$xmapFlt)=@_;
  open IN,"grep -v '#' $xmap|sort -k3,3n -k6,6n -k7,7nr|";
  my ($ref,$start,$end);
  my %insideAln;
  my $maxHang=5;
  my %BNGcmap=%{$OMcmap};
  my %cmap=%{$seqCmap};
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
      	$insideAln{$t[2].$t[5].$t[6]}=1;
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
      	$insideAln{$t[2].$t[5].$t[6]}=1;
      	next;
      }else {
      	($ref,$start,$end)=@t[1,3,4];
      } 	         
    }
  } 
  close IN;
  system("rm tmp");
  
  open IN,"grep -v '#' $xmap|sort -k3,3n -k6,6n |";
  open OUT,">$xmapFlt.tmp";
  my %insideCtg;
  my %scafNum;
  while (<IN>) {
  	chomp;
  	my @t=split /\t/;
  	if ($insideAln{$t[2].$t[5].$t[6]}) {
  	  my @aln;
  	  while ($t[13]=~m/\((\d+)\,(\d+)\)/g) {
        my ($rsite,$qsite)=($1,$2);	
        push @aln,[($rsite,$qsite)];  
      } 
      my $BNGleftHang=$aln[0][0]-1;
  	  my $BNGrightHang=$BNGcmap{$t[2]}{"ns"}-$aln[$#aln][0];
  	  
  	  my $NGSleftHang=$aln[0][1]-1;
  	  my $NGSrightHang=$cmap{$t[1]}{"ns"}-$aln[$#aln][1];  	   
      if ($t[7] eq "-") {      	
      	$NGSleftHang=$cmap{$t[1]}{"ns"}-$aln[0][1];
  	    $NGSrightHang=$aln[$#aln][1]-1;  	    
      }
      if ( ($NGSleftHang<=2) and ($NGSrightHang<=2) ) {
      	$insideCtg{$t[1]}=1;
      } 
      next;
  	}
  	if (!$ref) {
      ($ref,$start,$end)=@t[2,5,6]; 
      print OUT "$_\n"; 
      $scafNum{$t[2]}++;    
    }else {
      if ( ($t[2] eq $ref) and ($t[5]>=$start) and ($t[6]<=$end) ) {
      	my $siteEnd=$cmap{$t[1]}{"ns"};
      	$siteEnd=$cmap{$t[1]}{$siteEnd};
      	my $siteStart=$cmap{$t[1]}{"1"};      	
      	if (($t[6]-$t[5])>0.8*($siteEnd-$siteStart) ) {
      	  $insideCtg{$t[1]}=$_;
      	  next;	
      	}        	
      }else {
      	($ref,$start,$end)=@t[2,5,6];
      } 	    
      my @t=split /\t/;
      my @aln;
  	  while ($t[13]=~m/\((\d+)\,(\d+)\)/g) {
        my ($rsite,$qsite)=($1,$2);	
        push @aln,[($rsite,$qsite)];  
      } 
      my $BNGleftHang=$aln[0][0]-1;
  	  my $BNGrightHang=$BNGcmap{$t[2]}{"ns"}-$aln[$#aln][0];
  	  
  	  my $NGSleftHang=$aln[0][1]-1;
  	  my $NGSrightHang=$cmap{$t[1]}{"ns"}-$aln[$#aln][1];  	   
      if ($t[7] eq "-") {      	
      	$NGSleftHang=$cmap{$t[1]}{"ns"}-$aln[0][1];
  	    $NGSrightHang=$aln[$#aln][1]-1;  	    
      }     
      if (($NGSrightHang>$maxHang) and ($NGSleftHang>$maxHang)) {
     	print "NGS-R,NGS-L filter: $_\n";
      }elsif (($NGSrightHang>=$maxHang) and ($BNGrightHang>=$maxHang)) {
     	print "NGS-R,BNG-R filter: $_\n";
      }elsif (($NGSleftHang>=$maxHang) and ($BNGleftHang>=$maxHang)) {
     	print "NGS-L,BNG-L filter: $_\n";
      }else {
        $scafNum{$t[2]}++;
        print OUT "$_\n";	
      }      
    }
  }
  close OUT;close IN;
  
  open IN,"grep -v '#' $xmapFlt.tmp|sort -k3,3n -k6,6n -k7,7nr|";
  open OUT,">$xmapFlt";
  while (<IN>) {
  	chomp;
  	my @t=split /\t/;
  	next if ($insideAln{$t[2].$t[5].$t[6]});
  	if ($scafNum{$t[2]}>1) {
  	  print OUT "$_\n";
  	}
  }
  close OUT;close IN;
  
  
  return \%insideCtg;
}


sub nickScafGaps {
  my ($refHash,$start,$end,$tail1,$tail2,$nick)=@_;
  my %refHash=%{$refHash};
  my $nick1=$refHash{"pos"}{$start};
  my $nick2=$refHash{"pos"}{$end};
  my $gapSeq;
  my $dist1=$refHash{"site"}{$nick1+1}-$start;  	
  my $dist2=$end-$refHash{"site"}{$nick2-1};
  #print "tail: $tail1\t$tail2\t dist:$dist1\t$dist2\t Nick\t$nick1\t$nick2\n";
  if ($nick1==$nick2-2) {  	
    if ($dist1>=$tail1) {
      my $gap=$dist1-$tail1;
      $gapSeq="N" x  $gap;
      $gapSeq=$gapSeq.$nick;
    }else {
      $gapSeq=$nick;
    }
    
    if ($dist2>=$tail2) {
      my $gap=$dist2-$tail2;
      my $nseq="N" x $gap;
      $gapSeq=$gapSeq.$nseq;
    }else {    	      
    } 
    
  }elsif ($nick1<$nick2-2) {
    if ($dist1>=$tail1) {
      my $gap=$dist1-$tail1;
      $gapSeq="N" x  $gap;
      $gapSeq=$gapSeq.$nick;
    }else {
      $gapSeq=$nick;
    }
    #print "nick: $nick1\t$nick2\n";
    foreach my $k ($nick1+2..$nick2-1) {
      my $dist=$refHash{"site"}{$k}-$refHash{"site"}{$k-1}-5;
      my $nseq="N" x $dist;
      $gapSeq=$gapSeq.$nseq.$nick; 
      #print "K:   $k\t$dist\n";           
    }
    if ($dist2>=$tail2) {
      my $gap=$dist2-$tail2;
      my $nseq="N" x $gap;
      $gapSeq=$gapSeq.$nseq;
    }else {
    	      
    }  	
  }else {
  	###weird
  	
  }
  #print "GapSeq: $gapSeq\n";
  #exit;
  return $gapSeq;
  
}