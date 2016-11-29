#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: OM.breaks.ctg.pl
#
#        USAGE: ./OM.breaks.ctg.pl  
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
#      CREATED: 08/04/2016 11:29:09 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;

#my ($seqBreaks,$seqFasta,$OMbreaks,$OMcmap,$outdir)=@ARGV;

my ($seqBreaks,$seqFasta,$OMbreaks,$OMcmap,$outdir,$help);
GetOptions (
  'seqBreaks|s=s' => \$seqBreaks, 
  'seqFasta|f=s' => \$seqFasta,
  'OMbreaks|m=s' => \$OMbreaks,
  'OMcmap|c=s' => \$OMcmap,
  'outDir|o=s'    => \$outdir,
  'help|h+'      => \$help,
);

my $usage="
Usage: $0 seqBreaks seqFasta OMbreaks OMcmap outdir 
Options:
  -seqBreaks|-s	File of breaking positions of assembly sequences
  -seqFasta|-f	File of assembly sequences
  -OMbreaks|-m	File of breaking positions of optical consensus maps
  -OMcmap|-c	File of optical consensus maps
  -outDir|-o output directory
  -help|-h print this help
";
die $usage if (defined $help);
die $usage if ( (!-e $seqBreaks) or (!-e $seqFasta) or (!-e $OMbreaks) or (!-e $OMcmap));


system("mkdir -p $outdir");

my $OMcmaphash=hashCmap($OMcmap);
my $OMbrokenCmap="$outdir/OM.broken.cmap";
breakOMcmap($OMbreaks,$OMcmap,$OMcmaphash,$OMbrokenCmap,1);

my $seq=hashFastaSeq($seqFasta);
my $out="$outdir/seq.broken.fasta";
breakCtgSeq2($seqBreaks,$seq,$out);

sub hashFastaSeq {
  my ($in)=@_;
  my $seqin=new Bio::SeqIO(-format=>"fasta",-file=>"$in");
  my %seq;
  while (my $seqobj=$seqin->next_seq()) {
  	my $id=$seqobj->id();
  	my $seq=$seqobj->seq();
  	$seq{$id}=$seq;
  }
  return \%seq;
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


sub breakCtgSeq2 {
  my ($break,$seq,$out)=@_;
  
  my %seq=%{$seq};
  my ($cmapID,$pos,$i);
  
  open OUT,">$out"; 
  my %br;
  open IN,$break;
  while (<IN>) {
    chomp;
    next if (/^#/);
    my @t=split /\t/;
    $br{$t[0]}=1;

    if (!$cmapID) {
      if ($t[3] eq "HRbreaks") {        
        $cmapID=$t[0];
        $i=0;
        my $s=substr($seq{$t[0]},0,$t[1]);
        print OUT ">$cmapID\_0\n$s\n";
        $pos=$t[2];
        next;
      }
      $cmapID=$t[0];
      $i=0;	    	      
      if ($t[6]=~m/\,/) {
      	my @tt=split /\,/,$t[7];
        my $s=substr($seq{$t[0]},0,int($tt[0])+10); ###add 10bp to include the last nick site        
        print OUT ">$cmapID\_0\n$s\n";
        
        $pos=int($tt[$#tt]);
        next if ($tt[0]==$tt[$#tt]);
        $s=substr($seq{$t[0]},int($tt[0]),int($tt[$#tt]-$tt[0]));
        print OUT ">$cmapID\_01\n$s\n";
                                
      }else {
        if ($t[7]<$t[9]) {
          my $s=substr($seq{$t[0]},0,int($t[7])+10);          
          print OUT ">$cmapID\_0\n$s\n";
          
          $s=substr($seq{$t[0]},int($t[7]),int($t[9]-$t[7]));
          print OUT ">$cmapID\_01\n$s\n";
          $pos=int($t[9]);	
        }else {
          my $s=substr($seq{$t[0]},0,int($t[9])+10);          
          print OUT ">$cmapID\_0\n$s\n";
          
          $s=substr($seq{$t[0]},int($t[9]),int($t[7]-$t[9]));
          print OUT ">$cmapID\_01\n$s\n";
          $pos=int($t[7]);	
        }               
      }                       
    }elsif ($cmapID ne $t[0]) {
      my $s=substr($seq{$cmapID},$pos,length($seq{$cmapID})-$pos);
      $i++;
      print OUT ">$cmapID\_$i\n$s\n";
      
      if ($t[3] eq "HRbreaks") {        
        $cmapID=$t[0];
        $i=0;
        my $s=substr($seq{$t[0]},0,$t[1]);
        print OUT ">$cmapID\_0\n$s\n";
        $pos=$t[2];
        next;
      }
      $cmapID=$t[0];
      $i=0;
      if ($t[6]=~m/\,/) {
      	my @tt=split /\,/,$t[7];
        my $s=substr($seq{$t[0]},0,int($tt[0])+10); ###add 10bp to include the last nick site        
        print OUT ">$cmapID\_0\n$s\n";
        
        $pos=int($tt[$#tt]); 
        next if ($tt[0]==$tt[$#tt]);
        $s=substr($seq{$t[0]},int($tt[0]),int($tt[$#tt]-$tt[0]));
        print OUT ">$cmapID\_01\n$s\n";                       
      }else {
        if ($t[7]<$t[9]) {
          my $s=substr($seq{$t[0]},0,int($t[7])+10);          
          print OUT ">$cmapID\_0\n$s\n";
          
          $s=substr($seq{$t[0]},int($t[7]),int($t[9]-$t[7]));
          print OUT ">$cmapID\_01\n$s\n";
          $pos=int($t[9]);	
        }else {
          my $s=substr($seq{$t[0]},0,int($t[9])+10);          
          print OUT ">$cmapID\_0\n$s\n";
          
          $s=substr($seq{$t[0]},int($t[9]),int($t[7]-$t[9]));
          print OUT ">$cmapID\_01\n$s\n";
          $pos=int($t[7]);	
        } 	
      }          	      	
    }else {
      $i++;
      if ($t[3] eq "HRbreaks") {                
        my $s=substr($seq{$t[0]},0,$t[1]);
        print OUT ">$cmapID\_$i\n$s\n";
        $pos=$t[2];
        next;
      }
      if ($t[6]=~m/\,/) {
      	my @tt=split /\,/,$t[7];
        my $s=substr($seq{$t[0]},$pos-10,int($tt[0])-$pos+20); ###add 10bp to include the last nick site        
        print OUT ">$cmapID\_$i\n$s\n";
        
        $pos=int($tt[$#tt]); 
        next if ($tt[0]==$tt[$#tt]);
        my $j=$i+1;
        $s=substr($seq{$t[0]},int($tt[0]),int($tt[$#tt]-$tt[0]));
        print OUT ">$cmapID\_$i$j\n$s\n";                       
      }else {
      	my $j=$i+1;
        if ($t[7]<$t[9]) {
          my $s=substr($seq{$t[0]},$pos-10,int($t[7])-$pos+20);          
          print OUT ">$cmapID\_$i\n$s\n";                              
          $s=substr($seq{$t[0]},int($t[7]),int($t[9]-$t[7]));
          print OUT ">$cmapID\_$i$j\n$s\n";

          $pos=int($t[9]);	
        }else {
          my $s=substr($seq{$t[0]},$pos-10,int($t[9])-$pos+20);          
          print OUT ">$cmapID\_$i\n$s\n";
          
          $s=substr($seq{$t[0]},int($t[9]),int($t[7]-$t[9]));
          print OUT ">$cmapID\_$i$j\n$s\n";

          $pos=int($t[7]);	
        } 	
      } 
    }
  }
  my $s=substr($seq{$cmapID},$pos,length($seq{$cmapID})-$pos);
  $i++;
  print OUT ">$cmapID\_$i\n$s\n";
  
  foreach my $id (sort keys %seq) {
  	if (!$br{$id}) {
  	  print OUT ">$id\n$seq{$id}\n";  	
  	}
  }
  close IN;close OUT;
}



sub breakOMcmap {
  my ($break,$cmap,$hashCmap,$brokenCmap,$isBNG)=@_;
  my ($cmapID,$br,$bp);  
  my %cmap=%{$hashCmap};
  open IN,$break;
  my @breaks;my %breaks;
  while (<IN>) {
    chomp;
    next if (/^#/);
    my @t=split /\t/;
    shift @t if (!/^\d/);
    if (!$cmapID) {
      ($cmapID,$br,$bp)=@t[0,6,7];
      if ($t[6]=~m/\,/) {
      	my @tt=split /\,/,$t[6];
      	$br=($tt[0]+$tt[$#tt])/2;
      	$bp=int($cmap{$t[0]}{$tt[0]}+$cmap{$t[0]}{$tt[$#tt]})/2;
      }else {
      	$br=($t[6]+$t[8])/2;
      	$bp=int($cmap{$t[0]}{$t[6]}+$cmap{$t[0]}{$t[6]})/2;
      }
      push @breaks,[($br,$bp)];
    }elsif ($cmapID ne $t[0]) {
      $breaks{$cmapID}=[@breaks];
      ($cmapID,$br,$bp)=@t[0,6,7];
      @breaks=();
      if ($t[6]=~m/\,/) {
      	my @tt=split /\,/,$t[6];
      	$br=($tt[0]+$tt[$#tt])/2;
      	$bp=int($cmap{$t[0]}{$tt[0]}+$cmap{$t[0]}{$tt[$#tt]})/2;
      }else {
      	$br=($t[6]+$t[8])/2;
      	$bp=int($cmap{$t[0]}{$t[6]}+$cmap{$t[0]}{$t[6]})/2;
      }
      push @breaks,[($br,$bp)];	
    }else {
      ($cmapID,$br,$bp)=@t[0,6,7];
      if ($t[6]=~m/\,/) {
      	my @tt=split /\,/,$t[6];
      	$br=($tt[0]+$tt[$#tt])/2;
      	$bp=int($cmap{$t[0]}{$tt[0]}+$cmap{$t[0]}{$tt[$#tt]})/2;
      }else {
      	$br=($t[6]+$t[8])/2;
      	$bp=int($cmap{$t[0]}{$t[6]}+$cmap{$t[0]}{$t[6]})/2;
      }
      push @breaks,[($br,$bp)];
    }
  }
  $breaks{$cmapID}=[@breaks];
  @breaks=();
  
  close IN;
  open IN,$cmap;
  open OUT,">$brokenCmap";
  $cmapID="";
  my $bID="";
  my $i=0; my $p;my $flag=0;
  while (<IN>) {
    chomp;
    if (!/^\d/) {
      print OUT "$_\n";
    }else {
      my @t=split /\t/;      
      if (!$breaks{$t[0]}) {
      	if ($p) {
      	  print OUT "$p\n";
      	  $p="";
      	}      	
      	print OUT "$_\n";
      	next;
      }                  
      if (!$cmapID) {
      	$cmapID=$t[0];
      	my @breaks=@{$breaks{$t[0]}};
      	$t[0]=$t[0]."00000";
      	#$t[1]=sprintf("%.1f",$breaks[0][1]);      	
      	$t[2]=int($breaks[0][0]);
      	$t[1]=int( ($cmap{$cmapID}{$t[2]}+$cmap{$cmapID}{$t[2]+1})/2 );
      	$t[1]=sprintf("%.1f",$t[1]); 
      	$t[5]=sprintf("%.1f",$t[5]);
      	$p=join("\t",@t);
      	#print "$p\n";exit;
      	$i=0;
      	      	
      }elsif ($cmapID ne $t[0]) {
      	if ($p) {
      	  print OUT "$p\n";
      	  my @tmp=(split /\t/,$p);
      	  $tmp[4]=0;
      	  $tmp[5]=$tmp[1];
      	  $tmp[6]=sprintf("%.1f",0.0);
      	  $tmp[7]=sprintf("%.1f",1.0);
      	  $tmp[8]=sprintf("%.1f",1.0);
      	  if ($isBNG) {
      	    $tmp[9]=sprintf("%.3f",0.000);
      	    $tmp[10]=sprintf("%.3f",0.000);
      	    $tmp[11]=0;	
      	  }      	  
      	  $p=join("\t",@tmp);
      	  print OUT "$p\n";		
      	}else {      	  
      	}
      	$cmapID=$t[0];      	
      	my @breaks=@{$breaks{$t[0]}};
      	$t[0]=$t[0]."00000";
      	#$t[1]=sprintf("%.1f",$breaks[0][1]);      	
      	$t[2]=int($breaks[0][0]);
      	$t[1]=int( ($cmap{$cmapID}{$t[2]}+$cmap{$cmapID}{$t[2]+1})/2 );
      	$t[1]=sprintf("%.1f",$t[1]);
      	$t[5]=sprintf("%.1f",$t[5]);
      	$p=join("\t",@t);
      	$i=0;
      	
      }else {
      	my @breaks=@{$breaks{$t[0]}};
      	
      	while ($i<=$#breaks) {
      	  if ($t[3]<=$breaks[$i][0]) {
      	  	
      	  	print OUT "$p\n";
      	  	$t[0]=$t[0]."0000".$i;
      	  	if ($i==0) {
      	      #$t[1]=sprintf("%.1f",$breaks[$i][1]);   
      	      $t[2]=int($breaks[$i][0]);  
      	      $t[1]=int( ($cmap{$cmapID}{$t[2]}+$cmap{$cmapID}{$t[2]+1})/2 );
      	      $t[1]=sprintf("%.1f",$t[1]);
      	      #print "$0000\t$t[2]\n";exit;    	      
      	      $t[3]=$t[3];
      	      $t[5]=sprintf("%.1f",$t[5]);
      	    }else {
      	      #$t[1]=sprintf("%.1f",$breaks[$i][1]-$breaks[$i-1][1]);
      	      $t[2]=int($breaks[$i][0])-int($breaks[$i-1][0]);
      	      my $br1=int($breaks[$i-1][0]);
      	      my $br2=int($breaks[$i][0]);
      	      my $pos1=int( ($cmap{$cmapID}{$br1}+$cmap{$cmapID}{$br1+1})/2 );
      	      my $pos2=int( ($cmap{$cmapID}{$br2}+$cmap{$cmapID}{$br2+1})/2 );
      	      $t[1]=sprintf("%.1f",$pos2-$pos1);      	            	      
      	      $t[3]=$t[3]-int($breaks[$i-1][0]);
      	      $t[5]=sprintf("%.1f",$t[5]-$pos1);
      	    } 
      	    $p=join("\t",@t);
      	    last;      	          	  	     	  	
      	  }else {
      	  	print OUT "$p\n";
      	  	my @tmp=(split /\t/,$p);
      	  	$tmp[3]+=1;
      	    $tmp[4]=0;
      	    $tmp[5]=$tmp[1];
      	    $tmp[6]=sprintf("%.1f",0.0);
      	    $tmp[7]=sprintf("%.1f",1.0);
      	    $tmp[8]=sprintf("%.1f",1.0);
      	    if ($isBNG) {
      	      $tmp[9]=sprintf("%.3f",0.000);
      	      $tmp[10]=sprintf("%.3f",0.000);
      	      $tmp[11]=0;	
      	    } 
      	    $p=join("\t",@tmp);     	  	
      	  	$i++;
      	  	
      	  	if ($i>$#breaks) {
      	  	  last;
      	  	}
      	  	print OUT "$p\n";
      	  	$t[0]=$t[0]."0000".$i;
      	  	#print "$_\n";print "$t[0]\n";exit;
      	  	if ($i==0) {
      	      #$t[1]=sprintf("%.1f",$breaks[$i][1]);
      	      $t[2]=int($breaks[$i][0]);
      	      $t[1]=int( ($cmap{$cmapID}{$t[2]}+$cmap{$cmapID}{$t[2]+1})/2 );
      	      $t[1]=sprintf("%.1f",$t[1]);      	      
      	      $t[3]=$t[3];
      	      $t[5]=sprintf("%.1f",$t[5]);
      	    }else {
      	      #$t[1]=sprintf("%.1f",$breaks[$i][1]-$breaks[$i-1][1]);
      	      
      	      $t[2]=int($breaks[$i][0])-int($breaks[$i-1][0]); 
      	      my $br1=int($breaks[$i-1][0]);
      	      my $br2=int($breaks[$i][0]);
      	      my $pos1=int( ($cmap{$cmapID}{$br1}+$cmap{$cmapID}{$br1+1})/2 );
      	      my $pos2=int( ($cmap{$cmapID}{$br2}+$cmap{$cmapID}{$br2+1})/2 );
      	      $t[1]=sprintf("%.1f",$pos2-$pos1);      	      
      	      $t[3]=$t[3]-int($breaks[$i-1][0]);
      	      $t[5]=sprintf("%.1f",$t[5]-$pos1);
      	    } 
      	    $p=join("\t",@t); 
      	    last;    	  	
      	  }
      	}
      	if ($i==$#breaks+1) {
      	  print OUT "$p\n";      	  
      	  $t[0]=$t[0]."0000$i";
      	  #$t[1]=sprintf("%.1f",$t[1]-$breaks[$#breaks][1]);      	  
      	  $t[2]=$t[2]-int($breaks[$#breaks][0]);
      	  my $br1=int($breaks[$#breaks][0]);
      	  my $prepos=int( ($cmap{$cmapID}{$br1}+$cmap{$cmapID}{$br1+1})/2 );
      	  $t[1]=sprintf("%.1f",$cmap{$cmapID}{"leng"}-$prepos);
      	  
      	  $t[3]=$t[3]-int($breaks[$#breaks][0]);
      	  $t[5]=sprintf("%.1f",$t[5]-$prepos);
      	  $p=join("\t",@t);
      	  if ($t[4]==0) {
      	  	print OUT "$p\n";
      	  	$p="";
      	  }      	        		
      	}      	      	
      }      
    }
  }
  print OUT "$p\n" if ($p);
  close OUT;close IN;
}






