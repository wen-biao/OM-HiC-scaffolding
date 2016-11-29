#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: pacbio.break.ctg.pl
#
#        USAGE: ./pacbio.break.ctg.pl  
#
#  DESCRIPTION: break assembly sequences and consensus maps with input of break information without filtering steps 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wen-Biao Jiao (), 
# ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
#      VERSION: 1.0
#      CREATED: 05/10/2016 11:09:11 AM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;


my ($seqKey,$seqFasta,$seqbreaks,$BNGbreaks,$outputdir,$outdir,$help);

GetOptions (  
  'key|k=s'    => \$seqKey,
  'seq|f=s' => \$seqFasta,
  'seqMis|s=s' =>\$seqbreaks,
  'OMmis|m=s' =>\$BNGbreaks,    
  'aln|o=s'    => \$outputdir,
  'outdir|d=s'	=>\$outdir,
  'help|h+'      => \$help,
);

my $usage="
Usage: $0 seqKeyFile seqFastaFile seqMisInfor omMisInfor inital_alignment_dir output_dir 
Options:
  -key|-k In silico CMAP key file of sequence assembly
  -seq|-f assembly fasta file
  -seqMis|-s misassembly information file found by script OM.find.conflict-regions or OM.merge.conflict-regions.pl
  -OMmis|-m  optical consensus maps misassembly information file found by script OM.find.conflict-regions or OM.merge.conflict-regions.pl
  -aln|-o output directory of initial alignment between assembly sequences and optical consensus maps 
  -outDir|-d output directory of running
  -help|-h print this help
";


die $usage if (defined $help);
die $usage if ((!$seqKey) or (!$seqFasta) or (!$seqbreaks) or (!$BNGbreaks) or (!$outputdir));
die "can't find the file" if ((!-e $seqKey) or (! -e $seqFasta) or (! -e $seqbreaks) or (! -e $BNGbreaks) );

system("mkdir -p $outdir");
my $seqcmap="$outputdir/align1/align1_r.cmap";  ###Don't use the original CMAP from $seqfa, some close two nick sites were merged into one site after running RefAligner  
my $seqcmaphash=hashCmap($seqcmap);
my $seqcmapTrans=hashNGScmapKey($seqKey);

my $BNGcmap="$outputdir/align1/align1_q.cmap";
my $BNGcmaphash=hashCmap($BNGcmap);

my $seqbrokenCmap="$outdir/seq.assembly.conflicts.broken_BspQI_0Kb_0labels.cmap";
breakCmap($seqbreaks,$seqcmap,$seqcmaphash,$seqbrokenCmap,0);
my $BNGbrokenCmap="$outdir/BNG.assembly.conflicts.broken.cmap";
breakCmap($BNGbreaks,$BNGcmap,$BNGcmaphash,$BNGbrokenCmap,1);

my $seqBrokenFa="$outdir/seq.assembly.conflicts.broken.fasta";
my $newSeqKey="$outdir/seq.broken.new.key.txt";
my $seq=hashFastaSeq($seqFasta);
breakCtgSeq($seqbreaks,$seq,$seqcmapTrans,$seqBrokenFa,$newSeqKey);


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


sub hashFastaSeq {
  my ($in)=@_;
  my $seqin=new Bio::SeqIO(-format=>"fasta",-file=>"$in");
  my %seq;
  while (my $seqobj=$seqin->next_seq()) {
  	my $id=$seqobj->id();
  	my $seq=$seqobj->seq();
  	$seq{$id}=$seq;
  }
  #print scalar keys %seq,"\n";
  return \%seq;
}




sub breakCmap {
  my ($break,$cmap,$hashCmap,$brokenCmap,$isBNG)=@_;
  my ($cmapID,$br,$bp);  
  open IN,$break;
  my @breaks;my %breaks;
  while (<IN>) {
    chomp;
    next if (/^#/);
    my @t=split /\t/;
    shift @t if (!/^\d/);
    if (!$cmapID) {
      ($cmapID,$br,$bp)=@t[0,4,5];
      push @breaks,[($br,$bp)];
    }elsif ($cmapID ne $t[0]) {      
      $breaks{$cmapID}=[@breaks];
      ($cmapID,$br,$bp)=@t[0,4,5];
      @breaks=();
      push @breaks,[($br,$bp)];	
    }else {
      ($cmapID,$br,$bp)=@t[0,4,5];
      push @breaks,[($br,$bp)];
    }
  }
  $breaks{$cmapID}=[@breaks];
  @breaks=();
  
  my %cmap=%{$hashCmap};
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
  
}


sub breakCtgSeq {
  my ($break,$seq,$oldKey,$out,$newKey)=@_;
  
  my %seq=%{$seq};
  my $cmapID; my $pos; my $i;
  my %oldKey=%{$oldKey};
  my %transKey;
  foreach my $id (sort keys %oldKey) {
  	$transKey{$oldKey{$id}}=$id;
  }
  
  open OUT,">$out"; 
  open OUT2,">$newKey";
  my %br;
  open IN,$break;
  while (<IN>) {
    chomp;
    next if (/#/);
    my @t=split /\t/;
    $br{$t[1]}=1;
    if (!$cmapID) {
      $cmapID=$t[1];
      $pos=int($t[5]);
      $i=0;
      my $s=substr($seq{$t[1]},0,int($t[5]));
      print OUT ">$cmapID\_0\n$s\n";
      my $leng=length($s);
      print OUT2 "$transKey{$cmapID}00000\t$cmapID\_0\t$leng\n";
      
    }elsif ($cmapID ne $t[1]) {
      my $s=substr($seq{$cmapID},$pos,length($seq{$cmapID})-$pos);
      $i++;
      my $leng=length($s);
      print OUT ">$cmapID\_$i\n$s\n";
      print OUT2 "$transKey{$cmapID}0000$i\t$cmapID\_$i\t$leng\n";
      $cmapID=$t[1];
      $i=0;
      $pos=int($t[5]);
      $s=substr($seq{$cmapID},0,int($t[5]));
      print OUT ">$cmapID\_0\n$s\n";
      $leng=length($s);
      print OUT2 "$transKey{$cmapID}00000\t$cmapID\_0\t$leng\n";	      	
    }else {
      my $s=substr($seq{$cmapID},$pos,int($t[5])-$pos);
      $i++;
      print OUT ">$cmapID\_$i\n$s\n";
      my $leng=length($s);
      print OUT2 "$transKey{$cmapID}0000$i\t$cmapID\_$i\t$leng\n";
      $pos=int($t[5]);
    }
  }
  my $s=substr($seq{$cmapID},$pos,length($seq{$cmapID})-$pos);
  $i++;
  print OUT ">$cmapID\_$i\n$s\n";
  my $leng=length($s);
  print OUT2 "$transKey{$cmapID}0000$i\t$cmapID\_$i\t$leng\n";
  
  foreach my $id (sort keys %seq) {
  	if (!$br{$id}) {
  	  print OUT ">$id\n$seq{$id}\n";  	
  	  my $leng=length($seq{$id});  
  	  print OUT2 "$transKey{$id}\t$id\t$leng\n";
  	}
  }
  close IN;close OUT; close OUT2;
}
