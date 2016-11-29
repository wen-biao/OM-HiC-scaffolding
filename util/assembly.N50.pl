#!/usr/bin/perl
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;


my ($contig,$stats,$est,$help)=@ARGV;

GetOptions (  
  'seq|f=s' => \$seqFasta,
  'out|o=s'	=>\$stats,
  'gsize|s=s'	=>\$est,
  'help|h+'      => \$help,
);

my $usage="
Usage: $0  assembly N50_output estimated_genome_size 
Options:
  -seq|-f  assembly fasta file
  -out|-o  assembly N50 metrics output
  -gsize|s estimated genome size
  -help|-h print this help
";


die $usage if (defined $help);
die $usage if ((!$contig) or (!$est) or (!$stats));

open OUT,">$stats";
my $seqin=new Bio::SeqIO('-format'=>"fasta",'-file'=>"$contig");
while (my $seqobj=$seqin->next_seq() ){
    $id=$seqobj->id;
    $length=$seqobj->length;
    if ($length>=200) {
       $sum+=$length;
       push @len,$length;
       $n++;
    }
}

print OUT "TotalNumber\t$n\n";
print OUT "TotalLenth\t$sum\n";

$len25=$sum*0.25;
$len50=$sum*0.5;
$len75=$sum*0.75;
$len90=$sum*0.9;

$est25=$est*0.25;
$est50=$est*0.5;
$est75=$est*0.75;
$est90=$est*0.9;


my $l=0;
@so=sort {$b<=>$a} @len;
$con=0;
my $k=0;
my ($n25,$n50,$n75,$n90)=($so[0],$so[0],$so[0],$so[0]);
my ($l25,$l50,$l75,$l90)=(0,0,0,0);



while ($con<$len90) {
  while ($con<$len75) {
    while ($con<$len50) {
  	  while ($con<$len25) {
  	    $con+=$so[$k];
  	    $n25=$so[$k];
  	    $l25=$k;
  	    $k++;
      }
      $con+=$so[$k];
      $n50=$so[$k];
      $l50=$k;
      $k++;
    }	
    $con+=$so[$k];
    $n75=$so[$k];
    $l75=$k;
    $k++;
  }
  $con+=$so[$k];
  $n90=$so[$k];
  $l90=$k;
  $k++;
}

$l25++;$l50++;$l75++;$l90++;

print OUT "N25\t$n25\n";
print OUT "N50\t$n50\n";
print OUT "N75\t$n75\n";
print OUT "N90\t$n90\n";


print OUT "L25\t$l25\n";
print OUT "L50\t$l50\n";
print OUT "L75\t$l75\n";
print OUT "L90\t$l90\n";

$k=0;
$con=0;
my ($ng25,$ng50,$ng75,$ng90)=($so[0],$so[0],$so[0],$so[0]);
my ($lg25,$lg50,$lg75,$lg90)=(0,0,0,0);
while (($con<$est90) and ($k<=$#so))  {
  while ($con<$est75) {
    while ($con<$est50) {
  	  while ($con<$est25) {
        if ($sum<$est25) {
          ($lg25,$lg50,$lg75,$lg90)=(-1,-1,-1,-1);	
          ($ng25,$ng50,$ng75,$ng90)=(-1,-1,-1,-1);
          last;		
        }
  	    $con+=$so[$k];
  	    $ng25=$so[$k];
  	    $lg25=$k;
  	    $k++;
      }
      if ($sum<$est50) {
      	($lg50,$lg75,$lg90)=(-1,-1,-1);	
        ($ng50,$ng75,$ng90)=(-1,-1,-1);
        last;      	
      }
      $con+=$so[$k];
      $ng50=$so[$k];
      $lg50=$k;
      $k++
    }
    if ($sum<$est75) {
      ($lg75,$lg90)=(-1,-1);	
      ($ng75,$ng90)=(-1,-1);
      last;    	
    }
    $con+=$so[$k];
    $ng75=$so[$k];
    $lg75=$k;
    $k++;
  }
  if ($sum<$est90){
    $lg90=-1;	
   	$ng90=-1;
   	last;
  }
  $con+=$so[$k];
  $ng90=$so[$k];
  $lg90=$k;
  $k++;
}

$lg25++ if ($lg25>0);
$lg50++ if ($lg50>0);
$lg75++ if ($lg75>0);
$lg90++ if ($lg90>0);

print OUT "estSize\t$est\n";

print OUT "NG25\t$ng25\n";
print OUT "NG50\t$ng50\n";
print OUT "NG75\t$ng75\n";
print OUT "NG90\t$ng90\n";

print OUT "LG25\t$lg25\n";
print OUT "LG50\t$lg50\n";
print OUT "LG75\t$lg75\n";
print OUT "LG90\t$lg90\n";



my ($max,$min)=($so[0],$so[$#so]);
print OUT "Min\t$min\n";
print OUT "Max\t$max\n";

my $small=10000;
$k=$#so;
my $nSma=0;
my $nSSize=0;
while ($so[$k]<$small) {
  $nSma++;
  $nSSize+=$so[$k];
  $k--;	
}

print OUT "TotalSmallContigs(<$small)\t$nSma\n";
print OUT "SmallContigLength(<$small)\t$nSSize\n";




