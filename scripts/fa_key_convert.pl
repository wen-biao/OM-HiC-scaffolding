#!/usr/bin/perl -w
use strict;
use warnings;
use File::Basename;
use Cwd qw(abs_path);

# This adds "${CURRENT_SCRIPT_PATH}/perl5/" direcory at run time to the @INC array
# This script sould sit one level above the additional Perl modules directory.
BEGIN {
	my $script_path = abs_path(dirname($0));
	my $module_path = abs_path($script_path . "/perl5");
	my $lib = glob("~/scripts/HybridScaffold/scripts/perl5");
	unshift @INC, $module_path;
	unshift @INC, $lib;
	my $lib2; my $lib3;
	if ($] >= 5.010000 && $] <= 5.011000) {
		$module_path = $module_path."/5.10.1"; 
		$lib = $lib."/5.10.1";
		$lib2 = $lib."/x86_64-linux-thread-multi";
		$lib3 = $module_path."/x86_64-linux-thread-multi";}
	elsif ($] >= 5.014000 && $] <= 5.015000) {
		$module_path = $module_path."/5.14.4"; 
		$lib = $lib."/5.14.4";
		$lib2 = $lib."/x86_64-linux-thread-multi";
		$lib3 = $module_path."/x86_64-linux-thread-multi";}
	else {
		print "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X\n";
		die "ERROR: Unsupported Perl version found: $]. Supported Perl versions include 5.10.X and 5.14.X\n"; 
		exit; }
	unshift @INC, $module_path;
	unshift @INC, $lib;
	unshift @INC, $lib2;
	unshift @INC, $lib3;
	#print "$0 library paths:\n\t"; print join("\n\t", @INC); print "\n";
}
use BNG::Utility;

my ($IN, $KEY, $OUT);
my ($filename, $line, $id);
my %hash;  # FastaHeader -> CMapID
my $count = 0;

open($KEY, $ARGV[1]) || dieLog ("ERROR: Can't open $ARGV[1]: $!\n");
while($line = <$KEY>){
	chomp $line;
	$line =~ s/\r//g;
	next if($line =~ /^#/ || $line =~ /^C/);
	
	my @x = split("\t", $line);
	dieLog ("ERROR: Key file format error!\n") unless (@x == 3);
	
	$hash{$x[1]} = $x[0];
}
close($KEY);

#$filename = $ARGV[0];
#$filename =~ s/(\S+)\.(\w+)$/$1_CmapIdHeaders.fa/;
$filename = $ARGV[1];
#$filename =~ s/(\S+)\.(\w+)$/$1_CmapIdHeaders.fa/;
$filename =~ s/(\S+)\.(\w+)$/$1/;
$filename =~ s/(\S+)\_(\w+)$/$1_CmapIdHeaders.fa/;

open($IN, $ARGV[0]) || dieLog ("ERROR: Can't open $ARGV[0]: $!\n");
dieLog ("ERROR: Output file: $filename alredy exists!\n") if(-e $filename);

open($OUT, ">".$filename) || dieLog ("ERROR: Can't open $filename: $!\n");
while($line = <$IN>){
	chomp $line;
	$line =~ s/\r//g;
	
	if($line =~ /^>/){
		if(!defined( $hash{substr($line, 1)} )){
			$count++;
			$id = substr($line, 1);
			print "WARNING: ", substr($line, 1), " has no corresponding CMAP ID in the KEY file.\n";
		}
		else{
			$id = $hash{substr($line, 1)};
		}
		print $OUT(">$id\n");
	}
	else{
		print $OUT("$line\n");
	}
}
close($IN);
close($OUT);

print "File conversion finished successfully with $count warning(s)!\n";

__END__

./fa_key_convert.pl test.fa test_key.txt

