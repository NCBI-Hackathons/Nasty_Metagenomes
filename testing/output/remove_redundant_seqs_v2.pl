#!usr/bin/perl -w
use strict;

my $header = shift @ARGV or die;
my $fasta = shift @ARGV or die;
my $out = shift @ARGV or die;

my %hash = ();
open (my $H, $header) or die "can't open $header\n";
while (my $line = <$H>) {
    chomp $line;
    if (defined $hash{$line}) {
	$hash{$line} ++;
    }
    else {
	$hash{$line} = 1;
    }
}
close $H;

open (my $FAS, $fasta) or die "can't open $fasta\n";
open (my $OUT, ">$out") or die "can't open $out\n";
my $line_counter = 0;
my $def;
my $finder;
while (my $f_line = <$FAS>) {
    chomp $f_line;
    $line_counter ++;
	
	if ($f_line =~ m/^>(\S+)/){
		$finder = $f_line;
	} 
	
    if ($f_line =~ m/^>(\S+)/ && $hash{$finder} == 1) {
	print $OUT "$f_line\n";
	next;
    }

    if ($f_line =~ m/^>(\S+)/ && $hash{$finder} > 1 && $hash{$finder} < 10000) {
	print $OUT "$f_line\n";
	$hash{$finder} = 10000;
	next;
    }

    if ($f_line =~ m/^>(\S+)/ && $hash{$finder} == 10000) {
	next;
    }

    if ($line_counter % 2 == 0){
	
		if ($hash{$finder} == 10000 | $hash{$finder} == 1){
			print $OUT "$f_line\n";
			$hash{$finder} = 20000
		} else{
			next;
		}
	}
}
close $OUT;
exit;