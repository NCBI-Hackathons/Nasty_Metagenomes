#!usr/bin/perl -w
use strict;

my $blast = shift @ARGV or die;
my $out = shift @ARGV or die;

my %blast = ();
open (my $BLAST, $blast) or die "can't open $blast\n";
my $ref;
my $query;
my $aln_len;
my $query_end;
my $pid;
my $bitscore;
while ( my $b_line = <$BLAST> ){
    chomp $b_line;
    my @blast_array = split /\t/, $b_line;
    $ref = $blast_array[1];
    $query = $blast_array[0];
    $pid = $blast_array[2];
    $aln_len = $blast_array[3];
    $bitscore = $blast_array[-1];
    if ($blast_array[6] > $blast_array[7]) {
	$query_end = $blast_array[6];
    }
    else {
	$query_end = $blast_array[7];
    }
    if ($pid >= 85 && $bitscore >= 90) {
	$blast{$query} -> {$query_end} -> {$bitscore} -> {$ref} = "$aln_len\t$pid";
    }
}
my $key_number = keys %blast;
print "number of entries in blast hash: $key_number\n";
close $BLAST;
#die;

open (my $OUT, ">$out") or die "can't open $out\n";
foreach my $key (sort keys %blast) {
    my $end = 1;
    foreach my $endkey (sort {$b<=>$a} keys %{$blast{$key}}) {
	if (abs($endkey - $end) >= 200) {
	    $end = $endkey;
	    my $counter1;
	    foreach my $scorekey (sort {$b<=>$a} keys %{$blast{$key}{$endkey}}) {
		$counter1++;
		if ($counter1 == 1) {
		    my $counter2;
		    foreach my $refkey (sort keys %{$blast{$key}{$endkey}{$scorekey}}) {
			$counter2++;
			if ($counter2 == 1) {
			    print $OUT "$key\t$refkey\t$blast{$key}{$endkey}{$scorekey}{$refkey}\n";
			}
		    }
		}
	    }
	}
    }
}

close $OUT;
exit;
