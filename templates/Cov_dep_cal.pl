#!/usr/local/bin/perl -w
use strict;

my \$dep = '$depth_file';
my \$fas = '$amr_reference';
my \$out = 'out.tab';

#read in the depth file and store as hash with reference number as the first key layer, coordinate as the second key layer and coverage as the value
my %depth = ();
open (my \$DEP, \$dep) or die "can't open \$dep\n";
while (my \$dep_line = <\$DEP>) {
    chomp \$dep_line;
    my @dep_array = split /\\s+/, \$dep_line;
    \$depth{\$dep_array[0]} -> {\$dep_array[1]} = \$dep_array[2];
}
close \$DEP;

#read in the reference fasta file and store as hash with reference number as the key and sequence length as the value, reference number should match with that in the depth file
my %fasta = ();
open ( my \$FAS, \$fas ) or die "can't open \$fas\n";
my \$def;
while ( my \$f_line = <\$FAS> ){
    chomp \$f_line;
    if (\$f_line =~ m/\\|(NG_\\d+)\\|(\\d)/) {
	\$def = "\$1"."_"."\$2";
	\$fasta{\$def} = 1;
    }
    else {
	\$fasta{\$def} += length(\$f_line);
    }
}
close \$FAS;

#loop through the depth hash to tally number of base pairs covered and total read coverage in each reference, and calculate the percent coverage and average coverage over covered portion for each reference
open (my \$OUT, ">\$out") or die "can't open \$out\n";
print \$OUT "ref\tpercent_cov\taverage_cov\n";
foreach my \$key (sort keys %depth) {
    my \$line_counter;
    my \$ref_len;
    my \$cov;
    my \$ref_cov;
    my \$ave_cov;
    if (defined \$fasta{\$key}) {
        \$ref_len = \$fasta{\$key};
    }
    foreach my \$coord (sort {\$a<=>\$b} keys %{\$depth{\$key}}) {
	\$line_counter ++;
	\$cov += \$depth{\$key}{\$coord};
    }
    \$ref_cov = \$line_counter/\$ref_len*100;
    \$ref_cov = sprintf("%.2f",\$ref_cov);
    \$ave_cov = \$cov/\$line_counter;
    \$ave_cov = sprintf("%.2f",\$ave_cov);
    print \$OUT "\$key\t\$ref_cov\t\$ave_cov\n";
}
close \$OUT;
exit;
