#!/usr/bin/env perl

######### FIX TVC VCF #######
# TVC outputs "AF", which is VAF calculated from flow-based counts
# this script calculates "FA", which is VAF calculated from non-flow-based counts
# The flow-based counts are supposed to be more accurate but there are too many false positives, so do both!
################


use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('hf:n:t:', \%opt);

my $usage = <<ENDL;
perl fix_gatk_vcf.pl <STDIN>
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

my $now = localtime;

while (my $l = <>) {
	chomp $l;
	if ($l =~ /\#\#FORMAT=\<ID=AF,/) {
		print $l."\n";
		print "\#\#FORMAT=\<ID=FA,Number=R,Type=Float,Description=\"Variant allele fraction based on raw counts\"\>\n";
	} elsif ($l =~ /\#/) {
		print $l."\n"; 
	} else {
		my @line = split /\t/, $l;
		my $ao = my $dp = ""; #these are the non-flow counts
		my @format = split /:/, $line[8];
		for (my $i = 0; $i < scalar @format; $i++) {
			if ($format[$i] eq "AO") { $ao = $i;
			} elsif ($format[$i] eq "DP") { $dp = $i;
			}
		}
		if ($ao ne "" && $dp ne "") {
			push @format, "FA";
			$line[8] = join ':', @format;
			for (my $i = 9; $i <= (scalar @line) -1; $i++) {
				my @fields = split /:/, $line[$i];
				my @aos = split /,/, $fields[$ao];
				my $dps = $fields[$dp];
				my @fas = ();
				for (my $j = 0; $j < scalar @aos; $j++) {
					if ($dps > 0) {
						push @fas, $aos[$j]/$dps;
					} else { push @fas, "."; }
				}
				push @fields, (join ',',@fas);
				$line[$i] = join ':', @fields;
			}
		}
		print join "\t", @line; print "\n";
	}
}




