#!/usr/bin/env perl

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
		my $fao = my $fdp = ""; #these are the flow counts
		my @format = split /:/, $line[8];
		for (my $i = 0; $i < scalar @format; $i++) {
			if ($format[$i] eq "FAO") { $fao = $i;
			} elsif ($format[$i] eq "FDP") { $fdp = $i;
			}
		}
		if ($fao ne "" && $fdp ne "") {
			push @format, "FA";
			$line[8] = join ':', @format;
			for (my $i = 9; $i <= (scalar @line) -1; $i++) {
				my @fields = split /:/, $line[$i];
				my @faos = split /,/, $fields[$fao];
				my $fdps = $fields[$fdp];
				my @fas = ();
				for (my $j = 0; $j < scalar @faos; $j++) {
					if ($fdps > 0) {
						push @fas, $faos[$j]/$fdps;
					} else { push @fas, "."; }
				}
				push @fields, (join ',',@fas);
				$line[$i] = join ':', @fields;
			}
		}
		print join "\t", @line; print "\n";
	}
}




