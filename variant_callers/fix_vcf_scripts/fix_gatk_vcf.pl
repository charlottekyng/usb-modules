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
	if ($l =~ /\#\#FORMAT=\<ID=AD,/) {
		print $l."\n";
		print "\#\#FORMAT=\<ID=FA,Number=A,Type=Float,Description=\"Allele fraction of the alternate allele with regard to reference\"\>\n"; 
	} elsif ($l =~ /\#\#FORMAT=\<ID=AF,/) {
		$l =~ s/AF,/FA,/;
		print $l."\n";
	} elsif ($l =~ /\#/) {
		print $l."\n"; 
	} else {
		my @line = split /\t/, $l;
		my $ad = my $dp = my $af = "";
		my @format = split /:/, $line[8];
		for (my $i = 0; $i < scalar @format; $i++) {
			if ($format[$i] eq "AD") {
				$format[$i] = "AD:FA";
				$ad = $i;
			} elsif ($format[$i] eq "DP") {
				$dp = $i;
#			} elsif ($format[$i] eq "AF") {
#				$format[$i] = "FA";
#				$af = $i;
			}
		}
		$line[8] = join ':', @format;
		if ($ad ne "") {
		for (my $i = 9; $i < scalar @line; $i++) {
			next if ($line[$i] eq "./.");
			my @fields = split /:/, $line[$i];
			my @ads = split /,/, $fields[$ad];
			my $dps = $fields[$dp];
			my @fas = ();
			for (my $j = 1; $j < scalar @ads; $j++) {
				my $val = ".";
				if ($dps > 0) { 
					$val = $ads[$j]/$dps;
				}
				push @fas, $val;
			}
			$fields[$ad] = $fields[$ad].":".join ',',@fas;
			$line[$i] = join ':', @fields;
		}}
		print join "\t", @line; print "\n";
	}
}




