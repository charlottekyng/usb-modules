#!/usr/bin/env perl

####### FIX VARSCAN VCF #######
# change FREQ header and from % to numeric
# change header from NORMAL/TUMOR to the actual sample names
# change AD implementation to match that of GATK
# change empty ALT to <NON_REF> to allow GATK processing
########

use strict;
use warnings;

use Getopt::Std;
my %opt;
getopts('hf:n:t:', \%opt);

my $usage = <<ENDL;
perl varscanToVcf.pl -f [ref fasta] -t [tumor sample name] -n [normal sample name]
ENDL

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

my $now = localtime;
my $normal = $opt{n};
my $tumor = $opt{t};


while (my $l = <>) {
	chomp $l;
	if ($l =~ /\#\#/) {
		if ($l =~ /\#\#FORMAT\=\<ID\=FREQ,/) {
			$l = "\#\#FORMAT=\<ID=FA,Number=A,Type=Float,Description=\"Allele fraction of the alternate allele with regard to reference\"\>";
		} elsif ($l =~ /\#\#FORMAT\=\<ID\=AD,/) {
			$l = "\#\#FORMAT\=\<ID\=AD,Number\=A,Type\=Float,Description\=\"Allelic depths for the ref and alt alleles in the order listed\"\>";
		}
		print $l."\n"; 
	} elsif ( $l =~ /\#CHROM/) { 
		$l =~ s/TUMOR/$tumor/;
		$l =~ s/NORMAL/$normal/;
		print $l."\n";
	} else {
		$l =~ s/str10/PASS/;
		my @line = split /\t/, $l;
		my $freq = my $ad = my $rd = "";
		if ($line[4] eq "") { $line[4] = "\<NON_REF\>";}

		my @format = split /:/, $line[8];
		for (my $i = 0; $i < scalar @format; $i++) {
			if ($format[$i] eq "FREQ") {
				$format[$i] = "FA";
				$freq = $i;
			} elsif ($format[$i] eq "AD") {
				$ad = $i;
			} elsif ($format[$i] eq "RD") {
				$rd = $i;
			}
		}
		$line[8] = join ':', @format;
		for (my $i = 9; $i <= (scalar @line - 1); $i++) {
			my @fields = split /:/, $line[$i];
			$fields[$freq] =~ s/%//;
			$fields[$freq] = $fields[$freq]/100;

			$fields[$ad] = $fields[$rd].",".$fields[$ad];

			$line[$i] = join ':', @fields;
		}
		print join "\t", @line; print "\n";
	}
}



