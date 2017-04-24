#!/usr/bin/env perl

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
		$l =~ s/FREQ/FA/; 
		if ($l =~ /\#\#FORMAT\=\<ID\=AD,/) {
			$l = "\#\#FORMAT\=\<ID\=AD,Number\=R,Type\=Integer,Description\=\"Allelic depths for the ref and alt alleles in the order listed\"\>";
		}
		print $l."\n"; 
	} elsif ( $l =~ /\#CHROM/) { 
		$l =~ s/TUMOR/$tumor/;
		$l =~ s/NORMAL/$normal/;
		print $l."\n";
	} else {
		$l =~ s/str10/PASS/;
		my @line = split /\t/, $l;
		my $freq = my $ad = my $dp4 = "";
		my @format = split /:/, $line[8];
		for (my $i = 0; $i < scalar @format; $i++) {
			if ($format[$i] eq "FREQ") {
				$format[$i] = "FA";
				$freq = $i;
			} elsif ($format[$i] eq "AD") {
				$ad = $i;
			} elsif ($format[$i] eq "RD") {
				$dp4 = $i;
			}
		}
		$line[8] = join ':', @format;
		for (my $i = 9; $i <= (scalar @line - 1); $i++) {
			my @fields = split /:/, $line[$i];
			$fields[$freq] =~ s/%//;
			$fields[$freq] = $fields[$freq]/100;

			my @dp4 = split /,/, $fields[$dp4];
			$fields[$ad] = ($dp4[0]+$dp4[1]).",".$fields[$ad];

			$line[$i] = join ':', @fields;
		}
		print join "\t", @line; print "\n";
	}
}



