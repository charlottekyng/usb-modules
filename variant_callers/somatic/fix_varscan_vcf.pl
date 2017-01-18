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
	if ($l =~ /\#\#/) { print $l."\n"; 
	} elsif ( $l =~ /\#CHROM/) { 
		$l =~ s/TUMOR/$tumor/;
		$l =~ s/NORMAL/$normal/;
		print $l."\n";
	} else {
		$l =~ s/str10/PASS/;
		my @line = split /\t/, $l;
		my $freq = "";
		my @format = split /:/, $line[8];
		for (my $i = 0; $i < scalar @format; $i++) {
			if ($format[$i] eq "FREQ") {
				$format[$i] = "AF";
				$freq = $i;
			}
		}
		$line[8] = join ':', @format;
		for (my $i = 9; $i <= 10; $i++) {
			my @fields = split /:/, $line[$i];
			$fields[$freq] =~ s/%//;
			$fields[$freq] = $fields[$freq]/100;
			$line[$i] = join ':', @fields;
		}
		print join "\t", @line; print "\n";
	}
}



