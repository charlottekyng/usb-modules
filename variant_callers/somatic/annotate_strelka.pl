use strict;

while (my $line = <>) {
	chomp $line;
	if ($line =~ /^\#/) { 
		print $line; print "\n"; 
		if ($line =~ /^\#\#FORMAT\=\<ID\=TOR/) { 
			print "\#\#FORMAT\=\<ID\=AD,Number=2,Type=Float,Description=\"AD computed from tier 1\"\>\n";
			print "\#\#FORMAT\=\<ID\=AF,Number=2,Type=Float,Description=\"AF computed from tier 1\"\>\n";
		}
	} else {
		my @arr = split /\t/, $line;
		my @format = split /:/, $arr[8];
		my $n_dp = my $n_tir = 0;
		for (my $i=0; $i < scalar @format; $i++) {
			if ($format[$i] eq "DP") {
				$n_dp = $i; #last;
			} elsif ($format[$i] eq "TIR") {
				$n_tir = $i; #last;
			}
		}

		my @normal = split /:/, $arr[9];
		my @tumour = split /:/, $arr[10];

		$normal[$n_tir] =~ /^(\d+),/;
		my $normal_ad = $1.",".($normal[$n_dp]-$1);
		my $normal_af = $1/$normal[$n_dp];

		$tumour[$n_tir] =~ /^(\d+),/;
		my $tumour_ad = $1.",".($tumour[$n_dp]-$1);
		my $tumour_af = $1/$tumour[$n_dp];

		push @format, "AD", "AF";
		push @normal, $normal_ad, $normal_af;
		push @tumour, $tumour_ad, $tumour_af;

		$arr[8] = join ":", @format;
		$arr[9] = join ":", @normal;
		$arr[10] = join ":", @tumour;

		print join "\t", @arr; print "\n";
	}
}
		
	
