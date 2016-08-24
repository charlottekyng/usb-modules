use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;
use List::MoreUtils qw(uniq);

# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# x    define needed variables    x
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# command line arguments normal and tumor bam and output file
my $normalbam  = $ARGV[0];
my $tumorbam   = $ARGV[1];
my $outfile    = $ARGV[2];
my $reference  = $ARGV[3];
my $bedfile    = $ARGV[4];
my $min_ndepth = $ARGV[5];
my $max_ndepth = $ARGV[6];
my $snplocprefix = $ARGV[7];

# prefix for the pileup command
#my $cmdprefix = 'source /etc/profile.d/lmod.sh; module purge; source $HOME/.bashrc; module load SAMtools/1.3-goolf-1.7.20; ';
my $cmdprefix = 'module load SAMtools/1.3-goolf-1.7.20; ';
$cmdprefix .= "samtools mpileup -q15 -Q20 -B -f $reference -l $bedfile -r ";

# chromosomes
my @chroms = ();
open(BED, $bedfile) or die "Could not open $bedfile";
while (my $blah=<BED>) {
	my @chrom = split(/\t/, $blah);
	push @chroms, $chrom[0];
}
@chroms = uniq @chroms;
print STDERR join ',',@chroms;

# file prefix for snp location data
open(REF, $reference) or die "Could not open $reference";
my $line = <REF>;
if ($line =~ /^>chr/) {
	$cmdprefix .= "chr";
}
close REF;

# open the output file to print the pileup
open(my $ofh, "| /bin/gzip -c > $outfile") or die "Could not open file '$outfile' $!";

for my $chrom (@chroms) {
    # read in all snp locations
    my $snplocfile = $snplocprefix . "$chrom";
    # print "$snplocfile\n";
    open(my $ifh, '<', "$snplocfile") or die "Could not open snp locations file";
    my @allsnplocs = <$ifh>;
    close $ifh;
    # number of snps
    my $nsnps = $#allsnplocs + 1;
    # initialize variables
    my $currentsnpid = 0;
    my $currentloc = 0;
    my $currentsnp = $allsnplocs[$currentsnpid];
    my $acgtref =0;
    my @nacgt =();
    my @tacgt =();
    # create the samtools mpileup command for the data stream
    my $cmd = $cmdprefix . "$chrom $normalbam $tumorbam |";
    
    # open the pileup pipe
    open(PILEUP, "$cmd") || die "Open failed \n"; 
    # process until the pipe stops
    while (my $blah=<PILEUP>) {
	# split the input into fields (pileup uses tab not space)
	my @ipile = split(/\t/, $blah);
	
	# check normal read count is at least min_ndepth but no more than max_ndepth
	if ($ipile[3] >= $min_ndepth && $ipile[3] <= $max_ndepth) {
	    # current location
	    $currentloc = $ipile[1];
	    # check if current snp is not before current loc
	    while ($currentloc > $currentsnp) {
		$currentsnpid += 1;
		if ($currentsnpid < $nsnps) {
		    $currentsnp = $allsnplocs[$currentsnpid];
		} else {
		    # no more snps; set snp location to 2^28 (> chrom lengths)
		    $currentsnp = 268435456;
		}
	    }
	    if ($currentloc == $currentsnp){
		# print $ofh "$blah";
		# counts for the normal
		$acgtref = read2counts($ipile[2], $ipile[4]);
		@nacgt = @$acgtref;
		# counts for the tumor (provided tumor read count > 0)
		if ($ipile[6] > 0) {
		    $acgtref = read2counts($ipile[2], $ipile[7]);
		    @tacgt = @$acgtref;
		} else {
		    @tacgt = (0,0,0,0,0,0,0,0,0,0);
		}
		$ipile[0] =~ s/chr//;
		print $ofh "@ipile[0..1] $ipile[3] $nacgt[0] $ipile[6] $tacgt[0]\n";
	    }
	}
    }
}
# close the output file
close $ofh;

sub read2counts {
    # hash table to handle match to the ref allele
    my %base_hash = ('A',0, 'C', 2, 'G', 4, 'T', 6);
    # input arguments ref allele and read string
    my ($refallele, $readstring) = ($_[0], $_[1]);
    # counter for ACGT counts
    my @acgt = ();

    #remove mapping quality information
    $readstring =~ s/\^.//g;
    #remove insertions and deltions annotations
    while ($readstring =~ m/([\+|-])(\d+)/){
	my $indel_length = $2; 
	$readstring =~ s/([\+|-])(\d+)(\w{$indel_length})/@/;
    }
    #substitute tail information
    $readstring =~ s/\$//g;

    # start counting ref reads
    my $REF = ($readstring =~ tr/\.//);
    my $ref = ($readstring =~ tr/\,//);

    # add the two refs
    $REF += $ref;

    # put them in a vector
    @acgt = ($REF);

    # return result
    return(\@acgt);
}
