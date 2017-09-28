#!/scicore/soft/apps/Perl/5.22.2-goolf-1.7.20/bin/perl
# qsub wrapper script

use strict;
use warnings;
use lib "/scicore/home/terracci/GROUP/usr_nobackup/local/perl5/lib/perl5/x86_64-linux-thread-multi/";

use Schedule::DRMAAc qw/ :all /;
use File::Temp();
use Cwd qw/ realpath /;
use Cwd;
use Getopt::Std;
use Scalar::Util qw(looks_like_number);
use POSIX qw(strftime);
use Sys::Hostname;

my %opt;
getopts('hco:', \%opt);

my $usage = <<ENDL;
Usage: perl qsub.pl -h -- [qsub args]
    -o [file]: check file for non-zero size
    -c: check file across nodes for same file size
ENDL

my $host = hostname;

sub check_file {
    my ($cwd, $file) = @_;
#   my @nodes = qw/lbi01 lbi02 lbi03 lbi04 lbi05 lbi06 lbi07 lbi08 lbi09 lhi01 lhi02 lhi03 lhi04 lhi05 lhi06 lhi07 lhi08 lhi09 lhi10 lhi11 lhi12 lhi13 lhi14 lii01 lii02 lii03 lii04 lii05 lii06 lii07 lii08 lii09 lii10 lii11 lii12 lii13 lii14 lii15 lii16 lii17 lii18 lii19 lii20 lii21 lii22 lii23 lii24 lii25 lii26 lii27 lii28 lmi01 lmi02 lmi03 lmi04 lmi05 lmi06 lmi07 lmi08 login10 login12 login13 login14 login15 login17 login18 sgi01 smi01 usi101 usi102 usi103 usi104 usi105 usi106 usi107 usi108 usi109 usi110 usi111 usi112 usi113 usi114 usi115 usi116 usi117 usi118 usi119 usi120 usi121 usi122 usi123 usi124 usi125 usi126 usi127 usi128 usi129 usi130 usi131 usi132 usi133 usi134 usi135 usi136 usi137 usi138 usi139 usi140/;
    my @nodes = qw/coi01 lbi01 lbi02 lbi03 lbi04 lbi05 lbi06 lbi07 lbi08 lbi09 lhi01 lhi02 lhi03 lhi04 lhi05 lhi06 lhi07 lhi08 lhi09 lhi10 lhi11 lhi12 lhi13 lhi14 lii01 lii02 lii03 lii04 lii05 lii06 lii07 lii08 lii09 lii10 lii11 lii12 lii13 lii14 lii15 lii16 lii17 lii18 lii19 lii20 lii21 lii22 lii23 lii24 lii25 lii26 lii27 lii28 lmi01 lmi02 lmi03 lmi04 lmi05 lmi06 lmi07 lmi08 login12 login13 login14 login15 login17 login18 sbe01 sbe02 sbe03 sbe04 sbe05 sbe06 sbe07 sbe08 sbe09 sbe10 sbe11 sgi01 sgi21 sgi22 sgi23 sgi24 sgi25 sgi26 smi01 uni01 uni02 uni03 uni04 uni05 uni06 uni07 uni08 uni09 uni10 uni11 uni12 uni13 uni14 uni15 uni16 uni17 uni18 uni19 uni20 uni21 uni22 uni23 uni24 uni25 uni26 uni27 uni28 uni29 uni30 uni31 uni32 uni33 uni34 uni35 uni36 uni37 uni38 uni39 uni40 usi01 usi02 usi04 usi05 usi06 usi07 usi08 usi09 usi10 usi101 usi102 usi103 usi104 usi105 usi106 usi107 usi109 usi11 usi110 usi111 usi112 usi113 usi114 usi115 usi116 usi117 usi118 usi119 usi12 usi120 usi121 usi122 usi123 usi124 usi125 usi126 usi127 usi128 usi129 usi13 usi130 usi131 usi132 usi133 usi134 usi135 usi136 usi137 usi138 usi139 usi14 usi140 usi15 usi16 usi17 usi18 usi19 usi20 usi21 usi22 usi23 usi24 usi25 usi26 usi27 usi28 usi29 usi30 usi31 usi32 usi33 usi34 usi35 usi36 usi37 usi38 usi39 usi40 usi41 usi42 usi43 usi44 usi45 usi46 usi47 usi48 usi49 usi50 usi51 usi52 usi53 usi54 usi55 usi56 usi57 usi58 usi59 usi60 usi91/;
    my $maxConnectionFails = 1;
    my $fileSize = `stat -c\%s $cwd/$file`;
    my $nowString = strftime "%a %b %e %H:%M:%S %Y", localtime;
    print "[$host $nowString] failed local file size check\n" and return 0 if $?;
    chomp $fileSize;
    #print "checking $cwd/$file on nodes ($fileSize)\n";
    my $connectionFails = 0;
    for my $node (@nodes) {
        my $nodeFileSize = `ssh -x $node stat -c\%s $cwd/$file 2>&1`;
        chomp $nodeFileSize;
        if ($? || $nodeFileSize =~ /^stat/ || $nodeFileSize =~ /^ssh/ || !looks_like_number($nodeFileSize)) {
            print "[$node $nowString] failed remote file size check: $nodeFileSize\n";
            $connectionFails++;
        } elsif ($fileSize != $nodeFileSize) {
            print "[$node $nowString] file size does not match: $fileSize != $nodeFileSize\n";
            return 0;
        }
    }
    #print "$cwd/$file: all file sizes match\n";
    if ($connectionFails > $maxConnectionFails) {
        print "[$host $nowString] too many connection failures\n";
        return 0;
    } else {
        return 1;
    }
}

sub HELP_MESSAGE {
   print STDERR $usage;
   exit(1);
}

HELP_MESSAGE if $opt{h};

my $scriptFile = File::Temp->new(TEMPLATE => 'tempXXXXX', DIR => '/scicore/home/terracci/GROUP/tmp_nobackup', SUFFIX => '.sge');

my $args = join " ", @ARGV;
while (<STDIN>) {
    print $scriptFile $_;
}
close $scriptFile;

my ($error, $diagnosis) = drmaa_init(undef);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, my $jt, $diagnosis) = drmaa_allocate_job_template();
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_REMOTE_COMMAND, $scriptFile->filename);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_NATIVE_SPECIFICATION, $args);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_set_attribute($jt, $DRMAA_WD, getcwd());
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, my $jobid, $diagnosis) = drmaa_run_job($jt);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_delete_job_template($jt);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

sub signalHandler {
    my ($error, $diagnosis) = drmaa_control($jobid, $DRMAA_CONTROL_TERMINATE);
    die drmaa_strerror($error) . "\n" . $diagnosis if $error;
    die "Received interrupt: terminating job\n";
}

$SIG{INT} = \&signalHandler;
$SIG{TERM} = \&signalHandler;

# loop to give a chance to receive sigint/sigterms
my $stat;
do {
    ($error, my $jobidOut, $stat, $diagnosis) = drmaa_wait($jobid, 10);
} until ($error != $DRMAA_ERRNO_EXIT_TIMEOUT);

# pull all exit-related codes
($error, my $exitStatus, $diagnosis) = drmaa_wexitstatus($stat);
die drmaa_strerror($error) . "\n" . $diagnosis if $error;
($error, my $aborted, $diagnosis) = drmaa_wifaborted( $stat );
die drmaa_strerror($error) . "\n" . $diagnosis if $error;
($error, my $signaled, $diagnosis ) = drmaa_wifsignaled( $stat );
die drmaa_strerror($error) . "\n" . $diagnosis if $error;
($error, my $coreDumped, $diagnosis ) = drmaa_wcoredump( $stat );
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

($error, $diagnosis) = drmaa_exit();
die drmaa_strerror($error) . "\n" . $diagnosis if $error;

my $exitCodes = $exitStatus + $aborted + $signaled + $coreDumped;
my $fileStatus = 0;
if ($exitCodes == 0) {
    sleep 60; # wait for file sync
#    if ($opt{o} && $opt{o} ne "NULL") {
#        my $i = 0;
#        while (!check_file(getcwd(), $opt{o})) {
#            if ($i++ > 50) {
#                if (!-e $opt{o}) {
#                    $fileStatus = 66;
#                    my $nowString = strftime "%a %b %e %H:%M:%S %Y", localtime;
#                    print STDERR "[$host $nowString] ERROR $opt{o}: file does not exist\n";
#                } else {
#                    $fileStatus = 77;
#                    my $nowString = strftime "%a %b %e %H:%M:%S %Y", localtime;
#                    print STDERR "[$host $nowString] ERROR $opt{o}: file sizes do not match across nodes\n";
#                    system("rm $opt{o}");
#                }
#                last;
#            }
#            sleep 20;
#        }
#    }
}
#$exitCodes += $fileStatus;

# check for zero-size output file and remove it
if ($exitCodes == 0) {
    my $i = 0;
    while ($opt{c} && $opt{o} && $opt{o} ne "NULL" && !-s $opt{o}) {
        if ($i++ > 30) { 
            $fileStatus = 99;
            print STDERR "ERROR $opt{o}: file is size 0\n";
            system("rm $opt{o}");
            last;
        }
        sleep 20;
    }

}
$exitCodes += $fileStatus;

exit $exitCodes;
