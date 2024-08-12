use strict;
use warnings;
use Getopt::Std;
#use PerlIO::gzip;
use vars qw($opt_i $opt_o $opt_l);
getopts('i:o:l:');
use File::Path qw(make_path);


if ($opt_i =~ m/\.bam$/) {
    open IN, "samtools view -h $opt_i |" or die "Could not open BAM: $!";
} else {
    open IN, "$opt_i" or die;
}


$opt_o =~ s/\/$//;
$opt_o = "${opt_o}/";

unless (-d $opt_o) {
    make_path($opt_o) or die "Failed to create path: $opt_o\n";
}


my %sams;
my $annolines = '';
while (<IN>) {
    if ($_ =~ m/^@/) {
        $annolines .= $_;
        next;  
    }
    chomp;
    my ($readname,$flag,$mirname) = (split(/\t/, $_))[0,1,2];
    if ($flag & 4) {
        next;
    }
    my $readname2 = (split(/\|:_:\|/,$readname))[1];
    my ($xx, $yy) = (split(/_/,$readname2))[0,1];
    $sams{"${xx}x$yy"} .= "$_\n";
}

my @allpos = keys(%sams);
for my $var (@allpos) {
    open(OUT, "| samtools view -Sb - > $opt_o${var}.bam") or die "Cannot open samtools pipe: $!";
    print OUT "$annolines";
    print OUT "$sams{$var}";
    close OUT;
    system "samtools index $opt_o${var}.bam";
}


close IN;

