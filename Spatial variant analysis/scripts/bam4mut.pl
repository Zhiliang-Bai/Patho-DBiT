use strict;
use warnings;

use Getopt::Std;
use vars qw($opt_i $opt_o $opt_p);
getopts('i:o:p:');




open(my $bam, "-|", "samtools view -h $opt_i") or die "Cannot process input BAM file: $opt_i";
open(my $modified_bam, "|-", "samtools view -b -o $opt_o -") or die "Cannot write to output BAM file: $opt_o";

my %dedup;
while (my $line = <$bam>) {
    unless ($line =~ /^\@/) {
        my @fields = split(/\t/, $line);
        if ($fields[1] & 4) {
            next;
        }
        if (exists($dedup{$fields[0]})) {
            $dedup{$fields[0]} = 2;
        }else{
            $dedup{$fields[0]} = 1;
        }
    }
}
close $bam;

open(my $bam2, "-|", "samtools view -h $opt_i") or die "Cannot process input BAM file: $opt_i";
while (my $line = <$bam2>) {
    if ($line =~ /^\@/) {
        print $modified_bam "$line";
    } else {
        my @fields = split(/\t/, $line);
        if ($fields[1] & 4) {
            next;
        }
        if ($dedup{$fields[0]} == 2) {
            next;
        }
        print $modified_bam join("\t", @fields);
    }
}


close $bam2;
close $modified_bam;