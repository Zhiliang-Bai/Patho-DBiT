use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_o $opt_r $opt_i);
getopts('r:o:i:');
#perl compareDNA_RNA.pl -r cache/DNA_germline.pos.txt -o 123 -i cache/RNA_germline.pos.txt

open IN, "$opt_i" or die;
open REF, "$opt_r" or die;

my %ref;
while (<REF>) {
        chomp;
        my ($chr, $location, $old, $type) = split(/:/, $_);
        if ($type =~ m/,/) {
            my @types = split(/,/, $type);
            for my $var (@types) {
                $ref{"$chr:$location:$old:$var"} = 1;
            }
        }else{
            $ref{"$chr:$location:$old:$type"} = 1;
        }
}
close REF;

my $count = keys %ref;
print "$count\n";

open OUT, '>', "$opt_o" or die;
open OUTT, '>', "$opt_o.unmap" or die;
while (<IN>) {
        chomp;
        my ($chr, $location, $old, $type) = split(/:/, $_);
        if ($type =~ m/,/) {
            my @types = split(/,/, $type);
            for my $var (@types) {
                my $entry = "$chr:$location:$old:$var";
                if (exists($ref{$entry})) {
                    print OUT "$entry\n";
                }else{
                    print OUTT "$entry\n";
                }
            }
        }else{
            my $entry = "$chr:$location:$old:$type";
            if (exists($ref{$entry})) {
                print OUT "$entry\n";
            }else{
                print OUTT "$entry\n";
            }
        }
}
close IN;
close OUT;
close OUTT;



