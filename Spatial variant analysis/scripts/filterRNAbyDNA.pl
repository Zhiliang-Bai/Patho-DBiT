use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_o $opt_r $opt_i);
getopts('r:o:i:');
#perl filterRNAbyDNA.pl -r cache/DNA_germline.pos.txt -o 111 -i result/Strelka/RNA/germline/results/variants/variants.vcf.gz
#

open IN, "gzip -dc $opt_i|" or die;
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
    next if $_ =~ m/^#/;
    chomp;
    my @inputarray = split(/\t/, $_);
    my ($chr, $location, $old, $type) = @inputarray[0,1,3,4];
    next unless $inputarray[6] =~ m/PASS/;
    if ($type =~ m/,/) {
        my @types = split(/,/, $type);
        for my $var (@types) {
            my $entry = "$chr:$location:$old:$var";
            $inputarray[4] = $var;
            my $outputline = join("\t",@inputarray);
            if (exists($ref{$entry})) {
                print OUT "$outputline\n";
            }else{
                print OUTT "$outputline\n";
            }
        }
    }else{
        my $entry = "$chr:$location:$old:$type";
        $inputarray[4] = $type;
        my $outputline = join("\t",@inputarray);
        if (exists($ref{$entry})) {
            print OUT "$outputline\n";
        }else{
            print OUTT "$outputline\n";
        }
    }
}
close IN;
close OUT;
close OUTT;



