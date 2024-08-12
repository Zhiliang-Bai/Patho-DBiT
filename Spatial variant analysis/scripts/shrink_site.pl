use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_i $opt_o);
getopts('i:o:');
#perl shrink_site.pl -o ~/gibbs/share/Zhiliang/LM0623-SNV/ChisqFiltDNAvar/SiteAnno_brief.txt -i result/VarMAT/Pixel/ChisqFiltDNAvar/SiteAnno.txt

open OUT, '>', "$opt_o" or die;
open IN, "$opt_i" or die;

my %gene2mut;
my %variant2gene;

while(<IN>) {
    next if $_ =~ m/^#/;
    chomp;
    my ($variant, $mutation, $gene) = (split(/\t/, $_))[0,6,13];
    if ($gene =~ m/SYMBOL=([^=;]+);/) {
        $gene = $1;
    }else{
        next;
    }
    if (exists($variant2gene{$variant})) {
        $variant2gene{$variant} = "$variant2gene{$variant}".",$gene";
    }else{
        $variant2gene{$variant} = $gene;
    }
    my $vg = "$variant---$gene";
    if (exists($gene2mut{$vg})) {
        $gene2mut{$vg} = "$gene2mut{$vg}".",$mutation";
    }else{
        $gene2mut{$vg} = $mutation;
    }
}


my @variant2gene = sort {$a cmp $b} keys(%variant2gene);

for my $variant (@variant2gene) {
    my $genes = $variant2gene{$variant};
    my @genes = sort { $a cmp $b } keys %{ { map { $_ => 1 } split(/,/, $genes) } };
    for my $genei (@genes) {
        my $vg = "$variant---$genei";
        my @muts = sort { $a cmp $b } keys %{ { map { $_ => 1 } split(/,/, $gene2mut{$vg}) } };
        my $muts = join(",",@muts);
        print OUT "$variant\t$genei\t$muts\n";
    }

}

close IN;
close OUT;