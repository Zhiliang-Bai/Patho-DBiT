use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_i $opt_o $opt_d);
getopts('i:o:d:');


my @locations;
open DD, "$opt_d" or die;
while (<DD>) {
	chomp;
	my ($x, $y) = (split(/\t/, $_))[1,2];
    push(@locations, "${x}_$y");
}
close DD;


open IN, "$opt_i" or die;
open OUT, '>', "$opt_o" or die;


my %matrix;
my %mirnames;

while (<IN>) {
	chomp;
	my ($x, $y, $mirname) = (split(/\t/, $_))[0,1,2];
	$matrix{"${x}_${y}_$mirname"} += 1;
	$mirnames{$mirname} += 1;
}
close IN;




my @mirnames = sort {$a cmp $b} keys(%mirnames);
print OUT "";
for my $var (@mirnames) {
	print OUT "\t$var";
}
print OUT "\n";
 

for my $var (@locations) {
	my $varout = $var;
	$varout =~ s/_/x/;
	print OUT "$varout";	
	for my $var1 (@mirnames) {
		my $count;
		if (exists($matrix{"${var}_$var1"})) {
			$count = $matrix{"${var}_$var1"};
		}else{
			$count = 0;
		}
		print OUT "\t$count";
	}
	print OUT "\n";	 
	
}

close OUT;