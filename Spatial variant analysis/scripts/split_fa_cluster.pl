use strict;
use warnings;
use Getopt::Std;
#use PerlIO::gzip;
use vars qw($opt_i $opt_o $opt_r $opt_n $opt_t);
getopts('i:o:r:n:t:');




print "Use : to divide normal clusters!\n";


open IN, "$opt_r" or die;
my %allcluster;
my $line1 = <IN>;
my %loc2clu;
while(<IN>) {
    chomp;
    s/"//g;
    my ($loc, $cluster) = (split(/,/, $_))[0,1];
    $loc2clu{$loc} = $cluster;
    $allcluster{$cluster} = 1;
}
close IN;
my @allcluster = keys(%allcluster);

my @normalcluster;
my @tumorcluster;

###
###
###
my %tumornormalhash;
if (defined $opt_n) {
    my @normalcluster = split(/[^0-9]+/,$opt_n);
    my @tumorcluster = split(/[^0-9]+/,$opt_t);
    %tumornormalhash = map { $_ => 'normal' } @normalcluster;
    @tumornormalhash{@tumorcluster} = ('tumor') x @tumorcluster;
    for my $loc (keys(%loc2clu)) {
        if (exists($tumornormalhash{$loc2clu{$loc}})) {
            $loc2clu{$loc} = $tumornormalhash{$loc2clu{$loc}};
        }else{
            delete $loc2clu{$loc};
        }
        
    }
    @allcluster = qw(normal tumor);
}

###
###
###


if ($opt_i =~ m/\.bam$/) {
    open IN, "samtools view -h $opt_i |" or die "Could not open BAM: $!";
} else {
    open IN, "$opt_i" or die;
}

my %file_handles;

for my $var (@allcluster) {
    open my $fh, '>', "${opt_o}${var}.sam" or die $!;
    $file_handles{$var} = $fh;
}


my %sams;
my $annolines = '';
my $annoout = 0;
while (<IN>) {
    if ($_ =~ m/^@/) {
        $annolines .= $_;
        next;  
    }elsif($annoout == 0) {
        $annoout = 1;
        for my $fh (values %file_handles) {
            print { $fh } "$annolines";
        }
    }
    chomp;
    #0_19_10:TCTATGCACT
    my ($readname,$flag,$mirname) = (split(/\t/, $_))[0,1,2];
    next if $flag == 4;
    my $readname2 = (split(/:/,$readname))[0];
    my ($xx, $yy)  = (split(/_/,$readname2))[1,2];
    my $position = "${xx}x${yy}";
    if (exists($loc2clu{$position})) {
        my $cluster = $loc2clu{$position};
        print { $file_handles{$cluster} } "$_\n";
    }
    
}

for my $fh (values %file_handles) {
    close $fh;
}
for my $var (@allcluster) {
    system "samtools view -bS ${opt_o}${var}.sam > ${opt_o}${var}.bam && rm ${opt_o}${var}.sam";
    system "samtools index ${opt_o}${var}.bam";
}



close IN;

