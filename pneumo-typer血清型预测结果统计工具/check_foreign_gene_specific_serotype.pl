#perl t.pl tbl_part Serotype_sta.out gene.txt

use strict;
use warnings;

opendir (F, $ARGV[0]);
my @f = readdir F;
@f = grep ($_ !~ /^\./, @f);
closedir F;

my %g;
open (G, $ARGV[2]);
while(<G>){
    chomp;
    my @a = split /\t/, $_;
    push @{$g{$a[0]}}, $a[1];    # serotype vs gene
}
close G;

my %s;
my %w;
open (S, $ARGV[1]);
while(<S>){
    chomp;
    my @s = split /\t/, $_;
    my $serotype = $s[2];
    $serotype =~ s/^0//;
    $s{$s[0]}=$serotype if $serotype ne "15B|15C";    #genome vs serotype
    $w{$s[0]} = $s[3];  #genome vs right or wrong
}
close S;


foreach my $f (@f){
    open (T, "$ARGV[0]/$f");
    while(<T>){
        chomp;
        my @g = split /\t/, $_;
        if (defined $s{$f}){
            if (($w{$f} eq "Right") && ($g[1] ne "intron")){
                print "$f\t$s{$f}\t$w{$f}\t$g[1]\n" if !grep($g[1] eq $_, @{$g{$s{$f}}});
            }
        }
    }
    close T;
}
