#perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/script/serotype_stastatic.pl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/Pneumo-Typer_workplace/Serotype.out /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/genomes_metadata.txt /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/Pneumo-Typer_workplace/Serotype.out_statistics.out

use strict;
use warnings;

my $meta = $ARGV[1];
open(META, $meta);
my %meta;
while(<META>){
    chomp;
    my @arr = split /\t/, $_;
    my $serotype = $arr[4];
    next unless $serotype !~ /NA|blank|not|NT|Not/;
    $serotype =~ s/^0//g if $serotype =~ /^0/;
    $meta{$arr[0]} = $serotype;
}
close META;

open (OUT, ">$ARGV[2]");
my $correct_sero = $ARGV[0];
open(MATRIX, $correct_sero);
while(<MATRIX>){
    chomp;
    next unless $_ !~ /^Strain	Correct serotype/;
    my @array = split /\t/, $_;
    next unless $meta{$array[0]};
    print  OUT "$array[0]\tpredict: $array[1]\tCorrect: \tTure: $meta{$array[0]}\n" if $array[1] eq $meta{$array[0]};
    print  OUT "$array[0]\tpredict: $array[1]\tWrong: \tTure: $meta{$array[0]}\n" if $array[1] ne $meta{$array[0]};
}
close MATRIX;
close OUT;
