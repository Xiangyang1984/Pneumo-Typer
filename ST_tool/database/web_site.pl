# perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/database/web_site.pl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/database/core_genome.list


use strict;
use warnings;

open (WEB, $ARGV[0]);
while(<WEB>){
    chomp;

    system("wget -c --no-check-certificate $_");
    my $name = $_;

    #'https://pubmlst.org/bigsdb?db=pubmlst_spneumoniae_seqdef&page=downloadAlleles&locus=SPNE00001'

    #$name =~ s/\'//g;
    $name =~ s/https:\/\/pubmlst.org\///;
    $name =~ /.*\=(.*)\'/;
    #print "$name\t$1\n";
    system("mv $name $1.fas");
}


