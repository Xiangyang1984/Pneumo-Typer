#perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/mcr.pl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/Pneumo-Typer_workplace/CPS_cluster_mapping.result /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/Pneumo-Typer_workplace/Gcluster_workplace/blast_homologs_cluster/all_orthomcl.out

use strict;
use warnings;

open (CPS, $ARGV[0]);

my %hash;
my %taxa;
while (<CPS>){
    chomp;
    my @array = split /\t/, $_;
    my $gene = $array[0];
    $gene =~ s/#.*?;(.*)//g;
    my $class = $array[4];
    $class =~ s/\-.*//g;
    push (@{$hash{$class}}, "$class;$gene");
    push (@{$taxa{$class}}, $1);
   
}
close CPS;

open (OUT, ">$ARGV[1]");
my $count=-1;
foreach my $gc (keys %hash){
    $count++;
    my %derepeat;
    foreach my $t (@{$taxa{$gc}}){
        $derepeat{$t}=$t;
    }
    print OUT "homologous_gene_cluster_", $count, "(", scalar @{$hash{$gc}}, " genes,", scalar keys %derepeat, " taxa): ", join (" ", @{$hash{$gc}}), "\n";
}
close OUT;
