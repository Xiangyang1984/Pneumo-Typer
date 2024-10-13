use strict;
use warnings;

#usage: perl $home_directory/script/mcr.pl $workplace/CPS_cluster_mapping.result $workplace/cps_cluster_workplace/blast_homologs_cluster/all_orthomcl.out

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
foreach my $gc (sort keys %hash){
    $count++;
    my %derepeat;
    foreach my $t (sort @{$taxa{$gc}}){
        $derepeat{$t}=$t;
    }
    print OUT "homologous_gene_cluster_", $count, "(", scalar @{$hash{$gc}}, " genes,", scalar keys %derepeat, " taxa): ", join (" ", @{$hash{$gc}}), "\n";
}
close OUT;
