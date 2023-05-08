#perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/interested_gene.pl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/Pneumo-Typer_workplace/result_statistics/tbl_part /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/Pneumo-Typer_workplace/Gcluster_workplace/interested_gene.txt /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/Pneumo-Typer_workplace/Gcluster_workplace/Sub_TFT 20000

use strict;
use warnings;

my $parttbl_dir = $ARGV[0];
opendir (PAR, $parttbl_dir);
my @par = readdir PAR;
@par = grep($_ !~ /^\./, @par);
close PAR;

open (OUT, ">$ARGV[1]");
my $Gcluster_tft_part= $ARGV[2];
mkdir $Gcluster_tft_part;
foreach my $file (@par){
    my %hash;
    my %tft;
    my $eachfile = "$parttbl_dir/$file";
    my $outfile = "$Gcluster_tft_part/$file.tbl_XXX.part";
    my @position;
    open (EA, $eachfile);
    while(<EA>){
        chomp;
        my @arr = split /\t/, $_;
        my $gene_locus = $arr[5];
        $gene_locus =~ s/.*;//g;

        push @{$tft{$arr[7]}}, "$arr[2]\t$arr[3]\tCDS\t$gene_locus\t$arr[6]\t$arr[7]";
        push @position, $arr[2];
    }
    close EA;
    my $median = median_value(\@position);
    my $c_arr;
    my $max=0;
    foreach (keys %tft){
        if (scalar @{$tft{$_}} > $max){
            $max = scalar @{$tft{$_}};
            $c_arr = $_; 
        }
    }


    open(TFT_PART, ">$outfile");
    foreach (@{$tft{$c_arr}}){
        my @array = split /\t/, $_; 
        print TFT_PART "$_\n" if abs($array[0]-$median) <= $ARGV[3];   
        push @{$hash{$c_arr}}, "$array[3]\t#$array[4];$file" if abs($array[0]-$median) <= $ARGV[3];   
    }
    print OUT ${$hash{$c_arr}}[0], "\n";
}



close OUT;
close TFT_PART;


sub median_value {

    my $arr_ref = shift;
    my $median;
    my $mid = int @$arr_ref/2;
    my @sorted_values = sort {$a<=>$b} @$arr_ref;
    if (@$arr_ref % 2) {
        $median = $sorted_values[ $mid ];
    } else {
        $median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
    } 

    return $median;
}


