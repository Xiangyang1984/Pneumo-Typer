use strict;
use warnings;

# perl $home_directory/script/interested_gene.pl $home_directory/DATABASE/gene.txt $workplace/Serotype.out $workplace/result_statistics/tbl_part $gcluster_workplace/interested_gene.txt $subtft 20000 1
my $gene_order = $ARGV[0];
my $serotype_out = $ARGV[1];
my $parttbl_dir = $ARGV[2];
my $interestedgene_out = $ARGV[3];
my $Gcluster_tft_part= $ARGV[4];
mkdir $Gcluster_tft_part;

my $dis_len = $ARGV[5];
my $gene_num = $ARGV[6];

open(GO, $gene_order);
my %g; #serotype."_".gene vs "order number"
while(<GO>){
    chomp;
    my ($s, $g, $o) = split /\t/, $_;
    $g{$s."_".$g} = $o;
}
close GO;

open(S, $serotype_out);
my %gs; #genome vs firstserotype
my %gsfull; #genome vs fullserotype
while(<S>){
    chomp;
    my ($gg, $ss) = split /\t/, $_;
    $gsfull{$gg} = $ss;
    $ss =~ s/\|.*|^0//g;  #using the first serotype when genome was predicted to belong to multiple serotypes, e.g A|B
    $gs{$gg} = $ss;
}
close S;


opendir (PAR, $parttbl_dir);
my @par = readdir PAR;
@par = grep($_ !~ /^\./, @par);
close PAR;

open (OUT, ">$interestedgene_out");

foreach my $file (@par){
    my %go; # "gene locus" vs "gene order"
    my %tft; # "contig" vs "tbl inf", store tbl information for each contig
    my $eachfile = "$parttbl_dir/$file";
    my $outfile = "$Gcluster_tft_part/$file.tbl_XXX.part";  #add the serotype inf $gs{$file} in the file name
    my %position;
    open (EA, $eachfile);
    while(<EA>){
        chomp;
        my @arr = split /\t/, $_;
        my $gene_locus = $arr[5];
        $gene_locus =~ s/.*;//g;
        $arr[2] =~ s/^\>//;
        $arr[2] =~ s/^\<//;
        push @{$tft{$arr[7]}}, "$arr[2]\t$arr[3]\tCDS\t$gene_locus\t$arr[6]\t$arr[7]";
        push @{$position{$arr[7]}}, $arr[2];
        $go{$gene_locus} = $g{$gs{$file}."_".$arr[1]} if defined $g{$gs{$file}."_".$arr[1]}; # serotype."_".gene may be not defined in %g, as some predicted genes are not contained in type serotype (93 serotypes gbk files, e.g. 19F did not contains glf gene)
    }
    close EA;

    my $c_arr; #contig where cps gene cluster located
    my $max=0;
    foreach (keys %tft){
        if (scalar @{$tft{$_}} > $max){
            $max = scalar @{$tft{$_}};
            $c_arr = $_; 
        }
    }

    my $median = median_value(\@{$position{$c_arr}});
    
    my %hash;
    open(TFT_PART, ">$outfile") if scalar @{$tft{$c_arr}} > $gene_num; #only genome having the number of cps gene >$gene_num is shown
    foreach (@{$tft{$c_arr}}){
        next unless scalar @{$tft{$c_arr}} > $gene_num;
        my @array = split /\t/, $_; 
        print TFT_PART "$_\n" if abs($array[0]-$median) <= $dis_len;   
        push @{$hash{$c_arr}}, "$array[3]\t#$array[4];$file:$array[5]" if abs($array[0]-$median) <= $dis_len;   
    }

    my @go;
    my %gf; # "gene order" vs "#gene locus\tproduct;genome:contig"
    foreach (@{$hash{$c_arr}}){
        my $gf = $_;
        $gf =~ s/\t.*//g;       
        push @go, $go{$gf} if defined $go{$gf}; 
        $gf{$go{$gf}} = $_ if defined $go{$gf};
        #print "111\t$gf{$gf}\t$gf\t$go{$gf}\n";
    }

    @go = sort{$a<=>$b} @go;
   
    print OUT $gf{$go[0]}, "\n" if ( (scalar @{$tft{$c_arr}} > $gene_num) && (scalar @go > 0) ); 
    print OUT ${$tft{$c_arr}}[0], "\n" if ( (scalar @{$tft{$c_arr}} > $gene_num) && (scalar @go eq 0) ); 


}



close OUT;
close TFT_PART;


sub median_value {

    my $arr_ref = shift;
    my $median;
    my $mid = int @$arr_ref/2;
    my @sorted_values = sort {$a<=>$b} @$arr_ref;
    if (@$arr_ref % 2) {
        $median = $sorted_values[$mid];
    } else {
        $median = ($sorted_values[$mid-1] + $sorted_values[$mid])/2;
    } 

    return $median;
}


