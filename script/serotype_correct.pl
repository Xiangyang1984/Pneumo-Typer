#perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/script/serotype_correct.pl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/serotype_matrix.txt /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/Pneumo-Typer_workplace/result_statistics/Statistics_OUT/Antibiotic_LocusTag_statistics /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/Pneumo-Typer_workplace/result_statistics/Statistics_OUT/Antibiotic_LocusTag_statistics.out

use strict;
use warnings;

my $correct_matrix = $ARGV[0];
open(MATRIX, $correct_matrix);
my %hash;
while(<MATRIX>){
    chomp;
    next unless $_ !~ /^\s*$/;
    next unless $_ !~ /^#/;

    my @array = split /\t/, $_;
    my @arr = split /, /, $array[1];
    my @ele;
    foreach my $ele (@arr){
        $ele =~ s/: .*//g if $array[2] =~ /^#presence_absence|^#frameshift mutation/;
        push @ele, (split /\|/, $ele) if $ele =~ /\|/;
        push @ele, $ele if $ele !~ /\|/;       
    }
    my $ref = join ("|",@ele);
    $hash{$ref} = $array[0];
    #print "$ref\n";
}
close MATRIX;

my $primary = $ARGV[1];
my $out = $ARGV[2];
open (OUT, ">$out");

open(PRI, $primary);
while(<PRI>){
    chomp;    
    print  OUT "Strain\tCorrect serotype\n" if $_ =~ /^Strain/;
    next unless $_ !~ /^Strain/;

    my $full = $_;
    my @pri = split /\t/, $_;
    my $name = shift @pri;
    my $pre_serotype = pop @pri;

    $pre_serotype =~ s/\/.*//g;
    $pre_serotype =~ s/^Serotype //;
    $pre_serotype =~ s/SPC//g;
    my @ser;
    push @ser, $pre_serotype if $pre_serotype !~ /\|/;
    @ser = split /\|/, $pre_serotype if $pre_serotype =~ /\|/;
    #my $st = pop @pri;
    #my $true_sp = pop @pri;
    #my $full_predict = pop @pri;
    #next unless $true_sp !~ /NA|blank|not|NT|Not/;
    #$true_sp  =~ s/.* //g;
    #$true_sp =~ s/^0//g;

# to find unique gene based on predicted serotype 
############################################
    my $unique_gene="XXX";
    foreach my $key (sort keys %hash){
        my @a = split /\|/, $key;
        my @tmp = @a;
        push @tmp, @ser;
        my %num;
        foreach (@tmp){
            $num{$_}=$_;
        }
        if (scalar @a == scalar keys %num){
            $unique_gene = $hash{$key};
            last;
        }
    }
############################################

    my $correct_ser = $pre_serotype;
    my $string = join ("", @pri);
    if ($string =~ /\|\Q$unique_gene\E\-SPC(.*?)\|/){
        $correct_ser = $1; 
    }
    elsif ($string =~ /\|\Q$unique_gene\E\|/){
        $correct_ser = $unique_gene;
        $correct_ser =~ s/.*\-SPC//g;
        $correct_ser =~ s/$/A/ if $unique_gene =~ /22$/;  #22A and 22F
        #print "$name\t$correct_ser\t$pre_serotype\n";
    }
    elsif ($correct_ser =~ /22/){
        $correct_ser = "22A" if $string =~ /\|\Q$unique_gene\E\|/;
        $correct_ser = "22F" if $string !~ /\|\Q$unique_gene\E\|/;
    }
    $correct_ser =~ s/^0//g if $correct_ser =~ /^0/;
    #print OUT "$name\t$correct_ser\t$pre_serotype\n";
    print  OUT "$name\t$correct_ser\n";
}

close PRI;
close OUT;

sub getMax {
    my (%temp, $max, @ret);
    foreach ( @_ ) {
        my $count = $temp{$_}++;
        $max = $count if $count > $max; 
    }
    foreach ( keys %temp ) {
        push @ret, $_ if $temp{$_} == $max
    }    
    return @ret;
}
