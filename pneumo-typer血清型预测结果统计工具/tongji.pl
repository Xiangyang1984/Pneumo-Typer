# perl /home/xiangyang/Desktop/tool_dev/pneumo-typer-v1.0.2/pneumo-typer血清型预测结果统计工具/tongji.pl  Serotype_sta.out> /home/xiangyang/Desktop/tool_dev/pneumo-typer-v1.0.2/pneumo-typer血清型预测结果统计工具/accrute.txt

use strict;
use warnings;

open (Y, "$ARGV[0]");
my %r_hash;
my %w_hash;
my %winf_hash;
my %t_hash;
while(<Y>){
    chomp;
    my ($a, $b, $c, $d) = split /\t/, $_; #GCA_001098885.1	1	1	Right
    $c =~ s/^0//g;
    $c =~ s/\|0/\|/g;
    $b =~ s/^0//g;
    $b =~ s/\|0/\|/g;
    $r_hash{$c}++ if $d eq "Right";
    $w_hash{$c}++ if $d eq "Wrong";
    push @{$winf_hash{$c}}, $b if $d eq "Wrong";
    $t_hash{$c}++;
#print "$c\t$t_hash{$c}\n";

}
close Y;

my %w_inf;
foreach my $s (keys %w_hash){
    my %ref_hash = &arr_num(\@{$winf_hash{$s}});
    my @str;
    foreach my $e (keys %ref_hash){
        push @str, $e.",".$ref_hash{$e};
    }
    $w_inf{$s} = join(";", @str);
}

foreach (sort keys %t_hash){
    $r_hash{$_} = 0 if !defined $r_hash{$_};
    $w_hash{$_} = 0 if !defined $w_hash{$_};
    #$t_hash{$_} = 0 if !defined $t_hash{$_};
    my $p = ($r_hash{$_}/$t_hash{$_})*100;
    $p = sprintf("%.2f", $p);
    my $seq = $_;
    $seq =~ s/[A-Z]//g;
    print "$seq\t$_\t$t_hash{$_}\t$p\t$r_hash{$_}\t$w_hash{$_}\t$w_inf{$_}\n" if defined $w_inf{$_};
    print "$seq\t$_\t$t_hash{$_}\t$p\t$r_hash{$_}\t$w_hash{$_}\n" if !defined $w_inf{$_};

}

sub arr_num {
    my $ref_arr = shift;

    my $key;
    my $value;
    my %n_hash;
    foreach (@$ref_arr){
        ++$n_hash{$_};
    }

    return %n_hash;
}
