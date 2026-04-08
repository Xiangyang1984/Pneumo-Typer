use strict;
use warnings;

# perl capsule_genotype.pl CPS_LocusTag_statistics capsule_genotype.out
my %h;
my $CPS_LocusTag_statistics = $ARGV[0];
my $capsule_genotype_out = $ARGV[1];
open(CLS, $CPS_LocusTag_statistics);
<CLS>; #remove header
while(<CLS>){
    chomp;
    $_ =~ s/\n|\r//g;
    my @a = split /\t/, $_;
    my $h = $a[0];
    my @d;
    foreach my $e (@a){
        next unless $e =~ /\|/;        
        if ($e =~ /\:/){
            my @t0;            
            my @b0 = split /\:/, $e;
            foreach(@b0){
                my @c0 = split /\|/, $_;
                push @t0, $c0[1];
            }
            my $str0 = join(":", @t0);
            push @d, $str0;
        }else{
            my @e = split /\|/, $e;
            push @d, $e[1];
        }       
    }
    my $str = join("|", @d);
   
    my $pr = $a[-1]; #Serotype SPC22F|SPC22A/18/88.89%
    $pr =~ s/(.*) //;
    $pr =~ s/\/(.*)/\t\/$1/;
    $pr =~ s/SPC//g;
    #print "$pr\n";
    my @pr = split /\t/, $pr;
    my @pr1 = split /\|/, $pr[0] if $pr[0] =~ /\|/;
    #print "$pr[1]\n";
    my $pr1;
    $pr1 = join("|SPC", sort @pr1) if $pr[0] =~ /\|/;
    $pr1 = $pr[0] if $pr[0] !~ /\|/;
    $pr1 = "SPC".$pr1;
    my $pr2 = "$pr1"."$pr[1]";
    #print "$pr2\n";
    push @{$h{$pr2}}, "$h\t$str";
    
}
close CLS;

open(CAP, ">$capsule_genotype_out");

print CAP "Genome\tCapsule genotype (Serotype/gene number/frequency)\tCapsule gene details\tSubtype number of the capsule genotype\tSubtype index of the capsule genotype\tGenomes number for each subtype of the capsule genetype\tUnique genes in Ref subtype of the capsule genotype\tUnique genes or different count in none-Ref subtype of the capsule genotype\n";
foreach my $k (sort keys %h){
    my @g = sort @{$h{$k}};
    my %a;
    foreach(@g){
        my ($genome, $profile) = split /\t/, $_; #remove genome
        $a{$profile}++;
    }
    my %t;
    my $type=0;
    foreach(sort keys %a){
        $type++;
        $t{$_} = $type;
    }
    my $first_arr = $g[0];
    
    $first_arr =~ s/.*\t//; #remove genome 
    my @ref = split /\|/, $first_arr;
    my %ref_counts = count_elements(\@ref); # 计算参照数组的元素出现次数 
    foreach (@g){
        my ($genome, $profile) = split /\t/, $_; #remove genome
        #$_ =~ s/(.*)\t//; #remove genome
        my @vs = split /\|/, $profile;
        my @genome = split /\#/, $genome;
        my $capsulegenotype; 
        if (scalar keys %a > 1){
            $capsulegenotype = $k."-".$t{$profile} if scalar keys %a > 1;
        } else {
            $capsulegenotype = $k;
        }
        print CAP "$genome[0]\t$capsulegenotype\t$profile\t", scalar keys %a, "\t$t{$profile}\t$a{$profile}\t";
        print CAP &compare_arrays_with_ref(\%ref_counts, \@vs, "array"); # 比较并打印每个数组与参照数组之间的差异
        print CAP "\n";         
    }
}
close CAP;
 
  
# 函数来计算元素在数组中的出现次数  
sub count_elements {  
    my $array = shift;  
    my %counts;  
    foreach my $elem (@$array) {  
        $counts{$elem}++;  
    }  
    return %counts;  
}  
  
# 函数来比较数组与参照数组的差异  
sub compare_arrays_with_ref {  
    my ($ref_counts, $comp_array, $array_name) = @_;

   
    my %comp_counts = count_elements(\@$comp_array); 

    my @diff_in_ref = ();  
    my @diff_not_in_ref = ();  
  
    # 查找在参照数组中但出现次数多于比较数组的元素  
    foreach my $elem (keys %$ref_counts) {  
        if (!exists $comp_counts{$elem} || $ref_counts->{$elem} > $comp_counts{$elem}) {  
            push @diff_in_ref, "$elem, $ref_counts->{$elem} in ref" if !exists $comp_counts{$elem};  
            push @diff_in_ref, "$elem, $ref_counts->{$elem} in ref vs $comp_counts{$elem} in $array_name" if exists $comp_counts{$elem};
        }  
    }  
  
    # 查找在比较数组中但不在参照数组中的元素，或者出现次数多于参照数组的元素  
    foreach my $elem (keys %comp_counts) {  
        if (!exists $ref_counts->{$elem} || $comp_counts{$elem} > $ref_counts->{$elem}) {  
            push @diff_not_in_ref, "$elem, $comp_counts{$elem} in $array_name" if !exists $ref_counts->{$elem};  
            push @diff_not_in_ref, "$elem, $comp_counts{$elem} in $array_name vs $ref_counts->{$elem} in ref" if exists $ref_counts->{$elem};
        }  
    }  
  
    # 打印结果  
    my @result;
    if (@diff_in_ref) {  
        #print "In ref_array but not the same in $array_name: ";  
        push  @result, join("|", @diff_in_ref);  
    }  
  
    if (@diff_not_in_ref) {  
        #print "In $array_name but not in ref_array or different count: ";  
        push  @result, join("|", @diff_not_in_ref);  
    }  

    return join("\t", @result);
}
