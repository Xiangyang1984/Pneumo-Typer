#perl /home/xiangyang/Desktop/tool_dev/pneumo-typer-v1.0.2/script/serotype_correct.pl /home/xiangyang/Desktop/tool_dev/pneumo-typer-v1.0.2/DATABASE/serotype_matrix.txt /home/xiangyang/Desktop/tool_dev/pneumo-typer-v1.0.2/CPS_LocusTag_statistics /home/xiangyang/Desktop/tool_dev/pneumo-typer-v1.0.2/Serotype_sta.out

use strict;
use warnings;

my $correct_matrix = $ARGV[0];
open(MATRIX, $correct_matrix);
my %hash;
my %absence; # record gene absence information that is used to determine serotype keys: gene, values: serotypes were joined with "|"
my %presence; # record gene presence information that is used to determine serotype
while(<MATRIX>){
    chomp;
    next unless $_ !~ /^\s*$/;
    next unless $_ !~ /^#/;

    my @array = split /\t/, $_;    
    my @arr = split /, /, $array[1]; #wciG	35D: 0, 35B: 1
    
    my @ele; #record serotype
    my @temp0; #record serotype in which gene is absence
    my @temp1; #record serotype in which gene is presence
    #my $serogroup;
    foreach my $ele (@arr){
        if ($array[2] =~ /^#absence /){
            push @temp0, $1 if $ele =~ /(.*?): 0/;
            push @temp1, $1 if $ele =~ /(.*?): 1/;
            #$serogroup = $1;
            #$serogroup =~ s/[A-Z]$//;
        }
        $ele =~ s/: .*//g if $array[2] =~ /^#presence |^#absence|^#frameshift mutation/; 
        push @ele, (split /\|/, $ele) if $ele =~ /\|/;
        push @ele, $ele if $ele !~ /\|/;       
    }
    my $ref = join ("|",@ele);
    $hash{$ref} = $array[0];  #ref: all serotypes were joined with "|" for each row in serotype_matrix.txt   $array[0]: gene
    $absence{$array[0]} = join ("|", @temp0) if scalar @temp0 > 0; #record the gene belonging to absence 
    $presence{$array[0]} = join ("|", @temp1) if scalar @temp1 > 0; #record the gene belonging to absence

}
close MATRIX;

my $primary = $ARGV[1];
my $out = $ARGV[2];
open (OUT, ">$out");

open(PRI, $primary);
my @header_gene;
while(<PRI>){
    chomp;    
    print  OUT "Strain\tSerotype\n" if $_ =~ /^Strain/;
    @header_gene = split /\t/, $_ if $_ =~ /^Strain/;
    shift @header_gene if $_ =~ /^Strain/; #remove "Strain"
    pop @header_gene if $_ =~ /^Strain/;   #remove "Serotype/number/frequency"
    next unless $_ !~ /^Strain/;

    my $full = $_;
    my @pri = split /\t/, $_;    #cps gene distributed in genome
    my $name = shift @pri;       #genome name
    my $pre_serotype = pop @pri;

    $pre_serotype =~ s/\/.*//g;  #Serotype SPC03/8/100.00%
    $pre_serotype =~ s/^Serotype //;
    $pre_serotype =~ s/SPC//g;
    $pre_serotype =~ s/06EA/06A/;
    $pre_serotype =~ s/06EB/06B/;
    $pre_serotype =~ s/23B1/23B/;

    my @ser; #store serotype that is the first predict result
    push @ser, $pre_serotype if $pre_serotype !~ /\|/;
    @ser = split /\|/, $pre_serotype if $pre_serotype =~ /\|/;

# obtain the unique gene if @ser is the subset of @a
# the last is the prior unique gene
############################################
    my $unique_gene="XXX";
    my @unique_gene; #stroe gene for crossed serotype
    #push @unique_gene, $unique_gene;
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
            push @unique_gene, $unique_gene."#".$key;
            #last;
        }
    }
    #print scalar @header_gene, "\n";

    #put unique gene into @header_gene if unique gene was not found in the first row of CPS_LocusTag_statistics
    push @header_gene, $unique_gene if ( grep ($_ ne $unique_gene, @header_gene) && defined $absence{$unique_gene} );
    push @pri, "" if ( grep ($_ ne $unique_gene, @header_gene) && defined $absence{$unique_gene} );
    #print scalar @header_gene, "\n";
############################################

#obtain the gene index for gene absence in serotype_matrix.txt
############################################
    my $index = -1;
    my @rindex;
    foreach (@header_gene){
        $index++;
        push @rindex, $index if ( (defined $absence{$header_gene[$index]}) && ($header_gene[$index] eq $unique_gene) );

    }
############################################

    my $correct_ser = join ("|", sort @ser);
    my $string = join ("", @pri);
    $string =~ s/06EA/06A/g;
    $string =~ s/06EB/06B/g;
    $string =~ s/23B1/23B/g;
    if ( $unique_gene ne "XXX" ){

        if (scalar @unique_gene == 1){

            if ($string =~ /\|\Q$unique_gene\E\-SPC(.*?)\|/){
                $correct_ser = $1; #gene-SPC          
                if (scalar @rindex > 0) { #gene to determine the serotype
                    foreach my $rindex(@rindex){
                        $correct_ser = $presence{$header_gene[$rindex]} if ( (defined $absence{$header_gene[$rindex]}) && ($rindex != -1) && ($unique_gene eq $header_gene[$rindex]) );
                    }
                }
                $correct_ser =~ s/#SPC/\|/g if $correct_ser  =~ /#SPC/;
            } 
       
            elsif ($string =~ /\|\Q$unique_gene\E\|/){
                $correct_ser = $unique_gene; #gene-SPC
                $correct_ser =~ s/.*\-SPC//g;
            }

            elsif(scalar @rindex > 0){
                foreach my $rindex(@rindex){
                    $correct_ser = $absence{$header_gene[$rindex]} if ( (defined $absence{$header_gene[$rindex]}) && ($rindex != -1) && ($unique_gene eq $header_gene[$rindex]) && ($pri[$rindex] eq "") );
                }
            }
        }
        elsif (scalar @unique_gene > 1){
            my $k=0;
            foreach my $u (@unique_gene){
                $k++;
                my ($g, $s) = split /\#/, $u;
                my @temp_serotype = split /\|/, $s;

                if ($k==1){
                    if ($string =~ /\|\Q$g\E/){
                        $correct_ser = $temp_serotype[0];
                    }else{
                        $correct_ser = $temp_serotype[1];
                    }
                }
                
                if ($k==2){
                    print "$u\t$g\n";
                    if ($string =~ /\|\Q$g\E\|/){
                        $correct_ser = $g; #gene-SPC
                        $correct_ser =~ s/.*\-SPC//g;
                    }

                    elsif ($string =~ /\|\Q$g\E\-SPC(.*?)\|/){
                        $correct_ser = $temp_serotype[0]; #gene-SPC          
                        if (scalar @rindex > 0) { #gene to determine the serotype
                            foreach my $rindex(@rindex){
                                $correct_ser = $presence{$header_gene[$rindex]} if ( (defined $absence{$header_gene[$rindex]}) && ($rindex != -1) && ($g eq $header_gene[$rindex]) );
                            }
                        }
                        $correct_ser =~ s/#SPC/\|/g if $correct_ser  =~ /#SPC/;

                    }
                    else{
                        if(scalar @rindex > 0){
                            foreach my $rindex(@rindex){
                                $correct_ser = $absence{$header_gene[$rindex]} if ( (defined $absence{$header_gene[$rindex]}) && ($rindex != -1) && ($g eq $header_gene[$rindex]) && ($pri[$rindex] eq "") );
                            }
                        }
                    }
                }

            }
        }
        
    }

    print OUT "$name\t$correct_ser\n";
    pop @header_gene if ( grep ($_ ne $unique_gene, @header_gene) && defined $absence{$unique_gene} ); #related to line 90
    pop @pri if ( grep ($_ ne $unique_gene, @header_gene) && defined $absence{$unique_gene} ); #related to line 91
    
}


close PRI;
close OUT;
