#!/usr/bin/perl -w

use strict;
use threads;
use Bio::SeqIO;
use Getopt::Long;
use FindBin;

# perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/script/ser_join_ST.pl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/Pneumo-Typer_workplace/Serotype.out /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/Pneumo-Typer_workplace/ST_out.txt /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/Pneumo-Typer_workplace/cgST_out.txt /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/Pneumo-Typer_workplace/Serotype_ST.out /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/genomes_metadata.txt

my $pred_serotype = $ARGV[0];
my $st = $ARGV[1];
my $cgst = $ARGV[2];
my $final = $ARGV[3];

#######
my %meta;
if ($ARGV[4]){  #metedata extracted from gbk files
    open(META, $ARGV[4]);
    while(<META>){
        chomp;
        $meta{(split /\t/, $_)[0]} = (split /\t/, $_)[4];
    }
    close META;
}


#######
my %ST;
if (-e $st){
    open(ST, $st);
    while(<ST>){
        chomp;
        my @genome_array = split "\t", $_;
        my $assembly = $genome_array[0];
        #$assembly =~ s/\-/_/g;
        $ST{$assembly} = $genome_array[1];
    }
    close ST;
}

#######
my %cgST;
if (-e $cgst){
    open(CGST, $cgst);
    while(<CGST>){
        chomp;
        my @cggenome_array = split "\t", $_;
        my $cgassembly = $cggenome_array[0];
        #$assembly =~ s/\-/_/g;
        $cgST{$cgassembly} = $cggenome_array[1];
    }
    close CGST;
}

#######
open(LOCUS_FILE, $pred_serotype);
open (FINAL, ">$final");
my @array = <LOCUS_FILE>;
my $header = shift @array;
chomp $header;
if (defined $meta{(split /\t/, $array[0])[0]}){
    print FINAL "$header\tTrue serotype\tST\tcgST\n" if (-e $st && -e $cgst); 
    print FINAL "$header\tTrue serotype\tcgST\n" if ( !(-e $st) && -e $cgst); 
    print FINAL "$header\tTrue serotype\tST\n" if (-e $st && !(-e $cgst) );
    print FINAL "$header\n" if ( !(-e $st) && !(-e $cgst) ); 
}else{
    print FINAL "$header\tST\tcgST\n" if (-e $st && -e $cgst); 
    print FINAL "$header\tcgST\n" if ( !(-e $st) && -e $cgst); 
    print FINAL "$header\tST\n" if (-e $st && !(-e $cgst) );
    print FINAL "$header\n" if ( !(-e $st) && !(-e $cgst) );
}

foreach(@array){
    chomp;
    if (defined $meta{(split /\t/, $_)[0]}){
        print FINAL "$_\t", $meta{(split /\t/, $_)[0]}, "\t", $ST{(split /\t/, $_)[0]}, "\t", $cgST{(split /\t/, $_)[0]}, "\n" if (-e $st && -e $cgst); 
        print FINAL "$_\t", $meta{(split /\t/, $_)[0]}, "\t", $ST{(split /\t/, $_)[0]}, "\n" if (-e $st && !(-e $cgst) );
        print FINAL "$_\t", $meta{(split /\t/, $_)[0]}, "\t", $cgST{(split /\t/, $_)[0]}, "\n" if ( !(-e $st) && -e $cgst); 
        print FINAL "$_\t", $meta{(split /\t/, $_)[0]}, "\n" if ( !(-e $st) && !(-e $cgst) );
    }else{
        print FINAL "$_\t", $ST{(split /\t/, $_)[0]}, "\t", $cgST{(split /\t/, $_)[0]}, "\n" if (-e $st && -e $cgst);
        print FINAL "$_\t", $ST{(split /\t/, $_)[0]}, "\n" if (-e $st && !(-e $cgst) );
        print FINAL "$_\t", $cgST{(split /\t/, $_)[0]}, "\n" if ( !(-e $st) && -e $cgst); 
        print FINAL "$_\n" if ( !(-e $st) && !(-e $cgst) );
    }
}

close LOCUS_FILE;
close FINAL;

