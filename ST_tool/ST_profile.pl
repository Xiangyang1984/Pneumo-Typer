use strict;
use warnings;
use Getopt::Long;
use FindBin;
use File::Basename qw<basename dirname>;
### perl /home/xiangyang/MG_Pipeline_v1.01/ST_tool/ST_profile.pl /home/xiangyang/MG_Pipeline_v1.01/Kp_genomes/prokka_workplace/prokka_gbk_folder $path_gene /home/xiangyang/MG_Pipeline_v1.01/ST_tool/Klebsiella_pneumoniae 190 /home/xiangyang/MG_Pipeline_v1.01/Kp_genomes/ST_workplace

my $path_genbank = $ARGV[0];
my $path_gene = $ARGV[1];
my $mlsa_data_folder = $ARGV[2];
my $thread_number = $ARGV[3];
my $ST_workplace = $ARGV[4];
mkdir $ST_workplace;

my $mlsa_fasta;
my $mlsa_profile;
foreach my $file (glob "$mlsa_data_folder/*"){
    next unless $file !~ /^\./;
    if (-e "$file/"){
        $mlsa_fasta = $file;
    }else{
        $mlsa_profile = $file;
    }

}
#print "$mlsa_fasta\n$mlsa_profile\n";

my $home_directory = $FindBin::Bin;   # obtaining the home directory where ST_profile.pl located

#############################################################################
#### extract gene and protein sequences from gbks under the same folder, and obtain path_gene and path_protein folder
#############################################################################
#system("perl $home_directory/gene_bacth_extract.pl $path_genbank $thread_number $ST_workplace");   



#############################################################################
#### do blast, and obtain blastn_out.txt file
#############################################################################
system("perl  $home_directory/bacth_blast.pl $path_gene $mlsa_fasta $thread_number $ST_workplace"); 


my $blastn_out = "$ST_workplace/blastn_out.txt";


my %hash_profile;

open (PROFILE, $mlsa_profile);
my @mlsa_profile = <PROFILE>;
my $first_line = shift @mlsa_profile;
chomp ($first_line);
my @first_line = split "\t", $first_line;

# obtain mlst gene name as array, e.g.: ("gapA", "infB", "mdh","pgi", "phoE", "rpoB", "tonB");
my @standard_ST = ($first_line[1], $first_line[2], $first_line[3], $first_line[4], $first_line[5], $first_line[6], $first_line[7]);
#print join ("-", @standard_ST), "\n";
foreach(@mlsa_profile){
    chomp;
    #next unless $_ !~ /^ST/;
    my @profile = split "\t", $_;
    #$hash_profile{$profile[1]."-".$profile[2]."-".$profile[3]."-".$profile[4]."-".$profile[5]."-".$profile[6]."-".$profile[7]}=$profile[0]; #gene number profile
    $hash_profile{"$standard_ST[0]_".$profile[1]."-"."$standard_ST[1]_".$profile[2]."-"."$standard_ST[2]_".$profile[3]."-"."$standard_ST[3]_".$profile[4]."-"."$standard_ST[4]_".$profile[5]."-"."$standard_ST[5]_".$profile[6]."-"."$standard_ST[6]_".$profile[7]}=$profile[0]; #gene_number profile
}

close PROFILE;

my %input_hash;
open(IN, $blastn_out) or die $!;

while (<IN>){
    
    chomp;
    $_ =~ s/ //g;
    my @array = split /\t/, $_;
    next unless ( ($array[2] == 100) && ($array[4] eq $array[5]) );
    $array[0] =~ s/.*;//g;  # obtain strain name or assembly name
    my $gene = $array[1];
    $gene =~ s/_\d+//g;  #obtain gene e.g. "Klebsiella_pneumoniae_MLST_gapA"
    
    #$input_hash{$array[0]."_".$array_gene[0]} = $array_gene[1];  #key: assembly_gene   value: gene number
    $input_hash{$array[0]."_".$gene} = $array[1];  #key: assembly_gene   value: gene_number
    #print $array[0]."_".$gene,"\t", $array[1],"\n";
}

close IN;

chdir $path_genbank;
my $strain_name =  "$ST_workplace/strain_name.txt";
system ("ls > $strain_name");

open(STRAIN, $strain_name) or die $!;

my $ST_out = dirname($ST_workplace)."/ST_out.txt";
open (ST_OUT, ">$ST_out");
print ST_OUT "Strain_name\tST\tProfile\n";
while (<STRAIN>) {
    chomp $_;
    my $isolate = $_;
    #print "$isolate\t";
    my $seq_type;
    for (my $i=0; $i<7; $i++) {

        if( defined $input_hash{$isolate."_".$standard_ST[$i]} ) {
            $seq_type .=$input_hash{$isolate."_".$standard_ST[$i]}."-";    
        }else { 
            $seq_type .="#"."-";
        }
      
    }
    $seq_type =~ s/\-$//g;
    #print "$seq_type\n";
    if (defined $hash_profile{$seq_type}){
        print ST_OUT "$isolate\t$hash_profile{$seq_type}\t$seq_type\n";
    }else{        
        if ($seq_type !~ /#/){
            print ST_OUT "$isolate\tNovel\t$seq_type\n";
        }else{
            print ST_OUT "$isolate\tUnknown\t$seq_type\n";
        }
    }
}

close STRAIN;
