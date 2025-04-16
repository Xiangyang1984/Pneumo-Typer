use strict;
use warnings;
use Getopt::Long;
use FindBin;
use File::Basename qw<basename dirname>;
use Bio::SeqIO;
use threads;
#use vars qw(@match_str);

# perl ST_profile_cg.pl $path_genbank $path_fa database/cgmlst $thread_number $workplace;

### perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/ST_profile_cg.pl /home/xiangyang/MG_Pipeline_v1.01/Kp_genomes/prokka_workplace/prokka_gbk_folder $path_gene /home/xiangyang/MG_Pipeline_v1.01/ST_tool/Klebsiella_pneumoniae 190 /home/xiangyang/MG_Pipeline_v1.01/Kp_genomes/ST_workplace

# perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/ST_profile_cg.pl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/test /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/Pneumo-Typer_workplace/Whole_gene /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/database/cgmlst 3 /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/Pneumo-Typer_workplace/ST_workplace > /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/Pneumo-Typer_workplace/ST_workplace/profile.txt


# perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/ST_tool/ST_profile_cg2.pl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/test1 /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/ST_tool/database/cgmlst 5 /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/ST_tool/TEST_OUT
#our %match_str;
my $path_genbank = $ARGV[0];
my $path_fa = $ARGV[1];
my $mlsa_data_dir = $ARGV[2];
my $thread_number = $ARGV[3];
my $ST_workplace = $ARGV[4];
mkdir $ST_workplace;

my $mlsa_fasta="$mlsa_data_dir/cgMLSA_loci";
my $mlsa_profile="$mlsa_data_dir/cgMLST_profiles";
my $mlsa_genename="$mlsa_data_dir/cgMLSA_loci.name";


my $home_directory = $FindBin::Bin;   # obtaining the home directory where ST_profile.pl located

#############################################################################
#### do blast, and obtain blastn_out.txt file
#############################################################################

system("perl $home_directory/bacth_blat.pl $path_fa $mlsa_fasta $mlsa_genename $thread_number $ST_workplace"); 


my %hash_profile;

open (PROFILE, $mlsa_profile);
my @mlsa_profile = <PROFILE>;
my $first_line = shift @mlsa_profile;
chomp ($first_line);
my @standard_ST = split "\t", $first_line;

foreach(@mlsa_profile){
    chomp;
    my @profile = split "\t", $_;
    my $key;
    for (my $k=1; $k<1223; $k++){
        $profile[$k] =~ s/N/N\|\\d\+/;
        $key .= "(".$profile[$k].")"."-" if $profile[$k] =~ "N";
        $key .= $profile[$k]."-" if $profile[$k] !~ "N";
    }
    $key =~ s/\-$//g;
    $hash_profile{$key}=$profile[0]; #gene_number profile
    
}

close PROFILE;


my $blastn_out_dir = "$ST_workplace/cgblat_parse_dir";

chdir $path_genbank;
my $strain_name =  "$ST_workplace/strain_name.txt";
system ("ls > $strain_name");

my $temp_cgST_dir = "$ST_workplace/temp_cgST_dir";
mkdir $temp_cgST_dir;

my $contain_hash;
my $ST_out = dirname ($ST_workplace). "/cgST_out.txt";
#open ST_out, ">$ST_out";
#print ST_out "Strain_name\tcgST\tProfile\n";

&bacth_macth($strain_name, \%hash_profile, \@standard_ST, $thread_number, $blastn_out_dir, $temp_cgST_dir);

chdir $temp_cgST_dir;
system ("find . -type f -print | xargs cat > $ST_out");

#system("rm -rf $temp_cgST_dir");
#system("rm -rf $ST_workplace/all.txt");


sub bacth_macth{

    my ($strain_name, $ref_hash_profile, $ref_standard_ST, $thread_number, $blastn_out_dir, $temp_cgST_dir) = @_;

 
    open(STRAIN, $strain_name) or die $!;
    my @strain = <STRAIN>;

    my $sub_file_number = scalar @strain;  
  
    my $thread;
    my @threads;
    my $job_number=0;
    #print "sub_file_number\t$sub_file_number\n";
    while(){ 
        last if ($job_number eq $sub_file_number);  
        while(scalar(threads->list())<$thread_number) {     #set thread number
            $job_number++;  
	    #print "$job_number\n";
            my $isolate = $strain[$job_number-1];
            chomp $isolate;
            my $ef = "$blastn_out_dir/$isolate.cgblastn_out";
            my $out = "$temp_cgST_dir/$isolate";
            $threads[$job_number-1]=threads->new(\&match_array, $isolate, $ef, $ref_hash_profile, $ref_standard_ST, $out); 
            last if ($job_number eq $sub_file_number);  
        }
        foreach $thread(threads->list(threads::all)){
            if($thread->is_joinable()) {              
                $thread->join();                     
            }
        }
    }

    foreach my $th (threads->list(threads::all)) {
        $th->join();
    }

}


sub match_array {
    my ($isolate, $ef, $ref_hash_profile, $ref_standard_ST, $out) = @_;
    
    chomp $isolate;
    my %input_hash;
    open EF, $ef;
    while(<EF>){
        chomp;
        $_ =~ s/ //g;
        my @array = split /\t/, $_;
        my $mark = shift @array;

        $array[0] =~ s/.*;//g;  # obtain strain name or assembly name
        my $gene = $array[1];
        $gene =~ s/_(\d+)//g;  #obtain gene e.g. "Klebsiella_pneumoniae_MLST_gapA"
        $input_hash{$array[0]."_".$gene} = $1 if $mark eq "T";  #key: assembly_gene   value: gene_number
        $input_hash{$array[0]."_".$gene} = "N" if $mark eq "F";  #key: assembly_gene   value: gene_number
    }
    close EF;

    my %hash_profile = %$ref_hash_profile;
    my @standard_ST = @$ref_standard_ST;

    my $seq_type;
    for (my $i=1; $i<1223; $i++) {
        if( defined $input_hash{$isolate."_".$standard_ST[$i]} ) {
            $seq_type .=$input_hash{$isolate."_".$standard_ST[$i]}."-";    
        }else { 
            $seq_type .="#"."-";
        }     
    }
    $seq_type =~ s/\-$//g;
 
    my @temp_arr;
    open (OUT, ">$out");
    if (my @temp_arr = grep{$seq_type =~ /^$_$/} keys %hash_profile){
        print OUT "$isolate\t$hash_profile{$temp_arr[0]}\t$seq_type\n"; 
    }else{
        if ($seq_type !~ /N/){
            print OUT "$isolate\tNovel\t$seq_type\n";
        }else{
            print OUT "$isolate\tUnknown\t$seq_type\n";
        }
    }
    close OUT;
}

