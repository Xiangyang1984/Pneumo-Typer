use strict;
use warnings;
use Getopt::Long;
use FindBin;
use File::Basename qw<basename dirname>;
use Bio::SeqIO;
use threads;

### perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/ST_profile_cg.pl /home/xiangyang/MG_Pipeline_v1.01/Kp_genomes/prokka_workplace/prokka_gbk_folder $path_gene /home/xiangyang/MG_Pipeline_v1.01/ST_tool/Klebsiella_pneumoniae 190 /home/xiangyang/MG_Pipeline_v1.01/Kp_genomes/ST_workplace

# perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/ST_profile_cg.pl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/test /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/Pneumo-Typer_workplace/Whole_gene /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/database/cgmlst 3 /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/Pneumo-Typer_workplace/ST_workplace > /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/Pneumo-Typer_workplace/ST_workplace/profile.txt


# perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/ST_tool/ST_profile_cg.pl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/test1 /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/ST_tool/database/cgmlst 5 /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/ST_tool/TEST_OUT

my $path_genbank = $ARGV[0];
my $mlsa_data_folder = $ARGV[1];
my $thread_number = $ARGV[2];
my $ST_workplace = $ARGV[3];
mkdir $ST_workplace;

my $path_fa = "$ST_workplace/whole_fasta";
mkdir $path_fa;

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
#### do blast, and obtain blastn_out.txt file
#############################################################################

&batch_genbank_sequence_extract ($path_genbank, $path_fa, $thread_number);
opendir PATH_FA, $path_fa or die "could not open $path_fa";
my @PATH_FA = readdir PATH_FA; 
@PATH_FA =grep ($_!~/^\./ ,@PATH_FA);  #delete hidden file . ..
closedir PATH_FA;

my $number = scalar @PATH_FA;

my $count=0;
print "    Blat_percent: ";
foreach (@PATH_FA) {
    my $progress_record = int (($count/$number)*100);
    $count++;
    my $Blat_perform_progress = int (($count/$number)*100); 
    print "$Blat_perform_progress%","..." if ($count == 1 or ( ($Blat_perform_progress%10 ==0) && ($progress_record <$Blat_perform_progress)) );

    my $fa_file = "$path_fa/$_";
    system("perl $home_directory/bacth_blat.pl $mlsa_fasta $fa_file $thread_number $ST_workplace"); 
}
print "done\n";

my $blastn_out = "$ST_workplace/cgblastn_out.txt";


my %hash_profile;

open (PROFILE, $mlsa_profile);
my @mlsa_profile = <PROFILE>;
my $first_line = shift @mlsa_profile;
chomp ($first_line);
my @standard_ST = split "\t", $first_line;
#shift @standard_ST; # remove the first element ST or cgST
# obtain mlst gene name as array, e.g.: ("gapA", "infB", "mdh","pgi", "phoE", "rpoB", "tonB");
#my @standard_ST = ($first_line[1], $first_line[2], $first_line[3], $first_line[4], $first_line[5], $first_line[6], $first_line[7]);
#print join ("-", @standard_ST), "\n";
foreach(@mlsa_profile){
    chomp;
    #next unless $_ !~ /^ST/;
    my @profile = split "\t", $_;
    #$hash_profile{$profile[1]."-".$profile[2]."-".$profile[3]."-".$profile[4]."-".$profile[5]."-".$profile[6]."-".$profile[7]}=$profile[0]; #gene number profile
    my $key;
    for (my $k=1; $k<1223; $k++){

        $key .= "$standard_ST[$k]_".$profile[$k]."-";
    }
    $key =~ s/\-$//g;
    $hash_profile{$key}=$profile[0]; #gene_number profile
    #print "$profile[0]\t$key\n";
    
}

close PROFILE;

my %input_hash;
open(IN, $blastn_out) or die $!;

while (<IN>){
    
    chomp;
    $_ =~ s/ //g;
    my @array = split /\t/, $_;
    my $mark = shift @array;
    #next unless ( ($array[2] == 100) && ($array[4] eq $array[5]) );

    #next unless ( ($array[2] == 100) && ($array[3] eq ($array[9]-$array[8]+1)) );
    #print "$array[1]\t$array[0]\t$array[2]\t$array[3]\t$array[8]\t$array[9]\t$array[10]\n";
    $array[0] =~ s/.*;//g;  # obtain strain name or assembly name
    my $gene = $array[1];
    $gene =~ s/_\d+//g;  #obtain gene e.g. "Klebsiella_pneumoniae_MLST_gapA"
    $input_hash{$array[0]."_".$gene} = $array[1] if $mark eq "T";  #key: assembly_gene   value: gene_number
    $input_hash{$array[0]."_".$gene} = $gene."_N" if $mark eq "F";  #key: assembly_gene   value: gene_number
}

close IN;

chdir $path_genbank;
my $strain_name =  "$ST_workplace/strain_name.txt";
system ("ls > $strain_name");

open(STRAIN, $strain_name) or die $!;

my $ST_out = dirname($ST_workplace)."/cgST_out.txt";
open (ST_OUT, ">$ST_out");
print ST_OUT "Strain_name\tST\tProfile\n";
while (<STRAIN>) {
    chomp $_;
    my $isolate = $_;
    #print "$isolate\t";
    my $seq_type;
    for (my $i=1; $i<1223; $i++) {

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
        if ($seq_type !~ /N/){
            print ST_OUT "$isolate\tNovel\t$seq_type\n";
        }else{
            print ST_OUT "$isolate\tUnknown\t$seq_type\n";
        }
    }
}

close STRAIN;



####################################################
####################################################
sub batch_genbank_sequence_extract {

    my ($path, $fasta, $thread_number)=@_;

    opendir DIR, $path or die $!;
    my @dir = readdir DIR; 
    @dir = grep ($_!~/^\./, @dir);
    closedir DIR;

    my @input;
    my @outfile;

    foreach my $file(@dir){

        my $input="$path/$file";
        push (@input,$input);  
        my $output="$fasta/$file";
        push (@outfile,$output);
    }

    my $sub_file_number = scalar @input;
    my $thread;
    my @threads;
    my $job_number=0;

    while(){ 

        last if ($job_number eq $sub_file_number);  

        while(scalar(threads->list())<$thread_number) {     #set thread number
            $job_number++;    
            my $input_file = $input[$job_number-1];
            my $output_file = $outfile[$job_number-1];
            #print "Genbank_extraction: $job_number\n";
            $threads[$job_number-1]=threads->new(\&genbank_fasta_sequence_extract, $input_file, $output_file); 

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


####################################################
####################################################
sub genbank_fasta_sequence_extract {

    my ($input, $output)=@_;
    open (TEMP_IN, $input);
    my @temp_in = <TEMP_IN>;
    my $first_line = shift @temp_in;
    my $po_seq = join ("", @temp_in);
    open(OUTPUT, ">$output");
    my $filename = basename $input;
    if ($first_line =~ /^>/){
        print OUTPUT ">$filename\n$po_seq";      
    }else{
        my $DDD;
        my $in = new Bio::SeqIO( -file => $input, -format => 'genbank' );
        my $out = new Bio::SeqIO( -file => $output, -format => 'fasta' );
        while (my $seqObject = $in->next_seq()){
            my $acc = $seqObject->accession;
            my $desc = $seqObject->desc;
            my $dna = $seqObject->seq();
            my $id = $seqObject->display_id();
            $DDD .=$dna;
       
        }
        $DDD=print_sequence_into_file($DDD, 60);  # print_sequence_into_file($dna, 70);  # 60 characters per line
        print OUTPUT ">$filename\n$DDD";
    }
    close OUTPUT;

}


# subroutine print_sequence_into_file
# print $sequence into $FILE with $length character per line
sub print_sequence_into_file{   
    my($sequence, $length) = @_;
    my $DNA;
    my $string;
    # Print sequence in lines of $length
    for ( my $pos = 0 ; $pos < length($sequence) ; $pos += $length ){   
        $DNA=substr($sequence, $pos, $length)."\n";
        $string.=$DNA;
        
    }
    return $string;
}
