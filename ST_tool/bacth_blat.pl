#perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/bacth_blat.pl $ST_workplace/path_gene $mlsa_folder $thread_number $ST_workplace

#perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/bacth_blat.pl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/database/cgmlst/cgMLSA_loci /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/TEST_OUT/whole_fasta/GCA_900056415.2_13414_2_92_genomic.gbff 5 /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/TEST_OUT

use strict;
use warnings;
use threads;
use File::Basename qw<basename dirname>;
#my $now_time = localtime;
#print "\n$now_time: start...\n\n";

my $temp_protein_folder = $ARGV[0];
my $mlsa_database = $ARGV[1];
my $thread_number = $ARGV[2];
my $out_folder = $ARGV[3];
mkdir $out_folder;
#############################################################################
#### do blast
#############################################################################

opendir PATH_PROTEIN, $temp_protein_folder or die "could not open $temp_protein_folder";
my @PATH_PROTEIN = readdir PATH_PROTEIN; 
@PATH_PROTEIN =grep ($_!~/^\./ ,@PATH_PROTEIN);  #delete hidden file . ..
closedir PATH_PROTEIN;

my $sub_file_number = scalar @PATH_PROTEIN;
if ($sub_file_number <= $thread_number) {
    $thread_number = $sub_file_number;
}

my @input;
my @outfile;
my $blast_fold = "$out_folder/cgblast_folder";
mkdir $blast_fold;
chdir $blast_fold;
foreach my $file_protein(sort @PATH_PROTEIN){
    
    my $input="$temp_protein_folder/$file_protein"; 
    push (@input, $input);  

    my $output="$blast_fold/$file_protein.blastout";
    push (@outfile, $output);

}

my $e_value = 0;
my $blat = `which blat`;
$blat =~ s/\n//g;

#my $mlsa_database = "$out_folder/mlsa_database.fasta";
#system("cat $mlsa_folder/* > $mlsa_database");   # database

my $blat_out = "$out_folder/cgblastn_out.txt";
my @sub_function_parameters = (\@input, $mlsa_database, \@outfile, $e_value, $blat, $blat_out);###edit


&bacth_blat($sub_file_number, $thread_number, "do_blat", \@sub_function_parameters);


system ("rm -rf $blast_fold");
#chdir $blast_fold;
#system ("find . -name '*.blastout' | xargs cat > $blat_out");  # For large number of files, directly using cat may causing error
#my $finish_time = localtime;
#print "\n$finish_time: done!\n\n";
#########################################################
#########################################################
sub bacth_blat {

    my ($work_number, $thread_number, $subfunction, $sub_function_parameters) = @_;
    my ($input_file, $db_file, $output_file, $e_value, $blat, $blatparse_out) = @$sub_function_parameters;
 
    my @input = @$input_file;
    my @outfile = @$output_file;
    my $thread;
    my @threads;
    my $job_number=0;
    #print "    Blat_percent: ";
    while(){ 
        last if ($job_number>=$work_number);                         
        while(scalar(threads->list())<$thread_number) {     #set threadnumber；
            #my $progress_record = int (($job_number/$work_number)*100);
            $job_number++;                                 
            my $input_file = $input[$job_number-1];
            my $output_file = $outfile[$job_number-1];;

            #my $Blast_perform_progress = int (($job_number/$work_number)*100); 
            #print "$Blast_perform_progress%","..." if ($job_number == 1 or ( ($Blast_perform_progress%10 ==0) && ($progress_record <$Blast_perform_progress)) );

            $threads[$job_number-1]=threads->new(\&$subfunction, $input_file, $db_file, $output_file, $e_value, $blat, $blatparse_out);   
            last if ($job_number eq $work_number);    
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
    #print "done\n";
}



#########################################################
#########################################################
sub do_blat {
    my ($input, $db_file, $outfile, $e_value, $blat, $blatparse_out) =@_;
    system ("$blat $db_file $input $outfile -noHead > temp.txt");  
    open (OT, $outfile);
    open (PAR_OUT, ">>$blatparse_out");
    my %hash;
    my $check=0;
    while (<OT>){    
        chomp;
        $_ =~ s/ //g;
        my @array = split /\t/, $_;
        my $sum = $array[1]+$array[2]+$array[3]+$array[4]+$array[5]+$array[6]+$array[7];
        if ( ($array[0] eq $array[10]) && ($sum == 0) ){
            print PAR_OUT "T\t$array[13]\t$array[9]\t$array[15]\t$array[16]\t$array[10]\t$array[0]\n";
            $check=1;
        }else{
            my $iden = $array[0]/$array[10];
            $hash{$iden} = "$array[13]\t$array[9]\t$array[15]\t$array[16]\t$array[10]\t$array[0]";
        }
    }
    close OT;
    if ($check==0){
        my @a = sort{$a<=>$b} keys %hash;
        print PAR_OUT "F\t", $hash{$a[-1]}, "\n" if defined $a[-1]; 
        my $gene = basename $input;
        my $strain = basename $db_file;
        print PAR_OUT "F\t$strain\t$gene\tNo subject\n" if !defined $a[-1]; 
    } 
    close PAR_OUT;
}


