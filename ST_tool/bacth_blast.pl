#perl /home/xiangyang/Desktop/7588_genome/MLSA-sequence/bacth_blast.pl $ST_workplace/path_gene $mlsa_folder $thread_number $ST_workplace

use strict;
use warnings;
use threads;

my $temp_protein_folder = $ARGV[0];
my $mlsa_folder = $ARGV[1];
my $thread_number = $ARGV[2];
my $out_folder = $ARGV[3];

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
my $blast_fold = "$out_folder/blast_folder";
mkdir $blast_fold;

foreach my $file_protein(@PATH_PROTEIN){
    
    my $input="$temp_protein_folder/$file_protein"; 
    push (@input, $input);  

    my $output="$blast_fold/$file_protein.blastout";
    push (@outfile, $output);

}

my $e_value = "0.0000000001";
my $formatdb = `which formatdb`;
$formatdb =~ s/\n//g;
my $blastn = `which blastn`;
$blastn =~ s/\n//g;

my $mlsa_database = "$out_folder/mlsa_database.fasta";
system("cat $mlsa_folder/* > $mlsa_database");   # database
system ("$formatdb -i $out_folder/mlsa_database.fasta -p F");

my @sub_function_parameters = (\@input, $mlsa_database, \@outfile, $e_value, $blastn);###edit

&bacth_blast($sub_file_number, $thread_number, "do_blastn", \@sub_function_parameters);

my $blastn_out = "$out_folder/blastn_out.txt";
#system ("cat $blast_fold/* > $blastn_out");
chdir $blast_fold;
system ("find . -name '*.blastout' | xargs cat > $blastn_out");  # For large number of files, directly using cat may causing error

#########################################################
#########################################################
sub bacth_blast {

    my ($work_number, $thread_number, $subfunction, $sub_function_parameters) = @_;
    my ($input_file, $db_file, $output_file, $e_value, $blastn) = @$sub_function_parameters;
 
    my @input = @$input_file;
    my @outfile = @$output_file;
    my $thread;
    my @threads;
    my $job_number=0;
    print "    Blastn_percent: ";
    while(){ 
        last if ($job_number>=$work_number);                         
        while(scalar(threads->list())<$thread_number) {     #set threadnumber；
            my $progress_record = int (($job_number/$work_number)*100);
            $job_number++;                                 
            my $input_file = $input[$job_number-1];
            my $output_file = $outfile[$job_number-1];;

            my $Blast_perform_progress = int (($job_number/$work_number)*100); 
            print "$Blast_perform_progress%","..." if ($job_number == 1 or ( ($Blast_perform_progress%10 ==0) && ($progress_record <$Blast_perform_progress)) );

            $threads[$job_number-1]=threads->new(\&$subfunction, $input_file, $db_file, $output_file, $e_value, $blastn);   
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
    print "done\n";
}



#########################################################
#########################################################
sub do_blastn {
    my ($input, $db_file, $outfile, $e_value, $blastn) =@_;
    system ("$blastn -outfmt '6 qseqid sseqid pident qlen slen length qcovhsp qcovs bitscore' -num_alignments 1  -evalue $e_value -num_threads 1 -query $input -db $db_file -out $outfile");  

}



#########################################################
#########################################################
sub do_blastP {
    my ($input, $db_file, $outfile, $e_value, $blastp) =@_;
    system("$blastp -outfmt '6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore' -evalue $e_value -max_hsps 1 -max_target_seqs 1 -num_threads 1 -query $input -db $db_file -out $outfile");  #-max_hsps -num_descriptions -num_alignments

}
