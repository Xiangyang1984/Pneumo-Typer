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
my $makeblastdb = `which makeblastdb`;
$makeblastdb =~ s/\n//g;
my $blastn = `which blastn`;
$blastn =~ s/\n//g;

my $mlsa_database = "$out_folder/mlsa_database.fasta";
system("cat $mlsa_folder/* > $mlsa_database");   # database
system ("$makeblastdb -in $out_folder/mlsa_database.fasta -dbtype nucl > $out_folder/temp.txt");

&bacth_blast_run(\@input, $mlsa_database, \@outfile, $e_value, $blastn, $thread_number);

my $blastn_out = "$out_folder/blastn_out.txt";
#system ("cat $blast_fold/* > $blastn_out");
chdir $blast_fold;
system ("find . -type f -name '*.blastout' | xargs cat > $blastn_out");  # For large number of files, directly using cat may causing error


#########################################################
#########################################################
sub bacth_blast_run {

    my ($input_file, $db_file, $output_file, $e_value, $blastn, $thread_number) = @_;
 
    my @input = @$input_file;
    my @output = @$output_file;

    my $sub_file_number =  scalar @input;
    my $job_number=0;
    print "    Blastn_percent: ";

    use Parallel::Runner;
    my $runner = Parallel::Runner->new($thread_number);

    foreach (@input){
        my $progress_record = int (($job_number/$sub_file_number)*100);
        $job_number++;
        my $Blast_perform_progress = int (($job_number/$sub_file_number)*100);
        print "$Blast_perform_progress%","..." if ($job_number == 1 or ( ($Blast_perform_progress%10 ==0) && ($progress_record <$Blast_perform_progress)) );

        my $inputfile = $input[$job_number-1];
        my $outputfile = $output[$job_number-1];  

        $runner->run(
            sub{
                my $cmd = "$blastn -outfmt '6 qseqid sseqid pident qlen slen length qcovhsp qcovs bitscore' -num_alignments 1  -evalue $e_value -num_threads 1 -query $inputfile -db $db_file -out $outputfile";
                system ("$cmd 2>/dev/null");
            }
        )#run end
    }
    $runner->finish;
    print "done\n";
}

