
#perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/ST_tool/bacth_blat2.pl $path_fa $mlsa_database $mlsa_genename $thread_number $out_folder

use strict;
use warnings;
use threads;
use File::Basename qw<basename dirname>;
use vars qw(@loci_name);
#my $now_time = localtime;
#print "\n$now_time: start...\n\n";

my $path_fa = $ARGV[0];
my $mlsa_database = $ARGV[1];
my $mlsa_genename = $ARGV[2];
my $thread_number = $ARGV[3];
my $out_folder = $ARGV[4];
mkdir $out_folder;

opendir PATH_FA, $path_fa or die "could not open $path_fa";
my @PATH_FA = readdir PATH_FA; 
@PATH_FA =grep ($_!~/^\./ ,@PATH_FA);  #delete hidden file . ..
closedir PATH_FA;

my @loci_name;
open (NA, $mlsa_genename);
while(<NA>){
    chomp;
    @loci_name = split /\t/, $_;
}
close NA;
#############################################################################
#### do blat
#############################################################################



my $sub_file_number = scalar @PATH_FA;
if ($sub_file_number <= $thread_number) {
    $thread_number = $sub_file_number;
}

my $blast_fold = "$out_folder/cgblast_folder";
mkdir $blast_fold;
chdir $blast_fold;

my $blat_parse_dir = "$out_folder/cgblat_parse_dir";
mkdir $blat_parse_dir;


my $blat = `which blat`;
$blat =~ s/\n//g;


&bacth_blat_run($thread_number, $path_fa, \@PATH_FA, $mlsa_database, $blat, $blast_fold, $blat_parse_dir);


#########################################################
#########################################################
sub bacth_blat_run {

    my ($thread_number, $path_fa, $input_dir, $db_dir, $blat, $blat_fold, $blat_parse_fold) = @_;
 
    my @input = @$input_dir;
    my $sub_file_number =  scalar @input;
    my $job_number=0;
    print "    Blat_percent: ";
    use Parallel::Runner;
    my $runner = Parallel::Runner->new($thread_number);
    foreach (@input){
        my $input = "$path_fa/$_";  
        my @a = ("SPNE00213.fas", "cgMLSA_loci.fas");
        my $progress_record = int (($job_number/$sub_file_number)*100);
        $job_number++;
        my $Blat_percent = int (($job_number/$sub_file_number)*100); 
        print "$Blat_percent%","..." if ($job_number == 1 or ( ($Blat_percent%10 ==0) && ($progress_record <$Blat_percent)) );

        $runner->run(
            sub{
                    my $blatout_temp_dir = $blat_fold."/". (basename $input);
                    mkdir $blatout_temp_dir;

                    my $blatparse_out = $blat_parse_fold."/".basename $input.".cgblastn_out";
                    open (PAR_OUT, ">$blatparse_out");

                    
                    foreach my $db (@a){
                        my $db_file = "$db_dir/$db";
                        my $output = $blatout_temp_dir."/$db.blastout";
                        system ("$blat $input $db_file $output -fastMap -noHead > temp.txt") if $db ne "SPNE00213.fas"; 
                        system ("$blat $input $db_file $output -noHead > temp.txt") if $db eq "SPNE00213.fas"; #-fastMap could not be used due to the sequence length for "SPNE00213.fas" is >5000
                    } 
		   
                    my $outputfinal = $blatout_temp_dir."/final.blastout";
                    system ("cat $blatout_temp_dir/* > $outputfinal");

                    #################################    
                    open (OT, $outputfinal);
                    my @gene_loci;
                    while (<OT>){    
                        chomp;
                        $_ =~ s/ //g;
                        my @array = split /\t/, $_;
                        my $sum = $array[1]+$array[2]+$array[3]+$array[4]+$array[5]+$array[6]+$array[7];
                        if ( ($array[0] eq $array[10]) && ($sum == 0) ){
                            print PAR_OUT "T\t$array[13]\t$array[9]\t$array[15]\t$array[16]\t$array[10]\t$array[0]\n";
                            my $na = $array[9];
                            $na =~ s/_.*//g; 
                            push @gene_loci, $na;
                
                        }
                    }
                    close OT;

                    my $strain = basename $input;
                    foreach my $e(@loci_name){             
                        print PAR_OUT "F\t$strain\t$e\tNo hit\n" if !grep{$e eq $_}@gene_loci; 
                    } 
                    close PAR_OUT;
                    #################################
  
                    system ("rm -rf $blatout_temp_dir");
                }
        )#run end
    }
    $runner->finish;
    print "done\n";
}


