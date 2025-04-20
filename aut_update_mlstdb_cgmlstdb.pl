use strict;
use warnings;
use FindBin;
use Getopt::Long;
use lib "$FindBin::Bin/lib";
use POSIX;

my $usage = <<USAGE; 

# Usage: perl update_mlstdb_cgmlstdb.pl [-m T] [-c T] [-t parallel_number]

####################
=ARGUMENTS
=======================                          
    OPTIONAL ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~ 
       -m, --mlst
           Update mlst dataset (Default: T). 
       -c, --cgmlst
           Update cgmlst dataset (Default: T).
       -t, --multiple_threads
           download cgmlst loci using multiple threads (Default: 4)
       -h, --help
           Show this message.
=======================
=AUTHOR

Dr. Xiangyang Li (E-mail: lixiangyang\@fudan.edu.cn), Kaili University; Bacterial Genome Data mining & Bioinformatic Analysis (https://www.microbialgenomic.cn/).

=COPYRIGHT

Copyright 202e, Xiangyang Li. All Rights Reserved.

USAGE


my %options = (
    'mlst'                       => "T",   
    'cgmlst'                     => "T",  
    'multiple_threads'           => "4",
    'help'                       => undef
);

GetOptions(
    'm|mlst=s'                   => \$options{mlst},    
    'c|cgmlst=s'                 => \$options{cgmlst},
    't|multiple_threads=i'       => \$options{multiple_threads},
    'h|help'                     => \$options{help}
);

if ( defined( $options{help} ) ) {
    print $usage;
    exit(0);
}


my $home_dir = $FindBin::Bin;
my $log = "$home_dir/ST_tool/database/Update_time.txt";
open (LOG, ">$log");

my $mlst_list = "$home_dir/ST_tool/database/mlst_list.txt";
my $cgmlst_list = "$home_dir/ST_tool/database/cgmlst_list.txt";

# ***RESTful_API_download_pubmlst.pl needs token which is located in "$home_dir/ST_tool"
chdir "$home_dir/ST_tool";

my $nowtime;
### mlst dataset
if ($options{mlst} eq "T"){
    $nowtime = &get_log_time;
    print "[$nowtime] Update mlst dataset: \n";
    # download MLST gene sequences (7 genes)
    $nowtime = &get_log_time;
    print "[$nowtime] start to download mlst gene sequences\n";
    system("rm $home_dir/ST_tool/database/mlst/MLSA_senven_loci/*.fas");
    open(MLST, $mlst_list);
    while(<MLST>){
        chomp;
	#my $downsite = "https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/loci/".$_."/alleles_fasta";
	#system("wget -q --no-check-certificate -O $home_dir/ST_tool/database/mlst/MLSA_senven_loci/$_.fas -c $downsite");
	
	#OAuth credentials for download full database from pubMLST (https://bigsdb.readthedocs.io/en/latest/rest.html#api-oauth)
	my $downsite_loci = "loci/".$_."/alleles_fasta";
	system("perl $home_dir/ST_tool/RESTful_API_download_pubmlst.pl --route $downsite_loci --output $home_dir/ST_tool/database/mlst/MLSA_senven_loci/$_.fas");
        
    }
    
    # download MLST profile
    $nowtime = &get_log_time;
    print "[$nowtime] strat to download mlst profile\n";
    system("rm $home_dir/ST_tool/database/mlst/MLST_profiles");
    #system("wget -q --no-check-certificate -O $home_dir/ST_tool/database/mlst/MLST_profiles -c https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/schemes/1/profiles_csv");
    
    #OAuth credentials for download full database from pubMLST (https://bigsdb.readthedocs.io/en/latest/rest.html#api-oauth)
    my $website_profiles = "schemes/1/profiles_csv";
    system("perl $home_dir/ST_tool/RESTful_API_download_pubmlst.pl --route $website_profiles --output $home_dir/ST_tool/database/mlst/MLST_profiles");
}
$nowtime = &get_log_time;
print "[$nowtime] Finish mlst dataset update\n\n";
print LOG "mlst dataset update at $nowtime\n";
close MLST;



### cgmlst dataset
if ($options{cgmlst} eq "T"){
    $nowtime = &get_log_time;
    print "[$nowtime] Update cgmlst dataset: \n";
    # download cgMLST gene sequences (1222 genes)
    $nowtime = &get_log_time;
    print "[$nowtime] start to download cgmlst gene sequences using multiple threads\n";
    my $tmp_cgmlst = "$home_dir/ST_tool/database/tmp_cgmlst";
    mkdir $tmp_cgmlst;

    open(CGMLST, $cgmlst_list);
    my @cgmlst = <CGMLST>;
    my $job_number=0;
    my $sub_file_number = scalar @cgmlst;
    
    use Para::Runner;
    my $runner = Para::Runner->new($options{multiple_threads});
    print "    ";
    foreach(@cgmlst){
        chomp;
        my $progress_record = int (($job_number/$sub_file_number)*100);
        $job_number++;
        my $cgmlst_update_progress = int (($job_number/$sub_file_number)*100); 
        print "$cgmlst_update_progress%","..." if ($job_number == 1 or ( ($cgmlst_update_progress%10 ==0) && ($progress_record <$cgmlst_update_progress)) );
        $runner->run(
	    sub{
	        #my $downsite = "https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/loci/".$_."/alleles_fasta";
		#system("wget -q --no-check-certificate -O $tmp_cgmlst/$_.fas -c $downsite");
		
		#OAuth credentials for download full database from pubMLST (https://bigsdb.readthedocs.io/en/latest/rest.html#api-oauth)
                my $downsite_cgloci = "loci/".$_."/alleles_fasta";
                system("perl $home_dir/ST_tool/RESTful_API_download_pubmlst.pl --route $downsite_cgloci --output $tmp_cgmlst/$_.fas");
            }
        )#run end    
    }
    $runner->finish;
    print "done\n";
    
    system("mv $tmp_cgmlst/SPNE00213.fas $home_dir/ST_tool/database/cgmlst/cgMLSA_loci");
    chdir $tmp_cgmlst;
    system("find . -type f -name '*.fas' | xargs cat > $home_dir/ST_tool/database/cgmlst/cgMLSA_loci/cgMLSA_loci.fas");
    system("rm -r $tmp_cgmlst");

    
    # download cgMLST profile
    $nowtime = &get_log_time;
    print "[$nowtime] strat to download cgmlst profile\n";
    system("rm $home_dir/ST_tool/database/cgmlst/cgMLST_profiles");
    #system("wget -q --no-check-certificate -O $home_dir/ST_tool/database/cgmlst/cgMLST_profiles -c https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/schemes/2/profiles_csv");
    #OAuth credentials for download full database from pubMLST (https://bigsdb.readthedocs.io/en/latest/rest.html#api-oauth)
    my $website_cgprofiles = "schemes/2/profiles_csv";
    system("perl $home_dir/ST_tool/RESTful_API_download_pubmlst.pl --route $website_cgprofiles --output $home_dir/ST_tool/database/cgmlst/cgMLST_profiles");

    $nowtime = &get_log_time;
    print "[$nowtime] Finish cgmlst dataset update\n\n";
    print LOG "cgmlst dataset update at $nowtime\n";
}
close MLST;



sub get_log_time {
    my $time = strftime("%Y-%m-%d %H:%M:%S", localtime);
    return ($time);
}
