#!/usr/bin/perl -w

use strict;
use warnings;
use threads;
use FindBin;
use lib "$FindBin::Bin/lib";
use Getopt::Long;
use File::Spec;
use PTyper;
#use vars qw(%options);
my $usage = <<USAGE; 

=NAME

pneumo-typer.pl

=DESCRIPTION

    A comprehensive prediction and visualization of serotype and sequence type for streptococcus pneumoniae using assembled genomes. 

=USAGE

pneumo-typer.pl -d genbank_file_directory [options]

FOR EXAMPLE: 

perl [absolute path to] pneumo-typer.pl -d Test_data -t 10  -p T -srf /Users/zilinyang/Desktop/桌面材料/李向阳-论文投稿汇总/Serotype_ST_manuscript/bioconda/pneumo-typer-v1.0.1/order_sFig.txt -e /Users/zilinyang/Desktop/桌面材料/李向阳-论文投稿汇总/Serotype_ST_manuscript/bioconda/pneumo-typer-v1.0.1/expored_list.txt


#######################################################################################################################################
=ARGUMENTS
=======================
    REQUIRED ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -d, --genbank_file_directory
           A directory containing files as GenBank format, FASTA format, or a combination of both.                           
    OPTIONAL ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -o, --output_directory
           An output directory holding all the generated files by pneumo-typer.pl. if this option is not set,  a directory named "ipneumo-pyper_workplace" will be created in the bin directory from where pneumo-typer.pl was invoked.
       -m, --multiple_threads
           Set thread number (Default: 1)
       -jb, --start_at_blast 
           Jump to a local blastn analysis, and Skips sequencing extraction (Default: F).  
       -hgc, --homologous_gene_cutoff
           Set E-value, Identify, Coverage (Query and Subject), Match_length (alignment length) cutoff in Blastn analysis (default: E-value=1e-5, Identify=70, Coverage=95, Match_length=100).
       -js, --jump_serotype
           After blastn analysis, users can jump the process of doing blastn file to perform Serotype analysis again.
       -p, --prodigal_annotation
           Annotate all genomes using prodigal. 
       -m, --mlst
           Perform mlst analysis (Default: T). 
       -c, --cgmlst
           Perform cgmlst analysis. It need >10 mins for one genome (Default: F).
       -Ts, --test
           Run pneumo-typer using Test_data as input to check whether pneumo-typer is installed successfully (Default: F).
       -V, --version
           The version of pneumo-typer.
       -h, --help
           Show this message.

=AUTHOR

Dr. Xiangyang Li (E-mail: lixiangyang\@fudan.edu.cn), Kaili University; Bacterial Genome Data mining & Bioinformatic Analysis (www.microbialgenomic.cn/).

=COPYRIGHT

Copyright 2023, Xiangyang Li. All Rights Reserved.

USAGE


my %options = (
    'genbank_file_directory'                       => undef,   
    'workplace_directory'                          => undef,  
    'multiple_threads'                             => "3",
    'start_at_blast'                               => "F",
    'homologous_gene_cutoff'                       => "1e-5,70,95,100", #array to set blast parse cutoff: E-value, Identify, Coverage, Match_length
    'jump_serotype'                                => "F",
    'prodigal_annotation'                          => "F",
    'mlst'                                         => "T",
    'cgmlst'                                       => "F",
    'test'                                         => "F",
    'version'                                      => undef,
    'help'                                         => undef

);
 

GetOptions(
    'd|genbank_file_directory=s'                     => \$options{genbank_file_directory},    
    'o|workplace_directory=s'                        => \$options{workplace_directory},
    't|multiple_threads=i'                           => \$options{multiple_threads}, 
    'jb|start_at_blast=s'                            => \$options{start_at_blast},       
    'hgc|homologous_gene_cutoff=s'                   => \$options{homologous_gene_cutoff}, 
    'js|jump_serotype=s'                             => \$options{jump_serotype},
    'p|prodigal_annotation=s'                        => \$options{prodigal_annotation},
    'm|mlst=s'                                       => \$options{mlst},
    'c|cgmlst=s'                                     => \$options{cgmlst},
    'Ts|test=s'                                      => \$options{test},
    'V|version'                                      => \$options{version},
    'h|help'                                         => \$options{help}

);

if ( defined( $options{help} ) ) {
    print $usage;
    exit(0);
}

if (defined $options{version}) {
    print "\nVersion: Pneumo-Typer v1.0.2\n\n";
    exit(1);  
}

my $now_time = localtime;
print "\n$now_time: pneumo-typer.pl start...\n\n";

my $home_directory = $FindBin::Bin;   # obtaining the home directory where Pneumo-Typer.pl located
######################################################################################################
# Set the absolute path for four programs
my $blastn        = `which blastn`;  # "/usr/bin/blastp";
$blastn =~ s/\n//g; 

my $makeblastdb   = `which makeblastdb`; #"/usr/bin/makeblastdb";
$makeblastdb =~ s/\n//g; 
######################################################################################################


#check for Pneumo-Typer.pl workplace options
if ($options{test} eq "T") {
    $options{genbank_file_directory} = "$home_directory/Test_data";
    &PTyper::check_tool($0, $options{genbank_file_directory});
    exit(2);  
}

#check for required options
if ( !( defined($options{genbank_file_directory}) ) ) {
    print $usage;
    exit(3);
}

my $path_genbank = File::Spec->rel2abs($options{genbank_file_directory});
my $db_file = "$home_directory/DATABASE/CPS_workplace/seq_file.fasta";
my $mappingList = "$home_directory/DATABASE/CPS_workplace/id_list";
my $thread_number = $options{multiple_threads};
my $prodigal_annotation = $options{prodigal_annotation};

my $workplace;
if ( defined( $options{workplace_directory} ) ) {
    $workplace = File::Spec->rel2abs($options{workplace_directory});
    mkdir $workplace;
}else {
    $workplace = "$home_directory/pneumo-pyper_workplace";
    $workplace =~ s/\/\//\//g;
    mkdir $workplace;
}

my $path_fa1 = $workplace."/Whole_fasta1"; #fasta sequences are extracted as contigs for each genome
mkdir $path_fa1;

my $path_fa2 = $workplace."/Whole_fasta2"; #fasta sequences are concatenated for each genome
mkdir $path_fa2;

my $path_gene = $workplace."/Whole_gene";
mkdir $path_gene;

my $directory_TFT = $workplace."/Whole_TFT"; 
mkdir $directory_TFT;


if ($options{start_at_blast} eq "F"){
    print "STEP-1: Dealing with genomes <extract genome sequence, gene seqeunces and gene feature table (TFT); annnotate genome which has no annotation infotmation using prodigal>\n";

    &PTyper::batch_genomenucleotide_extract_run ($path_genbank, $path_fa1, $path_fa2, $thread_number);

    &PTyper::batch_genenucleotide_TFT_extract_run ($path_genbank, $directory_TFT, $path_gene, $thread_number, $prodigal_annotation, $path_fa1) if $prodigal_annotation eq "F";

    &PTyper::prodigal_bacth_run ($path_fa1, $path_gene, $directory_TFT, $thread_number, $prodigal_annotation) if $prodigal_annotation eq "T";
}


print "\nSTEP-2: Determining the Sequence Type (ST)\n" if ( ($options{mlst} eq "T") or ($options{cgmlst} eq "T") );

if ($options{mlst} eq "T") {
    print "  \nSTEP-2.1: MLST analysis\n";
    system ("perl $home_directory/ST_tool/ST_profile.pl $path_genbank $path_gene $home_directory/ST_tool/database/mlst $thread_number $workplace/ST_workplace");
}


if ($options{cgmlst} eq "T") {
    print "  \nSTEP-2.2: cgMLST analysis\n";
    system ("perl $home_directory/ST_tool/ST_profile_cg.pl $path_genbank $path_fa2 $home_directory/ST_tool/database/cgmlst $thread_number $workplace/ST_workplace");
}

 
print "\nSTEP-3: Predicting serotype\n";
my ($e_value, $identify, $coverage, $match_length) = split /,/, $options{homologous_gene_cutoff};
my $best_output = "$workplace/BEST.BLASTOUT";
my $blastn_out_dir = "$workplace/blastn_out_dir";
mkdir $blastn_out_dir;

&PTyper::bacth_blast_best_run ($workplace, $path_gene, $blastn_out_dir, $db_file, $e_value, $blastn, $makeblastdb, $thread_number) if $options{jump_serotype} eq "F";
chdir $blastn_out_dir;
system ("find . -type f -name '*.blastout' | xargs cat > $best_output");  # For large number of files, directly using cat may causing error

my $all_vs_all_cluster = PTyper::blast_filter ($best_output, $e_value, $identify, $coverage, $match_length);


############################################################
############################################################
my $map_table = "$workplace/CPS_Mapping_table.txt";
&PTyper::SARG_map_transformation($db_file, $mappingList, $map_table);

my $cluster_mapping;
$cluster_mapping = "$workplace/CPS_cluster_mapping.result";
&PTyper::mapping ($map_table, $all_vs_all_cluster, $cluster_mapping);


############################################################
############################################################

&PTyper::result_statistics($workplace, $directory_TFT, $cluster_mapping);

system ("perl $home_directory/script/serotype_correct.pl $home_directory/DATABASE/serotype_matrix.txt $workplace/result_statistics/Statistics_OUT/CPS_LocusTag_statistics $workplace/Serotype.out");


print "\nSTEP-4: Output sequence type and serotype results\n";

system ("perl $home_directory/script/ser_join_ST.pl $workplace/Serotype.out $workplace/ST_out.txt $workplace/cgST_out.txt $workplace/Serotype_ST.out");

#system ("rm -rf $new_db_file*"); #delete all serotype database files


### Visualizing the result
##########################

print "\nSTEP-5: Heatmaping the cps gene distribution in genomes\n";
my $cmd_1 ="perl $home_directory/script/heatmap.pl -dir $workplace/result_statistics/tbl_heatmap_class -left 20 -scale 4 -label T -dis 9 -w 4 -l 0 -right 50 -cf $workplace/result_statistics/Statistics_OUT/classification_CPS -e $workplace/Serotype_ST.out -o $workplace";
system ("$cmd_1"); 

my $cmd_2 ="perl $home_directory/script/heatmap.pl -dir $workplace/result_statistics/tbl_heatmap_gene -left 20 -scale 4 -label T -dis 9 -w 4 -l 0 -right 50 -cf $workplace/result_statistics/Statistics_OUT/classification_gene -e $workplace/Serotype_ST.out -o $workplace";
system ("$cmd_2"); 

print "\nSTEP-6: Visualizing the cps gene cluster in each genome\n";

my $gcluster_workplace = "$workplace/cps_cluster_workplace";
mkdir $gcluster_workplace;
my $subtft = "$gcluster_workplace/Sub_TFT";
mkdir $subtft;
#obtain interested gene and sub TFT files, in which gene distance from (20 kb) cps genes were removed; 
#only genome having the number of cps gene > 1 is shown 
system ("perl $home_directory/script/interested_gene.pl $home_directory/DATABASE/gene.txt $workplace/Serotype.out $workplace/result_statistics/tbl_part $gcluster_workplace/interested_gene.txt $subtft 20000 1");

#obtain homologs cluster
my $homologs_cluster = "$gcluster_workplace/blast_homologs_cluster";
mkdir $homologs_cluster;
system ("perl $home_directory/script/mcr.pl $workplace/CPS_cluster_mapping.result $homologs_cluster/all_orthomcl.out");

my $cmd_3 = "perl $home_directory/script/cps_cluster.pl -dir $path_genbank -gene $gcluster_workplace/interested_gene.txt -m $thread_number -map T -o $gcluster_workplace -SVG T -n 40 -e $workplace/Serotype_ST.out";
system ("$cmd_3"); 

my $finish_time = localtime;
print "\n$finish_time: done!\n\n";


