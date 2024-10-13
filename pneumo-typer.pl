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

Pneumo-Typer

=DESCRIPTION

    A comprehensive prediction and visualization of serotype and sequence type for streptococcus pneumoniae using assembled genomes. 

=USAGE

pneumo-typer.pl -d genbank_file_directory [options]

FOR EXAMPLE: 

perl pneumo-typer.pl -d Test_data -t 10  -p T


#######################################################################################################################################
=ARGUMENTS
=======================
    REQUIRED ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -d, --genbank_file_directory
           A directory containing files in GenBank format, FASTA format, or a combination of both.                           
    OPTIONAL ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -o, --output_directory
           An output directory holding all the generated files by pneumo-typer.pl. if this option is not set,  a directory named "pneumo-typer_workplace" will be created in the bin directory from where pneumo-typer.pl was invoked.
       -t, --multiple_threads
           Set thread number (Default: 1)
       -Ss, --skip_sequence_processing 
           Skip the process of sequence processing (STEP-1) (Default: F).  
       -hgc, --homologous_gene_cutoff
           Set E-value, Identify, Coverage (Query and Subject), Match_length (alignment length) cutoff in Blastn analysis (default: E-value=1e-5, Identify=70, Coverage=95, Match_length=100).
       -Sb, --skip_blastn
           Skip the process of doing blastn during serotype analysis.
       -p, --prodigal_annotation
           Annotate all genomes using prodigal (Default: T). 
       -m, --mlst
           Perform mlst analysis (Default: T). 
       -c, --cgmlst
           Perform cgmlst analysis. It needs about 3 mins for one genome (Default: F).
       -Rh, --recreate_heatmap                             
           Re-create the heatmap of cps gene distribution in genomes (Default: F). At this step, users can add a parameter "phylogenetic_tree" or "strain_reorder_file". 
       -Rf, --recreate_figure
           Re-create the figure of the genetic organization of cps gene cluster for genomes (Default: F). At this step, users can add a parameter "phylogenetic_tree" or "strain_reorder_file".
       -tree, --phylogenetic_tree
           A Newick format tree file is used by Pneumo-Typer to automatically associate the genomes with their phylogeny. Meanwhile, Pneumo-Typer will output a file named "temp_strain_reorder_file-svg.txt", which contains the order information of genomes in the tree from up to down. It should be noted that all node names in the provided tree must completely match the input file names of all genomes.
       -srf, --strain_reorder_file
           A two-column tab-delimited text file is used to sort genomes from up to down according to users' requirements. Each row must consist of a strain name followed by the numerical order that is used for sorting genomes. It should be noted that all strain names must completely match the input file names of all genomes.
       -Ts, --test
           Run pneumo-typer using Test_data as input to check whether Pneumo-Typer is installed successfully (Default: F).
       -V, --version
           The version of Pneumo-Typer.
       -h, --help
           Show this message.

=AUTHOR

Dr. Xiangyang Li (E-mail: lixiangyang\@fudan.edu.cn), Kaili University; Bacterial Genome Data mining & Bioinformatic Analysis (www.microbialgenomic.cn/).

=COPYRIGHT

Copyright 2024, Xiangyang Li. All Rights Reserved.

USAGE


my %options = (
    'genbank_file_directory'                       => undef,   
    'workplace_directory'                          => undef,  
    'multiple_threads'                             => "3",
    'skip_sequence_processing'                     => "F",
    'homologous_gene_cutoff'                       => "1e-5,70,95,30", #array to set blast parse cutoff: E-value, Identify, Coverage, Match_length
    'skip_blastn'                                  => "F",
    'prodigal_annotation'                          => "T",
    'mlst'                                         => "T",
    'cgmlst'                                       => "F",
    'recreate_heatmap'                             => "F",
    'recreate_figure'                              => "F",
    'phylogenetic_tree'                            => undef,
    'strain_reorder_file'                          => undef,
    'test'                                         => "F",
    'version'                                      => undef,
    'help'                                         => undef

);
 

GetOptions(
    'd|genbank_file_directory=s'                     => \$options{genbank_file_directory},    
    'o|workplace_directory=s'                        => \$options{workplace_directory},
    't|multiple_threads=i'                           => \$options{multiple_threads}, 
    'Ss|skip_sequence_processing=s'                  => \$options{skip_sequence_processing},       
    'hgc|homologous_gene_cutoff=s'                   => \$options{homologous_gene_cutoff}, 
    'Sb|skip_blastn=s'                               => \$options{skip_blastn},
    'p|prodigal_annotation=s'                        => \$options{prodigal_annotation},
    'm|mlst=s'                                       => \$options{mlst},
    'c|cgmlst=s'                                     => \$options{cgmlst},
    'Rh|recreate_heatmap=s'                          => \$options{recreate_heatmap},
    'Rf|recreate_figure=s'                           => \$options{recreate_figure},
    'tree|phylogenetic_tree=s'                       => \$options{phylogenetic_tree},
    'srf|strain_reorder_file=s'                      => \$options{strain_reorder_file},
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

my $tblastn        = `which tblastn`;  # "/usr/bin/tblastn"; mmseqs
$tblastn =~ s/\n//g;

my $blastn        = `which blastn`;  # "/usr/bin/blastn";
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
    $workplace = "$home_directory/pneumo-typer_workplace";
    $workplace =~ s/\/\//\//g;
    mkdir $workplace;
}

my $path_fa1 = $workplace."/Whole_fasta1"; #fasta sequences are extracted as contigs for each genome
mkdir $path_fa1;

my $path_fa2 = $workplace."/Whole_fasta2"; #fasta sequences are concatenated for each genome
mkdir $path_fa2;

my $path_protein = $workplace."/Whole_protein";
mkdir $path_protein;

my $path_gene = $workplace."/Whole_gene";
mkdir $path_gene;

my $directory_TFT = $workplace."/Whole_TFT"; 
mkdir $directory_TFT;

$options{skip_sequence_processing} = "T" if $options{skip_blastn} eq "T";

my $map_cmd = $workplace."/map_cmd.txt"; #store command to creat heatmap and figure

if ( ($options{recreate_heatmap} eq "F")  && ($options{recreate_figure} eq "F") ){

    system("rm -rf $workplace/result_statistics") if -d "$workplace/result_statistics"; #clear folder "result_statistics"

    if ($options{skip_sequence_processing} eq "F"){
        print "STEP-1: Dealing with genomes extract genome sequence, gene sequences and gene feature table (TFT); annotate genome which has no annotation information using prodigal>\n";

        &PTyper::batch_genomenucleotide_extract_run ($path_genbank, $path_fa1, $path_fa2, $thread_number);

        &PTyper::batch_genenucleotide_TFT_extract_run ($path_genbank, $directory_TFT, $path_protein, $path_gene, $thread_number, $prodigal_annotation, $path_fa1) if $prodigal_annotation eq "F";

        &PTyper::prodigal_bacth_run ($path_fa1, $path_protein, $path_gene, $directory_TFT, $thread_number, $prodigal_annotation) if $prodigal_annotation eq "T";
    }


    print "\nSTEP-2: Determining the sequence type (ST/cgST)\n" if ( ($options{mlst} eq "T") or ($options{cgmlst} eq "T") );

    if ($options{mlst} eq "T") {
        print "  \nSTEP-2.1: MLST analysis\n";
        system ("perl $home_directory/ST_tool/ST_profile.pl $path_genbank $path_gene $home_directory/ST_tool/database/mlst $thread_number $workplace/ST_workplace");
    }


    if ($options{cgmlst} eq "T") {
        print "  \nSTEP-2.2: cgMLST analysis\n";
    system ("perl $home_directory/ST_tool/ST_profile_cg.pl $path_genbank $path_fa2 $home_directory/ST_tool/database/cgmlst $thread_number $workplace/ST_workplace");
    }

 
    print "\nSTEP-3: Predicting serotype\n";
    
    
    my $blastp_out_dir = "$workplace/blastp_out_dir";
    mkdir $blastp_out_dir;
    
    
    my $blastn_out_dir = "$workplace/blastn_out_dir";
    mkdir $blastn_out_dir;   

    &PTyper::bacth_blast_best_run ($workplace, $path_protein, $blastp_out_dir, $path_gene, $blastn_out_dir, $db_file, $tblastn, $blastn, $makeblastdb, $thread_number, $options{homologous_gene_cutoff}) if $options{skip_blastn} eq "F";
#sed -i 's/.*\tGCA_001131025.1_112\t.*/10\t1000\tCDS\tGCA_001131025.1_112\tCKZO01000001/g' file
    print "    Process data and obtain serotype...";

    my $all_vs_all_cluster = "$workplace/all_vs_all.cluster";
    chdir $blastn_out_dir;
    system ("find . -type f -name '*.parsedout' | xargs cat > $all_vs_all_cluster") if $options{skip_blastn} eq "F";  # For large number of files, directly using cat may causing error

    my $filter_pseudogene = "$home_directory/DATABASE/filter_pseudogene.txt";
    &PTyper::recheck_pseudogene_len($all_vs_all_cluster, $filter_pseudogene, $options{homologous_gene_cutoff});

    ############################################################
    ############################################################
    my $map_table = "$workplace/CPS_Mapping_table.txt";
    &PTyper::SARG_map_transformation($db_file, $mappingList, $map_table);

    my $cluster_mapping;
    $cluster_mapping = "$workplace/CPS_cluster_mapping.result";
    my $serotype_negtive_strain = "$workplace/Serotype.no.detected.out";
    &PTyper::mapping ($map_table, $all_vs_all_cluster, $cluster_mapping, $path_genbank, $serotype_negtive_strain);
    ### if $cluster_mapping is blank, quit due to no cps gene is detected in any analyzed genomes
    my @feature = stat ($cluster_mapping);
    die "\nWarning: no cps gene is detected in any analyzed genomes!!!\n" if $feature[7] == 0;
    ############################################################
    ############################################################

    &PTyper::result_statistics($workplace, $directory_TFT, $cluster_mapping);

    system ("perl $home_directory/script/serotype_correct.pl $home_directory/DATABASE/serotype_matrix.txt $workplace/result_statistics/Statistics_OUT/CPS_LocusTag_statistics $workplace/Serotype.out");

    print "done\n"; #serotype prediction finished

    print "\nSTEP-4: Output sequence type and serotype results\n";

    system ("perl $home_directory/script/ser_join_ST.pl $workplace/Serotype.out $workplace/ST_out.txt $workplace/cgST_out.txt $workplace/Serotype_ST.out");



    ### Visualizing the result
    ##########################

    open (MAP_CMD, ">$map_cmd");

    print "\nSTEP-5: Heatmaping the cps gene distribution in genomes\n";
    my $cmd_1 ="perl $home_directory/script/heatmap.pl -dir $workplace/result_statistics/tbl_heatmap_class -left 20 -scale 4 -label T -dis 9 -w 4 -l 0 -right 50 -cf $workplace/result_statistics/Statistics_OUT/classification_CPS -e $workplace/Serotype_ST.out -o $workplace";
    system ("$cmd_1"); 

    my $cmd_2 ="perl $home_directory/script/heatmap.pl -dir $workplace/result_statistics/tbl_heatmap_gene -left 20 -scale 4 -label T -dis 9 -w 4 -l 0 -right 50 -cf $workplace/result_statistics/Statistics_OUT/classification_gene -e $workplace/Serotype_ST.out -o $workplace";
    system ("$cmd_2"); 

    print "\nSTEP-6: Visualizing the cps gene cluster for genomes\n";

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

    print MAP_CMD "heatmap_class:\n$cmd_1\n\nheatmap_gene:\n$cmd_2\n\ncps gene cluster:\n$cmd_3\n";
    close MAP_CMD;

}


my ($remap_cmd1, $remap_cmd2, $remap_cmd3) = &PTyper::obtain_map_cmd ($map_cmd, \%options);

#remap heatmap
if ($options{recreate_heatmap} eq "T"){

    print "\n    Re-create the heatmap of cps gene distribution in genomes (STEP-5)\n\n";
    system ("$remap_cmd1");
    system ("$remap_cmd2");

}

#remap figure
if ($options{recreate_figure} eq "T"){

    print "\n    Re-create the figure of the genetic organlization of cps gene cluster for genomes (STEP-6)\n\n";
    system ("$remap_cmd3");

}



my $finish_time = localtime;
print "\n$finish_time: done!\n\n";


