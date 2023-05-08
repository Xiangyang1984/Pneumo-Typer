#!/usr/bin/perl -w

use strict;
use warnings;
use threads;
use FindBin;
use lib "$FindBin::Bin/lib";
use Getopt::Long;
use File::Basename qw<basename dirname>;
use File::Spec;
use Bio::SeqIO;

my $usage = <<USAGE; 

=NAME

Pneumo-Typer.pl

=DESCRIPTION

    Run this command to enble users to obtain predict serotype and determine sequence type (ST/cgST) for Streptococcus pneumoniae, depending on obtaining the potential cps genes by a local blastn analysis using multiple threads. 

=USAGE

Pneumo-Typer.pl -dir genbank_file_directory [options]

FOR EXAMPLE: 

perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.03/Pneumo-Typer.pl -d /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer_v1.02/test -t 10 -s T -c T -p T -b T

nohup perl /home/xiangyang/Mazhongrui/Streptococcus_pneumoniae/Pneumo-Typer_v1.03/Pneumo-Typer.pl -d /home/xiangyang/Mazhongrui/Streptococcus_pneumoniae/prokka_workplace/prokka_gbk_folder_plus72 -t 190 -s T -c T -m F -b T

#######################################################################################################################################
=ARGUMENTS
=======================
    REQUIRED ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -d, --genbank_file_directory
           A directory containing annotated genomes as Genbank format file. To avoid a mistake, genome names cannot use special character,
           such as space, equal. For large number of genomes, users are recommended to download using Aspera, a high-speed file transfer
           tool (https://downloads.asperasoft.com/).                           
    OPTIONAL ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -o, --output_directory
           An output directory holding all the generated files by Pneumo-Typer.pl. if this option is not set, interested_gene_generation.pl will create a directory named "interested_gene_workplace" in the bin directory from where AMRG_Anno.pl was invoked.
       -m, --multiple_threads
           Set thread number (Default: 1)
       -b, --start_at_blast 
           Jump to a local blastn analysis, and Skips sequencing extraction (Default: F).  
       -hgc, --homologous_gene_cutoff
           Set E-value, Identify, Coverage (Query and Subject), Match_length (alignment length) cutoff in Blastn analysis (default: E-value=1e-5, Identify=70, Coverage=90, Match_length=100).
       -s, --split_file
           Split the query file into serveral sub-query files for parallal works, and the spilt number equie to the used threads number.
       -p, --prodigal_annotation
           Annotate all genomes using prodigal. 
       -m, --mlst
           Perform mlst analysis (Default: T). 
       -c, --cgmlst
           Perform cgmlst analysis. It need >10 mins for one genome (Default: F).
       -h, --help
           Show this message.

=AUTHOR

Dr. Xiangyang Li (E-mail: lixiangyang\@fudan.edu.cn), Kaili University; Bacterial Genome Data mining & Bioinformatic Analysis (www.microbialgenomic.cn/).

=COPYRIGHT

Copyright 2023, Xiangyang Li. All Rights Reserved.

USAGE



my $home_directory = $FindBin::Bin;   # obtaining the home directory where Pneumo-Typer.pl located
######################################################################################################
# Set the absolute path for four programs
my $blastn        = `which blastn`;  # "/usr/bin/blastp";
$blastn =~ s/\n//g; 

my $makeblastdb   = `which formatdb`; #"/usr/bin/makeblastdb";
$makeblastdb =~ s/\n//g; 

my $interested_gene_name = "TEST";
######################################################################################################

my %options = (
    'genbank_file_directory'                       => undef,   
    'workplace_directory'                          => undef,  
    'multiple_threads'                             => "3",
    'start_at_blast'                               => "F",
    'homologous_gene_cutoff'                       => "1e-5,70,90,100", #array to set blast parse cutoff: E-value, Identify, Coverage, Match_length
    'split_file'                                   => "F",
    'prodigal_annotation'                          => "F",
    'mlst'                                         => "T",
    'cgmlst'                                       => "F",
    'help'                                         => undef

);
 

GetOptions(
    'd|genbank_file_directory=s'                     => \$options{genbank_file_directory},    
    'o|workplace_directory=s'                        => \$options{workplace_directory},
    't|multiple_threads=i'                           => \$options{multiple_threads}, 
    'b|start_at_blast=s'                             => \$options{start_at_blast},       
    'hgc|homologous_gene_cutoff=s'                   => \$options{homologous_gene_cutoff}, 
    's|split_file=s'                                 => \$options{split_file},
    'p|prodigal_annotation=s'                        => \$options{prodigal_annotation},
    'm|mlst=s'                                       => \$options{mlst},
    'c|cgmlst=s'                                     => \$options{cgmlst},
    'h|help'                                         => \$options{help}

);


my $now_time = localtime;
print "\n$now_time: Pneumo-Typer.pl start...\n\n";

my $path_genbank = File::Spec->rel2abs($options{genbank_file_directory});
my $db_file = "$home_directory/DATABASE/CPS_workplace/seq_file.fasta";
my $mappingList = "$home_directory/DATABASE/CPS_workplace/id_list";
my $thread_number = $options{multiple_threads};
my $output_title = "Qseqid\tSseqid\tBitscore\tE-value\tPidenty\tQ_coverage\tS_coverage\tMacth_length";

#check for Pneumo-Typer.pl workplace options
my $workplace;
if ( defined( $options{workplace_directory} ) ) {
    $workplace = File::Spec->rel2abs($options{workplace_directory});
    mkdir $workplace;
}else {

    $workplace = "$home_directory/Pneumo-Typer_workplace";
    $workplace =~ s/\/\//\//g;
    mkdir $workplace;
}

my $path_gene = $workplace."/Whole_gene";
mkdir $path_gene;

my $directory_TFT = $workplace."/Whole_TFT"; 
mkdir $directory_TFT;


my $new_db_file = $workplace."/".basename($db_file);
system ("cp $db_file $new_db_file");
system ("$makeblastdb -i $new_db_file -p F");


if ($options{start_at_blast} eq "F"){
    print "STEP-1: Dealing with genomes <extract gene seqeunces and gene feature table; annnotate genome which has no annotation infotmation using prodigal>\n";
    &batch_genbank_sequence_TFT_extract($path_genbank, $directory_TFT, $path_gene, $thread_number) if $options{prodigal_annotation} eq "F";
    &prodigal_bacth($path_genbank, $path_gene, $directory_TFT, $thread_number)  if $options{prodigal_annotation} eq "T";
}

print "\nSTEP-2: Determining the Sequence Type (ST)\n" if ( ($options{mlst} eq "T") or ($options{cgmlst} eq "T") );

if ($options{mlst} eq "T") {
    print "  \nSTEP-2.1: MLST analysis\n";
    system ("perl $home_directory/ST_tool/ST_profile.pl $path_genbank $path_gene $home_directory/ST_tool/database/mlst $thread_number $workplace/ST_workplace");
}

if ($options{cgmlst} eq "T") {
    print "  \nSTEP-2.2: cgMLST analysis\n";
    system ("perl $home_directory/ST_tool/ST_profile_cg.pl $path_genbank $home_directory/ST_tool/database/cgmlst $thread_number $workplace/ST_workplace");
}

 
print "\nSTEP-3: Predicting serotype (ST)\n";
my ($e_value, $identify, $coverage, $match_length) = split /,/, $options{homologous_gene_cutoff};
chdir $path_gene;

#system("cat * > $workplace/all.all.fasta");
system ("find . -name '*.fasta' | xargs cat > $workplace/all.all.fasta");  # For large number of files, directly using cat may causing error
#system("cat $path_gene/* > $workplace/all.all.fasta");
my $total_gene = "$workplace/all.all.fasta";
my $best_output = "$workplace/BEST.BLASTOUT";

if ($options{split_file} eq "T"){

    my $split_protein_path = "$workplace/split_diretory";
    mkdir $split_protein_path;
    &split_to_subfiles ($total_gene, $thread_number, $split_protein_path);
    system("rm -f $total_gene");


    opendir SPLIT_PROTEIN_PATH, $split_protein_path or die "could not open $split_protein_path";
    my @SPLIT_PROTEIN_PATH = readdir SPLIT_PROTEIN_PATH; 
    @SPLIT_PROTEIN_PATH =grep ($_!~/^\./ ,@SPLIT_PROTEIN_PATH);  #delete hidden file . ..
    closedir SPLIT_PROTEIN_PATH;
    my (@input, @outfile);
    foreach my $splitprotein(@SPLIT_PROTEIN_PATH){
 
        my $input="$split_protein_path/$splitprotein"; 
        push (@input,$input);  
        #@input = sort @input;
        my $output="$split_protein_path/$splitprotein.blastout";
        push (@outfile,$output);
        #@outfile = sort @outfile;

    }

    my $file_number = scalar @input;
    if ($file_number <= $options{multiple_threads}) {
        $thread_number = $file_number;
    }

    my @sub_parameters = (\@input, $new_db_file, \@outfile, $e_value, $blastn);
    &bacth_blast_best($file_number, $thread_number, "blastn_best", \@sub_parameters);
    chdir $split_protein_path;
    system ("find . -name '*.blastout' | xargs cat > $best_output");  # For large number of files, directly using cat may causing error
    system ("rm -rf $split_protein_path"); #delete split_diretory used for split blastn analysis
}else {

    &blastn_best ($total_gene, $new_db_file, $best_output, $e_value, $blastn);               # A vs A blast

}



my $all_vs_all_cluster = blast_filter ($best_output, $e_value, $identify, $coverage, $match_length);


############################################################
############################################################
my $map_table = dirname ($new_db_file)."/CPS_Mapping_table.txt";
SARG_map_transformation($new_db_file, $mappingList, $map_table);

my $cluster_mapping;
$cluster_mapping = "$workplace/CPS_cluster_mapping.result";
&mapping ($map_table, $all_vs_all_cluster, $cluster_mapping);

############################################################
############################################################

&result_statistics($workplace, $directory_TFT, $cluster_mapping);

system ("perl $home_directory/script/serotype_correct.pl $home_directory/DATABASE/serotype_matrix.txt $workplace/result_statistics/Statistics_OUT/CPS_LocusTag_statistics $workplace/Serotype.out");

print "\nSTEP-4: Output sequence type and serotype results\n";

system ("perl $home_directory/script/ser_join_ST.pl $workplace/Serotype.out $workplace/ST_out.txt $workplace/cgST_out.txt $workplace/Serotype_ST.out");

system ("rm -rf $new_db_file*"); #delete all serotype database files

# Visualizing the result
print "\nSTEP-5: Heatmaping the cps gene distribution in genomes\n";
system ("perl $home_directory/script/heatmap.pl -dir $workplace/result_statistics/tbl_heatmap_antibiotic -left 20 -scale 5 -label T -dis 8 -w 4 -l 0 -cf $workplace/result_statistics/Statistics_OUT/classification_CPS -o $workplace"); 

print "\nSTEP-6: Visualizing the cps gene cluster in each genome\n";

my $gcluster_workplace = "$workplace/Gcluster_workplace";
mkdir $gcluster_workplace;
my $subtft = "$gcluster_workplace/Sub_TFT";
mkdir $subtft;
#obtain interested gene and sub TFT files, in which gene distance from (20 kb) cps genes were removed
system ("perl $home_directory/script/interested_gene.pl $workplace/result_statistics/tbl_part $gcluster_workplace/interested_gene.txt $subtft 20000");

#obtain homologs cluster
my $homologs_cluster = "$gcluster_workplace/blast_homologs_cluster";
mkdir $homologs_cluster;
system ("perl $home_directory/script/mcr.pl $workplace/CPS_cluster_mapping.result $homologs_cluster/all_orthomcl.out");


system ("perl $home_directory/script/Gcluster.pl -dir $path_genbank -gene $gcluster_workplace/interested_gene.txt -m $thread_number -map T -o $gcluster_workplace -SVG T -PNG F -n 20"); 


my $finish_time = localtime;
print "\n$finish_time: done!\n\n";







###########################################################
###########################################################
sub bacth_blast_best {

    my ($work_number, $thread_number, $subfunction, $sub_parameters) = @_;
    my ($input_file, $db_file, $output_file, $e_value, $blastn) = @$sub_parameters;
 
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
            my $output_file = $outfile[$job_number-1];
            
            my $Blast_perform_progress = int (($job_number/$work_number)*100); 
            print "$Blast_perform_progress%","..." if ($job_number == 1 or ( ($Blast_perform_progress%10 ==0) && ($progress_record <$Blast_perform_progress)) );

            $threads[$job_number-1]=threads->new(\&$subfunction, $input_file, $db_file, $output_file, $e_value, $blastn);   
            last if ($job_number>=$work_number);                         
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




#****** Warining: this subrouine is only used to find the best hit for each query gene against database, and a loose cutoff for evalue
sub blastn_best {
    my ($input, $db, $outfile, $e, $b) = @_;
    system("$b -outfmt '6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore' -max_target_seqs 2 -evalue $e -num_threads 1 -query $input -db $db -out $outfile"); #BLASTn -max_target_seqs -num_alignments 

}



#########################################################
#########################################################
sub result_statistics {

    my ($workplace, $tbl, $mapping_result) = @_;
    
    my $result_statistics_workplace = "$workplace/result_statistics";
    mkdir $result_statistics_workplace;

    my $tbl_folder = extract_tbl($result_statistics_workplace, $tbl, $mapping_result);  # obtain part_tbl folder

    my $output_statistics = "$result_statistics_workplace/Statistics_OUT";
    mkdir $output_statistics;
    #print "$result_statistics_workplace\n$tbl_folder\n$output_statistics\n";

    my $out1 = "$output_statistics/Gene_number_statistics";
    my $out2 = "$output_statistics/CPS_number_statistics";
    my $out3 = "$output_statistics/Gene_LocusTag_statistics";
    my $out4 = "$output_statistics/CPS_LocusTag_statistics";
    my $out_classs_gene = "$output_statistics/classification_gene";
    my $out_class_anti = "$output_statistics/classification_CPS";

    open (OUT_1, ">$out1") or die "could not open infile $out1\n";  
    open (OUT_2, ">$out2") or die "could not open infile $out2\n";
    open (OUT_3, ">$out3") or die "could not open infile $out3\n";  
    open (OUT_4, ">$out4") or die "could not open infile $out4\n";   

    my $temp_total_tbl = "$result_statistics_workplace/temp_total_tbl";
    chdir $tbl_folder;
 
    #system ("cat * > $temp_total_tbl");
    system ("find ! -name '.' | xargs cat > $temp_total_tbl");  # all files not be .; For large number of files, directly using cat may causing error


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    my %gene_hash;
    my %antibiotic_hash;
    open (TEMP_TBL, $temp_total_tbl);
    my @temp_tbl = <TEMP_TBL>;
    close TEMP_TBL;

    foreach (@temp_tbl){
        chomp;
        my @array = split "\t", $_;
        my $new_keys = "$array[1]|$array[0]"; # antibiotic|gene

        $gene_hash{$new_keys} = $array[1];
        $antibiotic_hash{$array[1]} = $array[0];
        #print "$array[0]\t$array[1]\n";
    }
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    print OUT_1 "Gene_type\t"; 
    print OUT_3 "Gene_type\t"; 

    my $type_check ="";
    foreach (sort keys %gene_hash) {
        $_ =~ s/\|.*//;
        $type_check = $_;
###################################################################################
#       if ($type_check ne $_){
#           $type_check = $_;
            print OUT_1 "$type_check\t";  # print Antibiotic gene name
            print OUT_3 "$type_check\t";  # print Antibiotic gene name

#        }else{
#           print OUT_1 "\t";  # print Antibiotic gene name
#           print OUT_3 "\t";  # print Antibiotic gene name

#        }
###################################################################################
    }

    print OUT_1 "\nStrain\t"; 
    print OUT_2 "Strain\t"; 
    print OUT_3 "\nStrain\t"; 
    print OUT_4 "Strain\t"; 

    foreach (sort keys %gene_hash) {
        $_ =~ s/(.*)\|//;
        print OUT_1 "$_\t";  # print Antibiotic gene name
        print OUT_3 "$_\t";  # print Antibiotic gene name
    }
    print OUT_1 "\n";
    print OUT_3 "\n";

    foreach (sort keys %antibiotic_hash) {
        print OUT_2 "$_\t";  # print Antibiotic name
        print OUT_4 "$_\t";  # print Antibiotic name
    }
    print OUT_2 "\n";
    print OUT_4 "Serotype/number/frequency\n"; # added by xiangyang 2022-02-13

###-------------------------------------------------------------------------------------------------------------------------
    opendir (TBL_FOLDER, $tbl_folder);
    my @total_protein = readdir TBL_FOLDER;
    @total_protein =grep ($_!~/^\./ ,@total_protein);  #delete hidden file . ..
    @total_protein =grep ($_!~/\~$/ ,@total_protein);  #delete temp file ~
    closedir TBL_FOLDER;
    
    my (@array_gene_num, @array_antbiotic_num);
    foreach my $eachfile (sort @total_protein) {
        open (EACHFILE, "$tbl_folder/$eachfile") or die "could not open infile $eachfile\n"; 
        print OUT_1 "$eachfile\t";
        print OUT_2 "$eachfile\t";
        print OUT_3 "$eachfile\t";
        print OUT_4 "$eachfile\t";
        my (@arr_temp_each1, @arr_temp_each2);
        my (%hash_LocusTag1, %hash_LocusTag2);
        while(<EACHFILE>) {
            chomp;
            my @sequencelines = split /\t/;
            push (@arr_temp_each1, $sequencelines[0]);  #Antibiotic_gene (gene)
            push (@arr_temp_each2, $sequencelines[1]);  #Antibiotic (class)
            $hash_LocusTag1{$sequencelines[5]} = $sequencelines[0]; #keys: LocusTag; values: Antibiotic_gene (gene)
            $hash_LocusTag2{$sequencelines[5]."|".$sequencelines[0]."|".$sequencelines[4]} = $sequencelines[1]; #keys: LocusTag; values: Antibiotic (class) #added by xiangyang 2022-02-14
        }
        close EACHFILE;

        #print join (": ", @arr_temp_each1);
        #print "\n";
        my @array_num_temp1;
        foreach my $keygene (sort keys %gene_hash) {
            #print "$keygene\n";
            $keygene =~ s/(.*)\|//;
            my $macth_num1 = grep ($keygene eq $_, @arr_temp_each1);
            print OUT_1 "$macth_num1\t";
            push @array_num_temp1, $macth_num1;   # recode gene number temply;

            my @macth_num1 = grep ($keygene eq $_, @arr_temp_each1);
            my @macth_LocusTag1;
            foreach my $key_LocusTag (keys %hash_LocusTag1) {
                push(@macth_LocusTag1, $key_LocusTag)  if grep ($hash_LocusTag1{$key_LocusTag} eq $_, @macth_num1);
            }
            print OUT_3 join (": ", @macth_LocusTag1), "\t";

        }
        print OUT_1 "\n";
        print OUT_3 "\n";
        @array_num_temp1 = sort{$b<=>$a}@array_num_temp1;
        push @array_gene_num, $array_num_temp1[0];  

        my @serotype; # added by xiangyang 2022-02-13
        my @array_num_temp2;
        foreach my $keyantibiotic (sort keys %antibiotic_hash) {
            my $macth_num2 = grep ($keyantibiotic eq $_, @arr_temp_each2);
            print OUT_2 "$macth_num2\t";
            push @array_num_temp2, $macth_num2;   # recode antibiotic number temply;

            my @macth_num2 = grep ($keyantibiotic eq $_, @arr_temp_each2);
            my @macth_LocusTag2;
            foreach my $key_LocusTag (keys %hash_LocusTag2) {
                push(@macth_LocusTag2, $key_LocusTag)  if grep ($hash_LocusTag2{$key_LocusTag} eq $_, @macth_num2);
                my $c_nu = $key_LocusTag;            # added by xiangyang 2022-02-13
                $c_nu =~ m/(.*?)\|(.*?)\|(.*?)/; # added by xiangyang 2022-02-13
                my $cc_nu = $2; # added by xiangyang 2022-02-13
                $cc_nu =~ s/.*\-//g; # added by xiangyang 2022-02-13
                push @serotype, $cc_nu if grep ($hash_LocusTag2{$key_LocusTag} eq $_, @macth_num2);# added by xiangyang 2022-02-13
            }
            print OUT_4 join (": ", @macth_LocusTag2), "\t";

        }
        print OUT_2 "\n";
########################################################## # added by xiangyang 2022-02-13
        my @f_arr;
        foreach my $c_count (@serotype){
             push @f_arr, (split /#/, $c_count) if $c_count =~ /#/;
             push @f_arr, $c_count if $c_count !~ /#/;
        }
        my %c_hash;
        foreach (@f_arr){
            $c_hash{$_}++;
        }


        my $key_max;
        my @temp_max = sort {$a<=>$b} values %c_hash;
        my $c_max= $temp_max[-1];
        for my $c_key (keys %c_hash) {
           if ($c_hash{$c_key} == $c_max){
               $key_max .= "$c_key|";   #multipe values have the same keys for max
           }
        }
        $key_max =~ s/\|$//; #delete "|" at the end of the string
        my $number = scalar @serotype;
        my $occur_pro = $c_max/(scalar @serotype)*100; # added by xiangyang 2022-02-13
        $occur_pro = sprintf "%.2f",$occur_pro; # added by xiangyang 2022-02-13
###########################################################
        print OUT_4 "Serotype $key_max/$number/$occur_pro","%\n"; # added by xiangyang 2022-02-13

        @array_num_temp2 = sort{$b<=>$a}@array_num_temp2;
        push @array_antbiotic_num, $array_num_temp2[0]; 
    }

    @array_gene_num = sort{$b<=>$a}@array_gene_num;
    my $max_gene_number = $array_gene_num[0];

    @array_antbiotic_num = sort{$b<=>$a}@array_antbiotic_num;
    my $max_antibiotic_number = $array_antbiotic_num[0];

    #print "max_gene_number: $max_gene_number\n";
    #print "max_antibiotic_number: $max_antibiotic_number\n";
    close OUT_1;
    close OUT_2;
    close OUT_3;
    close OUT_4;
    system("rm $temp_total_tbl");
    my $tbl_heatmap_dir_gene = tbl_heatmap_generation($out1, $result_statistics_workplace, "tbl_heatmap_gene", $out_classs_gene, $max_gene_number);
    my $tbl_heatmap_dir_antibiotic = tbl_heatmap_generation($out2, $result_statistics_workplace, "tbl_heatmap_antibiotic", $out_class_anti, $max_antibiotic_number);
###-------------------------------------------------------------------------------------------------------------------------

}



#########################################################
#########################################################
sub extract_tbl {

    my ($home_directory, $tbl_folder, $mapping_result) = @_;

    my $out_protein_folder = "$home_directory/tbl_part";

    mkdir $out_protein_folder;
    #system ("rm $out_protein_folder/*");

    my $temp_list_split_dir = split_list_tbl ($mapping_result, $home_directory);

    opendir (TBL_FOLDER, $tbl_folder);
    my @total_protein = readdir TBL_FOLDER;
    @total_protein =grep ($_!~/^\./ ,@total_protein);  #delete hidden file . ..
    closedir TBL_FOLDER;

    foreach my $eachfile (@total_protein) {
        my %seq_hash;
        open (EACHFILE, "$tbl_folder/$eachfile") or die "could not open infile $eachfile\n"; 

        while(<EACHFILE>) {
            chomp;
            next if($_ eq '');
            my @sequencelines = split /\t/ if $_ !~ /^>Contig/;
            $sequencelines[3] =~ s/.*;//g if $_ !~ /^>Contig/;
            $seq_hash{$sequencelines[3]} = $_ if $_ !~ /^>Contig/; #key: locus tag; value: all feature in one row in tbl file
        }
        close EACHFILE;
        my $matck_file = $eachfile; # added by xiangyang 2022-02-15
        $matck_file =~ s/\.tbl.*part/\.tbl/g; # added by xiangyang 2022-02-15
        if (-e "$temp_list_split_dir/$matck_file") { # replace $eachfile with $matck_file added by xiangyang 2022-02-15
            open (LIST_TP, "$temp_list_split_dir/$matck_file"); # replace $eachfile with $matck_file added by xiangyang 2022-02-15

            while (<LIST_TP>){
                chomp;

                my @array=split ("\t", $_);
                $array[0] =~ s/#/;/g;
                my @array_2 = split ";", $array[0];
##########################################################################################
                my $fhs;
                my $source = $array_2[2];
                #print $source, "\n";
                unless( $fhs->{$source} ) { #当hash键(第一列)不存在时，创建字符串句柄 
                    open my $fh, ">>$out_protein_folder/$source" or die '...';
                    $fhs->{$source} = $fh; 
                } 
                my $identity = $seq_hash{$array_2[0]} if defined $seq_hash{$array_2[0]}; # added by xiangyang 2022-02-13
                my $iden_cov = $array[8]."_".$array[9]."-".$array[10]; # added identity and qcoverage-scoverage for blastn from CPS_cluster_mapping.result file by xiangyang 2022-02-13
                $identity =~ s/\tCDS\t/\t\Q$iden_cov\E\t/  if defined $seq_hash{$array_2[0]}; # added by xiangyang 2022-02-13
                $identity =~ s/\\//g;
                $identity =~ s/\.00//g;
                print { $fhs->{$source} } "$array[2]\t$array[4]\t$identity\n" if defined $seq_hash{$array_2[0]}; 
##########################################################################################
            }
            close LIST_TP;
        }
   
    }
 
    system("rm -rf $temp_list_split_dir");
    return $out_protein_folder;
}


#########################################################
#########################################################
sub split_list_tbl {

    my ($list, $home_directory) = @_;

    my $temp_list_split_dir = "$home_directory/temp_tbl_split_dir";
    mkdir $temp_list_split_dir;
    open (LIST_PARSE, "$list");

    while (<LIST_PARSE>){
        chomp;
        $_ =~ s/\t(.*)//g;
        my $right_content = $1;
##########################################################################################
        my $fhs;
        my $source = $_;
        $source =~ s/.*;//g;
        unless( $fhs->{$source} ) { #当hash键(第一列)不存在时，创建字符串句柄 
            open my $fh, ">>$temp_list_split_dir/$source.tbl" or die '...';
            $fhs->{$source} = $fh; 
        } 
        print { $fhs->{$source} } "$_\t$right_content\n";
##########################################################################################
    }

    close LIST_PARSE;
    return $temp_list_split_dir;
}



#########################################################
#########################################################
sub tbl_heatmap_generation {

    my ($list, $home_directory, $directory_name, $out_classification, $max_number) = @_;

    my $tbl_heatmap_dir = "$home_directory/$directory_name";
    mkdir $tbl_heatmap_dir;
    open (LIST_PARSE, "$list");
    my @list_parse = <LIST_PARSE>;
    my @Antibiotic_name = split "\t", shift @list_parse;

    open (OUT_CLASSIFICATION, ">$out_classification") or die "could not open infile $out_classification\n";  
    print OUT_CLASSIFICATION "$max_number\n";
####### caculate the number of gene number belonging to specific antibiotic
    my %hash_number;
    #print "number: ",  scalar @Antibiotic_name, "\n";
    foreach my $element_array (@Antibiotic_name){
        #chomp;
        if ($element_array !~ /^\s|^Strain/){
            $hash_number{$element_array} = grep ($element_array eq $_, @Antibiotic_name);
            #push (@Antibiotic_name_new, $element_array);
        }
    }

    my %count;
    my @Antibiotic_name_new = grep { ++$count{ $_ } < 2; } @Antibiotic_name;
    foreach (@Antibiotic_name_new){
        print OUT_CLASSIFICATION "$_\t",  $hash_number{$_}, "\n" if $_ !~ /^\s|^Strain/;
    }

#######caculate the number of gene number belonging to specific antibiotic

    my @gene_name = split "\t", shift @list_parse if $list =~ /Gene_number_statistics/;
    foreach (@list_parse){
        chomp;
        my @array=split ("\t", $_);
##########################################################################################
        my $fhs;
        my $source = shift @array;
        unless( $fhs->{$source} ) { #当hash键(第一列)不存在时，创建字符串句柄 
            open my $fh, ">>$tbl_heatmap_dir/$source.tbl" or die '...';
            $fhs->{$source} = $fh; 
        } 
        my $length_count=0;
        foreach(@array){
            my $start = 1+50*$length_count;
            my $end = 50*($length_count+1);
            $length_count++;
            print { $fhs->{$source} } "$start\t$end\tCDS\t$gene_name[$length_count]\t$_\n" if $list =~ /Gene_number_statistics/;
            print { $fhs->{$source} } "$start\t$end\tCDS\t$Antibiotic_name[$length_count]\t$_\n" if $list !~ /Gene_number_statistics/;
        }
##########################################################################################
    }

    close LIST_PARSE;
    return $tbl_heatmap_dir;
}



#########################################################
#########################################################
sub mapping {

    my ($mappingFile, $filter_blast_cluster, $cluster_mapping)= @_;
    my %mapping;
 
    open (MAPPING_OUT, ">$cluster_mapping");
    open (MAPPING, $mappingFile);
    my $cout_test=0;

    while (my $line = <MAPPING>) {
        chomp($line);
        $line =~ s/\r//g;
        ## BacMet-Scan Version 2
        $cout_test++;
        my ($GI_number, $GenBank_ID, $Gene_name, $Organism, $Compound, $NCBI_annotation) = split("\t",$line);

        $mapping{$GI_number} = "$Gene_name\t$Organism\t$Compound\t$NCBI_annotation";
        $mapping{$GenBank_ID} = "$Gene_name\t$Organism\t$Compound\t$NCBI_annotation";

    }
    close MAPPING;

    open (ALL_ALL, $filter_blast_cluster);
    while (my $ele = <ALL_ALL>) {
        chomp($ele);
        next unless $ele !~ /Qseqid/; #???
        my @arr_temp = split('\t',$ele);

        print MAPPING_OUT "$arr_temp[0]\t$arr_temp[1]\t$mapping{$arr_temp[1]}\t$arr_temp[2]\t$arr_temp[3]\t$arr_temp[4]\t$arr_temp[5]\t$arr_temp[6]\t$arr_temp[7]\n";

    }
    close ALL_ALL;

    close MAPPING_OUT;

    return $cluster_mapping;

}



#########################################################
#########################################################
sub SARG_map_transformation {

    my ($database, $list, $map_table) = @_;

my %hash_id;
open (SEQ, $database);
while (<SEQ>){
    chomp;

    if (/^>(.+)/) {

    my $id = $1;
    $id =~ s/ /\t/;

    $id =~ s/\[/\t\[/  if $id =~ /\[/;
    #print "SEQ: $id\n";
    my @arr_id = split "\t", $id;
    $hash_id{$arr_id[0]} = "$arr_id[1]\t$arr_id[2]" if $id =~ /\[/;
    $hash_id{$arr_id[0]} = "$arr_id[1]\t[Organism_unavailable]" if $id !~ /\[/;

    #print "$count\t$arr_id[0]\t$hash_id{$arr_id[0]}\n";

    }

}
close SEQ;

open (IN, $list);
open (OUT, ">$map_table");

print OUT "GI_number\tCategories_in_database\tGene_name\tOrganism\tCompound\tNCBI_annotation\n";

my %text_hash;
while (<IN>){
    chomp $_;
    next unless $_ =~ /\_\_/;
    $_ =~ s/'//g;
    $_ =~ s/, /,/g;
    $_ =~ s/[\[\]]//g;
    my @array = split "\t", $_;
    
    my @array_id = split /,/, $array[1];

    foreach my $id(@array_id) {

        $text_hash{$id} = $array[0];

    }
    
}
close IN;

foreach my $key (keys %text_hash) {
    my $value = $text_hash{$key};
    my @value_split = split "\_\_", $value; 
    my ($NCBI_annotation_seq, $Organism_seq) = split "\t", $hash_id{$key} if defined $hash_id{$key};
    print OUT "$key\t$value\t$value_split[1]\t$Organism_seq\t$value_split[0]\t$NCBI_annotation_seq\n" if defined $hash_id{$key};
}
close OUT;
}



#########################################################
#########################################################
sub split_to_subfiles{

    my ($file,$split,$output_dir,$log_file);
    open(my $out_fh, '>-') or die "Could not open stdout for writing";

    if ( @_ && $_[0] eq __PACKAGE__)
    {
        GetOptions('i|input-file=s' => \$file,
                   'o|output-dir=s' => \$output_dir,
                   'n|split-number=i' => \$split,
                   'l|log-file=s' => \$log_file) or die "Invalid options\n";
    }
    else
    {
        ($file,$split,$output_dir,$log_file) = @_;
    }

       
    
    #this will find records with a quickness
    (my $total) = `grep -c ">" $file`;
    chomp $total;

    die "no fasta records identified, exiting.\n" unless $total;
    
    my $records_per_file = int ($total / $split);
    my $leftover_records = $total % $split;
    
    if ($split > $total) {
       warn "Total records in file ($total) less than splitnum ($split); adjusting toatal\n";
       $records_per_file = 1;
    }
    
    my $input_file_name = basename($file);
    my $output_base_path = (defined $output_dir) ? "$output_dir/$input_file_name" : $input_file_name;
    my $in = new Bio::SeqIO (-file=>$file, -format=>"fasta");
    my $x;
    my @outs;
    my $records;
    for my $x (1..$split) { 
      #adjust # of records per file for leftover records
      my $adjusted_records_per_file; 
      $adjusted_records_per_file = $records_per_file+1 if $x <= $leftover_records;
      push @outs, new Bio::SeqIO (-file=>">$output_base_path.$x", -format=>"fasta")
    }
    
    my $out = shift @outs;
    my $filecounter = 1;
    my $recordcounter =1;
    while (my $seq = $in->next_seq) {

      my $adjusted_records_per_file = $filecounter<=$leftover_records?$records_per_file+1:$records_per_file;
      $out->write_seq($seq);
      if (++$records>=$adjusted_records_per_file) {
        $out = shift @outs; $records =0; $filecounter++;   
      }
    }
}




#########################################################
#########################################################
sub blast_filter {  #blastn 2023-01-10 xiangyang li edited

    my ($infile, $e_value, $identify, $coverage, $match_length) = @_;
    open (INFILE, $infile);
    my $all_vs_all_cluster = dirname ($infile)."/all_vs_all.cluster";

    open (OUTFILE, ">$all_vs_all_cluster");
    print OUTFILE "$output_title\n";
    my @infile;
    while(<INFILE>){
        chomp;
        $_ =~ s/ //g;

        my @blast_result0=split '\t', $_;
        my $qcover0 =  $blast_result0[5]/$blast_result0[3] * 100;  #added by xiangyang
        $qcover0 = sprintf("%.2f", $qcover0);
        my $scover0 = $blast_result0[5]/$blast_result0[4] * 100;  #added by xiangyang
        $scover0 = sprintf("%.2f", $scover0);
        if (($blast_result0[10] <= $e_value) && ($blast_result0[2] >= $identify) && ($qcover0 >= $coverage) && ($scover0 >= $coverage) && ($blast_result0[5] >= $match_length)){ 
            push @infile, $_;
        }
    }
    close INFILE;
    my $add_row = "Addition_row	Addition_row	100	1000	1100	1000	1	1000	1	1100	0.0	 0";
    push @infile, $add_row;
    my $row=0;
    my %st_hash;
    for(my $i=0; $i<scalar @infile -1; $i++){

        my @blast_result = split '\t', $infile[$i];
        my @blast_result_k = split '\t', $infile[$i+1];
        my $qcover =  ($blast_result[7] - $blast_result[6] +1)/$blast_result[3] * 100;  #added by Yangzi lin
        my $scover = ($blast_result[9] - $blast_result[8] +1)/$blast_result[4] * 100;  #added by xiangyang
        $scover = sprintf("%.2f", $scover);

        my $qcover_k =  ($blast_result_k[7] - $blast_result_k[6] +1)/$blast_result_k[3] * 100;  #added by xiangyang
        $qcover_k = sprintf("%.2f", $qcover);
        my $scover_k = ($blast_result_k[9] - $blast_result_k[8] +1)/$blast_result_k[4] * 100;  #added by xiangyang
        $scover_k = sprintf("%.2f", $scover_k);
 
            #set cut-off of e-value, identity, coverage, match_length for homologous protein
            #if ($blast_result[5]>0) {$blast_result[11] = $blast_result[11]/$blast_result[5];} else {$blast_result[11] = $blast_result[11];}

            if ($blast_result[0] eq $blast_result_k[0]) {
                if ($blast_result[11] <= $blast_result_k[11]) {
                    if ($blast_result[3] == $blast_result[4]){
                        $st_hash{$blast_result[0]} = "$blast_result[0]\t$blast_result[1]\t$blast_result[11]\t$blast_result[10]\t$blast_result[2]\t$qcover\t$scover\t$blast_result[5]\n";
                        $i = $i+1;

                    }else{
                        $st_hash{$blast_result_k[0]} = "$blast_result_k[0]\t$blast_result_k[1]\t$blast_result_k[11]\t$blast_result_k[10]\t$blast_result_k[2]\t$qcover_k\t$scover_k\t$blast_result_k[5]\n";
                        $i = $i+1;

                    }
                    
                }else{
                    $st_hash{$blast_result[0]} = "$blast_result[0]\t$blast_result[1]\t$blast_result[11]\t$blast_result[10]\t$blast_result[2]\t$qcover\t$scover\t$blast_result[5]\n";

                    $i = $i+1;
                }

            }else{               
                $st_hash{$blast_result[0]} = "$blast_result[0]\t$blast_result[1]\t$blast_result[11]\t$blast_result[10]\t$blast_result[2]\t$qcover\t$scover\t$blast_result[5]\n";    
            }


    }


    foreach (sort keys %st_hash){
        print OUTFILE $st_hash{$_};
    }
    close OUTFILE;

    return $all_vs_all_cluster;

}



#########################################################
#########################################################
sub de_repeat_array {
    
    my @array = @_;
    my %hash; #定义一个空hash
    my @result=grep {++$hash{$_}<2} @array; #去除冗余元素
    #print join ("\n", @result), "\n";
    return @result;
}




##################################################################################################
###### Subrounting--do bacth work
###### Function:
###### do bacth work
##################################################################################################
sub batch_genbank_sequence_TFT_extract {
    my ($path_genbank, $path_TFT, $path_gene, $thread_number)=@_;

    opendir PATH_GENBANK, $path_genbank or die "could not open $path_genbank";
    my @path_genbank_temp = readdir PATH_GENBANK;
    foreach (@path_genbank_temp){
        if (/[ =]/) {
            my $new_name = $_;
            $new_name =~ s/ /_/g;   #repalce the space with "_"  for all GenBank files
            $new_name =~ s/=/_/g;   #repalce "=" with "_"  for all GenBank files
            $new_name =~ s/_+/_/g;  
            system ("mv '$path_genbank/$_' '$path_genbank/$new_name'");  
        }

    }
    closedir PATH_GENBANK;

    opendir PATH_GENBANK, $path_genbank or die "could not open $path_genbank";
    my @path_genbank = readdir PATH_GENBANK;
    @path_genbank =grep ($_!~/^\./ ,@path_genbank);  #delete hidden file . .. 
    closedir PATH_GENBANK;

    my $sub_file_number = scalar @path_genbank;

    my @input;
    my @outfile_1;
    my @outfile_2;
    my @outfile_3;
 
    foreach my $file_genbank(@path_genbank){ 
        $file_genbank =~ s/ /_/g;
        my $input="$path_genbank/$file_genbank"; 
        push (@input,$input);  
        #@input = sort @input;  

        #@outfile_1 = sort @outfile_1;
        my $output_2="$path_TFT/$file_genbank.tbl";
        push (@outfile_2,$output_2);
        #@outfile_2 = sort @outfile_2;
        my $output_3="$path_gene/$file_genbank.fasta";
        push (@outfile_3,$output_3);

    }

    my $thread;
    my @threads;
    my $job_number=0;
    print "    GenBank_extraction_percent: ";
    while(){ 

        last if ($job_number eq $sub_file_number);  

        while(scalar(threads->list())<$thread_number) {     #set thread number
            my $progress_record = int (($job_number/$sub_file_number)*100);
            $job_number++;    
            my $input_file = $input[$job_number-1];
            my $output_file_2 = $outfile_2[$job_number-1];
            my $output_file_3 = $outfile_3[$job_number-1];
            my $GenBank_extraction_progress = int (($job_number/$sub_file_number)*100); 
            print "$GenBank_extraction_progress%","..." if ($job_number == 1 or ( ($GenBank_extraction_progress%10 ==0) && ($progress_record <$GenBank_extraction_progress)) );
            $threads[$job_number-1]=threads->new(\&genbank_sequence_TFT_extract, $input_file, $output_file_2, $output_file_3); 

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
    print "done\n";
}
##################################################################################################
###### Subrounting--genbank_sequence_TFT_extract
###### Function:
###### extract protein sequence and tbl (cds/rRNA/tRNA information) from genbank file
##################################################################################################
sub genbank_sequence_TFT_extract {

    my ($input, $output_2, $output_3)=@_;
    my $file_name = basename($input);
    open(OUTPUT_2, ">$output_2");
    open(OUTPUT_3, ">$output_3");

        my $in = new Bio::SeqIO( -file => $input, -format => 'genbank' );

        while ( my $seqObject = $in->next_seq() ) {

            my $acc = $seqObject->accession;
            my @features = $seqObject->get_SeqFeatures();
            my $count = 0;       
            
            foreach (@features) {
                my $feat = $_;	        
                unless ($feat->primary_tag =~ /^CDS|^rRNA|^tRNA|^source/) {
                       next;
                }

                my $primary;
                ($primary) = $feat->primary_tag;
                if ($primary eq "source") {
                    print OUTPUT_2 ">Contig: $acc\n";
                }
                else {
                    $count++;
                    #print "genome_$count\n";
              
                    my ($start, $stop);
                    my $location  = $feat->location;
                    my $locString = $location->to_FTstring;
                
                    if ($locString =~ /,/) {
                        if ($locString =~ /complement/) {
                            my @loc = split( /,/, $locString);
                            if ($loc[0] =~ /(\d+)\.\.(\d+)/){ 
                                my $length_part = $2-$1+1;
                                $loc[1] =~ /(\d+)\.\.(\d+)/;
                                $start = $2;
                                $stop = $1-$length_part;
                             }
                             else{ 
                                 $loc[1] =~ /(\d+)\.\.(\d+)/;
                                 $start = $2;
                                 $stop = $1-1;
                             }                   
                         } 
                        else { 
                            my @loc = split( /,/, $locString );
                            if ($loc[0] =~ /(\d+)\.\.(\d+)/) { 
                                my $length_part = $2-$1+1;
                                $loc[1] =~ /(\d+)\.\.(\d+)/;
                                $stop = $2;
                                $start = $1-$length_part;
                            } 
                            else{ 
                                $loc[1] =~ /(\d+)\.\.(\d+)/;
                                $stop = $2;
                                $start = $1-1;
                            }
                        }
                    }
                    else {
                        if ($locString =~ /complement\((.+)\.\.(.+)\)/) {
                            $start = $2;
                            $stop = $1;
                         }
                         else {
                             $locString =~ /(.+)\.\.(.+)/;
                             $start = $1;
                             $stop = $2;
                         }
                    }
                    my $gene_name;
                    if ( $feat->has_tag('gene') ) {
                        ($gene_name) = $feat->get_tag_values('gene');
                    }

                    my $locus_tag;
                    if ( $feat->has_tag('locus_tag') ) {
                        ($locus_tag) = $feat->get_tag_values('locus_tag');
                    }else {$locus_tag = "CDS_$acc"."_".$count;}
               
                    my $product;
                    if ( $feat->has_tag('product') ) {
                        ($product) = $feat->get_tag_values('product');
                        $product =~ s/ /_/g;
                        $product =~ s/;/_/g;
                        $product =~ s/,/_/g;
                    }else {
                        $product = "No_product_tag";
                    }

                    my $pseudo;
                    if ( $feat->has_tag('pseudo') ) {
                        $pseudo = "pseudo";
                        if (defined $gene_name){
                            print OUTPUT_2 "$start\t$stop\t$primary\t$gene_name;$locus_tag\t$product\t$acc\t$pseudo\n";
                        } else {
                            print OUTPUT_2 "$start\t$stop\t$primary\t$locus_tag\t$product\t$acc\t$pseudo\n";
                          }
                    }
                    else { 
                        if (defined $gene_name){
                            print OUTPUT_2 "$start\t$stop\t$primary\t$gene_name;$locus_tag\t$product\t$acc\n";
                        } else {
                            print OUTPUT_2 "$start\t$stop\t$primary\t$locus_tag\t$product\t$acc\n";

                          }
                    }
 
                    my $dna;
                    if ( $feat->spliced_seq->seq ) {
                        $dna = $feat->spliced_seq->seq;
                        print OUTPUT_3 ">$locus_tag#$product;$file_name\n$dna\n"; 
                    }
                }
            }
        }
    close OUTPUT_2;
    close OUTPUT_3;
    my @feature = stat ($output_3);
    my $size = $feature[7];
    &prodigal_tool($input, $output_3, $output_2) if $size == 0;
 
}

#########################################
# bacth annnotate genome using prodigal
#########################################
sub prodigal_bacth {
    my ($input_dir, $prokka_workplace, $gff_dir, $thread_number) = @_;

    opendir (DIR_SEQ, $input_dir) or die "could not open $input_dir";
    my @input_dir = readdir DIR_SEQ; 
    @input_dir = grep ($_!~/^\./ ,@input_dir);
    closedir DIR_SEQ;


    my @input;
    my @output;
    my @outgff;
    foreach(sort @input_dir){
        push @input, "$input_dir/$_";
        push @output, "$prokka_workplace/$_";    
        push @outgff, "$gff_dir/$_";    
    }

    my $sub_file_number =  scalar @input;
    my $thread;
    my @threads;
    my $job_number=0;

    print "    Annotating genome using prodigal: ";
    while(){ 

        last if ($job_number eq $sub_file_number);  
        while(scalar(threads->list())<$thread_number) {     #set threadnumber；
        my $progress_record = int (($job_number/$sub_file_number)*100);
        $job_number++;    
        my $input_file = $input[$job_number-1];
        my $output_file = $output[$job_number-1];
        my $output_gff = $outgff[$job_number-1];

        my $genome_annotation_percent = int (($job_number/$sub_file_number)*100); 
        print "$genome_annotation_percent%","..." if ($job_number == 1 or ( ($genome_annotation_percent%10 ==0) && ($progress_record <$genome_annotation_percent)) );

        $threads[$job_number-1]=threads->new(\&prodigal_tool, $input_file, $output_file, $output_gff); 
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

    print "done\n";

}


#########################################
# gene and gff file obtain using prodigal tool
#########################################
sub prodigal_tool {
    my ($input, $outgene, $outgff) =@_;
    $outgene =~ s/.fasta$//;
    $outgff =~ s/.tbl$//;
    system ("prodigal -i $input -f gbk -g 11 -q -d $outgene -o $outgff");
    &gene_file($outgene);
    &tbl_file($outgff);
}


#########################################
# gene file conversion 
#########################################
sub gene_file {
    my ($out) = @_;

    my $genome_name = basename($out);
    my $id = substr($genome_name, 0, 15);
    my $count=0;

    my $term     = $/; # input record separator;
    my %seq_hash;
    open (FAS, $out);
    $/ = ">"; # input record separator;
    open (OUT, ">$out.fasta");
    while(<FAS>){
        chomp;
        next if($_ eq '');
        $count++;
        my ($fid, @sequencelines) = split /\n/;
        foreach my $line (@sequencelines) {
            $seq_hash{$fid} .= $line;
        }
        my $gene = $id."_".$count;
        print OUT ">$gene#XXX;$genome_name\n", $seq_hash{$fid}, "\n";
    }
    $/ = $term;
    close FAS;
    close OUT;
    system("rm -rf $out");


} # end of gene_file


#########################################
# gff file was transfomed to tbl file 
#########################################
sub tbl_file {
    my ($gff) = @_;

    my $genome_name = basename($gff);
    my $id = substr($genome_name, 0, 15);
    my $count =0;
    my %contig;
    open (IN, $gff);
    open (OUT, ">$gff.tbl");
    while(<IN>){
        chomp;
        if ($_ =~ /seqhdr/){
            my $contig = $_;
            $contig =~ s/.*?"//;
            $contig =~ s/".*//;
            print OUT ">Contig: $contig\n";
            $contig{$genome_name}=$contig;
        }
        if ($_ =~ /^     CDS/){
            $count++;
            my $cds = $_;
            $cds =~ s/^     CDS.*?(\S)/$1/;#complement(571..1734)
            $cds =~ s/\(|\)//g;

            my ($s,$e);
            if ($cds =~ /complement/){
                $cds =~ s/complement//; 
                $cds =~ /(.*?)\.\.(.*)/; 
                $s= $2;
                $e= $1;
            }else{
                $cds =~ /(.*?)\.\.(.*)/; 
                $s= $1;
                $e= $2;
            }
            my $gene = $id."_".$count;
            print OUT "$s\t$e\tCDS\t$gene\tXXX\t$contig{$genome_name}\n";
        }
    }
    close IN;
    close OUT;
    system("rm -rf $gff");

} # end of gene_file
