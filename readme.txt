my $usage = <<USAGE; 

=NAME

Pneumo-Typer.pl

=DESCRIPTION

    Run this command to enble users to obtain predict serotype and determine sequence type (ST/cgST) for Streptococcus pneumoniae, depending on obtaining the potential cps genes by a local blastn analysis using multiple threads. 

=USAGE

Pneumo-Typer.pl -dir genbank_file_directory [options]

FOR EXAMPLE: 

perl /Users/zilinyang/Desktop/桌面材料/李向阳-论文投稿汇总/Serotype_ST_manuscript/pneumo-typer/pneumo-typer/pneumo-typer.pl -d /Users/zilinyang/Desktop/桌面材料/李向阳-论文投稿汇总/Serotype_ST_manuscript/Pneumo-Typer_v1.03/final_93_gbk -t 10 -s T  -p T -o /Users/zilinyang/Desktop/桌面材料/李向阳-论文投稿汇总/Serotype_ST_manuscript/Pneumo-Typer_v1.03/final_93_gbk.OUT

perl /home/xiangyang/Mazhongrui/Streptococcus_pneumoniae/Pneumo-Typer_v1.03/Pneumo-Typer.pl -d /home/xiangyang/Mazhongrui/Streptococcus_pneumoniae/FIG_data -t 190 -s T -c T -m T -o /home/xiangyang/Mazhongrui/Streptococcus_pneumoniae/FIG_data_OUT

nohup perl /home/xiangyang/Mazhongrui/Streptococcus_pneumoniae/Pneumo-Typer_v1.03/Pneumo-Typer.pl -d /home/xiangyang/Mazhongrui/Streptococcus_pneumoniae/prokka_workplace/prokka_gbk_folder_plus72 -t 190 -b T -jsp T -jse T

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
       -jsp, --jump_split
           After splitting file, users can jump the process of spiltting file to perform blastn analysis again.
       -jse, --jump_serotype
           After blastn analysis, users can jump the process of doing blastn file to perform Serotype analysis again.
       -p, --prodigal_annotation
           Annotate all genomes using prodigal. 
       -m, --mlst
           Perform mlst analysis (Default: T). 
       -c, --cgmlst
           Perform cgmlst analysis. It need >10 mins for one genome (Default: F).
       -Ts, --test
           Run pneumo-typer using Test_data as input to check whether pneumo-typer is installed successfully (Default: F).
       -h, --help
           Show this message.

=AUTHOR

Dr. Xiangyang Li (E-mail: lixiangyang\@fudan.edu.cn), Kaili University; Bacterial Genome Data mining & Bioinformatic Analysis (www.microbialgenomic.cn/).

=COPYRIGHT

Copyright 2023, Xiangyang Li. All Rights Reserved.

USAGE