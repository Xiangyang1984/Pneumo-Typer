use strict;
use warnings;
use threads;
use Bio::SeqIO;
use File::Basename qw<basename dirname>;

# perl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/ST_tool/gene_bacth_extract.pl /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/ST_tool/67_gbk 3 /home/xiangyang/Mazhongrui/Serotype_ST/Pneumo-Typer/ST_tool/out_folder

my $path_genbank = $ARGV[0];
my $thread_number = $ARGV[1];
my $out_folder = $ARGV[2];
mkdir $out_folder;
my $path_protein = "$out_folder/path_protein";
mkdir $path_protein;

my $path_gene = "$out_folder/path_gene";
mkdir $path_gene;

print "\nExtract the Amino acid sequence and Nucleotide sequence of genes for each query genome:\n";
&batch_genbank_sequence_TFT_extract($path_genbank, $path_protein, $path_gene, $thread_number);

sub batch_genbank_sequence_TFT_extract {

    my ($path_genbank, $path_protein, $path_gene, $thread_number)=@_;

    system ("rename 's/ /_/g' $path_genbank/*");  #repalce the space with "_"  for all GenBank files
    system ("rename 's/:/_/g' $path_genbank/*");  #repalce the colon with "_"  for all GenBank files
    system ("rename 's/_=_/_/g' $path_genbank/*");  #repalce "=" with "_"  for all GenBank files
    opendir PATH_GENBANK, $path_genbank or die "could not open $path_genbank";

    my @path_genbank = readdir PATH_GENBANK;
    @path_genbank =grep ($_!~/^\./ ,@path_genbank);  #delete hidden file . .. 
    #@path_genbank =grep ($_!~/$\~/ ,@path_genbank);  #delete hidden file . .. 
    closedir PATH_GENBANK;

    my $sub_file_number = scalar @path_genbank;
    $thread_number = $sub_file_number if $thread_number >= $sub_file_number;
    my @input;
    my @outfile_1;
    my @outfile_2;
 
    foreach my $file_genbank(@path_genbank){ 

        my $input="$path_genbank/$file_genbank"; 
        push (@input,$input);  
        @input = sort @input;  

        my $output_1="$path_protein/$file_genbank";
        push (@outfile_1,$output_1);
        @outfile_1 = sort @outfile_1;

        my $output_2="$path_gene/$file_genbank";
        push (@outfile_2,$output_2);
        @outfile_2 = sort @outfile_2;

    }

    my $thread;
    my @threads;
    my $job_number=0;
    my @jindu;
    print "\nGenBank_extraction_percent (Number of GenBank files: $sub_file_number): ";
    while(){ 

        last if ($job_number eq $sub_file_number);  

        while(scalar(threads->list())<$thread_number) {     #set thread number
            my $progress_record = int (($job_number/$sub_file_number)*100);
            $job_number++;    
            my $input_file = $input[$job_number-1];
            my $output_file_1 = $outfile_1[$job_number-1];
            my $output_file_2 = $outfile_2[$job_number-1];

            my $GenBank_extraction_progress = int (($job_number/$sub_file_number)*100); 
            print "$GenBank_extraction_progress%","..." if ($job_number == 1 or ( ($GenBank_extraction_progress%10 ==0) && ($progress_record <$GenBank_extraction_progress)) );
            $threads[$job_number-1]=threads->new(\&genbank_protein_sequence_extract, $input_file, $output_file_1, $output_file_2); 

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


sub genbank_protein_sequence_extract {
####### Part of script in this subroutine comes from CGView Comparison Tool. Reference: Grant JR, Arantes AS, Stothard P (2012) Comparing thousands of circular genomes using the CGView Comparison Tool. BMC Genomics 13:202.

    my ($input, $output_1, $output_2)=@_;
    my $file_name = basename($input);
    open(OUTPUT_1, ">$output_1");
    open(OUTPUT_2, ">$output_2");

        my $in = new Bio::SeqIO( -file => $input, -format => 'genbank' );

        while ( my $seqObject = $in->next_seq() ) {
            my $acc = $seqObject->accession;
            my @features = $seqObject->get_SeqFeatures();
            my $count = 0;       
            
            foreach (@features) {
                my $feat = $_;	        
                unless ($feat->primary_tag =~ /^CDS/) {
                    next;
                }
                    if ( $feat->has_tag('pseudo') ) {

                    }else {$count++;}
                    #print "genome_$count\n";

                    my $locus_tag;
                    if ( $feat->has_tag('locus_tag') ) {
                        ($locus_tag) = $feat->get_tag_values('locus_tag');

                    }else {$locus_tag = "CDS_$acc"."_".$count;}
               
                    my $product;
                    if ( $feat->has_tag('product') ) {
                        ($product) = $feat->get_tag_values('product');
                        $product =~ s/ /_/g;
                        $product =~ s/;/|/g;
                    }
 
                    my $protein_seqeunce;
                    if ( $feat->has_tag('translation')) {
                        ($protein_seqeunce) = $feat->get_tag_values('translation');
                      
                        print OUTPUT_1 ">$locus_tag;$product;$file_name\n$protein_seqeunce\n";
                        
                    }

                    my $dna;
                    if ( $feat->spliced_seq->seq ) {
                        $dna = $feat->spliced_seq->seq;
                        print OUTPUT_2 ">$locus_tag;$product;$file_name\n$dna\n"; 
                    }
            }
        }

    close OUTPUT_1;
    close OUTPUT_2;
    #my @feature = stat ($output_1);
    #my $size = $feature[7];
    #print "$output_1\t$size\n" if $size == 0;

}

