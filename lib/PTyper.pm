package PTyper;


use strict;
use warnings;
use threads;
use Bio::SeqIO;
use File::Basename qw<basename dirname>;
use FindBin;
use lib "$FindBin::Bin/lib";

####################################################
####################################################
sub batch_genomenucleotide_extract {

    my ($path, $fasta1, $fasta2, $thread_number)=@_;

    opendir DIR, $path or die $!;
    ################################start rename the file in following case
    my @path_temp = readdir DIR;
    foreach (@path_temp){
        if (/[ =]/) {
            my $new_name = $_;
            $new_name =~ s/ /_/g;   #repalce the space with "_"  for all GenBank files
            $new_name =~ s/=/_/g;   #repalce "=" with "_"  for all GenBank files
            $new_name =~ s/_+/_/g;  
            system ("mv '$path/$_' '$path/$new_name'");  
        }

    }
    closedir DIR;
    #################################end rename the file in following case
    opendir DIR, $path or die $!;
    my @dir = readdir DIR; 
    @dir = grep ($_!~/^\./, @dir);
    closedir DIR;

    my @input;
    my @outfile1;
    my @outfile2;

    foreach my $file(@dir){

        my $input="$path/$file";
        push (@input,$input);  

        my $output1="$fasta1/$file";
        push (@outfile1,$output1);

        my $output2="$fasta2/$file";
        push (@outfile2,$output2);
    }

    my $sub_file_number = scalar @input;
    my $thread;
    my @threads;
    my $job_number=0;

    while(){ 

        last if ($job_number eq $sub_file_number);  

        while(scalar(threads->list())<$thread_number) {     #set thread number
            $job_number++;    
            my $input_file = $input[$job_number-1];
            my $output_file1 = $outfile1[$job_number-1];
            my $output_file2 = $outfile2[$job_number-1];
            #print "Genbank_extraction: $job_number\n";
            $threads[$job_number-1]=threads->new(\&genomenucleotide_extract, $input_file, $output_file1, $output_file2); 

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

}


####################################################
####################################################
sub genomenucleotide_extract {

    my ($input, $output1, $output2)=@_;
    open (TEMP_IN, $input);
    my @temp_in = <TEMP_IN>;
    my $first_line = shift @temp_in;
    my $po_seq = join ("", @temp_in);
    open(OUTPUT_1, ">$output1");
    open(OUTPUT_2, ">$output2");
    my $filename = basename $input;
    if ($first_line =~ /^>/){                  #input is fasta format
        my %hash_seq = parse_fasta($input);
        my $jdna;
        foreach (keys %hash_seq){
            print OUTPUT_1 ">$_\n", print_sequence_into_file($hash_seq{$_}, 60), "\n";
            $jdna .=$hash_seq{$_};
        }
        print OUTPUT_2 ">$filename\n", print_sequence_into_file($jdna, 60); 
                 
    }else{                                     #input is genbank format
        my $DDD;
        my $in = new Bio::SeqIO( -file => $input, -format => 'genbank' );
        #my $out = new Bio::SeqIO( -file => $output1, -format => 'fasta' );
        while (my $seqObject = $in->next_seq()){
            my $acc = $seqObject->accession;
            my $desc = $seqObject->desc;
            my $dna = $seqObject->seq();
            my $id = $seqObject->display_id();
            print OUTPUT_1 ">$id\n", print_sequence_into_file($dna, 60), "\n";
            $DDD .=$dna;
       
        }
        $DDD=print_sequence_into_file($DDD, 60);  # print_sequence_into_file($dna, 70);  # 60 characters per line
        print OUTPUT_2 ">$filename\n$DDD";
    }
    close OUTPUT_1;
    close OUTPUT_2;

}


# subroutine print_sequence_into_file
# print $sequence into $FILE with $length character per line
sub print_sequence_into_file{   
    my($sequence, $length) = @_;
    my $DNA;
    my $string;
    # Print sequence in lines of $length
    for ( my $pos = 0 ; $pos < length($sequence) ; $pos += $length ){   
        $DNA=substr($sequence, $pos, $length)."\n";
        $string.=$DNA;
        
    }
    return $string;
}




##################################################################################################
###### Subrounting--do bacth work
###### Function:
###### do bacth work
##################################################################################################
sub batch_genenucleotide_TFT_extract {
    my ($path_genbank, $path_TFT, $path_gene, $thread_number, $prodigal_annotation, $path_fa)=@_;

    opendir PATH_GENBANK, $path_genbank or die "could not open $path_genbank";
    my @path_genbank = readdir PATH_GENBANK;
    @path_genbank =grep ($_!~/^\./ ,@path_genbank);  #delete hidden file . .. 
    closedir PATH_GENBANK;

    my $sub_file_number = scalar @path_genbank;

    my @input;
    my @outfile_2;
    my @outfile_3;
 
    foreach my $file_genbank(@path_genbank){ 
        $file_genbank =~ s/ /_/g;
        my $input = "$path_genbank/$file_genbank"; 
        push (@input,$input);  

        my $output_2 = "$path_TFT/$file_genbank.tbl";
        push (@outfile_2,$output_2);

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
            $threads[$job_number-1]=threads->new(\&genenucleotide_TFT_extract, $input_file, $output_file_2, $output_file_3, $prodigal_annotation, $path_fa); 

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
###### extract gene sequence and tbl (cds/rRNA/tRNA information) from genbank file
##################################################################################################
sub genenucleotide_TFT_extract {

    my ($input, $output_2, $output_3, $prodigal_annotation, $path_fa)=@_;
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
    my $sec_input = "$path_fa/$file_name";
    &prodigal_tool($sec_input, $output_3, $output_2, $prodigal_annotation) if $size == 0;
 
}

#########################################
# bacth annnotate genome using prodigal
#########################################
sub prodigal_bacth {
    my ($input_dir, $gene_dir, $tft_dir, $thread_number, $prodigal_annotation) = @_;
    opendir (DIR_SEQ, $input_dir) or die "could not open $input_dir";
    my @input_dir = readdir DIR_SEQ; 
    @input_dir = grep ($_!~/^\./ ,@input_dir);
    closedir DIR_SEQ;


    my @input;
    my @outgene;
    my @outtft;
    foreach(sort @input_dir){
        push @input, "$input_dir/$_";
        push @outgene, "$gene_dir/$_";    
        push @outtft, "$tft_dir/$_";    
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
        my $output_gene = $outgene[$job_number-1];
        my $output_tft = $outtft[$job_number-1];

        my $genome_annotation_percent = int (($job_number/$sub_file_number)*100); 
        print "$genome_annotation_percent%","..." if ($job_number == 1 or ( ($genome_annotation_percent%10 ==0) && ($progress_record <$genome_annotation_percent)) );

        $threads[$job_number-1]=threads->new(\&prodigal_tool, $input_file, $output_gene, $output_tft, $prodigal_annotation); 
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
    my ($input, $outgene, $outtft, $prodigal_annotation) =@_;
    $outgene =~ s/.fasta$// if $prodigal_annotation eq "F";
    $outtft =~ s/.tbl$//;
    my $prodigal = `which prodigal`;  # "/usr/bin/blastp";
    $prodigal =~ s/\n//g; 
    system ("$prodigal -i $input -f gbk -g 11 -q -p meta -d $outgene -o $outtft");
    &gene_file($outgene);
    &tbl_file($outtft);
}


#########################################
# gene file conversion 
#########################################
sub gene_file {
    my ($out) = @_;
    #$out =~ s/.fasta$//;
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


sub parse_fasta {

    my ($infile) = @_;
    my $term     = $/; # input record separator;
    my %seq_hash = (); # key:seqid, val:seq
    open my $INFILE, $infile or die "could not open infile '$infile' : $! \n"; 
    $/ = ">"; # input record separator;
    while(<$INFILE>) {
        chomp;
        next if($_ eq '');
        my ($id, @sequencelines) = split /\n/;
        foreach my $line (@sequencelines) {
            $seq_hash{$id} .= $line;
        }
    }
    $/ = $term;
    return(%seq_hash);

} # end of parse_fasta




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
    my $cmd = "$b -outfmt '6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore' -max_target_seqs 2 -evalue $e -num_threads 1 -query $input -db $db -out $outfile";
    system("$cmd 2>/dev/null"); #BLASTn -max_target_seqs -num_alignments 

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
    my $out_class = "$output_statistics/classification_CPS";

    open (OUT_1, ">$out1") or die "could not open infile $out1\n";  
    open (OUT_2, ">$out2") or die "could not open infile $out2\n";
    open (OUT_3, ">$out3") or die "could not open infile $out3\n";  
    open (OUT_4, ">$out4") or die "could not open infile $out4\n";   

    my $temp_total_tbl = "$result_statistics_workplace/temp_total_tbl";
    chdir $tbl_folder;
 
    #system ("cat * > $temp_total_tbl");
    system ("find . ! -name '.' | xargs cat > $temp_total_tbl");  # all files not be .; For large number of files, directly using cat may causing error


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    my %gene_hash;
    my %class_hash;
    open (TEMP_TBL, $temp_total_tbl);
    my @temp_tbl = <TEMP_TBL>;
    close TEMP_TBL;

    foreach (@temp_tbl){
        chomp;
        my @array = split "\t", $_;
        my $new_keys = "$array[1]|$array[0]"; # class|gene

        $gene_hash{$new_keys} = $array[1];
        $class_hash{$array[1]} = $array[0];
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
            print OUT_1 "$type_check\t";  # print class gene name
            print OUT_3 "$type_check\t";  # print class gene name

#        }else{
#           print OUT_1 "\t";  # print class gene name
#           print OUT_3 "\t";  # print class gene name

#        }
###################################################################################
    }

    print OUT_1 "\nStrain\t"; 
    print OUT_2 "Strain\t"; 
    print OUT_3 "\nStrain\t"; 
    print OUT_4 "Strain\t"; 

    foreach (sort keys %gene_hash) {
        $_ =~ s/(.*)\|//;
        print OUT_1 "$_\t";  # print class gene name
        print OUT_3 "$_\t";  # print class gene name
    }
    print OUT_1 "\n";
    print OUT_3 "\n";

    foreach (sort keys %class_hash) {
        print OUT_2 "$_\t";  # print class name
        print OUT_4 "$_\t";  # print class name
    }
    print OUT_2 "\n";
    print OUT_4 "Serotype/number/frequency\n"; # added by xiangyang 2022-02-13

###-------------------------------------------------------------------------------------------------------------------------
    opendir (TBL_FOLDER, $tbl_folder);
    my @total_protein = readdir TBL_FOLDER;
    @total_protein =grep ($_!~/^\./ ,@total_protein);  #delete hidden file . ..
    @total_protein =grep ($_!~/\~$/ ,@total_protein);  #delete temp file ~
    closedir TBL_FOLDER;
    
    my (@array_gene_num, @array_class_num);
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
            push (@arr_temp_each1, $sequencelines[0]);  #class_gene (gene)
            push (@arr_temp_each2, $sequencelines[1]);  #class (class)
            $hash_LocusTag1{$sequencelines[5]} = $sequencelines[0]; #keys: LocusTag; values: class_gene (gene)
            $hash_LocusTag2{$sequencelines[5]."|".$sequencelines[0]."|".$sequencelines[4]} = $sequencelines[1]; #keys: LocusTag; values: class (class) #added by xiangyang 2022-02-14
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
        foreach my $keyclass (sort keys %class_hash) {
            my $macth_num2 = grep ($keyclass eq $_, @arr_temp_each2);
            print OUT_2 "$macth_num2\t";
            push @array_num_temp2, $macth_num2;   # recode class number temply;

            my @macth_num2 = grep ($keyclass eq $_, @arr_temp_each2);
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
        push @array_class_num, $array_num_temp2[0]; 
    }

    @array_gene_num = sort{$b<=>$a}@array_gene_num;
    my $max_gene_number = $array_gene_num[0];

    @array_class_num = sort{$b<=>$a}@array_class_num;
    my $max_class_number = $array_class_num[0];

    #print "max_gene_number: $max_gene_number\n";
    #print "max_class_number: $max_class_number\n";
    close OUT_1;
    close OUT_2;
    close OUT_3;
    close OUT_4;
    system("rm $temp_total_tbl");
    my $tbl_heatmap_dir_gene = tbl_heatmap_generation($out1, $result_statistics_workplace, "tbl_heatmap_gene", $out_classs_gene, $max_gene_number);
    my $tbl_heatmap_dir_class = tbl_heatmap_generation($out2, $result_statistics_workplace, "tbl_heatmap_class", $out_class, $max_class_number);
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
                #$identity = sprintf "%.2f",$identity;
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
    my @Class_name = split "\t", shift @list_parse;

    open (OUT_CLASSIFICATION, ">$out_classification") or die "could not open infile $out_classification\n";  
    print OUT_CLASSIFICATION "$max_number\n";
####### caculate the number of gene number belonging to specific class
    my %hash_number;
    #print "number: ",  scalar @Class_name, "\n";
    foreach my $element_class (@Class_name){
        if ($element_class !~ /^\s|^Strain/){
            $hash_number{$element_class} = grep ($element_class eq $_, @Class_name);
        }
    }

    my %count;
    my @Class_name_new = grep { ++$count{ $_ } < 2; } @Class_name;
    
        foreach (@Class_name_new){
            print OUT_CLASSIFICATION "Class_type\t$_\t",  $hash_number{$_}, "\n" if $_ !~ /^\s|^Strain|^Gene_type/;
        }
    
#######caculate the number of gene number belonging to specific class

    my @gene_name = split "\t", shift @list_parse if $list =~ /Gene_number_statistics/;

    if ($list =~ /Gene_number_statistics/){
        foreach my $element_gene (@gene_name){
            if ($element_gene !~ /^\s|^Strain/){
                $hash_number{$element_gene} = grep ($element_gene eq $_, @gene_name);
            }
        }
        my $index_gene=0;
        foreach (@gene_name){
            print OUT_CLASSIFICATION "Gene_type\t$_\t",  $hash_number{$_}, "\n" if $_ !~ /^\s|^Strain/; 
            $index_gene++;
        }
    }
    close OUT_CLASSIFICATION;

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
            print { $fhs->{$source} } "$start\t$end\tCDS\t$Class_name[$length_count]\t$_\n" if $list !~ /Gene_number_statistics/;
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
    my $output_title = "Qseqid\tSseqid\tBitscore\tE-value\tPidenty\tQ_coverage\tS_coverage\tMacth_length";
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
        my $qcover =  ($blast_result[7] - $blast_result[6] +1)/$blast_result[3] * 100;  #added by xiangyang
        $qcover = sprintf("%.2f", $qcover);
        my $scover = ($blast_result[9] - $blast_result[8] +1)/$blast_result[4] * 100;  #added by xiangyang
        $scover = sprintf("%.2f", $scover);

        my $qcover_k =  ($blast_result_k[7] - $blast_result_k[6] +1)/$blast_result_k[3] * 100;  #added by xiangyang
        $qcover_k = sprintf("%.2f", $qcover);
        my $scover_k = ($blast_result_k[9] - $blast_result_k[8] +1)/$blast_result_k[4] * 100;  #added by xiangyang
        $scover_k = sprintf("%.2f", $scover_k);
 
            #set cut-off of e-value, identity, coverage, match_length for homologous protein
            #if ($blast_result[5]>0) {$blast_result[11] = $blast_result[11]/$blast_result[5];} else {$blast_result[11] = $blast_result[11];}
        my $iden_blast_result = $blast_result[2];
        $iden_blast_result = sprintf("%.2f", $iden_blast_result);
        my $iden_blast_result_k = $blast_result_k[2];
        $iden_blast_result_k = sprintf("%.2f", $iden_blast_result_k);

            if ($blast_result[0] eq $blast_result_k[0]) {
                if ($blast_result[11] <= $blast_result_k[11]) {
                    if ($blast_result[3] == $blast_result[4]){
                        $st_hash{$blast_result[0]} = "$blast_result[0]\t$blast_result[1]\t$blast_result[11]\t$blast_result[10]\t$iden_blast_result\t$qcover\t$scover\t$blast_result[5]\n";
                        $i = $i+1;

                    }else{
                        $st_hash{$blast_result_k[0]} = "$blast_result_k[0]\t$blast_result_k[1]\t$blast_result_k[11]\t$blast_result_k[10]\t$iden_blast_result_k\t$qcover_k\t$scover_k\t$blast_result_k[5]\n";
                        $i = $i+1;

                    }
                    
                }else{
                    $st_hash{$blast_result[0]} = "$blast_result[0]\t$blast_result[1]\t$blast_result[11]\t$blast_result[10]\t$iden_blast_result\t$qcover\t$scover\t$blast_result[5]\n";

                    $i = $i+1;
                }

            }else{               
                $st_hash{$blast_result[0]} = "$blast_result[0]\t$blast_result[1]\t$blast_result[11]\t$blast_result[10]\t$iden_blast_result\t$qcover\t$scover\t$blast_result[5]\n";    
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





sub check_tool {

    ###Test-step1: Checks for pneumo-typer dependencies

    my ($script, $dir) = @_;
    print "\n\nTest-step1: Checks for pneumo-typer dependencies...\n";
    print "################################################################\n";
    #Checks for perl modules
    &check_perl_modules;
    print "\n----------------------------------------------------------------\n";
    #Chechs for additional programs dependencies
    &check_programs;
    print "################################################################\n";

    ###Test-step2: Begin test pneumo-typer.pl with the test_data
    print "\n\n\n\n\nTest-step2: Begin test pneumo-typer.pl...\n";
    print "################################################################\n";
    my $check_g = system ("perl $script -d $dir -t 10 -m T -c T -p T");
    print "################################################################\n";
    if ($check_g eq 0){
        print "Ok, pneumo-typer works successfully!!\n\n";
    }else {
        print "Not Ok, pneumo-typer works with some errors!\n\n";
    }
}



#check Perl modules dependencies;
sub check_perl_modules {
 
    my @test_modules = ("GD", "GD::SVG", "SVG", "threads", "File::Basename", "FindBin", "lib", "Getopt::Long", "Math::BigFloat", "Storable", "vars", "File::Spec", "Bio::SeqIO", "Bio::Tree::NodeI", "Bio::TreeIO");

    my $check_number=0;
    foreach (@test_modules){
        my $cmd = "perl -M$_ -e 'print \"***$_ Version\\t \".$_->VERSION.\"\tok.\n\"'";
        #perl -MGD::SVG -e 'print GD::SVG->VERSION. "\n"'
        my $test_out = system ($cmd);

        print  "***Warning:\t$_ is not installed\n" if $test_out ne 0;
        $check_number++ if $test_out ne 0;
    }
    if ($check_number == 0){
        print "!!!Ok, all dependencies Perl modulers are installed*\n";
    }else{
        print "!!!Not ok, it appears that certain dependencies Perl modulers are installed*\n";
    }
}


# checks software dependencies
sub check_programs{

    print "Checking for makeblastdb ... ";
    my $makeblastdb_path = `which makeblastdb`;
    chomp $makeblastdb_path;
    if (not -e $makeblastdb_path){
      print "error: makeblastdb is not installed\n";
    }
    else{
      print "OK, makeblastdb is installed at: $makeblastdb_path\n";
    }
    
    print "Checking for blastn ... ";
    my $blastn_path = `which blastn`;
    chomp $blastn_path;
    if (not -e $blastn_path){
      print "error: blastn is not installed\n";
    }
    else{
      print "OK, blastn is installed at: $blastn_path\n";
    }
    
    print "Checking for prodigal ... ";
    my $prodigal_path = `which prodigal`;
    chomp $prodigal_path;
    if (not -e $prodigal_path){
      print "error: prodigal is not installed\n";
    }
    else{
      print "OK, prodigal is installed at: $prodigal_path\n";
    }


    print "Checking for blat ... ";
    my $blat_path = `which blat`;
    chomp $blat_path;
    if (not -e $blat_path){
        print "error: blat is not installed\n";
    }
    else{
        print "OK, blat is installed at: $blat_path\n";
    }

}


####################################################
####################################################
sub batch_genomenucleotide_extract_run {

    my ($path, $fasta1, $fasta2, $thread_number)=@_;

    opendir DIR, $path or die $!;
    ################################start rename the file in following case
    my @path_temp = readdir DIR;
    foreach (@path_temp){
        if (/[ =]/) {
            my $new_name = $_;
            $new_name =~ s/ /_/g;   #repalce the space with "_"  for all GenBank files
            $new_name =~ s/=/_/g;   #repalce "=" with "_"  for all GenBank files
            $new_name =~ s/_+/_/g;  
            system ("mv '$path/$_' '$path/$new_name'");  
        }

    }
    closedir DIR;
    #################################end rename the file in following case
    opendir DIR, $path or die $!;
    my @dir = readdir DIR; 
    @dir = grep ($_!~/^\./, @dir);
    closedir DIR;

    my $job_number=0;
    my $sub_file_number = scalar @dir;
    print "    Genomenucleotide_extract_percent: ";

    use Para::Runner;
    my $runner = Para::Runner->new($thread_number);

    foreach my $file(@dir){

        my $input="$path/$file"; 
        my $output1="$fasta1/$file";
        my $output2="$fasta2/$file";
        my $progress_record = int (($job_number/$sub_file_number)*100);
        $job_number++;
        my $Genomenucleotide_extract_progress = int (($job_number/$sub_file_number)*100); 
        print "$Genomenucleotide_extract_progress%","..." if ($job_number == 1 or ( ($Genomenucleotide_extract_progress%10 ==0) && ($progress_record <$Genomenucleotide_extract_progress)) );
        $runner->run(
	    sub{
                open (TEMP_IN, $input);
                my @temp_in = <TEMP_IN>;
                my $first_line = shift @temp_in;
                my $po_seq = join ("", @temp_in);
                open(OUTPUT_1, ">$output1");
                open(OUTPUT_2, ">$output2");
                my $filename = basename $input;
                if ($first_line =~ /^>/){                  #input is fasta format
                    my %hash_seq = parse_fasta($input);
                    my $jdna;
                    foreach (keys %hash_seq){
                        print OUTPUT_1 ">$_\n", print_sequence_into_file($hash_seq{$_}, 60), "\n";
                        $jdna .=$hash_seq{$_};
                    }
                    print OUTPUT_2 ">$filename\n", print_sequence_into_file($jdna, 60); 
                 
                }else{                                     #input is genbank format
                    my $DDD;
                    my $in = new Bio::SeqIO( -file => $input, -format => 'genbank' );
                    while (my $seqObject = $in->next_seq()){
                        my $acc = $seqObject->accession;
                        my $desc = $seqObject->desc;
                        my $dna = $seqObject->seq();
                        my $id = $seqObject->display_id();
                        print OUTPUT_1 ">$id\n", print_sequence_into_file($dna, 60), "\n";
                        $DDD .=$dna;       
                    }
                    $DDD=print_sequence_into_file($DDD, 60);  # print_sequence_into_file($dna, 70);  # 60 characters per line
                    print OUTPUT_2 ">$filename\n$DDD";
                }
                close OUTPUT_1;
                close OUTPUT_2;

            }
        )#run end

    }
    $runner->finish;
    print "done\n";
}



##################################################################################################
###### Subrounting--do bacth work
###### Function:
###### do bacth work
##################################################################################################
sub batch_genenucleotide_TFT_extract_run {
    my ($path_genbank, $path_TFT, $path_gene, $thread_number, $prodigal_annotation, $path_fa)=@_;

    opendir PATH_GENBANK, $path_genbank or die "could not open $path_genbank";
    my @path_genbank = readdir PATH_GENBANK;
    @path_genbank =grep ($_!~/^\./ ,@path_genbank);  #delete hidden file . .. 
    closedir PATH_GENBANK;

    my $job_number=0;
    my $sub_file_number = scalar @path_genbank;
    print "    Genenucleotide_TFT_extract_percent: ";
    use Para::Runner;
    my $runner = Para::Runner->new($thread_number);
 
    foreach my $file_genbank(@path_genbank){ 
        $file_genbank =~ s/ /_/g;
        my $input = "$path_genbank/$file_genbank";  
        my $output_2 = "$path_TFT/$file_genbank.tbl";
        my $output_3="$path_gene/$file_genbank.fasta";
        my $progress_record = int (($job_number/$sub_file_number)*100);
        $job_number++;
        my $Genenucleotide_TFT_extract_progress = int (($job_number/$sub_file_number)*100); 
        print "$Genenucleotide_TFT_extract_progress%","..." if ($job_number == 1 or ( ($Genenucleotide_TFT_extract_progress%10 ==0) && ($progress_record <$Genenucleotide_TFT_extract_progress)) );
        $runner->run(
	    sub{
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
                my $sec_input = "$path_fa/$file_name";
                &prodigal_tool($sec_input, $output_3, $output_2, $prodigal_annotation) if $size == 0;
 
            }

        )#run end

    }
    $runner->finish;
    print "done\n";
}



#########################################
# bacth annnotate genome using prodigal
#########################################
sub prodigal_bacth_run {
    my ($input_dir, $gene_dir, $tft_dir, $thread_number, $prodigal_annotation) = @_;
    my $prodigal = `which prodigal`;
    $prodigal =~ s/\n//g;

    opendir (DIR_SEQ, $input_dir) or die "could not open $input_dir";
    my @input_dir = readdir DIR_SEQ; 
    @input_dir = grep ($_!~/^\./ ,@input_dir);
    closedir DIR_SEQ;

    my $sub_file_number =  scalar @input_dir;
    my $job_number=0;

    print "    Annotating genome using prodigal: ";
    use Para::Runner;
    my $runner = Para::Runner->new($thread_number);

    foreach(sort @input_dir){
        my $input = "$input_dir/$_";
        my $outgene = "$gene_dir/$_";    
        my $outtft = "$tft_dir/$_"; 
        my $progress_record = int (($job_number/$sub_file_number)*100);
        $job_number++;
        my $genome_annotation_percent = int (($job_number/$sub_file_number)*100); 
        print "$genome_annotation_percent%","..." if ($job_number == 1 or ( ($genome_annotation_percent%10 ==0) && ($progress_record <$genome_annotation_percent)) );
        $runner->run(
	    sub{
                $outgene =~ s/.fasta$// if $prodigal_annotation eq "F";
                $outtft =~ s/.tbl$//; 
                system ("$prodigal -i $input -f gbk -g 11 -q -p meta -d $outgene -o $outtft");
                &gene_file($outgene);
                &tbl_file($outtft);
            }
        )#run end   
    }
    $runner->finish;
    print "done\n";

}


#########################################################
#########################################################
sub bacth_blast_best_run {

    my ($workplace, $path_gene, $blastn_out_dir, $db_file, $evalue, $blastn, $makeblastdb, $thread_number) = @_;

    use File::Basename qw<basename dirname>;
    my $db = $workplace."/".basename($db_file);
    system ("cp $db_file $db");
    system ("$makeblastdb -in $db -dbtype nucl > $workplace/temp.txt");

    opendir (GENE_FOLDER, $path_gene);
    my @path_gene = readdir GENE_FOLDER;
    @path_gene = grep ($_!~/^\./, @path_gene);
    closedir GENE_FOLDER;

    my $sub_file_number =  scalar @path_gene;
    my $job_number=0;
    print "    Blastn_percent: ";
    use Para::Runner;
    my $runner = Para::Runner->new($thread_number);

    foreach my $eachfile (@path_gene) {
        my $input = "$path_gene/$eachfile";  
        my $outfile = "$blastn_out_dir/$eachfile.blastout";
        my $progress_record = int (($job_number/$sub_file_number)*100);
        $job_number++;
        my $Blastn_percent = int (($job_number/$sub_file_number)*100); 
        print "$Blastn_percent%","..." if ($job_number == 1 or ( ($Blastn_percent%10 ==0) && ($progress_record <$Blastn_percent)) );
        $runner->run(
            sub{
                #****** Warining: this subrouine is only used to find the best hit for each query gene against database, and a loose cutoff for evalue
                my $cmd = "$blastn -outfmt '6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore' -max_target_seqs 2 -evalue $evalue -num_threads 1 -query $input -db $db -out $outfile";
                system("$cmd 2>/dev/null"); #not print warning on terminal
               }
        )#run end
    }
    $runner->finish;
    print "done\n";

}



1;               
__END__  
