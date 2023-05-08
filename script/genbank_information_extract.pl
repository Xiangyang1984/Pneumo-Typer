#!/usr/bin/perl -w

use strict;
use threads;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename qw<basename dirname>;
use FindBin;

# perl /home/xiangyang/Mazhongrui/Ma_genbank_information_extract.pl -dir /home/xiangyang/Mazhongrui/surppment_SP/GBK_folder -m 190 -g /home/xiangyang/Mazhongrui/surppment_SP/metadata_workplace



my %options = (
    'genbank_file_directory=s'                 => undef,   
    'multiple_threads=s'                       => "1",
    'gi_out_dir=s'                             => undef

);

GetOptions(
    'dir|genbank_file_directory=s'             => \$options{genbank_file_directory},    
    'm|multiple_threads=s'                     => \$options{multiple_threads},
    'g|gi_out_dir=s'                           => \$options{gi_out_dir}


);

my $path_genbank = $options{genbank_file_directory};
my $thread_number = $options{multiple_threads};
my $gi_out_dir = $options{gi_out_dir};
mkdir $gi_out_dir;

my $contain_temp = "$gi_out_dir/metadata_temp_workplace";
mkdir $contain_temp;



opendir PATH_GENBANK, $path_genbank or die "could not open $path_genbank";

my @path_genbank = readdir PATH_GENBANK;
@path_genbank =grep ($_!~/^\./ ,@path_genbank);  #delete hidden file . .. 
#@path_genbank =grep ($_!~/$\~/ ,@path_genbank);  #delete hidden file . .. 
closedir PATH_GENBANK;

my $sub_file_number = scalar @path_genbank;


my @input;
my @outfile;
 
foreach my $file_genbank(@path_genbank){ 
    my $input="$path_genbank/$file_genbank"; 
    push (@input,$input);  
    @input = sort @input;  

    my $output="$contain_temp/$file_genbank.list";
    push (@outfile,$output);
    @outfile = sort @outfile;

}

    my $thread;
    my @threads;
    my $job_number=0;

while(){ 

    last if ($job_number eq $sub_file_number);  

    while(scalar(threads->list())<$thread_number) {     #set threadnumber；
        $job_number++;    
        my $input_file = $input[$job_number-1];
        my $output_file = $outfile[$job_number-1];
        print "$job_number\n";
        $threads[$job_number-1]=threads->new(\&genbank_protein_sequence_extract, $input_file, $output_file); 
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
 

chdir $contain_temp;
system("cat *.list > all_information.list");

my %hash_derepeat;
open(LIST, "all_information.list");
while(<LIST>){
    chomp;
    my @temp_array = split ("\t", $_);
    $hash_derepeat{$temp_array[0]} =$_;
}
close LIST;


open(FINAL_LIST, ">$gi_out_dir/genomes_metadata.txt");
    print FINAL_LIST "Assembly\tAcc\tOrganism\tStrain\tSerovar\tIsolation_source\tCountry\tHost\tCollection_date\n";

foreach(keys %hash_derepeat){
    print FINAL_LIST $hash_derepeat{$_},"\n"; 
}
close FINAL_LIST;

system("rm -rf $contain_temp");

sub genbank_protein_sequence_extract {
    my ($input, $output_1)=@_;

    open(OUTPUT_1, ">$output_1");

        my $in = new Bio::SeqIO( -file => $input, -format => 'genbank' );

        while ( my $seqObject = $in->next_seq() ) {
            my $acc = $seqObject->accession;
            my @features = $seqObject->get_SeqFeatures();
	    my $org = "";
	    my $strain = "";
	    my $host = "";
	    my $country = "";
            my $isolation_source = "";
	    my $collection_date = "";
            my $count = 0;       
            my $serovar;

            foreach (@features) {
                my $feat = $_;
#####################################################################################	 
 
		if ($feat->primary_tag() eq "source"){
                    if ( $feat->has_tag('organism') ) {
				($org) = $feat->get_tag_values('organism');
                    }else {$org = "blank";}

                    if ( $feat->has_tag('strain') ) {
				($strain) = $feat->get_tag_values('strain');
                    }else {$strain = "blank";}

                    if ( $feat->has_tag('isolation_source') ) {
				($isolation_source) = $feat->get_tag_values('isolation_source');
                    }else {$isolation_source = "blank";}

                    if ( $feat->has_tag('country') ) {
				($country) = $feat->get_tag_values('country');
                    }else {$country = "blank";}

                    if ( $feat->has_tag('collection_date') ) {
				($collection_date) = $feat->get_tag_values('collection_date');
                    }else {$collection_date = "blank";}

                    if ( $feat->has_tag('host') ) {
				($host) = $feat->get_tag_values('host');
                    }else {$host = "blank";}
                    if ( $feat->has_tag('serovar') ) {
				($serovar) = $feat->get_tag_values('serovar');
                    }else {$serovar = "blank";}

				next;
		} 
#####################################################################################	

                unless ($feat->primary_tag =~ /^CDS/) {
                    next;
                }
                    if ( $feat->has_tag('pseudo') ) {

                    }else {$count++;}
                    
                    #print "genome_$count\n";

            }
                    my $strain_name_l = basename($input);
                    $strain_name_l =~ s/_genomic.gbff//g;      
                    $strain_name_l =~ s/\-/_/g;                    
                    $acc =~ s/0000.*// if $acc =~ /0000/;
                    print OUTPUT_1 "$strain_name_l\t$acc\t$org\t$strain\t$serovar\t$isolation_source\t$country\t$host\t'$collection_date'\n";
        }

    close OUTPUT_1;

}




