#!/usr/bin/perl -w

use strict;
use threads;
use Getopt::Long;
use File::Basename qw<basename dirname>;
use File::Spec;

# The function of sEach_pTyper_list.pl is to group genomes based on a list file and then run pneumo-typer.pl on each grouped set of genomes.
# list file, each row contain a genome and its group information (e.g. serotype) separated by tab
# usage: perl $home_directory/script/sEach_pTyper_list.pl -f $file_dir -d $final_dir -m htread_number -list list_file

my %options = (
    'file_dir'                               => undef,
    'final_directory'                        => undef,   
    'list'                                   => undef,
    'cp_mv'                                  => "cp",   
    'multiple_threads'                       => "3"

);

GetOptions(
    'f|file_dir=s'                             => \$options{file_dir},
    'd|final_directory=s'                      => \$options{final_directory},    
    'l|list=s'                                 => \$options{list}, 
    'cm|cp_mv=s'                               => \$options{cp_mv},   
    'm|multiple_threads=i'                     => \$options{multiple_threads}

);

my $script_path = File::Spec->rel2abs($0);
my $pneumotyper_dir = dirname(dirname $script_path);
print "pneumotyper_dir\t$pneumotyper_dir\n";
my $file_dir = File::Spec->rel2abs($options{file_dir});

my $final_directory = File::Spec->rel2abs($options{final_directory});
mkdir $final_directory;

#print "$final_directory\n";
mkdir "$final_directory/FIG";

my $thread_number = $options{multiple_threads};


my $list_file = File::Spec->rel2abs($options{list});
#print "$list_file\n";

chdir $file_dir;
open (LIST, $list_file);
my %input;
 
while (<LIST>){
    chomp $_;
    next unless $_ !~ /^\s/;
    #next unless $_ =~ /^10D/;
    my @tmparray = split "\t", $_;
    #next unless $tmparray[1]  =~ /^22F|^33F|^35D/;
    push (@{$input{$tmparray[1]}}, $tmparray[0]);  
   
}


foreach (keys %input) {

    print "Subgroup: $_\n";
    my @input = @{$input{$_}};

    my $subfinal_directory = "$final_directory/$_";
    mkdir $subfinal_directory;

    my $sub_file_number =  scalar @input;
    print "genome number: $sub_file_number";
    my $thread;
    my @threads;
    my $job_number=0;

    while(){ 

        last if ($job_number eq $sub_file_number);  

        while(scalar(threads->list())<$thread_number) {     #set threadnumber；
            $job_number++;    
            my $input_file = $input[$job_number-1];
 
            
            $threads[$job_number-1]=threads->new(\&copy_file, $input_file, $subfinal_directory, $options{cp_mv}); 
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
    print "\n\n";
    system ("perl $pneumotyper_dir/pneumo-typer.pl -d $subfinal_directory -t $thread_number -m T -o $subfinal_directory.OUT -p T");
    
    system("cp $subfinal_directory.OUT/CPS_operon.svg $final_directory/FIG/$_.n$sub_file_number.CPS_operon.svg");
    system("cp $subfinal_directory.OUT/heatmap_class.svg $final_directory/FIG/$_.n$sub_file_number.heatmap_class.svg");
    system("cp $subfinal_directory.OUT/heatmap_gene.svg $final_directory/FIG/$_.n$sub_file_number.heatmap_gene.svg");
}
    #chdir $final_directory;
    #system("tar czf FIG.tar.gz FIG");
    #system("perl /home/xiangyang/Desktop/autorun.pl FIG.tar.gz");

sub copy_file {
    my ($input, $final_dir, $cp_mv) =@_;
 
    system ( "cp  $input -t $final_dir") if $cp_mv eq "cp";
    system ( "mv  $input -t $final_dir") if $cp_mv eq "mv";

}

