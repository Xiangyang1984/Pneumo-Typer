package SVG_heatmap;

use strict;
use warnings;
use Bio::SeqIO;
use File::Basename qw<basename dirname>;
use GD::SVG;
use SVG;
use FindBin;
use lib "$FindBin::Bin/lib";
use Bio::TreeIO;
use Bio::Tree::NodeI;
use font_rotate;
use process;
use metrics;

use vars qw($image %hash_nif %color_hash %cluster_color_hash %cluster_color_hash_2 %direction_F_R %shift_distance_x $yellow $green $blue $black $black2 $red $white $gray $dgray $tree_width @array_image %cluster_GeneName_hash @rgb_color %tab_hash $width_tab_file $Gene_height);

#create a svg image 

sub create_image_SVG {

    my ($home_directory, $directory_part_TFT, $strain_reorder_file, $phylogenetic_file, $show_tree_topology, $workplace, $ref_options, $heatmap_name_suffix) = @_;
    
    my %options = %$ref_options;

##############################################
#color bar create
##############################################
    open (TYPE, $options{classification_file});
    my @type_classification = <TYPE>;
    my $max_number = shift @type_classification;
    $max_number =~ s/\n//g; 
    color_spectrum ($max_number); #$max_number

#start to adjust the height by modify the up_shift when add the classname, line and gene name
##############
    my (@Class_type_length, @Gene_type_length);    
    foreach (@type_classification){
        my($type_gene, $gene_name, $gene_number) = split "\t", $_;
        push @Class_type_length, metrics::string_width_svg($gene_name, "$options{font_family}:$options{font_style}", $options{label_font_size}) if $type_gene eq "Class_type";
        push @Gene_type_length, metrics::string_width_svg($gene_name, "$options{font_family}:$options{font_style}", $options{label_font_size}) if $type_gene eq "Gene_type";
    } 
    @Class_type_length = sort{$a<=>$b}@Class_type_length;
    @Gene_type_length = sort{$a<=>$b}@Gene_type_length if $options{classification_file} =~ /classification_gene/;
    my $constant =3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068/180;
    my $degree = "-$options{rotate_gene_label}";
    my $Class_height = $Class_type_length[-1]*abs(sin($degree*$constant));
    #$Gene_height = 0;
    $Gene_height = $Gene_type_length[-1]*abs(sin($degree*$constant)) if $options{classification_file} =~ /classification_gene/;
    my $add_height;
    $add_height = $Class_height + 10 + $Gene_height if $options{classification_file} =~ /classification_gene/;  #class name+line+gene name 
    $add_height = $Class_height + 10 if $options{classification_file} !~ /classification_gene/;  #class name+line
    $options{up_shift} = $options{up_shift} + $add_height;

##############
#end to adjust the height

    #caculate the height size
    opendir DIR_PART_TFT, $directory_part_TFT or die $!;
    my  @dir_part_TFT_0 = readdir DIR_PART_TFT; 
    closedir DIR_PART_TFT;

    foreach (@dir_part_TFT_0){ 
        if (/~/) {
            system ("rm $directory_part_TFT/$_");            
        }
    }

    opendir DIR_PART_TFT, $directory_part_TFT or die $!;  # reopen $directory_part_TFT to caculate the number of sub-TFT file
    my  @dir_part_TFT = readdir DIR_PART_TFT; 
    @dir_part_TFT =grep ($_!~/^\./ ,@dir_part_TFT);  #delete hidden file . ..
    closedir DIR_PART_TFT;

    my $image_height=$options{distance_between_two_genomes}*(scalar @dir_part_TFT)+$options{up_shift}+$options{down_shift}+50;  # image_height 2019-06-20

    #subroutine to retrieve the lengthest width from FTL table files
    @array_image = process::figure_size_width_heatmap ($home_directory, $directory_part_TFT, $options{font_family}, $options{font_style}, $options{strain_name_font_size});  
    #%direction_F_R = %{$array_image[2]}; 
    #%shift_distance_x = %{$array_image[3]}; 
    my $temp_strain_reorder_file; 
    my $tree_strain_number=0;
    if (defined $options{phylogenetic_file} && $options{show_tree_topology} eq "T") {
        ($tree_width, $tree_strain_number, $temp_strain_reorder_file) = tree_width($phylogenetic_file, $options{show_tree_branch}, $options{left_shift}, $options{x_step}, $workplace);# caclulate the width and height of the tree
        $image_height=$options{distance_between_two_genomes}*$tree_strain_number+$options{up_shift}+$options{down_shift}+50;  # image_height 2019-06-20
    }else{

        $tree_width = $options{left_shift};
    }
    
    
    my $strain_number = 0;
    if (defined $options{strain_reorder_file}) {
        open (STRAIN_REORDER, $strain_reorder_file) or die "Must supply strain_reorder_file\n";
        
        while (<STRAIN_REORDER>){
            $strain_number++;
        }
        close STRAIN_REORDER;
        $image_height=$options{distance_between_two_genomes}*$strain_number+$options{up_shift}+$options{down_shift}+50;  # image_height 2019-06-20
    } #$strain_reorder_file

    if (defined $options{expand_tab}){

        my %length_str;
        open (T, $options{expand_tab});
        <T>;
        while(<T>){
            chomp;
            my @tab = split /\t/, $_;
            #print "$tab[0]\t$tab[1]\n";
            my $expand_genome = shift @tab;
            
            my $ccnum=0;
            foreach (@tab){
                push @{$length_str{$ccnum}}, metrics::string_width_svg($_, "$options{font_family}:$options{font_style}", $options{Ename_font_size});

                $ccnum++;
            }

            $tab_hash{$expand_genome} = $_;
            #print "$expand_genome\t$tab_hash{$expand_genome}\n";
        }
        close T;
        
        my @length_str_max;
        foreach (sort keys %length_str){

            my @length_str_temp = @{$length_str{$_}};           

            @length_str_temp = sort{$a<=>$b} @length_str_temp;
            push @length_str_max, $length_str_temp[-1];
            #print "$_\t", $length_str_temp[-1], "\n";
        }           
        
        my $tab_width;
        foreach (@length_str_max){
            $tab_width += $_ + 5;
            #print "$_\n";
        }
        #$tab_width = $tab_width - 5;
        #print "$tab_width\n";
        foreach (keys %tab_hash){
            my @tt = split /\t/, $tab_hash{$_};
            my $expand_shown;
            my $cccnum=0;
            shift @tt; #the first element is the genome name
            foreach (@tt){
                my $substr_len = $length_str_max[$cccnum];
                my $substr_text = $tt[$cccnum];
                $expand_shown .= "$substr_text"."#"."$substr_len\t";
                #print "$expand_shown", length ($expand_shown), "\n";
                $cccnum++;
            }
            $expand_shown =~ s/\t$//g; # delete the last tab
            $tab_hash{$_} = $expand_shown;
        }
        #$width_tab_file = &string_svg_width(\%tab_hash, $options{Ename_font_size}, $options{font_family}, $options{font_style});
        $width_tab_file = $tab_width+5;
    }else{
        $width_tab_file = 0;
    }
    

    # $array_image[2]: genome width; $tree_width=tree width + left shift; $width_tab_file: expand context width
    my $image_width = $array_image[0]*$options{figure_Scale_up_multiple}/20+$options{right_shift}+$array_image[2]+$tree_width+$width_tab_file;  # image_height 2019-06-26

    #print "\nStep 5(-SVG): Be going to create a size of ", $image_width, "x", "$image_height SVG image\n";

    #print "\nStep 6(-SVG): SVG format was chosen to map\n";

    #create a width x height size image 
        $image=GD::SVG::Image->new($image_width, $image_height);
  

        $yellow = $image->colorAllocate(255,255,0);
        $green = $image->colorAllocate(0,100,0);  
        $blue = $image->colorAllocate(0,0,255);
        $black2= $image->colorAllocate (255,0,0);
        $black= $image->colorAllocate (0,0,0);

        $red = $image->colorAllocate(255,0,0);
        $white= $image->colorAllocate(255,255,255);
        $gray= $image->colorAllocate (212,212,212);
        $dgray= $image->colorAllocate (153,153,153);
        #$ddgray   = $image->colorAllocate(127,127,127);

        $image->setThickness($options{line_drawing_width}); #set the line drawing width
        $image->filledRectangle(0, 0, $image_width, $image_height, $white); #add a white background color in image

    #%cluster_GeneName_hash = process::cluster_GeneName_hash($directory_part_TFT, $cluster_result);

    #print "\nStep 7(-SVG): SVG Figure will be generated, please wait!!!\n";
    # each tbl file is used as input file to generate map using draw_map subroutine 

    if (defined $options{phylogenetic_file}) {
        &TreeDraw($workplace, $phylogenetic_file, $options{show_tree_topology}, $options{show_tree_branch}, $options{show_tree_bootstrap}, $options{x_step}, $options{left_shift}, $options{up_shift}, $options{distance_between_two_genomes}, $options{font_family}, $options{strain_name_font_size});  # display the tree
        
    }

    my $eachfile=0;
    if (defined $options{strain_reorder_file} || defined $options{phylogenetic_file}) {
        my %strain_reorder_hash;
        if (defined $options{phylogenetic_file}) {open (STRAIN_REORDER, $temp_strain_reorder_file) or die "Must supply strain_reorder_file\n";} #$temp_strain_reorder_file
        if (defined $options{strain_reorder_file}) {open (STRAIN_REORDER, $strain_reorder_file) or die "Must supply strain_reorder_file\n";} #$strain_reorder_file

        while (<STRAIN_REORDER>){
            chomp;
            my @strain_reorder = split '\t', $_;
            $strain_reorder_hash{$strain_reorder[1]} = $strain_reorder[0];

        }

        foreach my $part_reorder (sort {$a<=>$b} keys %strain_reorder_hash){ 
            

            my $check_value = $strain_reorder_hash{$part_reorder};
            
            foreach my $tbl_part_file(sort @dir_part_TFT){ 
                my $tbl_file_temp = $tbl_part_file;
                if(($tbl_file_temp =~ /$check_value/) and ($tbl_file_temp =~ /tbl$/)) { # delete temp file generated when part_TFT file were edited to modify the genename
                    $tbl_file_temp =~ s/(.tbl)//g;
 
                    if ($tbl_file_temp eq $check_value) {
                        #$tbl_part_file .=$1; 

                        $eachfile++;
                        my $Y_parameter=$options{up_shift}+$options{distance_between_two_genomes}*$eachfile;    # strart Y 
                        my $tbl_part_infile = "$directory_part_TFT/$tbl_part_file" or die "Must supply input filename\n";
                        
                        draw_map_SVG ($tbl_part_infile, $Y_parameter, $tree_width, $options{arrow_relative_Length}, $options{arrow_relative_Height}, $options{strain_name_shift_Y},$options{gene_label_shift_Y}, $options{font_family}, $options{font_style}, $options{label_font_size}, $options{label_font_color}, $options{interested_gene_label_font_color}, $options{cds_color_border}, $options{pseudo_color_border}, $options{RNA_color_border}, $options{gene_no_color_filled}, $options{figure_Scale_up_multiple}, $options{rotate_gene_label}, $options{strain_name_font_size}, $options{strain_name_font_color}, $options{show_label}, $options{unification_label}, $eachfile, $options{Ename_font_size}, $options{Ename_font_color});

                    }

                 }
    
             }

        }

    }

    if (!defined $options{phylogenetic_file} && !defined $options{strain_reorder_file}) {
        foreach my $tbl_part_file(sort @dir_part_TFT){
            
            $eachfile++;

            my $Y_parameter=$options{up_shift}+$options{distance_between_two_genomes}*$eachfile;     # left_shift
        
            my $tbl_part_infile = "$directory_part_TFT/$tbl_part_file" or die "Must supply input filename\n";
            

            draw_map_SVG ($tbl_part_infile, $Y_parameter, $options{left_shift}, $options{arrow_relative_Length}, $options{arrow_relative_Height}, $options{strain_name_shift_Y}, $options{gene_label_shift_Y}, $options{font_family}, $options{font_style}, $options{label_font_size}, $options{label_font_color}, $options{interested_gene_label_font_color}, $options{cds_color_border}, $options{pseudo_color_border}, $options{RNA_color_border}, $options{gene_no_color_filled}, $options{figure_Scale_up_multiple}, $options{rotate_gene_label}, $options{strain_name_font_size}, $options{strain_name_font_color}, $options{show_label}, $options{unification_label}, $eachfile, $options{Ename_font_size}, $options{Ename_font_color});
    
        }
    }

#add class type for "tbl_heatmap_gene"
###########################################################################################################
my $tbl_dirname = basename $directory_part_TFT;
if ($tbl_dirname eq "tbl_heatmap_gene"){
    my ($type_class, $class_name, $class_number);
    my $left_number = 0;
    $image->setThickness($options{line_drawing_width}); #set the line drawing width
    foreach(@type_classification){
        chomp;
        next unless $_ =~ /^Class_type/;
        ($type_class, $class_name, $class_number) = split "\t", $_;
        my $head_start = $left_number*50*$options{figure_Scale_up_multiple}/20+50*$options{figure_Scale_up_multiple}/20*0.1+$tree_width+15+$array_image[2]+$width_tab_file;

        $left_number += $class_number; 
        my $head_end = $left_number*50*$options{figure_Scale_up_multiple}/20-50*$options{figure_Scale_up_multiple}/20*0.1+$tree_width+15+$array_image[2]+$width_tab_file;
        $image->line($head_start,$options{up_shift}-$Gene_height-10,$head_end,$options{up_shift}-$Gene_height-10,$red);
        $image->line($head_start,$options{up_shift}-$Gene_height-5,$head_start,$options{up_shift}-$Gene_height-10,$red);
        $image->line($head_end,$options{up_shift}-$Gene_height-5,$head_end,$options{up_shift}-$Gene_height-10,$red);
        $image->stringrotate(($head_end+$head_start+1)/2, $options{up_shift}-$Gene_height-11, $options{font_family}, $options{font_style}, $options{label_font_size}, $class_name, $blue, -$options{rotate_gene_label});
    }

}
   

#generate color_bar;
###########################################################################################################
my $length_count=0;
        foreach(@rgb_color){
            my $start =$options{left_shift}+5+10*$length_count;
            my $end = $options{left_shift}+5+10*($length_count+1);
            my @final_rgb0 = split ",", $rgb_color[$length_count];
            my $temp_color = $image->colorAllocate($final_rgb0[0], $final_rgb0[1], $final_rgb0[2]);
            #unless ($length_count<15){ next; }
            $length_count++;
            #print $_, "\t$length_count\n";
            $image->filledPolygon(rec_arrow_F($start, $end, 0, $options{arrow_relative_Height}*1.5, $image_height-40),$temp_color);
            #$image->polygon(rec_arrow_F($start+750,$end+750,0,10, 200), $temp_color);
            $image->stringrotate(($start+$end)/2, $image_height-40+$options{arrow_relative_Height}*1.5+$options{label_font_size}, $options{font_family}, $options{font_style}, $options{label_font_size}, $length_count-1, $black, 0);
        }
       $image->stringrotate($options{left_shift}+5, $image_height-40+$options{arrow_relative_Height}*1.5+$options{label_font_size}*2, $options{font_family}, $options{font_style}, $options{label_font_size}, "Gene number", $red, 0);
###########################################################################################################

    open (SVG_MAP, ">$workplace/heatmap_$heatmap_name_suffix.svg") or die("Cannot open file for writing");

    binmode SVG_MAP;

    print SVG_MAP $image->svg();
        
}


##################################################################################################
###### Subrounting--color_product
###### Function:
###### translate color (rgb format) into hash 
##################################################################################################
sub color_product {

    my $color_rgb_code = shift;
    my $line=0;
    open (COLOR, $color_rgb_code) or die "Can't read file"; # limit the color number into 1488

    while(<COLOR>){
        chomp;
        if ($_ !~ /^#/) {
            $line++;
            my @color=split(" = ", $_);      
        
            my @gbR=split (",", $color[1]); 

            $color_hash{$line}=$image->colorAllocate($gbR[0],$gbR[1],$gbR[2]);# limit the color number into 1488
        }

    }

}


##################################################################################################
###### Subrounting--draw_map_SVG
###### Function:
###### main subrouting to draw gene cluster map
##################################################################################################
sub draw_map_SVG {
    my %given_color_hash = (
        'blue'                 => $blue, 
        'red'                  => $red,
        'black'                => $black,
        'white'                => $white,
        'gray'                 => $gray,
        'dgray'                => $dgray
       );

    my $inputfile = shift;
    my $Y_centre = shift;       # position in Y axis
    my $left_shift = shift;
    my $arrow_length = shift;    # arrow length
    my $arrow_high = shift;       # arrow width
    my $strain_name_shift_Y = shift;
    my $gene_label_shift_Y = shift;
    my $font_family = shift;
    my $font_style = shift;
    my $label_font_size = shift;
    my $label_font_color = shift;
    my $interested_gene_label_font_color = shift;
    my $cds_color_border = shift;
    my $pseudo_color_border = shift;
    my $RNA_color_border = shift;
    my $gene_no_color_filled = shift;
    my $figure_Scale_up_multiple = shift;  # figure size Scale up multiple
    my $rotate_gene_label = shift;
    my $strain_name_font_size = shift;
    my $strain_name_font_color = shift;
    my $show_label = shift;    
    my $unification_label = shift;
    my $eachfile = shift;                              #added on 2020-11-21
    my $filename = basename $inputfile;
    my $filename_shift = basename $inputfile;
    my $Ename_font_size = shift;
    my $Ename_font_color = shift;
    $filename =~ s/.tbl(.*)//g;
    #$filename =~ s/_/ /g;
    my $strain_name_height = metrics::string_height_svg("$font_family:$font_style", $strain_name_font_size);
    $image->stringrotate($tree_width+2.5, $Y_centre+$strain_name_height/4-$strain_name_shift_Y, $font_family, $font_style, $strain_name_font_size, $filename, $given_color_hash{$strain_name_font_color}, 0); # genomic names
#print "ttttt: stringrotate\n";

#my $head_start = $tree_width+15+$array_image[2]+$width_tab_file;
#my $head_end   = $tree_width+15+$array_image[2]+$width_tab_file;
#$image->stringrotate(($head_end+$head_start+1)/2, $options{up_shift}-40, $options{font_family}, $options{font_style}, $options{label_font_size}, $class_name, $blue, -$options{rotate_gene_label});




    #print "$inputfile\n";
    if (defined $tab_hash{$filename}){

        my @tab_text = split /\t/, $tab_hash{$filename};
	my @add_head;
        @add_head = ("Serotype") if scalar @tab_text ==1;
        @add_head = ("Serotype", "ST") if scalar @tab_text ==2;
        @add_head = ("Serotype", "ST", "cgST") if scalar @tab_text ==3;
        my $extand_rigth_distance = 0;
        my $ccccnum = 0;
        foreach (@tab_text){
            my ($text, $max_width) = split /\#/, $_;
            $image->stringrotate($tree_width+10+$array_image[2]+5+$extand_rigth_distance, $Y_centre+$strain_name_height/4-$strain_name_shift_Y, $font_family, $font_style, $Ename_font_size, $text, $given_color_hash{$Ename_font_color}, 0); # if defined $tab_hash{$filename}; # tab names
            $image->stringrotate($tree_width+10+$array_image[2]+5+$extand_rigth_distance+$max_width/2, $Y_centre-$arrow_high-$gene_label_shift_Y, $font_family, $font_style, $Ename_font_size, $add_head[$ccccnum], $given_color_hash{$Ename_font_color}, -$rotate_gene_label) if $eachfile == 1; # Serotype, ST and cgST
            $extand_rigth_distance += $max_width+5;
            $ccccnum++;
        }
    }

    open (ABL, $inputfile);

    my @min=();

    while (<ABL>){
        my @arr= split "\t", $_;
        $arr[0]=~ s/[><]//g;
        $arr[1]=~ s/[><]//g;
        push (@min,@arr[0,1]);
    }
    @min = sort{$a<=>$b} @min;

    seek (ABL,0,0);
    #my $highlight_gene;
    while (<ABL>){
        chomp;
        my @cds= split "\t", $_;
        $cds[0]=~ s/[><]//g;
        $cds[1]=~ s/[><]//g;

        my %geneName_change_hash;
        my $tag;
        if ($cds[3] =~ /;/) {
            my @geneName_array = split ";", $cds[3];
            $geneName_change_hash{$geneName_array[1]} = $geneName_array[0];
            $geneName_change_hash{$geneName_array[1]} = $cluster_GeneName_hash{$geneName_array[1]} if ( (defined $cluster_GeneName_hash{$geneName_array[1]}) && ($unification_label eq "T") ); #uniform gene name for homologs
            #$geneName_change_hash{$geneName_array[1]} = $geneName_array[0] if (!defined $cluster_GeneName_hash{$geneName_array[1]});
            $tag = $geneName_array[1]; 
        }else {
            $geneName_change_hash{$cds[3]} = $cds[3];
            $geneName_change_hash{$cds[3]} = $cluster_GeneName_hash{$cds[3]} if ( (defined $cluster_GeneName_hash{$cds[3]}) && ($unification_label eq "T") ); #uniform gene name for homologs            
            #$geneName_change_hash{$cds[3]} = $cds[3] if (!defined $cluster_GeneName_hash{$cds[3]});
            $tag = $cds[3];
        }
	
        my ($strart, $end);
 
        $strart=1+($cds[0]-$min[0])*$figure_Scale_up_multiple/20+$left_shift+15+$array_image[2]+$width_tab_file;
        $end  = $cds[1]*$figure_Scale_up_multiple/20+$left_shift+15+$array_image[2]+$width_tab_file;

        ### gene_label mark
        if ($show_label eq "T") {

                  
           #$tag =~ s/;.*//g; #delete locus_Tag
            if ($eachfile == 1){
                $image->stringrotate(($strart+($end-$strart)/2), $Y_centre-$arrow_high-$gene_label_shift_Y, $font_family, $font_style, $label_font_size, $geneName_change_hash{$tag}, $given_color_hash{$label_font_color}, -$rotate_gene_label); # Default color is dgray for gene label of the other genes
            }

        }

        # CDS: forward direction
        #$image->setThickness(0.5);
        if ($end>$strart){                                            
                my $index_rgb_color;
                #$index_rgb_color = int(log2($cds[4])) if $cds[4]>0;
                $index_rgb_color = $cds[4];
                #$index_rgb_color = 0 if $cds[4]==0;
                my @final_rgb = split ",", $rgb_color[$index_rgb_color];
                my $color = $image->colorAllocate($final_rgb[0], $final_rgb[1], $final_rgb[2]);
                my $s = $strart;
                my $e = $end;
                my $Y_c = $Y_centre;
                my $w_h = $arrow_high;
                #$image->filledPolygon(rec_arrow_F($strart,$end,$arrow_length,$arrow_high, $Y_centre),$color); # if $cds[4]>0;#$rgb_color[$cds[4]], number = 0 is not shown
                #$image->filledPolygon(rec_arrow_F($strart,$end,$arrow_length,$arrow_high, $Y_centre),$gray) if $cds[4]==0;#$rgb_color[$cds[4]], number = 0 is not shown
                
                #$image->line($e-($e-$s+1)/2, $Y_c-$w_h/4,$e-($e-$s+1)/2, $Y_c+$w_h/4,$color) if $cds[4]>0; #verity
                #$image->line($e-($e-$s+1)/2, $Y_c-$w_h/4,$e-($e-$s+1)/2, $Y_c+$w_h/4,$black) if $cds[4]==0; #verity
                #$image->line($e, $Y_c-$w_h/2,$s, $Y_c-$w_h/2,$color) if $cds[4]>0; #level
                #$image->line($e, $Y_c-$w_h/2,$s, $Y_c-$w_h/2,$white) if $cds[4]==0; #level

                #$rec_arrow->addPt($e-$a_l, $Y_c-$w_h/2);      # piont3
                #$rec_arrow->addPt($s, $Y_c-$w_h/2);           # piont4
                #$rec_arrow->addPt($s, $Y_c+$w_h/2);           # piont5
                #$rec_arrow->addPt($e-$a_l, $Y_c+$w_h/2);      # piont6

#print $filename,":\t$rgb_color[$index_rgb_color]\t",$image->colorAllocate($rgb_color[$index_rgb_color]);
#my $len = $end-$strart;
#print "$cds[3]\t$cds[4]\t$len\t$strart\n";
                $image->polygon(rec_arrow_F($strart, $end, $arrow_length, $arrow_high, $Y_centre), "$white") if $cds[2] eq "CDS"; #set color of border for CDS genes
                $image->polygon(rec_arrow_F($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$pseudo_color_border}) if scalar @cds > 5; #set color of border for Pseudo genes
                $image->polygon(rec_arrow_F($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$RNA_color_border}) if $cds[2] =~ /^rRNA|^tRNA/; #set color of border for RNA genes
                $image->filledPolygon(rec_arrow_F($strart,$end,$arrow_length,$arrow_high, $Y_centre),$color); # if $cds[4]>0;#$rgb_color[$cds[4]], number = 0 is not shown
        }

        # CDS: reverse direction
        elsif ($end<$strart){                           
           
                #$image->filledPolygon(rec_arrow_R($strart,$end,$arrow_length,$arrow_high, $Y_centre), 15);
                #$image->polygon(rec_arrow_R($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$cds_color_border}) if $cds[2] eq "CDS"; #set color of border for CDS genes
                #$image->polygon(rec_arrow_R($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$pseudo_color_border}) if scalar @cds > 5; #set color of border for Pseudo genes
                #$image->polygon(rec_arrow_R($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$RNA_color_border}) if $cds[2] =~ /^rRNA|^tRNA/; #set color of border for RNA genes

        }

    }

 return $image;

} 


##################################################################################################
###### Subrounting--rec_arrow_F
###### Function:
###### generate senven points to draw polygon (CDS: forward direction) 
##################################################################################################
sub rec_arrow_F {
    my $s=shift;
    my $e=shift;
    my $a_l=shift;
    my $w_h=shift;
    $w_h = $w_h*2;
    my $Y_c=shift;
    my $rec_arrow=GD::Polygon-> new();
    #$rec_arrow->addPt($e, $Y_c);                  # piont1 (medium)
    #$rec_arrow->addPt($e-$a_l, $Y_c-$w_h);        # piont2 
    $rec_arrow->addPt($e-$a_l, $Y_c-$w_h/2);      # piont3
    $rec_arrow->addPt($s, $Y_c-$w_h/2);           # piont4
    $rec_arrow->addPt($s, $Y_c+$w_h/2);           # piont5
    $rec_arrow->addPt($e-$a_l, $Y_c+$w_h/2);      # piont6
    #$rec_arrow->addPt($e-$a_l, $Y_c+$w_h);        # piont7
    #$image->filledPolygon(rec_arrow_F($strart,$end,$arrow_length,$arrow_high, $Y_centre),$color)
return $rec_arrow;              

}


##################################################################################################
###### Subrounting--rec_arrow_R
###### Function:
###### generate senven points to draw polygon (CDS: reverse direction) 
##################################################################################################
sub rec_arrow_R {
    my $s=shift;
    my $e=shift;
    my $a_l=shift;
    my $w_h=shift;
    $w_h = $w_h*2;
    my $Y_c=shift;
    my $rec_arrow=GD::Polygon-> new();
    $rec_arrow->addPt($e, $Y_c);                  # piont1
    $rec_arrow->addPt($e+$a_l, $Y_c-$w_h);        # piont2
    $rec_arrow->addPt($e+$a_l, $Y_c-$w_h/2);      # piont3
    $rec_arrow->addPt($s, $Y_c-$w_h/2);           # piont4
    $rec_arrow->addPt($s, $Y_c+$w_h/2);           # piont5
    $rec_arrow->addPt($e+$a_l, $Y_c+$w_h/2);      # piont6
    $rec_arrow->addPt($e+$a_l, $Y_c+$w_h);        # piont7

return $rec_arrow; 

}


sub TreeDraw {

my ($workplace, $phylogenetic_file, $show_tree_topology, $use_branch, $bootstrap, $xstep, $left_shift, $up_shift, $ystep, $font, $font_size) = @_;

#my $temp_strain_reorder_file = "$workplace/temp_strain_reorder_file-svg.txt";

my %xx;            # horizontal coordinate for each node
my %yy;            # vertical coordinate for each node
my $tree_object;   # first Bio::Tree::Tree object
#$xstep = 10;         # branch length in drawing      
my $tip = 5; # extra space between tip and label

##########################################################################################################
    my $treeio = Bio::TreeIO->new(-format => 'newick',
                                  -file   => $phylogenetic_file); 
    $tree_object = $treeio->next_tree;  
    my @taxa1 = $tree_object->get_leaf_nodes;
    my $root1 = $tree_object->get_root_node;

    #open (STRAIN_REORDER_FILE, ">$temp_strain_reorder_file");
    #my $count_number=0;
    my @strain_order = reverse $tree_object->get_nodes(-order => 'depth');
    pop @strain_order; # skip root

    my $y = $up_shift+$ystep;
    for my $node_leaf (@strain_order) {
        #print $node_leaf->id;
        if ($node_leaf->is_Leaf){
            #$count_number++;
            #print STRAIN_REORDER_FILE $node_leaf->id, "\t", "$count_number", "\n";

            #print "test:   ", $node_leaf->id, "\n";
            $xx{$node_leaf} = 0;                 # a temp value
            $yy{$node_leaf} = $y;                # leaf y-value
            $y += $ystep;
        }
    }
    #close STRAIN_REORDER_FILE;

if($show_tree_topology eq "T") {

########################################################################################################## 
#set width of the image###################################################################################

    my @stack;
    my @queue; # postorder traversal
    push @stack, $tree_object->get_root_node;
    while (@stack) {
        my $node = pop @stack;
        push @queue, $node;
        foreach my $child ($node->each_Descendent(-sortby => 'internal_id')) {
            push @stack, $child;
        }
     }

    
########################################################################################################## 

    if ($use_branch eq "F") { # ragged right, ignoring branch lengths

    @queue = reverse @queue;
    my @floor;
    for my $node (@queue) {
        if (!$node->is_Leaf) {
            my @children = $node->each_Descendent;
            my $child = shift @children;
            my $xmin = $xx{$child};    # a temp value
            foreach $child (@children) {
	        $xmin = $xx{$child} if $xx{$child} < $xmin;
            }
            $xx{$node} = $xmin - $xstep;
            push @floor,  $xx{$node};
        }
        
    }

        @floor = sort {$a<=>$b} @floor;
        my $floor_tree = ($floor[$#floor] - $floor[0])/$xstep; 
        my $x = $left_shift + $xstep * ($floor_tree+1) + $tip;

        for my $taxon (reverse @taxa1) { 
            $xx{$taxon} = $x - $tip;
        }        

        for my $node (@queue) {
            if (!$node->is_Leaf) {

                my @children = $node->each_Descendent;
                my $child = shift @children;
                my $xmin = $xx{$child};
                my $ymin = my $ymax = $yy{$child};
                    foreach $child (@children) {
	                $xmin = $xx{$child} if $xx{$child} < $xmin;
	                $ymax = $yy{$child} if $yy{$child} > $ymax;
	                $ymin = $yy{$child} if $yy{$child} < $ymin;
                    }

                $xx{$node} = $xmin - $xstep;

                $yy{$node} = ($ymin + $ymax)/2;    # no-leaf y-value
            }
        
        }

        my @preorder = $tree_object->get_nodes(-order => 'depth');

        my $root_strain = shift @preorder; # skip root

        for my $node (@preorder) {

              $xx{$node} = $xx{$node->ancestor} + $xstep;    # determinate the all nodes values in x-axis

        }


    } else { # set to aspect ratio and use branch lengths if available
        
        @queue = reverse @queue;
        for my $node (@queue) {
            if (!$node->is_Leaf) {
                my @children = $node->each_Descendent;
                my $child = shift @children;
                my $ymin = my $ymax = $yy{$child};

                foreach $child (@children) {
	            $ymax = $yy{$child} if $yy{$child} > $ymax;
	            $ymin = $yy{$child} if $yy{$child} < $ymin;
                }

                $yy{$node} = ($ymin + $ymax)/2;

            }
        
        }

        my @length_x;
        for ( $root1->get_all_Descendents ) {
            $xx{$_} = $left_shift + depth1($_)*2000;   # determinate the all nodes values in x-axis
        }        

    }

    my @strain_tree = reverse $tree_object->get_nodes(-order => 'depth');
    pop @strain_tree; # skip root

    for my $taxon (@strain_tree) {
        #print $node_leaf->id;
        if ($taxon->is_Leaf){
            #$image->stringrotate($xx{$taxon} + $tip, $yy{$taxon}+$font_size / 3, $font, "bold", $font_size/1.5, $taxon->id, $red, 0); # draw strain
            dottedline (5+$tree_width,$yy{$taxon},$xx{$taxon},$yy{$taxon}, $dgray);  #align right strain name using dotted line 
            
        }
    }

##############################################################################################################
##############################################################################################################
    for my $node ($tree_object->get_nodes) {

        if ($node->ancestor) {
        
            if (defined $xx{$node->ancestor}) {

                $image->line($xx{$node},$yy{$node},$xx{$node->ancestor},$yy{$node},$black);


                $image->line($xx{$node->ancestor},$yy{$node},$xx{$node->ancestor},$yy{$node->ancestor},$black);
            } else {

                $image->line($xx{$node},$yy{$node},$left_shift,$yy{$node},$black);
                $image->line($left_shift,$yy{$node},$left_shift,$yy{$node->ancestor},$black);
            }

                if ( $bootstrap eq "T" ) {               
                    if (defined $node->ancestor->id) {
                        my $bootstrap_value =  int (($node->ancestor->id)*100);
 
                        $image->stringrotate($xx{$node->ancestor}+ $font_size/10, $yy{$node->ancestor}+ ($font_size / 3), $font, 'normal', $font_size*0.8, $bootstrap_value, $black, 0); # draw bootstrap value
                    } 

                }

        }

    }
########################################################################################################## 
#set root value of the image##############################################################################
    my $ymin = $yy{$root1};
    my $ymax = $yy{$root1};
    foreach my $child ($root1->each_Descendent) {
        $ymax = $yy{$child} if $yy{$child} > $ymax;
        $ymin = $yy{$child} if $yy{$child} < $ymin;
    }

    my $zz = ($ymin + $ymax)/2;

    if (!defined  $xx{$root1}) {
        $xx{$root1} = $left_shift;
 
    }
 
    $image->line($xx{$root1},$zz,$xx{$root1}- $xstep, $zz,$black);   


    #return $image;
}
    #return $temp_strain_reorder_file;
}

sub dottedline {   # output dotted line

    my ($x1, $y1, $x2, $y2, $color) = @_;

    my $dotted_step = 2;

    for (my $change_x=0; $change_x<=$x1-$x2-$dotted_step; $change_x++) {
        $change_x += $dotted_step;
        $image->line($x2 + $dotted_step+$change_x, $y1, $x2+$change_x, $y2, $color);
        $change_x += $dotted_step;
    }

    return $image;

}


sub depth1 {
   
   my $depth = 0;
   my $node = shift;
   while( defined $node->ancestor ) { 

       my $branch_length_node;
       if (defined $node->branch_length) {
           $branch_length_node = $node->branch_length;
       }else{  
           $branch_length_node = 0;  #the length of some branches length may lose
       } 
      
       $depth += $branch_length_node;
       $node = $node->ancestor;

   }
   return $depth;
}


sub tree_width {

    my ($phylogenetic_file, $use_branch, $left_shift, $xstep, $workplace) = @_;
    my %xx;        # horizontal coordinate for each node
    my $width;     # total drawing width
    my $temp_strain_reorder_file = "$workplace/temp_strain_reorder_file-svg.txt";

    my $treeio = Bio::TreeIO->new(-format => 'newick',
                                  -file   => $phylogenetic_file); 
    my $tree_obj = $treeio->next_tree;  # first Bio::Tree::Tree object
    my @taxa = $tree_obj->get_leaf_nodes;
    my $root = $tree_obj->get_root_node;
    

    open (STRAIN_REORDER_FILE, ">$temp_strain_reorder_file");
    my $count_number=0;
    my @strain_order1 = reverse $tree_obj->get_nodes(-order => 'depth');
    pop @strain_order1; # skip root
    for my $node_leaf1 (@strain_order1) {
        #print $node_leaf->id;
        if ($node_leaf1->is_Leaf){
            $count_number++;
            print STRAIN_REORDER_FILE $node_leaf1->id, "\t", "$count_number", "\n";
        }
    }
    close STRAIN_REORDER_FILE;


    for my $taxon (reverse @taxa) {     #strain_name in Y-aixs    xiangyang Li 2019-07-17
        $xx{$taxon} = 0;  # a temp value

    }

    #$xstep = 10;  # branch length in drawing


    my @stack;
    my @queue; # postorder traversal
    push @stack, $root;
    while (@stack) {
        my $node = pop @stack;
        push @queue, $node;
        foreach my $child ($node->each_Descendent(-sortby => 'internal_id')) {
            push @stack, $child;
        }
     }

    if ($use_branch eq "F") { # ignoring branch lengths

    @queue = reverse @queue;
    my @floor;
    for my $node (@queue) {
        if (!$node->is_Leaf) {
            my @children = $node->each_Descendent;
            my $child = shift @children;
            my $xmin = $xx{$child};    # a temp value
            foreach $child (@children) {
	        $xmin = $xx{$child} if $xx{$child} < $xmin;
            }
            $xx{$node} = $xmin - $xstep;
            push @floor,  $xx{$node};
        }
        
    }

        @floor = sort {$a<=>$b} @floor;
        my $floor_tree = ($floor[$#floor] - $floor[0])/$xstep;
 
        $width = $left_shift + $xstep * ($floor_tree+1) ; # set the width of the treemap 

    } else { # set to aspect ratio and use branch lengths if available

        my @length_x;
        for ( $root->get_all_Descendents ) {
            push @length_x, depth1($_);           
            $xx{$_} = $left_shift + depth1($_)*2000;  
 
        }        
        @length_x = sort {$a<=>$b} @length_x;
        $width = $left_shift + $length_x[$#length_x]*2000; # set the width of the treemap 

    }

    return ($width, $count_number, $temp_strain_reorder_file);

}



sub color_spectrum {

    my ($color_num) = @_;
    @rgb_color = ();
    my @rgb_color_source = ("212,212,212", "255,0,0","255,170,51","255,255,0","141,204,1","0,128,0","8,143,143","0,0,255","59,67,192","128,0,128","149,53,83");
    foreach my $i_num (0..$color_num){
        push @rgb_color, $rgb_color_source[$i_num];
        #push @rgb_color, $rgb_color_source[-1] if $i_num > scalar @rgb_color_source; #avoid color_num > the number of rgb_color_source
    }

}


sub log2 { 
    my $n = shift; 
      
    # using pre-defined log function 
    return log($n) / log(2); 
}


sub string_svg_width{

    my ($ref_hash, $sfont_size, $tfont_family, $tfont_style) = @_;
    my @tname_width_svg;
    foreach (keys %$ref_hash){
        my $filename = $$ref_hash{$_};
     
        #svg string_width
        my $tname_width = metrics::string_width_svg($filename, "$tfont_family:$tfont_style", $sfont_size); 
        push (@tname_width_svg, $tname_width);  
  
    }
    @tname_width_svg = sort{$b<=>$a} @tname_width_svg;

    return $tname_width_svg[0];

}



1;

__END__  
