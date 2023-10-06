Pneumo-Typer
===

Pneumo-Typer is a comprehensive prediction and visualization of serotype and sequence type for streptococcus pneumoniae using assembled genomes. It is freely available at https://www.microbialgenomic.cn/Pneumo-Typer.html and https://github.com/Xiangyang1984/Pneumo-Typer under an open-source GPLv3 license.


* [Installation](#installation)

 	* [Installing the Pneumo-Typer via Conda](#installing-the-Pneumo-Typer-via-conda)
	
 	* [Installing the Pneumo-Typer from Source Code](#installing-the-Pneumo-Typer-from-source-code)
	
 	* [Test the Pneumo-Typer with Example Data](#test-the-Pneumo-Typer-with-example-data)
	
 * [Usage](#usage)
	
 	* [Running Pneumo-Typer](#running-Pneumo-Typer)
	
	* [Detailed Explanations for Arguments in Pneumo-Typer](#detailed-explanations-for-arguments-in-Pneumo-Typer)


# Installation

Pneumo-Typer is a Perl script that doesn't need compilation. But before running, Pneumo-Typer needs to pre-install prodigal, blat, blast, and several Perl modules. There are two ways to install the Pneumo-Typer.

## Installing the Pneumo-Typer via Conda
We have built a bioconda package for Pneumo-Typer. Users are recommended to install the [conda](https://www.anaconda.com), and then install this package simply with the following command:

	$ conda install -c bioconda pneumo-typer

Once the installation is finished, the absolute paths for blat, prodigal, blastn and makeblastdb have been auto-configured well for pneumo-typer.pl, so users should be able to run Pneumo-Typer.

## Installing the Pneumo-Typer from Source Code
Installation of Pneumo-Typer can be accomplished by downloading the code (at https://www.microbialgenomic.cn/Pneumo-Typer.html and https://github.com/xiangyang1984/Pneumo-Typer.git) and then following the steps below.
#### Step 1: Download source code
Download Pneumo-Typer, and put the Pneumo-Typer directory into your PATH with the following command：

	```	
	$ wget -c https://www.microbialgenomic.cn/gz/pneumo-typer-v1.0.1.tar.gz (Recommended to use)

	or

	$ git clone https://github.com/xiangyang1984/Pneumo-Typer.git

	***Pneumo-Typer contains two large-size files which need to be downloaded manually.
	* cgMLST_profiles
	$ wget -c -O path_to_Pneumo-Typer/ST_tool/database/cgmlst/cgMLST_profiles https://media.githubusercontent.com/media/Xiangyang1984/Pneumo-Typer/main/ST_tool/database/cgmlst/cgMLST_profiles?download=true
	* cgMLSA_loci.fas
	$ wget -c -O path_to_Pneumo-Typer/ST_tool/database/cgmlst/cgMLSA_loci/cgMLSA_loci.fas https://media.githubusercontent.com/media/Xiangyang1984/Pneumo-Typer/main/ST_tool/database/cgmlst/cgMLSA_loci/cgMLSA_loci.fas?download=true

	$ export PATH=/path/to/Pneumo-Typer/:$PATH
	```

#### Step 2: Perl modules installation
Pneumo-Typer requires Perl as well as Perl modules including GD; GD::SVG, SVG; threads, File::Basename, FindBin, File::Spec, lib, Getopt::Long, Math::BigFloat, Storable, vars, Bio::SeqIO, Bio::Tree::NodeI, Bio::TreeIO. 

These modules can be installed with cpan using:

	$ sudo cpan install GD GD::SVG SVG threads File::Basenamey FindBin lib Getopt::Long Math::BigFloat Storable vars BioPerl


#### Step 3: Programs installation
Additional software dependencies for the Pneumo-Typer are as follows:

* makeblastdb and blastp
  
Both of them come from NCBI BLAST+, available at https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

* prodigal

This software is available at http://micans.org/mcl/

* blat
  
This software is available at http://micans.org/mcl/


## Test the Pneumo-Typer with Example Data
Once Pneumo-Typer installation is finished, a small dataset in the **./test_data** directory can be used to test whether Pneumo-Typer (**pneumo-typer.pl**) can run on your system (**Linux/MacOS**) successfully or not using the command as below:

	$ perl pneumo-typer -Ts T
	
	Test-step1: Checks for pneumo-typer dependencies...
	################################################################
	***GD Version  2.71 ok.
	***GD::SVG Version   0.33 ok.
	***SVG Version   2.84 ok.
	***threads Version   2.15 ok.
	***File::Basename Version  2.85 ok.
	***FindBin Version   1.51 ok.
	***lib Version   0.63 ok.
	***Getopt::Long Version  2.49 ok.
	***Math::BigFloat Version  1.999806 ok.
	***Storable Version  3.15 ok.
	***vars Version  1.03 ok.
	***File::Spec Version  3.75 ok.
	***Bio::SeqIO Version  1.007002 ok.
	***Bio::Tree::NodeI Version   ok.
	***Bio::TreeIO Version   1.007002 ok.
	!!!Ok, all dependencies Perl modulers are installed*
	----------------------------------------------------------------
	Checking for makeblastdb ... OK, makeblastdb is installed at: /miniconda3/bin//makeblastdb
	Checking for blastn ... OK, blastn is installed at: /miniconda3/bin//blastn
	Checking for prodigal ... OK, prodigal is installed at: /miniconda3/bin//prodigal
	Checking for blat ... OK, blat is installed at: /miniconda3/bin//blat
	################################################################


	Test-step2: Begin test pneumo-typer.pl...
	################################################################

	Fri Jul 21 19:56:13 2023: pneumo-typer.pl start...

	STEP-1: Dealing with genomes extract genome sequence, gene sequences and gene feature table (TFT);
	annotate genome which has no annotation information using prodigal>
	Annotating genome using prodigal: 20%...40%...60%...80%...100%...done

	STEP-2: Determining the Sequence Type (ST)
  
	STEP-2.1: MLST analysis
	Blastn_percent: 20%...40%...60%...80%...100%...done
 
	STEP-2.2: cgMLST analysis
	Blat_percent: 20%...40%...60%...80%...100%...done

	STEP-3: Predicting serotype
	Blastn_percent: 10%...20%...30%...40%...50%...60%...70%...80%...90%...100%...done

	STEP-4: Output sequence type and serotype results

	STEP-5: Heatmaping the cps gene distribution in genomes

	STEP-6: Visualizing the cps gene cluster in each genome

	Fri Jul 21 20:07:30 2023: done!

	################################################################
	Ok, pneumo-typer works successfully!

 

# Usage

It is very simple to use Pneumo-Typer. First, prepare input data, at least containing a Genbank_file_directory (containing files in GenBank format, FASTA format, or a combination of both); then, run Pneumo-Typer like this "perl pneumo-typer.pl -d Genbank_file_directory". 


## Running Pneumo-Typer

Here, we used 18 genomes as an example to show how to use Pneumo-Typer. 18 genomes are under a directory named "18_genomes_dir"([downoload 18_genomes.tar.gz](https://www.microbialgenomic.cn/temp_dir/18_genomes.tar.gz)).

#### Example 1: Run Pneumo-Typer is an easy task by using the following command


*Perform serotype prediction, heatmap and figure creation, and ST analysis with 10 threads:
	
	$ perl path_to_pneumo-typer.pl -d path_to_18_genomes_dir -t 10 -m T
 
Pneumo-Typer will output results as follows:
* a.ST results
* b.predicted serotype results
* c.create three maps with the ST and predicted serotype and information showed 
	*heatmap_gene.svg: a heatmap of the distribution of cps gene at gene level;
  
	*heatmap_class.svg: a heatmap of the distribution of cps gene at class level;
  
	*cps_cluster.svg: a figure showing the genetic organization of cps gene cluster
		
Setting "-c" to "T" will perform cgST analysis which takes quite a long time, and the cgST information will also be shown on maps.

#### Example 2: A Newick format tree file is used by Pneumo-Typer to automatically associate the distribution of cps gene and genetic organization of cps cluster with their phylogeny.

[18_genome_tree.nwk](https://www.microbialgenomic.cn/temp_dir/18_genome_tree.nwk): a tree file, it was a nwk format phylogenetic tree of 18 genomes using RaxML. 
 

*Create two heatmaps according to the distribution of the cps gene at the gene and class level using the following command:

	$ perl path_to_pneumo-typer.pl -d path_to_18_genomes_dir -Rh T -tree path_to_18_genome_tree.nwk		

*Create a map of the genetic organization of cps gene using the following command:
	
	$ perl path_to_pneumo-typer.pl -d path_to_18_genomes_dir -Rf T -tree path_to_18_genome_tree.nwk

	
#### Example 3: A two-column tab-delimited text file is used to sort genomes from up to down according to users' requirement

Here, we provided a srf file named ["18_genome_order.txt"](https://www.microbialgenomic.cn/temp_dir/18_genome_order.txt) that orders the maps by serotypes. For example, users can reorder genomes in the genetic organization of the cps cluster.

	$ perl path_to_pneumo-typer.pl -d path_to_18_genomes_dir -Rf T -srf path_to_18_genome_order.txt
		  

 
## Detailed Explanations for Arguments in Pneumo-Typer

```
--------------
    REQUIRED ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -d, --genbank_file_directory
           A directory containing files in GenBank format, FASTA format, or a combination of both.                           
    OPTIONAL ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
	-o, --output_directory
           An output directory holding all the generated files by pneumo-typer.pl. if this option is not set,  a directory named "ipneumo-pyper_workplace" will be created in the bin directory from where pneumo-typer.pl was invoked.
       -t, --multiple_threads
           Set thread number (Default: 1)
       -Ss, --skip_sequence_processing 
           Skip the process of sequence processing (STEP-1) (Default: F).  
       -hgc, --homologous_gene_cutoff
           Set E-value, Identify, Coverage (Query and Subject), Match_length (alignment length) cutoff in Blastn analysis (default: E-value=1e-5, Identify=70, Coverage=95, Match_length=100).
       -Sb, --skip_blastn
           Skip the process of doing blastn during serotype analysis.
       -p, --prodigal_annotation
           Annotate all genomes using prodigal. 
       -m, --mlst
           Perform mlst analysis (Default: T). 
       -c, --cgmlst
           Perform cgmlst analysis. It need >10 mins for one genome (Default: F).
       -Rh, --recreate_heatmap                             
           Re-create the heatmap of cps gene distribution in genomes (Default: F). At this step, users can add a parameter "phylogenetic_tree" or "strain_reorder_file". 
       -Rf, --recreate_figure
           Re-create the figure of the genetic organlization of cps gene cluster for genomes (Default: F). At this step, users can add a parameter "phylogenetic_tree" or "strain_reorder_file".
       -tree, --phylogenetic_tree
           A Newick format tree file is used by Pneumo-Typer to automatically accociate the genomes with their phylogeny. Meanwhile, Pneumo-Typer will output a file named "temp_strain_reorder_file", which contains the order information of genomes in tree from up to down. It should be noted that all nodes name in provided tree must completely match with the genbank files name of all genomes.
       -srf, --strain_reorder_file
           A two-column tab-delimited text file is used to sort genomes from up to down accoding to users requirement. Each row must consist of a strain name followed by the numerical order that is used for sorting genomes. It should be noted that all strains name must completely match with the genbank files name of all genomes.
       -Ts, --test
           Run pneumo-typer using Test_data as input to check whether Pneumo-Typer is installed successfully (Default: F).
       -V, --version
           The version of Pneumo-Typer.
       -h, --help
           Show this message.

```

## COPYRIGHT

Dr. Xiangyang Li (E-mail: lixiangyang\@fudan.edu.cn, lixiangyang1984\@gmail.com), Kaili University; Bacterial Genome Data mining & Bioinformatic Analysis (http://www.microbialgenomic.cn/).

Copyright 2023, Xiangyang Li. All Rights Reserved.
