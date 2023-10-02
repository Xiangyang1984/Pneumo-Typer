# Pneumo-Typer
pneumo-typer: A fast and simple-to-use tool for predicting serotype and determining sequence type (ST/cgST) for Streptococcus pneumoniae

pneumo-typer is a fast and simple-to-use tool for predicting serotype and determining sequence type (ST/cgST) for Streptococcus pneumoniae. It is freely available at http://www.microbialgenomic.cn/pneumo-typer.html and https://github.com/Xiangyang1984/pneumo-typer under an open source GPLv3 license. It is a stand-alone Perl application, which requires blat, NCBI BLAST+ and several Perl Modules (GD, GD::SVG) to be installed before use.

When you use pneumo-typer, please cited:
Xiangyang Li, Zilin Yang, Zhongrui Ma. pneumo-typer: A fast and simple-to-use tool for predicting serotype and determining sequence type (ST/cgST) for Streptococcus pneumoniae, Submitted to Bioinformatics, 2023XXXX.

Software download: pneumo-typer_v1.03.tar.gz
pneumo-typer User guide

1. Installation
=================


pneumo-typer is a Perl script which doesn't need compilation. But before running, pneumo-typer needs to pre-install prodigal, blat, blast, and several Perl modules. There are two ways to install the pneumo-typer.

1.1-1 Installing the pneumo-typer via Conda

We have build a bioconda package for pneumo-typer. Users are recommended to install the [conda](https://www.anaconda.com), then to install this package simply with the following command:
$ conda install -c xiangyang1984 pneumo-typer
**if occured "The environment is inconsistent, please check the package plan carefully", try to using "conda update --all" before install pneumo-typer.

1.1-2 Installing the pneumo-typer from Source Code

pneumo-typer is available at https://github.com/xiangyang1984/pneumo-typer.git. Installation pneumo-typer can be accomplished by downloading the code and then following the steps below.
#### Step 1: Download source code Download pneumo-typer，and put the pneumo-typer directory into your PATH with the following command： ```
$ git clone https://github.com/xiangyang1984/pneumo-typer.git
$ export PATH=/path/to/pneumo-typer/:$PATH ```
#### Step 2: Perl modules installation
The pneumo-typer requires Perl as well as Perl modules including GD; GD::SVG, SVG; threads, File::Basename, FindBin, File::Spec, lib, Getopt::Long, Math::BigFloat, Storable, vars, Bio::SeqIO, Bio::Tree::NodeI, Bio::TreeIO.
These modules can be installed with cpan using:
$ sudo cpan install GD GD::SVG SVG threads File::Basenamey FindBin lib Getopt::Long Math::BigFloat Storable vars Bio::SeqIO Bio::Tree::NodeI Bio::TreeIO
#### Step 3: Programs installation
Additional software dependencies for the pneumo-typer are as follows:
* makeblastdb and blastn
Both of them come from NCBI BLAST+, available at https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
* prodigal
This software is available at http://micans.org/mcl/
* blat
This software is available at http://micans.org/mcl/

1.2 Test the pneumo-typer with Example Data

Once pneumo-typer installation finished, a small dataset in the **./test_data** directory can be used to test whether pneumo-typer (for **pneumo-typer.pl**) can run on your system (**Linux/MacOS**) successfully or not using the command as below:
$ pneumo-typer -Ts T


test as the follow:
Test-step1: Checks for pneumo-typer dependencies...
################################################################
***GD Version 2.71 ok.
***GD::SVG Version 0.33 ok.
***SVG Version 2.84 ok.
***threads Version 2.15 ok.
***File::Basename Version 2.85 ok.
***FindBin Version 1.51 ok.
***lib Version 0.63 ok.
***Getopt::Long Version 2.49 ok.
***Math::BigFloat Version 1.999806 ok.
***Storable Version 3.15 ok.
***vars Version 1.03 ok.
***File::Spec Version 3.75 ok.
***Bio::SeqIO Version 1.007002 ok.
***Bio::Tree::NodeI Version ok.
***Bio::TreeIO Version 1.007002 ok.
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

STEP-1: Dealing with genomes extract genome sequence, gene seqeunces and gene feature table (TFT);
annnotate genome which has no annotation infotmation using prodigal>
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
Ok, pneumo-typer works success!
# Pneumo-Typer
