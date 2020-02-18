# Gcluster_v1.01
Gcluster is a simple tool for visualizing and comparing genome contexts for massive genomes. It is freely available at
http://www.microbialgenomic.com/Gcluster_tool.html and https://github.com/Xiangyang1984/Gcluster_v1.01 under an open source 
GPLv3 license. It is a stand-alone Perl application, which requires MCL, NCBI BLAST+ and several Perl Modules (e.g. GD, GD::SVG) to
be installed before use.
## Installing the OrthoMCL Pipeline
1. Installation of required Perl modules and programs
======================================================
Gcluster is a perl script which doesn't need compilation. But before running, Gcluster needs to pre-install several Perl Modules
and three extra programs. In addition, the paths of those three programs in Gcluster.pl need to be set.

Step 1: Perl Modules Dependencies:
----------------------------------
* GD;
* GD::SVG;
* SVG;
* threads;
* File::Basename;
* FindBin;
* lib;
* Getopt::Long;
* Math::BigFloat;
* Storable;
* vars;
* File::Spec;
* Bio::SeqIO; (part of BioPerl, http://bioperl.org)
* Bio::Tree::NodeI; (part of BioPerl, http://bioperl.org)
* Bio::TreeIO; (part of BioPerl, http://bioperl.org)

These can be installed with cpan using:
$ sudo cpan install GD GD::SVG SVG threads File::Basenamey FindBin lib Getopt::Long Math::BigFloat Storable vars Bio::SeqIO Bio::Tree::NodeI Bio::TreeIO


Step 2: Programs Dependencies:
------------------------------
* makeblastdb and blastp in NCBI BLAST+, which is available from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/.
* mcl (Markov Clustering algorithm), is available at http://micans.org/mcl/
              
Please set the absolute path for three programs in "Gcluster/Gcluster.pl", as in the following example:
my $blastp        = "/usr/bin/blastp";
my $makeblastdb   = "/usr/bin/makeblastdb";
my $MCL           = "/usr/bin/mcl";

 These software dependencies can be checked and the configuration file created using the **scripts/setup.pl** script as below:





2. Running Gcluster
====================
 
perl Gcluster.pl -dir genbank_file_directory -gene interested_gene_file [options]

FOR EXAMPLE: 
perl /home/xiangyang/Gcluster_v1.01-master/Gcluster.pl -dir /home/xiangyang/Gcluster_v1.01-master/test_data/gbk -gene
/home/xiangyang/Gcluster_v1.01-master/test_data/interested_gene_name.txt -tree /home/xiangyang/Gcluster_v1.01-
master/test_data/full_NJ_rRNA_tree.nwk -m 3

Large test data is available at website (http://www.microbialgenomic.com/160_genomes_testdata.tar.gz). It contains 160 annotated 
genomes used as input files to create Fig. 1 in manuscript.


3. Gcluster parameters
=======================

    REQUIRED ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -dir, --genbank_file_directory
             A Directory containing annotated genomes as Genbank format file. To avoid a mistake, genome names cannot use
             special character, such as space, equal. For large number of genomes, users are recommended to download using
             Aspera, a high-speed file transfer tool (https://downloads.asperasoft.com/). If enough rgb colors provided or not
             to color homologous genes, there are no limitation to the number of genomes and the number of genes flanking gene
             of interest for Gcluster.
       -gene, --interested_gene_file
             A list of the interested gene, in which each line contains a locus tag of the interested gene for individual
             genome. Users are recommended to use "interested_gene_generation.pl" in Gcluster package for generation this file.
             In this situation, user needs to provide a blast database file in FASTA format, which contains at least one protein
             sequence homologous to the gene of interest. To map genome contexts to a given phylogeny or to order the genome
             contexts using a "strain_reorder_file", only one gene of interest is allowed to provide for each genome.
             For example:
                 AX2_RS10405	#arsenite_oxidase_large_subunit;Achromobacter_xylosoxidans_NBRC_15126_ATCC_27061
                 KUC_RS10495	#arsenite_oxidase_large_subunit;Halomonas_boliviensis_LC1
                 KYC_RS14580	#arsenite_oxidase_large_subunit;Achromobacter_arsenitoxydans_SY8
                 ...

    OPTIONAL ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -o, --Gcluster_output_directory
             An output directory holding all the generated files by Gcluster.pl. if this option is not set, Gcluster will create a 
             directory named "Gcluster_workplace" in the bin directory from where Gcluster.pl was invoked. 
       -tree, --phylogenetic_file
             A Newick format tree file is used by Gcluster to automatically accociate the genomes with their phylogeny.
             Meanwhile, Gcluster will output a file named "temp_strain_reorder_file", which contains the order information of
             genomes in tree from up to down. It should be noted that all nodes name in provided tree must completely match with
             the genbank files name of all genomes. 
             Gcluster provides a perlscript in "Gcluster/script" directory for batch extraction of 16S rRNA gene sequences, 
             which can be used to build a 16S rRNA tree using software like MEGA (https://www.megasoftware.net/). 
       -topology, --show_tree_topology
             Display the tree topology, which is obtained from the tree file (Default: T).
       -branch, --show_tree_branch
             Draw tree using tree branch length, which is obtained from the tree file (Default: F).
       --x_step
             Draw tree using xstep instead of tree branch length (Default: 10).
       -bootstrap, --show_tree_bootstrap
             Display the tree bootstrap value, which is obtained from the tree file (Default: F).
       -srf, --strain_reorder_file
             A two-column tab-delimited text file is used to sort genomes from up to down accoding to users requirement. Each
             row must consist of a strain name followed by the numerical order that is used for sorting genomes. It should be
             noted that all strains name must completely match with the genbank files name of all genomes. Gcluster needs a
             "strain_reorder_file" or a "phylogenetic_file", but not both at the same time. 
             For example:
                 Achromobacter_xylosoxidans_NCTC10807	1
                 Achromobacter_xylosoxidans_NBRC_15126_ATCC_27061	2
                 Achromobacter_xylosoxidans_NCTC10808	9
                 Achromobacter_sp._2789STDY5608623	5
                 Achromobacter_sp._2789STDY5608621	4
                 Alcaligenes_faecalis_subsp._faecalis_NCIB_8687	10
                 Achromobacter_xylosoxidans_DPB_1	7
                 Achromobacter_xylosoxidans_B_1	8
                 Achromobacter_piechaudii_HLE	3
                 Achromobacter_marplatensis_B2	6
                 ...	...

       -n, --flanking_gene_number
             Number of genes flanking gene of interest are set to show. If enough rgb colors provided or not to color homologous
             genes, there are no limitation to the number of genomes and the number of genes flanking gene of interest (Default:
             10).
       -color_f, --gene_color_filled
             Color was used to fill homologous gene clusters or gene families of interest (Default: T), if choose F, all of the
             genes were filled with the color customized by "gene_no_color_filled" parameter).
       -pso, --percent_strain_homologouscluster_color
             Only color certain homologous gene clusters, in which the holding number of different genomes exceeds the threshold
             number (Default: 0). This is measured using (X/Y)*100. In this formula, X denotes the number of different genomes
             in a set of homologous gene cluster, and Y denotes the total number of genomes. This parameter is useful when no
             enough rgb colors are provided in colors_configure_file under color_configure direcoty. Users could to reduce the 
             number of colors used by setting a high value.  
       -c_color_b, --cds_color_border
             To color the border of the CDS genes (Default: black), users can choose from blue, black, red, white, gray, dgray.
       -p_color_b, --pseudo_gene_color_border
             To color the border of the Pseudo genes (Default: dgray), users can choose from blue, black, red, white, gray, 
             dgray.
       -r_color_b, --RNA_gene_color_border
             To color the border of the RNA (tRNA, rRNA) genes (Default: red), users can choose from blue, black, red, white, 
             gray, dgray.
       -no_color_f, --gene_no_color_filled
             To fill uniqe genes (including RNA genes), pseudo genes, and homologous gene clusters not meeting the criteria set
             by "percent_strain_homologouscluster_color" parameter with a single color (Default: white), users can choose from 
             blue, black, red, white, gray, dgray.
       -l, --arrow_relative_Length    
             Set the relative length of the gene arrow (Default: 4).
       -w, --arrow_relative_Height
             Set the relative Height of the gene arrow (Default: 6).
       -scale, --figure_Scale_up_multiple
             Adjust gene length through zooming (Default: 0.5).
       -s_Y, --strain_name_shift_Y
             Set the offset along Y-axis for strain names (Default: 0).
       -g_Y, --gene_label_shift_Y
             Set the offset along Y-axis for gene labels (Default: 2).
       -dis, --distance_between_two_genomes
             Set the distance between two genome contexts in Y-axis (Default: 70).
       -up, --up_shift
             Set the top margin of image in pixels (Default: 10).
       -down, --down_shift
             Set the bottom margin of image in pixels (Default: 20).
       -left, --left_shift
             Set the left margin of image in pixels (Default: 10).
       -right, --right_shift
             Set the right margin of image in pixels (Default: 10).
       -label, --show_label
             Display the gene label (gene Locus Tag or genename) (Default: T).
       -family, --font_family
             Set font family for the genome name and the gene label, e.g. Times New Roman, Arial, Verdana and so on (Default:
             Times New Roman). Users are suggested to choose font family listed in metrcis module, or causing a miscalculation
             of string width for genome name in SVG-format map.
       -style, --font_style
             Set font style for the genome name and the gene label, e.g. Normal, Bold, Italic (Default: Normal).
       -size, --label_font_size
             Set font size for gene label (Default: 8).
       -color, --label_font_color
             Set font color for gene label (Default: dgray).
       -i_color, --interested_gene_label_font_color
             Customize gene label color for gene of interest (Default: red), users can choose from blue, black, red, white, 
             gray, dgray.
       -r, --rotate_gene_label (Default: 30)
             Rotate the angle of the gene label, e.g. 30, 45, 135 and so on.
       --strain_name_font_size
             Set font size for genome name (Default: 12).   
       --Strain_name_font_color
             set font color for genome name (Default: black). Users can choose from blue, black, red, white, gray, dgray.
       -Bst, --homologous_gene_cutoff
             Array to set blast parse cutoff: E-value, Identify, Coverage, Match_length (Default: E-value=1-e5, Identify=0, 
             Coverage=50%, Match_length=0).
       -m, --multiple_threads
             Numbers of thread to use (Default: 1).
       -SVG, --SVG_image
             Create SVG format figure (Default: T).
       -PNG, --PNG_image
             Create PNG format figure (Default: T).
       -sub_TFT, --start_at_sub_TFT (Default: F)
             Jump to generate a collection of sub-TFT tables and perform homologous gene analysis (Default: F). Skips sequences
             extraction and TFT file generation.  
       -map, --start_at_map
             Jump to map generation (Default: F). Generation of a collection of sub-TFT tables and homologous gene clusters has
             already been done. This parameter is very useful to customize the map quickly. It should be noted that there's no 
             sense to reset "flanking_gene_number" parameter if this parameter set to "T".
             Importantly, at this step, users can revise the gene label by directly edition of the locus_tag in sub_TFT file or
             all_orthomcl.out. In sub_TFT files and all_orthomcl.out file, there are two forms of gene locus tag, (1) 
             "Locus_Tag", in this case, no genename is defined for a gene; (2) "GeneName;Locus_Tag", in this case, genename is
             given for a gene. For the first form, user can revise gene label by addition of a genename followed by a semicolon
             in the front of the Locus_Tag. For the second first form, user can revise gene label by modification of the 
             genename.
       -h, --help
             Show this message.


Dr. Xiangyang Li (E-mail: lixiangyang@fudan.edu.cn, lixiangyang1984@gmail.com), Fudan university; Kaili University; Bacterial
Genome Data mining & Bioinformatic Analysis (http://www.microbialgenomic.com/).

Copyright 2019, Xiangyang Li. All Rights Reserved.
