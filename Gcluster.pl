#!usr/bin/perl -w 


use strict;
use warnings;
use threads;
use FindBin;
use lib "$FindBin::Bin/lib";
use Getopt::Long;
use File::Basename qw<basename dirname>;
use File::Spec;

######################################################################################################
# Please set the absolute path for three programs before using Gcluster 
my $blastp        = "/usr/bin/blastp";
my $makeblastdb   = "/usr/bin/makeblastdb";
my $MCL           = "/usr/bin/mcl";

# set the inflation for MCL
my $inflation     = 1.5;
######################################################################################################


my $usage = <<USAGE;

=NAME

Gcluster.pl

=DESCRIPTION

Gcluster is a simple tool for visualizing and comparing genome contexts for massive genomes. It is freely available at http://www.microbialgenomic.com/Gcluster_tool.html and https://github.com/Xiangyang1984/Gcluster_v1.01 under an open source GPLv3 license. It is a stand-alone Perl application, which requires MCL, NCBI BLAST+ and several Perl Modules (e.g. GD, GD::SVG) to be installed before use.

=USAGE

Gcluster.pl -dir genbank_file_directory -gene interested_gene_file [options]

FOR EXAMPLE: 
perl /home/xiangyang/Gcluster_v1.01-master/Gcluster.pl -dir /home/xiangyang/Gcluster_v1.01-master/test_data/gbk -gene /home/xiangyang/Gcluster_v1.01-master/test_data/interested_gene_name.txt -tree /home/xiangyang/Gcluster_v1.01-master/test_data/full_NJ_rRNA_tree.nwk -m 3

Large test data is available at website (http://www.microbialgenomic.com/160_genomes_testdata.tar.gz). It contains 160 annotated genomes used as input files to create Fig. 1 in manuscript.

=ARGUMENTS

    REQUIRED ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -dir, --genbank_file_directory
             A Directory containing annotated genomes as Genbank format file. To avoid a mistake, genome names cannot use special character, such as space, equal. For large number of genomes, users are recommended to download using Aspera, a high-speed file transfer tool (https://downloads.asperasoft.com/). If enough rgb colors provided or not to color homologous genes, there are no limitation to the number of genomes and the number of genes flanking gene of interest for Gcluster.
       -gene, --interested_gene_file
             A list of the interested gene, in which each line contains a locus tag of the interested gene for individual genome. Users are recommended to use "interested_gene_generation.pl" in Gcluster package for generation this file. In this situation, user needs to provide a blast database file in FASTA format, which contains at least one protein sequence homologous to the gene of interest. To map genome contexts to a given phylogeny or to order the genome contexts using a "strain_reorder_file", only one gene of interest is allowed to provide for each genome.
             For example:
                 AX2_RS10405	#arsenite_oxidase_large_subunit;Achromobacter_xylosoxidans_NBRC_15126_ATCC_27061
                 KUC_RS10495	#arsenite_oxidase_large_subunit;Halomonas_boliviensis_LC1
                 KYC_RS14580	#arsenite_oxidase_large_subunit;Achromobacter_arsenitoxydans_SY8
                 ...

    OPTIONAL ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -o, --Gcluster_output_directory
             An output directory holding all the generated files by Gcluster.pl. if this option is not set, Gcluster will create a directory named "Gcluster_workplace" in the bin directory from where Gcluster.pl was invoked.  
       -tree, --phylogenetic_file
             A Newick format tree file is used by Gcluster to automatically accociate the genomes with their phylogeny. Meanwhile,  Gcluster will output a file named "temp_strain_reorder_file", which contains the order information of genomes in tree from up to down. It should be noted that all nodes name in provided tree must completely match with the genbank files name of all genomes. 
             Gcluster provides a perlscript in "Gcluster/script" directory for batch extraction of 16S rRNA gene sequences, which can be used to build a 16S rRNA tree using software like MEGA (https://www.megasoftware.net/). 
       -topology, --show_tree_topology
             Display the tree topology, which is obtained from the tree file (Default: T).
       -branch, --show_tree_branch
             Draw tree using tree branch length, which is obtained from the tree file (Default: F).
       --x_step
             Draw tree using xstep instead of tree branch length (Default: 10).
       -bootstrap, --show_tree_bootstrap
             Display the tree bootstrap value, which is obtained from the tree file (Default: F).
       -srf, --strain_reorder_file
             A two-column tab-delimited text file is used to sort genomes from up to down accoding to users requirement. Each row must consist of a strain name followed by the numerical order that is used for sorting genomes. It should be noted that all strains name must completely match with the genbank files name of all genomes. Gcluster needs a "strain_reorder_file" or a "phylogenetic_file", but not both at the same time. 
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
             Number of genes flanking gene of interest are set to show. If enough rgb colors provided or not to color homologous genes, there are no limitation to the number of genomes and the number of genes flanking gene of interest (Default: 10).
       -color_f, --gene_color_filled
             Color was used to fill homologous gene clusters or gene families of interest (Default: T), if choose F, all of the genes were filled with the color customized by "gene_no_color_filled" parameter).
       -pso, --percent_strain_homologouscluster_color
             Only color certain homologous gene clusters, in which the holding number of different genomes exceeds the threshold number (Default: 0). This is measured using (X/Y)*100. In this formula, X denotes the number of different genomes in a set of homologous gene cluster, and Y denotes the total number of genomes. This parameter is useful when no enough rgb colors are provided in colors_configure_file under color_configure direcoty. Users could to reduce the number of colors used by setting a high value.  
       -c_color_b, --cds_color_border
             To color the border of the CDS genes (Default: black), users can choose from blue, black, red, white, gray, dgray.
       -p_color_b, --pseudo_gene_color_border
             To color the border of the Pseudo genes (Default: dgray), users can choose from blue, black, red, white, gray, dgray.
       -r_color_b, --RNA_gene_color_border
             To color the border of the RNA (tRNA, rRNA) genes (Default: red), users can choose from blue, black, red, white, gray, dgray.
       -no_color_f, --gene_no_color_filled
             To fill uniqe genes (including RNA genes), pseudo genes, and homologous gene clusters not meeting the criteria set by "percent_strain_homologouscluster_color" parameter with a single color (Default: white), users can choose from blue, black, red, white, gray, dgray.
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
             Set the right margin of image in pixels (Default: 20).
       -label, --show_label
             Display the gene label (gene Locus Tag or genename) (Default: T).
       -ul, --unification_label
             Unify gene lable for homologous gene cluster (Default: T). Among a set of homologous gene cluster, if a gene is annotated with a name X, all other genes will be labeled with X.
       -family, --font_family
             Set font family for the genome name and the gene label, e.g. Times New Roman, Arial, Verdana and so on (Default: Times New Roman). Users are suggested to choose font family listed in metrcis module, or causing a miscalculation of string width for genome name in SVG-format map.
       -style, --font_style
             Set font style for the genome name and the gene label, e.g. Normal, Bold, Italic (Default: Normal). it should be noted that the font style "Bold" does not work when using to cearte a PNG format figure in MacOS.
       -size, --label_font_size
             Set font size for gene label (Default: 6).
       -color, --label_font_color
             Set font color for gene label (Default: dgray).
       -i_color, --interested_gene_label_font_color
             Customize gene label color for gene of interest (Default: red), users can choose from blue, black, red, white, gray, dgray.
       -r, --rotate_gene_label (Default: 30)
             Rotate the angle of the gene label, e.g. 30, 45, 135 and so on.
       --strain_name_font_size
             Set font size for genome name (Default: 12).   
       --Strain_name_font_color
             set font color for genome name (Default: black). Users can choose from blue, black, red, white, gray, dgray.
       -Bst, --homologous_gene_cutoff
             Array to set blast parse cutoff: E-value, Identify, Coverage, Match_length (Default: E-value=1-e5, Identify=0, Coverage=50%, Match_length=0).
       -m, --multiple_threads
             Numbers of thread to use (Default: 1).
       -SVG, --SVG_image
             Create SVG format figure (Default: T).
       -PNG, --PNG_image
             Create PNG format figure (Default: T).
       -sub_TFT, --start_at_sub_TFT (Default: F)
             Jump to generate a collection of sub-TFT tables and perform homologous gene analysis (Default: F). Skips sequences extraction and TFT file generation.  
       -map, --start_at_map
             Jump to map generation (Default: F). Generation of a collection of sub-TFT tables and homologous gene clusters has already been done. This parameter is very useful to customize the map quickly. It should be noted that there's no sense to reset "flanking_gene_number" parameter if this parameter set to "T".
             Importantly, at this step, users can revise the gene label by directly edition of the locus_tag in sub_TFT file or all_orthomcl.out. In sub_TFT files and all_orthomcl.out file, there are two forms of gene locus tag, (1) "Locus_Tag", in this case, no genename is defined for a gene; (2) "GeneName;Locus_Tag", in this case, genename is given for a gene. For the first form, user can revise gene label by addition of a genename followed by a semicolon in the front of the Locus_Tag. For the second form, user can revise gene label by modification of the genename.
       -h, --help
             Show this message.


=AUTHOR

Dr. Xiangyang Li (E-mail: lixiangyang\@fudan.edu.cn, lixiangyang1984\@gmail.com), Fudan university; Kaili University; Bacterial Genome Data mining & Bioinformatic Analysis (http://www.microbialgenomic.com/).

=COPYRIGHT

Copyright 2019, Xiangyang Li. All Rights Reserved.

USAGE

my %options = (
    'genbank_file_directory'                       => undef,  
    'interested_gene_file'                         => undef, 
    'Gcluster_output_directory'                    => undef, 
    'phylogenetic_file'                            => undef,
    'show_tree_topology'                           => "T",
    'show_tree_branch'                             => "F",
    'x_step'                                       => 10,  # Draw tree using xstep instead of tree branch length
    'show_tree_bootstrap'                          => "F",
    'strain_reorder_file'                          => undef,  
    'flanking_gene_number'                         => 10,
    'gene_color_filled'                            => "T",    #To color homologous gene clusters or gene families of interest
    'percent_strain_homologouscluster_color'       => "0",    #This parameter is useful when no enough colors are provided in
                                                              #colors_configure_file under color_configure direcoty. Users could
                                                              #to reduce the number of colors needed by setting a high value.
    'cds_color_border'                             => "black", #choose from blue, black, red, white, gray, dgray
    'pseudo_color_border'                          => "dgray", #mark Pseudo genes, choose from blue, black, red, white, gray, dgray
    'RNA_color_border'                             => "red",   #mark RNA genes, choose from blue, black, red, white, gray, dgray
    'gene_no_color_filled'                         => "white", #choose from blue, black, red, white, gray, dgray
    'arrow_relative_Length'                        => 4,
    'arrow_relative_Height'                        => 6,
    'figure_Scale_up_multiple'                     => 0.5,
    'strain_name_shift_Y'                          => 0,
    'gene_label_shift_Y'                           => 2,
    'distance_between_two_genomes'                 => 70,
    'up_shift'                                     => 5,
    'down_shift'                                   => 20,
    'left_shift'                                   => 10,
    'right_shift'                                  => 20,
    'show_label'                                   => "T",
    'unification_label'                            => "T",
    'font_family'                                  => "Times New Roman",  #e.g. Times New Roman, Arial, Verdana and so on.
    'font_style'                                   => "Normal",  #Normal, Bold, Italic
    'label_font_size'                              => 6, 
    'label_font_color'                             => "dgray", #choose from blue, black, red, white, gray, dgray
    'interested_gene_label_font_color'             => "red", #choose from blue, black, red, white, gray, dgray
    'rotate_gene_label'                            => 30,
    'strain_name_font_size'                        => 12, 
    'strain_name_font_color'                       => "black", #choose from blue, black, red, white, gray, dgray
    'homologous_gene_cutoff'                       => "1e-5,0,50,0", #array to set blast parse cutoff: E-value, Identify, Coverage, Match_length
    'multiple_threads'                             => "1",
    'SVG_image'                                    => "T", 
    'PNG_image'                                    => "T", 
    'start_at_sub_TFT'                             => "F",
    'start_at_map'                                 => "F",
    'help'                                         => undef

);
 

GetOptions(
    'dir|genbank_file_directory=s'                   => \$options{genbank_file_directory},    
    'gene|interested_gene_file=s'                    => \$options{interested_gene_file}, 
    'o|Gcluster_output_directory=s'                  => \$options{Gcluster_output_directory},     
    'tree|phylogenetic_file=s'                       => \$options{phylogenetic_file},     
    'topology|show_tree_topology=s'                  => \$options{show_tree_topology}, 
    'branch|show_tree_branch'                        => \$options{show_tree_branch}, 
    'x_step=i'                                       => \$options{x_step}, 
    'bootstrap|show_tree_bootstrap'                  => \$options{show_tree_bootstrap}, 
    'srf|strain_reorder_file=s'                      => \$options{strain_reorder_file},     
    'n|flanking_gene_number=i'                       => \$options{flanking_gene_number},  
    'color_f|gene_color_filled=s'                    => \$options{gene_color_filled}, 
    'pso|percent_strain_homologouscluster_color=f'   => \$options{percent_strain_homologouscluster_color},
    'c_color_b|cds_color_border=s'                   => \$options{cds_color_border}, 
    'p_color_b|pseudo_color_border=s'                => \$options{pseudo_color_border},
    'r_color_b|RNA_color_border=s'                   => \$options{RNA_color_border},  
    'no_color_f|gene_no_color_filled=s'              => \$options{gene_no_color_filled}, 
    'l|arrow_relative_Length=f'                      => \$options{arrow_relative_Length}, 
    'w|arrow_relative_Height=f'                      => \$options{arrow_relative_Height}, 
    'scale|figure_Scale_up_multiple=f'               => \$options{figure_Scale_up_multiple},
    's_Y|strain_name_shift_Y=i'                      => \$options{strain_name_shift_Y}, 
    'g_Y|gene_label_shift_Y=i'                       => \$options{gene_label_shift_Y}, 
    'dis|distance_between_two_genomes=f'             => \$options{distance_between_two_genomes},
    'up|up_shift=f'                                  => \$options{up_shift},
    'down|down_shift=f'                              => \$options{down_shift},
    'left|left_shift=f'                              => \$options{left_shift},
    'right|right_shift=f'                            => \$options{right_shift},
    'label|show_label=s'                             => \$options{show_label},
    'ul|unification_label=s'                         => \$options{unification_label},
    'family|font_family=s'                           => \$options{font_family},         
    'style|font_style=s'                             => \$options{font_style},          
    'size|label_font_size=i'                         => \$options{label_font_size},    
    'color|label_font_color=s'                       => \$options{label_font_color},              
    'i_color|interested_gene_label_font_color=s'     => \$options{interested_gene_label_font_color}, 
    'r|rotate_gene_label=i'                          => \$options{rotate_gene_label}, 
    'strain_name_font_size=f'                        => \$options{strain_name_font_size}, 
    'strain_name_font_color=s'                       => \$options{strain_name_font_color}, 
    'Bst|homologous_gene_cutoff=s'                   => \$options{homologous_gene_cutoff}, 
    'm|multiple_threads=i'                           => \$options{multiple_threads},     
    'SVG|SVG_image=s'                                => \$options{SVG_image}, 
    'PNG|PNG_image=s'                                => \$options{PNG_image}, 
    'sub_TFT|start_at_sub_TFT=s'                     => \$options{start_at_sub_TFT},
    'map|start_at_map=s'                             => \$options{start_at_map},
    'h|help'                                         => \$options{help}

);


if ( defined( $options{help} ) ) {
    print $usage;
    exit(0);
}

#check for required options
if ( !( defined( $options{genbank_file_directory} ) ) ) {
    print $usage;
    exit(1);
}

#check for required options
if ( !( defined( $options{interested_gene_file} ) ) ) {
    print $usage;
    exit(2);
}

#check for coexisting options
if ( defined($options{phylogenetic_file}) && defined($options{strain_reorder_file}) ) {
    print "Warning: Gcluster needs a --strain_reorder_file-- or a --phylogenetic_file--, but not both at the same time.\n";
    exit(3);
}

my $now_time = localtime;
print "\n$now_time: Gcluster.pl start...\n\n";

my $home_directory = $FindBin::Bin;           # obtaining the home directory where Gcluster.pl located
my $image;
my $cluster_result;
# define global variables
my %cluster_color_hash;                    # using in subrouting &cluster_color_hash();
my %cluster_color_hash_2;                  # using in subrouting &cluster_color_hash();

#check for Gcluster.pl workplace options
my $workplace;
if ( defined( $options{Gcluster_output_directory} ) ) {
    $workplace = File::Spec->rel2abs($options{Gcluster_output_directory});
    mkdir $workplace;
}else {

    $workplace = "$home_directory/Gcluster_workplace";
    $workplace =~ s/\/\//\//g;
    mkdir $workplace;
}

my $genbank_file_directory = File::Spec->rel2abs($options{genbank_file_directory});
my $interested_gene_file = File::Spec->rel2abs($options{interested_gene_file});
my $interested_gene_name = basename $interested_gene_file;
my $strain_reorder_file = File::Spec->rel2abs($options{strain_reorder_file}) if defined $options{strain_reorder_file};
my $phylogenetic_file = File::Spec->rel2abs($options{phylogenetic_file}) if defined $options{phylogenetic_file};
my $show_tree_topology = $options{show_tree_topology};


######protein sequence file
my $directory_sequence = "$workplace/Whole_proteome.".basename $genbank_file_directory;            
    mkdir $directory_sequence;
    opendir (DIR_SEQ, $directory_sequence) or die "could not open $directory_sequence";
    my @directory_sequence = readdir DIR_SEQ; 
    @directory_sequence = grep ($_!~/^\./ ,@directory_sequence);
    closedir DIR_SEQ;

######tft file
my $directory_TFT = "$workplace/Whole_TFT.".basename $genbank_file_directory; 
    mkdir $directory_TFT;
    opendir (DIR_TFT, $directory_TFT) or die "could not open $directory_TFT";
    my @directory_TFT = readdir DIR_TFT; 
    @directory_TFT = grep ($_!~/^\./ ,@directory_TFT);
    closedir DIR_TFT;

######Part of the protein sequence file
my $directory_protein_part = "$workplace/Sub_proteome.".basename $genbank_file_directory; 
    mkdir $directory_protein_part;
    opendir (DIR_PRP, $directory_protein_part) or die "could not open $directory_protein_part";
    my @directory_protein_part = readdir DIR_PRP; 
    @directory_protein_part = grep ($_!~/^\./ ,@directory_protein_part);
    closedir DIR_PRP;

######Part of the tft file
my $directory_part_TFT = "$workplace/Sub_TFT.".basename $genbank_file_directory; 
    mkdir $directory_part_TFT;
    opendir (DIR_PTFT, $directory_part_TFT) or die "could not open $directory_part_TFT";
    my @directory_part_TFT = readdir DIR_PTFT; 
    @directory_part_TFT = grep ($_!~/^\./ ,@directory_part_TFT);
    closedir DIR_PTFT;

######Sub blastout file
my $directory_sub_blastout = "$workplace/Sub_blast_out.".basename $genbank_file_directory; 
    mkdir $directory_sub_blastout;
    opendir (DIR_SBO, $directory_sub_blastout) or die "could not open $directory_sub_blastout";
    my @directory_sub_blastout = readdir DIR_SBO; 
    @directory_sub_blastout = grep ($_!~/^\./ ,@directory_sub_blastout);
    closedir DIR_SBO;

#######Homologs cluster file
my $directory_homologs_cluster = "$workplace/blast_homologs_cluster.".basename $genbank_file_directory; 
    mkdir $directory_homologs_cluster;
    opendir (DIR_ORC, $directory_homologs_cluster) or die "could not open $directory_homologs_cluster";
    my @directory_homologs_cluster = readdir DIR_ORC; 
    @directory_homologs_cluster = grep ($_!~/^\./ ,@directory_homologs_cluster);
    closedir DIR_ORC;


    if ( ($options{start_at_map} eq "F") && ($options{start_at_sub_TFT} eq "F") ) {
        system ("rm $directory_sequence/*") if (scalar @directory_sequence > 0);
        system ("rm $directory_TFT/*") if (scalar @directory_TFT > 0);
        system ("rm $directory_protein_part/*") if scalar @directory_protein_part > 0;
        system ("rm $directory_part_TFT/*") if scalar @directory_part_TFT > 0;
        system ("rm $directory_sub_blastout/*") if scalar @directory_sub_blastout > 0;
        system ("rm $directory_homologs_cluster/*") if scalar @directory_homologs_cluster > 0;
    }
    elsif ( ($options{start_at_map} eq "F") && ($options{start_at_sub_TFT} eq "T") ) {
        system ("rm $directory_protein_part/*") if scalar @directory_protein_part > 0;
        system ("rm $directory_part_TFT/*") if scalar @directory_part_TFT > 0;       
        system ("rm $directory_sub_blastout/*") if scalar @directory_sub_blastout > 0;
        system ("rm $directory_homologs_cluster/*") if scalar @directory_homologs_cluster > 0;

    }


    #transform locus tag list of interested gene into a hash using subroutine; 
    my %hash_nif = process::interested_gene_hash($interested_gene_file);
   
    if ($options{start_at_map} eq "T") {
        print "Gcluster-tool: generate a PNG or SVG figure \"###start_at_map###\"\n";
        $options{start_at_sub_TFT} = "T";

    }

    if ($options{start_at_sub_TFT} eq "F" ) {

        #subroutine to extract protein sequence and the Feature Table List (FTL, contain cds/rRNA/tRNA/source) from all genebank files in project folder. 
        process::batch_genbank_sequence_TFT_extract($genbank_file_directory, $directory_sequence, $directory_TFT, $options{multiple_threads});
        print "Step 1-1: Extract all predicted proteomic sequences of each genebank file in folder, and transformat protein sequence into a ID-sequence hash\n";
        print "\nStep 1-2: Extract all feature table of each genebank file in folder\n";
    } 
    elsif ($options{start_at_sub_TFT} eq "T" ) {

        if ($options{start_at_map} eq "F") {
            print "Gcluster-tool: generate a PNG or SVG figure \"###start_at_sub_TFT###\"\n";
        }
    }


    if ($options{start_at_map} eq "F") {

        #subroutine to extract a subpartial lines of TFT table (called Sub-TFT file) around interested gene, the default lines of eachside is 15 if $options{flanking_gene_number} is not provided.
        process::extract_flanking_TFT($directory_TFT, $options{flanking_gene_number}, $directory_part_TFT);   
        print "\nStep 2: Generate a sub-TFT table file containing the gene information flanking the interested gene\n";

        #subroutine to extract protein seqeunce using a sub-TFT table file generated above and all protein sequence files in whole_protein folder.
        if ($options{gene_color_filled} eq  "T") {  # addied by xiangyang Li 2019-10-21
            process::extract_protein_part ($directory_sequence, $directory_part_TFT, $directory_protein_part);
            print "\nStep 3: Extract the interested gene related sub-proteomic sequences according the sub-TFT table files\n";

            #subroutine to do local blastp analysis and group proteins into homologous family using mcl
            if (scalar keys %hash_nif <= $options{multiple_threads}) {
                $options{multiple_threads} = scalar keys %hash_nif;
            }
            use OrthoMCL_analysis;
            #$cluster_result = process::blast_homologs_cluster($directory_sub_blastout, $directory_protein_part, $directory_homologs_cluster, $interested_gene_name, $options{homologous_gene_cutoff}, $options{multiple_threads});
            $cluster_result = OrthoMCL_analysis::orthoMCL($directory_sub_blastout, $directory_protein_part, $directory_homologs_cluster, $interested_gene_name, $options{homologous_gene_cutoff}, $options{multiple_threads}, $blastp, $makeblastdb, $MCL, $inflation);
            
            print "\nStep 4: Obtain the homologous gene clusters information\n";
        }

    }
    else {

        opendir CLUSTER_0, $directory_homologs_cluster or die $!;
        my  @cluster_0 = readdir CLUSTER_0; 
        closedir CLUSTER_0;
    
        foreach (@cluster_0) {

            if (/.out/) {            
                
                $cluster_result = "$directory_homologs_cluster/$_";
                #print "$cluster_result\n";

            }
        }
 
    }
    

#create a png format image 
    if ($options{PNG_image} eq "T") {

    use PNG_image;

    PNG_image::create_image_PNG($home_directory, $directory_part_TFT, $cluster_result, $strain_reorder_file, $phylogenetic_file, $show_tree_topology, \%hash_nif, $workplace, \%options);

    }

#create a SVG format image 
    if ($options{SVG_image} eq "T") {

    use SVG_image;

    SVG_image::create_image_SVG($home_directory, $directory_part_TFT, $cluster_result, $strain_reorder_file, $phylogenetic_file, $show_tree_topology, \%hash_nif, $workplace, \%options);
    }


$now_time = localtime;
print "$now_time: Finished!\n";
