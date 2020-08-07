package PNG_image;

use strict;
use warnings;
use Bio::SeqIO;
use File::Basename qw<basename dirname>;
use FindBin;
use lib "$FindBin::Bin/lib";
use GD;
use Bio::TreeIO;
use Bio::Tree::NodeI;
use process;
use metrics;

use vars qw($image $font_path $courier %hash_nif %color_hash %cluster_color_hash %cluster_color_hash_2 %direction_F_R %shift_distance_x %string_height_png $yellow $green $blue $black $red $white $gray $dgray $tree_width @array_image %cluster_GeneName_hash);


sub create_image_PNG {

    my ($home_directory, $directory_part_TFT, $cluster_result, $strain_reorder_file, $phylogenetic_file, $show_tree_topology, $ref_hash_nif, $workplace, $ref_options) = @_;

    %hash_nif = %$ref_hash_nif;
    my %options = %$ref_options;

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

    my $image_height=$options{distance_between_two_genomes}*(scalar @dir_part_TFT)+$options{up_shift}+$options{down_shift};  # image_height 2019-06-20


    #subroutine to retrieve the lengthest width from FTL table files
    @array_image = process::figure_size_width ($home_directory, $directory_part_TFT, $options{font_family}, $options{font_style}, $options{strain_name_font_size});   
    %direction_F_R = %{$array_image[2]}; 
    %shift_distance_x = %{$array_image[3]};  
    %string_height_png = %{$array_image[6]};

    if (defined $options{phylogenetic_file} && $options{show_tree_topology} eq "T") {
        $tree_width = tree_width($phylogenetic_file, $options{show_tree_branch}, $options{left_shift}, $options{x_step});# caclulate the width and height of the tree
    }else{

        $tree_width = $options{left_shift};
    }   

    my $image_width = ($array_image[0]+$array_image[1])*$options{figure_Scale_up_multiple}/20+$options{right_shift}+$array_image[4]+$tree_width;  # image_height 2019-06-26

    print "\nStep 5(-PNG): Be going to create a size of ", $image_width, "x", "$image_height PNG image\n";

    print "\nStep 6(-PNG): PNG format was chosen to map\n";

        #create a width x height size image  
        $image = GD::Image->new($image_width,$image_height, 1);

        $yellow = $image->colorAllocate(255,255,0);
        $green = $image->colorAllocate(0,100,0);  
        $blue = $image->colorAllocate(0,0,255);
        $black= $image->colorAllocate (0,0,0);
        $red = $image->colorAllocate(255,0,0);
        $white=$image->colorAllocate(255,255,255);
        $gray= $image->colorAllocate (212,212,212);
        $dgray= $image->colorAllocate (153,153,153);
        #$ddgray   = $image->colorAllocate(127,127,127);
        $image->setThickness($options{line_drawing_width}); #set the line drawing width
        $image->filledRectangle(0, 0, $image_width, $image_height, $white); #add a white background color in image
        
        $font_path = "$home_directory/font_configure/6x10.bdf.fnt";
        $courier = GD::Font->load($font_path) or die "Can't load font";
        $image->useFontConfig(1);

    %cluster_GeneName_hash = process::cluster_GeneName_hash($directory_part_TFT, $cluster_result);

    if ($options{gene_color_filled} eq  "T") {

        #subroutine to tranform a list of rgb color value into a color hash

        my $color_rgb_code = "$home_directory/color_configure/colors_configure_file";
        color_product($color_rgb_code);

        #subroutine to related color value with each set of homologus gene cluster
        my @ref_cluster_color_hash = process::cluster_color_hash($cluster_result, $color_rgb_code, \%color_hash, scalar @dir_part_TFT, $options{percent_strain_homologouscluster_color});
        %cluster_color_hash = %{$ref_cluster_color_hash[0]};
        %cluster_color_hash_2 = %{$ref_cluster_color_hash[1]};
        
    }

    print "\nStep 7(-PNG): PNG Figure will be generated, please wait!!!\n";

    # each tbl file is used as input file to generate map using draw_map subroutine 

    my $temp_strain_reorder_file;
    if (defined $options{phylogenetic_file}) {
        $temp_strain_reorder_file = TreeDraw($workplace, $phylogenetic_file, $options{show_tree_topology}, $options{show_tree_branch}, $options{show_tree_bootstrap}, $options{x_step}, $options{left_shift}, $options{up_shift}, $options{distance_between_two_genomes}, $options{font_family},$options{font_style}, $options{strain_name_font_size});  # display the tree
    }

    my $eachfile=0;
    if (defined $options{strain_reorder_file} || defined $options{phylogenetic_file}) {
        my %strain_reorder_hash;
        if (defined $options{phylogenetic_file}) {open (STRAIN_REORDER, $temp_strain_reorder_file) or die "Must supply strain_reorder_file\n";} #$temp_strain_reorder_file
        if (defined $options{strain_reorder_file}) {open (STRAIN_REORDER, $strain_reorder_file) or die "Must supply strain_reorder_file\n";} #$strain_reorder_file

    while (<STRAIN_REORDER>){
        my @strain_reorder = split '\t', $_;
        $strain_reorder_hash{$strain_reorder[1]} = $strain_reorder[0];

    }

    foreach my $part_reorder (sort {$a<=>$b} keys %strain_reorder_hash){ 


        my $check_value = $strain_reorder_hash{$part_reorder};
        
   
        foreach my $tbl_part_file(sort @dir_part_TFT){ 
            my $tbl_part_file_new = $tbl_part_file;
            $tbl_part_file_new =~ s/.tbl_.*//g;                
	    
            if(($tbl_part_file_new eq $check_value) and ($tbl_part_file =~ /part$/)) { 
                $tbl_part_file =~ s/(.tbl_.*)//g;
                if ($tbl_part_file eq $check_value) {
                    $tbl_part_file .=$1; 
                    #print "$tbl_part_file\n";  
              
                    $eachfile++;
                    my $Y_parameter=$options{up_shift}+$options{distance_between_two_genomes}*$eachfile;    # strart Y 
                    my $tbl_part_infile = "$directory_part_TFT/$tbl_part_file" or die "Must supply input filename\n";
        
                    draw_map_PNG ($home_directory, $tbl_part_infile, $Y_parameter, $tree_width, $options{arrow_relative_Length}, $options{arrow_relative_Height}, $options{strain_name_shift_Y},$options{gene_label_shift_Y}, $options{font_family}, $options{font_style}, $options{label_font_size}, $options{label_font_color}, $options{interested_gene_label_font_color}, $options{cds_color_border}, $options{pseudo_color_border}, $options{RNA_color_border}, $options{gene_no_color_filled}, $options{figure_Scale_up_multiple}, $options{rotate_gene_label}, $options{strain_name_font_size}, $options{strain_name_font_color}, $options{show_label}, $options{unification_label});

                }

            }
    
         }

    }

}

    if (!defined $options{phylogenetic_file} && !defined $options{strain_reorder_file}) {

        foreach my $tbl_part_file(sort @dir_part_TFT){
        unless ($tbl_part_file =~ /.part$/){ next;}
        $eachfile++;
        my $Y_parameter=$options{up_shift}+$options{distance_between_two_genomes}*$eachfile;     # left_shift
        my $tbl_part_infile = "$directory_part_TFT/$tbl_part_file" or die "Must supply input filename\n";
        
            draw_map_PNG ($home_directory, $tbl_part_infile, $Y_parameter, $options{left_shift}, $options{arrow_relative_Length}, $options{arrow_relative_Height}, $options{strain_name_shift_Y},$options{gene_label_shift_Y}, $options{font_family}, $options{font_style}, $options{label_font_size}, $options{label_font_color}, $options{interested_gene_label_font_color}, $options{cds_color_border}, $options{pseudo_color_border}, $options{RNA_color_border}, $options{gene_no_color_filled}, $options{figure_Scale_up_multiple}, $options{rotate_gene_label}, $options{strain_name_font_size}, $options{strain_name_font_color}, $options{show_label}, $options{unification_label});

        }

    }
    open (PNG_MAP, ">$workplace/figure.png") or die("Cannot open file for writing");

    binmode PNG_MAP;

    print PNG_MAP $image->png;
    close (PNG_MAP);
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
###### Subrounting--angle
###### Function:
###### transform pi to angle
##################################################################################################
sub angle {
    my $rotate_angle = shift;
    use Math::BigFloat;
    my $pi = Math::BigFloat->bpi(100);
    
    my $angle =($pi/180)*$rotate_angle;

return $angle;

}


##################################################################################################
###### Subrounting--draw_map_PNG
###### Function:
###### main subrouting to draw gene cluster map
##################################################################################################
sub draw_map_PNG {

    my %given_color_hash = (
        'blue'                 => $blue, 
        'red'                  => $red,
        'black'                => $black,
        'white'                => $white,
        'gray'                 => $gray,
        'dgray'                => $dgray
       );

    my $home_directory = shift;
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

    my $filename = basename $inputfile;
    my $filename_shift = basename $inputfile; 
    $filename =~ s/.tbl(.*)//g;
    $filename =~ s/_/ /g;

    $image->stringFT($given_color_hash{$strain_name_font_color}, "$font_family:$font_style", $strain_name_font_size, 0, 10+$left_shift, $Y_centre+$string_height_png{$filename}/3-$strain_name_shift_Y, $filename); # black color for genomic names $options{strain_name_font_size}

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

        if ($direction_F_R{$filename_shift} ge 1) {
 
            $strart=($cds[0]-$min[0]+$shift_distance_x{$filename_shift})*$figure_Scale_up_multiple/20+$left_shift+$array_image[4]+5;
            $end=($cds[1]-$min[0]+$shift_distance_x{$filename_shift})*$figure_Scale_up_multiple/20+$left_shift+$array_image[4]+5;
 
        }

        else {

            $strart=(abs($cds[0]-$min[$#min])+$shift_distance_x{$filename_shift})*$figure_Scale_up_multiple/20+$left_shift+$array_image[4]+5; 
            $end=(abs($cds[1]-$min[$#min])+$shift_distance_x{$filename_shift})*$figure_Scale_up_multiple/20+$left_shift+$array_image[4]+5;   
        }

        #my $left=($end+$strart)/2;          
        #my $right=$Y_centre-1.25*$arrow_high;
        #my $tag=$cds[3];

        ### gene_label mark
        if ($show_label eq "T") {

        if (defined ($hash_nif{$tag})) {

            #$highlight_gene = $hash_nif{$tag};

            $image->stringFT($given_color_hash{$interested_gene_label_font_color}, "$font_family:$font_style", $label_font_size, angle($rotate_gene_label), ($end+$strart)/2, $Y_centre-$arrow_high-$gene_label_shift_Y, $geneName_change_hash{$tag});  # highligth interested gene

        }
        elsif (!(defined ($hash_nif{$tag}))){           
            #$tag =~ s/;.*//g;
            $image->stringFT($given_color_hash{$label_font_color}, "$font_family:$font_style", $label_font_size, angle($rotate_gene_label), ($end+$strart)/2, $Y_centre-$arrow_high-$gene_label_shift_Y, $geneName_change_hash{$tag});   # Default color is dgray for gene label of the other genes
         
        }

        }

        # CDS: forward direction
        if ($end>$strart){                                            

            if(exists $cluster_color_hash_2{$tag}){      

                $image->filledPolygon(rec_arrow_F($strart, $end, $arrow_length, $arrow_high, $Y_centre), $cluster_color_hash{$tag});
                $image->polygon(rec_arrow_F($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$cds_color_border}); #border

            }
            else{ 

                $image->filledPolygon(rec_arrow_F($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$gene_no_color_filled});
                $image->polygon(rec_arrow_F($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$cds_color_border}) if $cds[2] eq "CDS"; #set color of border for CDS genes
                $image->polygon(rec_arrow_F($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$pseudo_color_border}) if scalar @cds > 5; #set color of border for Pseudo genes
                $image->polygon(rec_arrow_F($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$RNA_color_border}) if $cds[2] =~ /^rRNA|^tRNA/; #set color of border for RNA genes
 
            }                 # map
        }

        # CDS: reverse direction
        elsif ($end<$strart){                           
    
            if(exists $cluster_color_hash_2{$tag}) { 
           
                $image->filledPolygon(rec_arrow_R($strart, $end, $arrow_length, $arrow_high, $Y_centre), $cluster_color_hash{$tag});
                $image->polygon(rec_arrow_R($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$cds_color_border}); #border

            }
            else{ 

                $image->filledPolygon(rec_arrow_R($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$gene_no_color_filled});
                $image->polygon(rec_arrow_R($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$cds_color_border}) if $cds[2] eq "CDS"; #set color of border for CDS genes
                $image->polygon(rec_arrow_R($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$pseudo_color_border}) if scalar @cds > 5; #set color of border for Pseudo genes
                $image->polygon(rec_arrow_R($strart, $end, $arrow_length, $arrow_high, $Y_centre), $given_color_hash{$RNA_color_border}) if $cds[2] =~ /^rRNA|^tRNA/; #set color of border for RNA genes

            }                     # map
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
    $rec_arrow->addPt($e, $Y_c);                  # piont1
    $rec_arrow->addPt($e-$a_l, $Y_c-$w_h);        # piont2
    $rec_arrow->addPt($e-$a_l, $Y_c-$w_h/2);      # piont3
    $rec_arrow->addPt($s, $Y_c-$w_h/2);           # piont4
    $rec_arrow->addPt($s, $Y_c+$w_h/2);           # piont5
    $rec_arrow->addPt($e-$a_l, $Y_c+$w_h/2);      # piont6
    $rec_arrow->addPt($e-$a_l, $Y_c+$w_h);        # piont7

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

my ($workplace, $phylogenetic_file, $show_tree_topology, $use_branch, $bootstrap, $xstep, $left_shift, $up_shift, $ystep, $font, $font_style, $font_size) = @_;

my $temp_strain_reorder_file = "$workplace/temp_strain_reorder_file-png.txt";

my %xx;                 # horizontal coordinate for each node
my %yy;                 # vertical coordinate for each node
my $tree_object;        # first Bio::Tree::Tree object
#my $xstep = 20;         # branch length in drawing, which should keep according with the $xstep in sub Tree_wirth     
my $tip = 5;            # extra space between tip and label

##########################################################################################################
    my $treeio = Bio::TreeIO->new(-format => 'newick',
                                  -file   => $phylogenetic_file); 
    $tree_object = $treeio->next_tree;  
    my @taxa1 = $tree_object->get_leaf_nodes;
    my $root1 = $tree_object->get_root_node;

    open (STRAIN_REORDER_FILE, ">$temp_strain_reorder_file");
    my $count_number=0;
    my @strain_order = reverse $tree_object->get_nodes(-order => 'depth');
    pop @strain_order; # skip root

    my $y = $up_shift+$ystep;
    for my $node_leaf (@strain_order) {
        #print $node_leaf->id;
        if ($node_leaf->is_Leaf){
            $count_number++;
            print STRAIN_REORDER_FILE $node_leaf->id, "\t", "$count_number", "\n";

            #print "test:   ", $node_leaf->id, "\n";
            $xx{$node_leaf} = 0;                 # a temp value
            $yy{$node_leaf} = $y;                # leaf y-value
            $y += $ystep;
        }
    }
    close STRAIN_REORDER_FILE;

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

                $yy{$node} = ($ymin + $ymax)/2;
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
            #$image->stringFT($green, "$font:$font_style", $font_size/1.5, 0, $xx{$taxon} + $tip, $yy{$taxon}+$font_size / 3, $taxon->id); # draw strain 
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

                        $image->stringFT($black, "$font:$font_style", $font_size*0.8, 0, $xx{$node->ancestor}+ $font_size/10, $yy{$node->ancestor}+ ($font_size / 3), $bootstrap_value); # draw bootstrap value####################################
                        
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
    return $temp_strain_reorder_file;
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

my ($phylogenetic_file, $use_branch, $left_shift, $xstep) = @_;
my %xx;        # horizontal coordinate for each node
my $width;     # total drawing width

    my $treeio = Bio::TreeIO->new(-format => 'newick',
                                  -file   => $phylogenetic_file); 
    my $tree_obj = $treeio->next_tree;  # first Bio::Tree::Tree object
    my @taxa = $tree_obj->get_leaf_nodes;
    my $root = $tree_obj->get_root_node;

    for my $taxon (reverse @taxa) {     #strain_name in Y-aixs    xiangyang Li 2019-07-17
        $xx{$taxon} = 0;  # a temp value

    }

    #my $xstep = 20;  # branch length in drawing

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

    return $width;

}


1;

__END__  
