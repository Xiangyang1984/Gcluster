package process;

use strict;
use warnings;
use threads;
use GD;
use Bio::SeqIO;
use File::Basename qw<basename dirname>;
use FindBin;
use lib "$FindBin::Bin/lib";

use vars qw(%hash_nif);


##################################################################################################
###### Subrounting--interested_gene_hash
###### Function:
###### tranformat into hash
##################################################################################################
sub interested_gene_hash {
    my $interested_gene_file = shift;
#    print "$interested_gene_file\n";
    open (NIFH, $interested_gene_file);
    while (<NIFH>){
        chomp $_;
        $_ =~ s/\t//g;                 
        $_ =~ s/#.*//g;                 # delete the explanation of the interested gene following marker "#"
        $hash_nif{$_}=$_;
    }
    
    return %hash_nif;

}


##################################################################################################
###### Subrounting--do bacth work
###### Function:
###### do bacth work
##################################################################################################
sub batch_genbank_sequence_TFT_extract {
    my ($path_genbank, $path_protein, $path_TFT, $thread_number)=@_;

    opendir PATH_GENBANK, $path_genbank or die "could not open $path_genbank";
    my @path_genbank_temp = readdir PATH_GENBANK;
    foreach (@path_genbank_temp){
        if (/[ =]/) {
            my $new_name = $_;
            $new_name =~ s/ /_/g;   #repalce the space with "_"  for all GenBank files
            $new_name =~ s/=/_/g;   #repalce "=" with "_"  for all GenBank files
            $new_name =~ s/_+/_/g;  
            system ("mv '$path_genbank/$_' '$path_genbank/$new_name'");  
        }

    }
    closedir PATH_GENBANK;

    opendir PATH_GENBANK, $path_genbank or die "could not open $path_genbank";
    my @path_genbank = readdir PATH_GENBANK;
    @path_genbank =grep ($_!~/^\./ ,@path_genbank);  #delete hidden file . .. 
    #@path_genbank =grep ($_!~/$\~/ ,@path_genbank);  #delete hidden file . .. 
    closedir PATH_GENBANK;

    my $sub_file_number = scalar @path_genbank;

    my @input;
    my @outfile_1;
    my @outfile_2;
 
    foreach my $file_genbank(@path_genbank){ 
        $file_genbank =~ s/ /_/g;
        my $input="$path_genbank/$file_genbank"; 
        push (@input,$input);  
        @input = sort @input;  

        my $output_1="$path_protein/$file_genbank.fasta";
        push (@outfile_1,$output_1);
        @outfile_1 = sort @outfile_1;
        my $output_2="$path_TFT/$file_genbank.tbl";
        push (@outfile_2,$output_2);
        @outfile_2 = sort @outfile_2;

    }

    my $thread;
    my @threads;
    my $job_number=0;
    print "GenBank_extraction_percent: ";
    while(){ 

        last if ($job_number eq $sub_file_number);  

        while(scalar(threads->list())<$thread_number) {     #set thread number
            $job_number++;    
            my $input_file = $input[$job_number-1];
            my $output_file_1 = $outfile_1[$job_number-1];
            my $output_file_2 = $outfile_2[$job_number-1];
            my $GenBank_extraction_progress = int (($job_number/$sub_file_number)*100);
            print "$GenBank_extraction_progress","...";
            $threads[$job_number-1]=threads->new(\&genbank_sequence_TFT_extract, $input_file, $output_file_1, $output_file_2); 

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
###### extract protein sequence and tbl (cds/rRNA/tRNA information) from genbank file
##################################################################################################
sub genbank_sequence_TFT_extract {

    my ($input, $output_1, $output_2)=@_;
    
    #record contigID for accession number
    open(GINPUT, $input);
    my @accession_num;
    while(<GINPUT>){
        chomp;
        my $locus = $1 if $_ =~ /^LOCUS       (.*?) /;
        #print "$locus\n" if $_ =~ /^LOCUS       (.*?) /;
        push @accession_num, $locus if $_ =~ /^LOCUS       (.*?) /;
    }
    close GINPUT;
    
    open(OUTPUT_1, ">$output_1");
    open(OUTPUT_2, ">$output_2");
        my $filen = basename $input;
        my $in = new Bio::SeqIO( -file => $input, -format => 'genbank' );
        my $contig_index = -1;
        while ( my $seqObject = $in->next_seq() ) {
            $contig_index++;
            my $acc = $seqObject->accession;
            $acc = $accession_num[$contig_index] if $acc eq "unknown";
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
                    }else {$locus_tag = "$filen"."_".$count;}
               
                    my $product;
                    if ( $feat->has_tag('product') ) {
                        ($product) = $feat->get_tag_values('product');
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
 
                    my $protein_seqeunce;
                    if ( $feat->has_tag('translation')) {
                        ($protein_seqeunce) = $feat->get_tag_values('translation');
                        if (defined $gene_name){
                            print OUTPUT_1 ">$gene_name;$locus_tag\n$protein_seqeunce\n";
                        }else {
                            print OUTPUT_1 ">$locus_tag\n$protein_seqeunce\n";
                         }
                    }
                }
            }
        }
 #   }

}


##################################################################################################
###### Subrounting--extract_flanking_TFT
###### Function:
###### extract cds/rRNA/tRNA information flanking interested gene using tbl file
##################################################################################################
sub extract_flanking_TFT {
    my $path_TFT=shift;
    my $gene_number_show=shift;
    my $directory_part_TFT=shift;
    opendir PATH_TFT, $path_TFT or die $!;
    my @path_TFT = readdir PATH_TFT; 
    closedir PATH_TFT;   

    #reading each tbl file
    foreach my $file_TFT(@path_TFT){
        my $tbl_infile = "$path_TFT/$file_TFT" or die "Must supply input filename\n";
        open (IN, $tbl_infile) or die "Unable to open $tbl_infile\n";
       
        #output rang_line
        my $key=0;
        my $Feature=0;
        my @mark_feature=();
        my @mark_line=();
        while (<IN>){
            chomp;  
            $key++;
            $Feature++;
            if (/>Contig:.*/){
                push (@mark_feature, $Feature);
            }
            elsif (!(/>Contig:.*/)){
                my @tag= split "\t", $_;
                my $cds_locus_tag = $tag[3];
                $cds_locus_tag =~ s/.*;//g;   # addied by xiangyang Li on 2019-10-21
                if ((exists $hash_nif{$cds_locus_tag}) or (exists $hash_nif{$tag[3]})) {
                    push (@mark_line,$key);
                }     
            } 
        }      

    #output rang_TFT file
        foreach my $k (@mark_line){
            my  $Get_min_value;
            if (scalar @mark_feature >1){
                $Get_min_value=Get_min_array($k, @mark_feature);
            }
            else {
                $Get_min_value=1;
            }

            my $tbl_part = "$directory_part_TFT/$file_TFT\_$k.part";
            open (TBL_OUT_PART, ">$tbl_part")  or die "Unable to open feature table output file\n";
            my $line=1;
                if ($k-$Get_min_value>$gene_number_show){ 
                    seek(IN,0,0);
                    while (<IN>){
                        chomp;  
                        if ($line >=$k-$gene_number_show && $line<=$k+$gene_number_show) { 
                            print TBL_OUT_PART "$_\n";
                            last if ($_=~ />Contig:.*/) ;
                        }  
                        $line++;  
                }  
            }
            elsif ($k-$Get_min_value<=$gene_number_show){  
                my $jump=0;
                $line=1;
                seek(IN,0,0);
                while (<IN>){
                    chomp;  
                    if ($line >=$k-$gene_number_show && $line<=$k+$gene_number_show && $k-$gene_number_show> 0) { 
                        $jump++;
                        if ($jump>$gene_number_show+$Get_min_value-$k+1) {
                            print TBL_OUT_PART "$_\n";  
                            last if ($_=~ />Contig:.*/) ;
                        }
                    } 
                    elsif ($line >$Get_min_value && $line<=$k+$gene_number_show && $k-$gene_number_show <= 0) { 
                        
                            print TBL_OUT_PART "$_\n";  
                            last if ($_=~ />Contig:.*/) ;
                    }  
                    $line++;  
                }  
            }
        }
    }

    system ("perl -pi -e 's/^>Contig:.*//' $directory_part_TFT/*.part"); 
    system ("perl -pi -e 's/^\n//' $directory_part_TFT/*.part"); 

}


##################################################################################################
###### Subrounting--Get_min_array
###### Function:
###### obtain the element from a array， which has a minimal distance with a given numerical string
##################################################################################################
sub Get_min_array {
    my $mark = shift;
    my @array = @_;
    my $Get_min;
    @array=sort {$a<=>$b}@array;    
        if ($mark>=$array[$#array]){
            $Get_min=$array[$#array];
        } 
        elsif ($mark<$array[$#array]){
            for (my $i=0; $i<scalar @array; $i++){
                if($mark-$array[$i]<0) {
                    $Get_min=$array[$i-1];
                    last;
                }
            }
        }
 
   return $Get_min;

}


##################################################################################################
###### Subrounting--extract_protein_part
###### Function:
###### extract protein seqeunce according to part_TFT file
##################################################################################################
sub extract_protein_part{
  my ($path_protein, $path_part_TFT, $path_protein_part) = @_;
  my %hash_protein;

    # transform protein sequence files in folder into a ID-sequence hash (%hash_protein)
    opendir PATH_PROTEIN, $path_protein or die "could not open $path_protein";
    my @path_protein = readdir PATH_PROTEIN; 
    closedir PATH_PROTEIN;
    foreach my $file_protein(@path_protein) { 
        my $input_protein="$path_protein/$file_protein"; 
        open FASTA, $input_protein or die $!;
        my @sequence=<FASTA>;
        my $length= scalar @sequence;
        for (my $i=0; $i<$length; $i=$i+2) { 
            $sequence[$i]=~ s/[\r\n]$//;
            $hash_protein{$sequence[$i]}=$sequence[$i+1];
        }

    }

    #Extract sequence from the ID-sequence hash (%hash_protein) according to a patial FTL list
    opendir PATH_PART_TFT, $path_part_TFT or die $!;
    my @path_part_TFT = readdir PATH_PART_TFT; 
    @path_part_TFT =grep ($_!~/^\./ ,@path_part_TFT);  #delete hidden file . ..
    closedir PATH_PART_TFT;

    foreach my $file_part_TFT(@path_part_TFT) {
        my $input_TFT="$path_part_TFT/$file_part_TFT";
        my $output_fatsa="$path_protein_part/$file_part_TFT.fasta";
        open(OUTPUT_FASTA, ">$output_fatsa");
        open LIST, $input_TFT or die $!;
         
            while (<LIST>) {
                chomp;
                my @locus_tag=split ("\t", $_);
                $locus_tag[3]=">".$locus_tag[3];
               
                if (exists $hash_protein{$locus_tag[3]} ){
                    print OUTPUT_FASTA $locus_tag[3],"\n",$hash_protein{$locus_tag[3]};
                }
            }

    }

}


##################################################################################################
###### Subrounting--blast_homologs_cluster
###### Function:
###### Cluster analysis of homologs protein
##################################################################################################
sub blast_homologs_cluster {
    
    my ($directory_sub_blastout, $path_protein_part, $directory_homologs_cluster, $interested_gene_name, $homologous_gene_cutoff, $thread_number) = @_;

    my ($e_value, $identify, $coverage, $match_length) = split ",", $homologous_gene_cutoff;
    
    system ("cat $path_protein_part/*.fasta > $directory_sub_blastout/$interested_gene_name.fasta");
    system ("makeblastdb -in $directory_sub_blastout/$interested_gene_name.fasta -dbtype prot -logfile $directory_sub_blastout/$interested_gene_name.fasta.makeblast.log");

##################################################################################################

    chdir $path_protein_part;
    
    my @input;
    my @outfile;
    opendir PATH_PROTEIN_PART, $path_protein_part or die "could not open $path_protein_part";
    my @path_protein_part = readdir PATH_PROTEIN_PART; 
    @path_protein_part =grep ($_!~/^\./ ,@path_protein_part);  #delete hidden file . ..
    closedir PATH_PROTEIN_PART;

    foreach my $each_protein_file(@path_protein_part){
 
        my $input="$path_protein_part/$each_protein_file"; 
        push (@input,$input);  
        @input = sort @input;
        my $output="$directory_sub_blastout/$each_protein_file.blastout";
        push (@outfile,$output);
        @outfile = sort @outfile;

    }

    my $sub_file_number = scalar @input;
    my $sub_function = "do_blastP";
    my @sub_function_parameters = (\@input, "$directory_sub_blastout/$interested_gene_name.fasta", \@outfile, $e_value);

    bacth_blast($sub_file_number, $thread_number, "do_blastP", \@sub_function_parameters);

    system ("cat $directory_sub_blastout/*.blastout > $directory_sub_blastout/all_vs_all.blast");

    my $all_vs_all_cluster = blast_filter($directory_sub_blastout, "$directory_sub_blastout/all_vs_all.blast", $e_value, $identify, $coverage, $match_length);

    system("mcl $all_vs_all_cluster --abc -V all -I 1.5 -o $directory_homologs_cluster/$interested_gene_name.cluster");
    #system("rm  $interested_gene_name.fasta");#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    my $cluster_result="$directory_homologs_cluster/$interested_gene_name.cluster";

    return $cluster_result;
}

##################################################################################################
###### Subrounting--do_blastP
###### Function:
###### do blastp analysis
##################################################################################################
sub do_blastP {
    my ($input, $db_file, $outfile, $e_value) =@_;    

    system("blastp -outfmt '6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore' -evalue $e_value -num_threads 1 -query $input -db $db_file -out $outfile");
    #return $outfile;

}


##################################################################################################
###### Subrounting--do_batch_blastP
###### Function:
###### do bacth blastp analysis
##################################################################################################
sub bacth_blast {

    my ($work_number, $thread_number, $subfunction, $sub_function_parameters) = @_;
    my ($input_file, $db_file, $output_file, $e_value) = @$sub_function_parameters;
 
    my @input = @$input_file;
    my @outfile = @$output_file;
    my $thread;
    my @threads;
    my $job_number=0;

    while(){ 
        last if ($job_number>=$work_number);                         
        while(scalar(threads->list())<$thread_number) {     #set threadnumber；
            $job_number++;                                 
            my $input_file = $input[$job_number-1];
            my $output_file = $outfile[$job_number-1];
            print "Blast_perform: $job_number\n";
            $threads[$job_number-1]=threads->new(\&$subfunction, $input_file, $db_file, $output_file, $e_value);   
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


##################################################################################################
###### Subrounting--blast_filter
###### Function:
###### blast filter parse
##################################################################################################
sub blast_filter {

    my ($directory_sub_blastout, $all_vs_all_out, $e_value, $identify, $coverage, $match_length) = @_;
    open (INFILE, $all_vs_all_out);
    my $all_vs_all_cluster = "$directory_sub_blastout/all_vs_all.cluster";
    open (OUTFILE, ">$all_vs_all_cluster");

    while (<INFILE>){
        chomp;
        my @blast_result=split '\t', $_;

        my $qcover =  ($blast_result[7] - $blast_result[6] +1)/$blast_result[3] * 100;  #added by xiangyang
        $qcover = sprintf("%.2f", $qcover);
        my $scover = ($blast_result[9] - $blast_result[8] +1)/$blast_result[4] * 100;  #added by xiangyang
        $scover = sprintf("%.2f", $scover);

        if (($blast_result[10] <= $e_value) && ($blast_result[2] >= $identify) && ($qcover >= $coverage) && ($scover >= $coverage) && ($blast_result[5] >= $match_length)){  

            #set cut-off of e-value, identity, coverage, match_length for homologous protein
            if ($blast_result[5]>0) {$blast_result[11] = $blast_result[11]/$blast_result[5];} else {$blast_result[11] = $blast_result[11];}

            print OUTFILE "$blast_result[0]\t$blast_result[1]\t$blast_result[11]\t$blast_result[10]\t$blast_result[2]\t$qcover\t$qcover\t$blast_result[5]\n";

        }

    }

    close OUTFILE;
    return $all_vs_all_cluster;

}


##################################################################################################
###### Subrounting--genome_number_hgc
###### Function:
###### calculate genome number, which are included in a set of homologous gene cluster.
# ***protein files used for homologous gene analysis must come from Gcluster, or an error may occurs
##################################################################################################
sub genome_number_hgc {

    my $each_hgc = shift;

    chomp $each_hgc;
    $each_hgc =~ s/.*: //g;
    my @hgc = split (" ", $each_hgc );
    my @copy_hgc = @hgc;
    my %de_repeat;

    foreach (@hgc) {
        $_ =~ s/.*;//g;
        $_ =~ s/_.*//g;   
        $de_repeat{$_} = $_;
    }

    my $genome_number = keys %de_repeat;

    return (\@copy_hgc, $genome_number);

}



##################################################################################################
###### Subrounting--cluster_color_hash
###### Function:
###### related color value with each set of homologous gene cluster 
##################################################################################################
sub cluster_color_hash {
    my $cluster_result = shift;
    my $color_rgb_code = shift; 
    my $ref_color_hash = shift;
    my $total_strain_number = shift;
    my $percent_strain_homologouscluster_color = shift;
    my %color_hash = %$ref_color_hash;
    my %cluster_color_hash_2;
    my %cluster_color_hash;
    my $value=0;    

    open (CLUSTER, $cluster_result);
    while (<CLUSTER>){
        chomp;
        my ($hgc_ref, $genome_number) = genome_number_hgc($_); # to dealwith cluster result from orthoMCL_analysis_out 2019-11-25
        my @cluster = @$hgc_ref;
        my $ture_percent = $genome_number/$total_strain_number*100;

        #only color homologous gene cluster with gene number >=2
        if ((scalar @cluster > 1) && ($ture_percent >=$percent_strain_homologouscluster_color)) {   
            $value++;
            foreach (@cluster) {

                $_ =~ s/.*\|//g;  # delete "strain name|" for OrthoMCL with mysql
                $_ =~ s/.*;//g;   # delete the geneName from "geneName;Locus_tag"
                if (scalar keys %color_hash >= $value) {
                    $cluster_color_hash{$_}=$color_hash{$value};
                    $cluster_color_hash_2{$_}=$_;
                }
                
            }
        }
    
    }
    close CLUSTER;

    print "\nWarning: a total of $value homologous gene cluster need to be colored, but only ", scalar keys %color_hash," rgb colors are provided in \"$color_rgb_code\", Please provide enough rgb colors\n" if $value > scalar keys %color_hash;

    return \%cluster_color_hash, \%cluster_color_hash_2;

}



##################################################################################################
###### Subrounting--cluster_GeneName_hash
###### Function:
###### To uniform gene name for homologs, and to relate the modified genename in sub_TFT file 
###### with homologous gene cluster
##################################################################################################
sub cluster_GeneName_hash {
    my ($directory_part_TFT, $cluster_result) = @_;
    my %cluster_GeneName_hash;    
    my %TFT_GeneName_hash;    

    opendir DIR_PART_TFT, $directory_part_TFT or die $!;  # reopen $directory_part_TFT to caculate the number of sub-TFT file
    my @dir_part_TFT = readdir DIR_PART_TFT; 
    @dir_part_TFT = grep ($_!~/^\./ ,@dir_part_TFT);  #delete hidden file . ..
    closedir DIR_PART_TFT;

    foreach my $TFT_file (@dir_part_TFT){
        open (TFT_FILE_1, "$directory_part_TFT/$TFT_file");
        while (<TFT_FILE_1>){
            my @TFT_gene_genename = split "\t", $_;
            $TFT_GeneName_hash{$TFT_gene_genename[3]} = $1 if $TFT_gene_genename[3] =~ s/(.*);//g;
        }
    }

    open (CLUSTER, $cluster_result);
 
    while (<CLUSTER>){
        chomp;

        $_ =~ s/.*: //g;   # to dealwith cluster result from orthoMCL_analysis 2019-11-25
        $_ =~ s/.*\|//g;   # delete "strain name|" for OrthoMCL with mysql

        my @cluster=split(" ", $_);

        if (scalar @cluster > 1) {   #only color homologous gene cluster with genes number >=2
            my @tmp_array;
            foreach my $first_clycle (@cluster) {
                push (@tmp_array, $1) if $first_clycle =~ s/(.*);//g;   # delete the GeneName from "GeneName;Locus_Tag" 
            }

            my $TFT_genename; # modified genename in sub_TFT file, it is used to relate modified genename in sub_TFT file with homologous gene cluster

            if (scalar @tmp_array >0) {# to obtain the modified genename in sub_TFT file, when a set of homologous gene cluster has  genename               
                foreach my $fst_cycle_0 (@cluster) {
                    $fst_cycle_0 =~ s/.*;//g;   
                    if ( (defined $TFT_GeneName_hash{$fst_cycle_0}) && (!grep {$TFT_GeneName_hash{$fst_cycle_0} eq $_ }@tmp_array) ) {
                   
                        $TFT_genename =  $TFT_GeneName_hash{$fst_cycle_0};

                    }

                }

            }else{# to obtain the modified genename in sub_TFT file, when a set of homologous gene cluster has not genename 

                foreach my $fst_cycle_1 (@cluster) {

                    $fst_cycle_1 =~ s/.*;//g;  
                    if (defined $TFT_GeneName_hash{$fst_cycle_1}){
                   
                        $TFT_genename =  $TFT_GeneName_hash{$fst_cycle_1};

                    }

                }

            }


            if (scalar @tmp_array >0) { 
                my $final_geneName = max_occurence_ele(\@tmp_array); 
                foreach my $sec_cycle_0 (@cluster) {
                    $sec_cycle_0 =~ s/.*;//g;  
                    $cluster_GeneName_hash{$sec_cycle_0} = $TFT_genename if defined $TFT_genename; # modified genename in sub_TFT file instead of genename in a set of homologous gene cluster
                    $cluster_GeneName_hash{$sec_cycle_0} = $final_geneName if !defined $TFT_genename;  

                }
            }else{ # modified genename in sub_TFT file instead of locus tag in a set of homologous gene cluster, which has not genename
                foreach my $sec_cycle_1 (@cluster) {
                    $sec_cycle_1 =~ s/.*;//g;  
                    $cluster_GeneName_hash{$sec_cycle_1} = $TFT_genename if defined $TFT_genename; 


                }

            }

        }
    
    }
    
    close CLUSTER;
    return %cluster_GeneName_hash;

}



##################################################################################################
###### Subrounting--max_occurence_ele
###### Function:
###### getting the element, which have max occurence in an array
##################################################################################################
sub max_occurence_ele {
    my $array_occurence_0 = shift;
    my @array_occurence = @$array_occurence_0;
    my %hash_occurence;
    my $max_occurence_ele;

    foreach my $ele (@array_occurence){
        $hash_occurence{$ele} += 1;
    }

    my $max_ele = (sort{$b<=>$a} values %hash_occurence)[0];

    foreach (keys %hash_occurence) {

        $max_occurence_ele = $_ if $hash_occurence{$_} == $max_ele;
    }

    return $max_occurence_ele;

}


##################################################################################################
###### Subrounting--figure_size_width
###### Function:
###### retrieved the max_left_length and max_rigth_length around the start position of interested gene, and then caclutate the width of figure size
##################################################################################################
sub figure_size_width {

    my %start_right;
    my @max_start_position;
    my @max_width_right;
    my %shift_distance_temp;
    my %direction_F_R;
    my %shift_distance_x;
    my @interested_gene_length;
    my @text_name_width_svg;
    my @text_name_width_png;
    my %string_height_png;

    my ($home_directory, $path_part_TFT, $font_family, $font_style, $strain_name_font_size) = @_;

    chdir $path_part_TFT;

    opendir DIR_PART_TFT, $path_part_TFT or die "could not open $path_part_TFT";
    my  @dir_part_TFT = readdir DIR_PART_TFT; 
    closedir DIR_PART_TFT;

# define a fictitious image for caculation of string width for genome name in PNG-format map
    my $image_tmp = GD::Image->new(100,100);
    my $font_path = "$home_directory/font_configure/6x10.bdf.fnt";
    my $courier = GD::Font->load($font_path) or die "Can't load font";
    $image_tmp->useFontConfig(1);
    my $black_tmp = $image_tmp->colorAllocate(0,0,0);

    foreach my $each_TFT_file(sort @dir_part_TFT){
 
        unless ($each_TFT_file =~ /.part$/){
            next;
        }
        my $filename = $each_TFT_file;

        $filename =~ s/.tbl(.*)//g;
        $filename =~ s/_/ /g;

#png string_width
        my @bounds = GD::Image->stringFT($black_tmp, "$font_family:$font_style", $strain_name_font_size, 0, 10, 10, "$filename");
        my $length_tmp =  $bounds[2]-$bounds[0];
        push (@text_name_width_png, $bounds[2]-$bounds[0]);
        $string_height_png{$filename} = $bounds[1]-$bounds[7];

#svg string_width
        my $text_name_width = metrics::string_width_svg($filename, "$font_family:$font_style", $strain_name_font_size); 
        push (@text_name_width_svg, $text_name_width);

        open (EACH_TFT_FILE, "$each_TFT_file") or die "Can' t open file '$each_TFT_file'";
        my @each_max_width;
        my @start_position;
        my @chech_dir;

       
        while (<EACH_TFT_FILE>){
            my @each_width= split "\t", $_;
            $each_width[0]=~ s/[><]//g;
            $each_width[1]=~ s/[><]//g;
            push (@each_max_width,@each_width[0,1]);

            # nif_start position
            $each_width[3] =~ s/.*;//g;  # # delete the geneName from "geneName;Locus_tag" for interested_gene
            if ( defined ($hash_nif{$each_width[3]}) ) {
                my $each_key_gene_start;
                if($each_width[0]<$each_width[1]) {
                    $each_key_gene_start = $each_width[0];
                } else {
                    $each_key_gene_start = $each_width[1];
                    } 
                push (@start_position, $each_key_gene_start);
                push (@interested_gene_length,abs($each_width[0]-$each_width[1]+1));
                push (@chech_dir, $each_width[0]);

             }
             # nif_start position

         }    

         @each_max_width = sort{$a<=>$b} @each_max_width;

         @start_position = sort{$a<=>$b} @start_position;
         my $start_position = scalar @start_position;   

         # obtain the start distance       
         my $median_value; 

         # check gene encoding direction
         my $check_number=0;
         foreach (@chech_dir) {
             if ($_ eq $start_position[int(($start_position-1)/2)]) {
                 $check_number++;
             }else { 

                }         
         }
         
         $direction_F_R{$each_TFT_file} = $check_number;
         if ($direction_F_R{$each_TFT_file} ge 1){                
             $median_value = $start_position[int(($start_position-1)/2)]-$each_max_width[0];
         }else{
             $median_value = $each_max_width[$#each_max_width]-$chech_dir[int(($start_position-1)/2)];           
             }

         # obtain the key-value of %start_right
         $start_right {$median_value} = $each_max_width[$#each_max_width]-$each_max_width[0];
         $shift_distance_temp{$each_TFT_file} = $median_value; 
  
    }
         @text_name_width_svg = sort{$b<=>$a} @text_name_width_svg;
         @text_name_width_png = sort{$b<=>$a} @text_name_width_png;

         @max_start_position = keys %start_right;
         @max_start_position = sort{$a<=>$b} @max_start_position;
         my $max_start = $max_start_position[$#max_start_position];

         foreach (@max_start_position) {          
             push (@max_width_right, $max_start-$_+$start_right{$_});
         }

         @max_width_right = sort{$a<=>$b} @max_width_right;

         foreach (keys %shift_distance_temp) {
             $shift_distance_x{$_} = $max_start - $shift_distance_temp{$_};
         }

         @interested_gene_length = sort{$a<=>$b} @interested_gene_length;
        my $interested_gene_lengthest = $interested_gene_length[$#interested_gene_length];

    return ($max_width_right[$#max_width_right], $interested_gene_lengthest, \%direction_F_R, \%shift_distance_x, $text_name_width_png[0], $text_name_width_svg[0], \%string_height_png);

}


sub geneName_blank_site_count {
    my $genome_file_name = shift;
    my $offset=0;
    my @blank_site_array;
    push (@blank_site_array, 0);
    my $result_split = index($genome_file_name, " ", $offset); #

    while ($result_split != -1) {
            push (@blank_site_array, $result_split);

            $offset = $result_split + 1;
            $result_split = index($genome_file_name, " ", $offset);
    }

            return @blank_site_array;


}


1;               
__END__  
