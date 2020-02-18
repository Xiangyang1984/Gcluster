#!/usr/bin/perl -w

use strict;
use warnings;
use threads;
use Bio::SeqIO;
use Getopt::Long;
use File::Basename qw<basename dirname>;
use FindBin;
use File::Spec;

#########################################################################################################################################
# Before using, interested_gene_generation.pl needs to pre-install several Perl Modules and NCBI BLAST+ programs. In addition, the paths of makeblastdb and blastp need to be set.
# Required Perl Modules: threads File::Basename Getopt::Long FindBin File::Spec Bio::SeqIO; (part of BioPerl, http://bioperl.org)
# Required Programs: makeblastdb and blastp in NCBI BLAST+, which is available from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/.

#Please set the absolute path for two programs in "Gcluster/interested_gene_generation.pl", as in the following example:

my $blastp       = "/miniconda3/bin/blastp";
my $makeblastdb  = "/miniconda3/bin/makeblastdb";
#########################################################################################################################################

my $usage = <<USAGE; 

=NAME

interested_gene_generation.pl

=DESCRIPTION

    Run this command to enble users to obtain a list of the interested gene (a two-column tab-delimited text file) by a local blastP analysis using multiple threads. 

=USAGE

interested_gene_generation.pl -dir genbank_file_directory -db database [options]

FOR EXAMPLE: 
# perl /home/xiangyang/Gcluster_v1.01-master/interested_gene_generation.pl -dir /home/xiangyang/Gcluster_v1.01-master/test_data/gbk -db /home/xiangyang/Gcluster_v1.01-master/test_data/aioB.fasta -m 3

=ARGUMENTS
=======================
    REQUIRED ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -dir, --genbank_file_directory
           A directory containing annotated genomes as Genbank format file. To avoid a mistake, genome names cannot use special character,
           such as space, equal. For large number of genomes, users are recommended to download using Aspera, a high-speed file transfer
           tool (https://downloads.asperasoft.com/).                           
       -db, --database
           A protein database in FASTA format, which contains at least one protein sequence homologous to the gene of interest.
 
    OPTIONAL ARGUMENTS:
    ~~~~~~~~~~~~~~~~~~~
       -o, --output_directory
           An output directory holding all the generated files by interested_gene_generation.pl. if this option is not set, interested_gene_generation.pl will create a directory named "interested_gene_workplace" in the bin directory from where interested_gene_generation.pl was invoked.
       -m, --multiple_threads
           set thread number (Default: 1)
       -b, --start_at_blast 
           Jump to a local blastp analysis, and Skips sequencing extraction (Default: T).  
       -e, --e_value
           set E-value cutoff in Blast analysi (default: 1e-5)
       -i, --identify
           set percent identity cutoff in Blast analysis (default: 30)
       -c, --coverage
           set percent coverage (Query and Subject) cutoff in Blast analysis (default: 50)
       -l, --match_length
           set alignment length cutoff in Blast analysis (default: 30) 
       -h, --help
           Show this message.

=AUTHOR

Dr. Xiangyang Li (E-mail: lixiangyang\@fudan.edu.cn), Fudan university; Kaili University; Bacterial Genome Data mining & Bioinformatic Analysis (www.microbialgenomic.com/).

=COPYRIGHT

Copyright 2019, Xiangyang Li. All Rights Reserved.

USAGE


my %options = (
    'genbank_file_directory'                   => undef,   
    'database'                                 => undef, 
    'output_directory'                         => undef,  
    'multiple_threads'                         => "1",
    'start_at_blast'                           => "F",
    'e_value'                                  => "1e-5",    #E-value cutoff, default = 1e-5
    'identify'                                 => "30",      #Percent identity cutoff, default = 30
    'coverage'                                 => "50",      #Percent coverage cutoff, default = 50
    'match_length'                             => "30",      #Length cutoff, default = 30
    'help'                                     => undef

);

GetOptions(
    'dir|genbank_file_directory=s'             => \$options{genbank_file_directory},    
    'db|database=s'                            => \$options{database},  
    'o|output_directory=s'                     => \$options{output_directory},     
    'm|multiple_threads=s'                     => \$options{multiple_threads},     
    'b|start_at_blast=s'                       => \$options{start_at_blast},
    'e|e_value=s'                              => \$options{e_value},
    'i|identify=i'                             => \$options{identify},
    'c|coverage=i'                             => \$options{coverage},
    'l|match_length=i'                         => \$options{match_length},
    'h|help'                                   => \$options{help}

);

if (defined $options{help}) {
    print "$usage\n";
    exit(0);
}


my $now_time = localtime;
print "\n$now_time: interested_gene_generation.pl start...\n\n";

my $e_value = $options{e_value};
my $identify = $options{identify};
my $coverage = $options{coverage};
my $match_length = $options{match_length};
my $output_title = "Qseqid\tSseqid\tBitscore\tE-value\tPidenty\tQ_coverage\tS_coverage\tMacth_length";

my $path_genbank = File::Spec->rel2abs($options{genbank_file_directory});
my $db_file = File::Spec->rel2abs($options{database});
my $thread_number = $options{multiple_threads};

my $home_directory = $FindBin::Bin;           # obtaining the home directory where interested_gene_generation.pl located
print $home_directory,"\n";
#check for interested_gene_generation.pl workplace options
my $workplace;
if ( defined( $options{output_directory} ) ) {
    $workplace = File::Spec->rel2abs($options{output_directory});
    mkdir $workplace;
}else {

    $workplace = "$home_directory/interested_gene_workplace";
    $workplace =~ s/\/\//\//g;
    mkdir $workplace;
}

my $path_protein = $workplace."/Whole_proteome";
mkdir $path_protein;

my $blast_fold = "$workplace/blast_out";
system ("rm -r $blast_fold")   if (-e $blast_fold);
mkdir $blast_fold;


if ($options{start_at_blast} eq "F"){
     
     batch_genbank_sequence_TFT_extract($path_genbank, $path_protein, $thread_number);
}

my $new_db_file = $workplace."/".basename($db_file);
system ("cp $db_file $new_db_file");
format_sequence($makeblastdb, $new_db_file);

#chdir $path_protein;
my @input;
my @outfile;
opendir PATH_PROTEIN, $path_protein or die "could not open $path_protein";
my @PATH_PROTEIN = readdir PATH_PROTEIN; 
@PATH_PROTEIN =grep ($_!~/^\./ ,@PATH_PROTEIN);  #delete hidden file . ..
closedir PATH_PROTEIN;

foreach my $file_protein(@PATH_PROTEIN){
 
    my $input="$path_protein/$file_protein"; 
    push (@input,$input);  
    @input = sort @input;
    my $output="$blast_fold/$file_protein.blastout";
    push (@outfile,$output);
    @outfile = sort @outfile;

}

my $sub_file_number = scalar @input;
if ($sub_file_number <= $options{multiple_threads}) {
    $thread_number = $sub_file_number;
}

my @sub_function_parameters = (\@input, $new_db_file, \@outfile, $e_value, $blastp);
bacth_blast($sub_file_number, $thread_number, "do_blastP", \@sub_function_parameters);


#chdir $blast_fold;
system ("cat $blast_fold/* > $blast_fold/all_vs_all.blast");

#do local blastp analysis and group proteins into homologs family using mcl
my $all_vs_all_cluster = blast_filter ("$blast_fold/all_vs_all.blast", $e_value, $identify, $coverage, $match_length);

my %delete_repeat;
my %delete_repeat_2;
my $all_vs_all_cluster_tophit = "$blast_fold/all_vs_all_cluster_tophit";
my $all_vs_all_cluster_nontop = "$blast_fold/all_vs_all_cluster_nontop";
open (DELETE_REPEAT, $all_vs_all_cluster);
open (NEW_TOP_CLUSTER, ">$all_vs_all_cluster_tophit");
open (NEW_NONTOP_CLUSTER, ">$all_vs_all_cluster_nontop");

print NEW_TOP_CLUSTER "$output_title\n";
print NEW_NONTOP_CLUSTER "$output_title\n";

while (<DELETE_REPEAT>){
    chomp;
    if ($_ ne $output_title) {
        $_ =~ s/;/#/ if ($_ !~ /#/);
        my @delete_repeat = split "\t", $_;
        my @subject_id =  split ";", $delete_repeat[0]; 
        $delete_repeat{$subject_id[1]} = $_;       #key--strain Vs value--gene
        $delete_repeat_2{$_} = $subject_id[1];     #key--gene Vs value--strain

    }
}
 
foreach my $strain (keys %delete_repeat) {
    my @strain_number;
    foreach my $eachgene (keys %delete_repeat_2) {

        if ( $delete_repeat_2{$eachgene} eq $strain ) {
   
            push (@strain_number, $eachgene);

        }

    }

    if (scalar @strain_number >= 2) {
        my %score;
        foreach (@strain_number){
    
            my @delete_repeat = split "\t", $_;
            $score{$_} = $delete_repeat[2];

        }
        my $maxscore = 0;
        my $tophit_gene;

################################################################## Record top hit query
        foreach my $key(keys %score) {

            if ($maxscore < $score{$key}){

                $maxscore = $score{$key}; 
                $tophit_gene = $key;
            }
        }
        print NEW_TOP_CLUSTER $tophit_gene,"\n";

################################################################## Record none top hit query
        foreach my $key_nontophit(keys %score) {

            if ($key_nontophit ne $tophit_gene){

                print NEW_NONTOP_CLUSTER $key_nontophit,"\n";        
            }
        }


    }else {print NEW_TOP_CLUSTER $delete_repeat{$strain},"\n"; } # Record top hit query
            
}

close NEW_NONTOP_CLUSTER;
close NEW_TOP_CLUSTER;

my %nontop_gene;
open (GENE_NONTOP_TXT, $all_vs_all_cluster_nontop);

while (<GENE_NONTOP_TXT>){
    chomp;
    if ($_ ne $output_title) {
        my @gene_nontop_txt = split ("\t", $_);
        my @nontophit_array = split ";", $gene_nontop_txt[0];

        $nontop_gene{$nontophit_array[1]} .= "\t".$gene_nontop_txt[0];
    }   
}

open (GENE_LIST, $all_vs_all_cluster_tophit);
open (GENE_TXT, ">$workplace/interested_gene_name.txt");
while (<GENE_LIST>){
    chomp;
    if ($_ ne $output_title) {
        my @gene_top = split ("\t", $_);

        if ($gene_top[0] =~ /#/) {    
            $gene_top[0] =~ s/#/\t#/;
        }
        else {
            $gene_top[0] =~ s/;/\t#/;
        }

        my @tophit_array = split ";", $gene_top[0];

        if (defined $nontop_gene{$tophit_array[1]}) {

            print GENE_TXT $gene_top[0], "\t|none_tophit_gene:", $nontop_gene{$tophit_array[1]}, "\n";
        }
        else{

            print GENE_TXT "$gene_top[0]\n";
        }
    }
}

  $now_time = localtime;
  print "$now_time: Finished!\n";

sub batch_genbank_sequence_TFT_extract {

    my ($path_genbank, $path_protein, $thread_number)=@_;

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
    #@path_genbank =grep ($_!~/^\~/ ,@path_genbank);  #delete temp file 
    closedir PATH_GENBANK;

    my $sub_file_number = scalar @path_genbank;

    my @input;
    my @outfile_1;
 
    foreach my $file_genbank(@path_genbank){ 

        my $input="$path_genbank/$file_genbank"; 
        push (@input,$input);  
        @input = sort @input;  

        my $output_1="$path_protein/$file_genbank.fasta";
        push (@outfile_1,$output_1);
        @outfile_1 = sort @outfile_1;

    }

    my $thread;
    my @threads;
    my $job_number=0;

    while(){ 

        last if ($job_number eq $sub_file_number);  

        while(scalar(threads->list())<$thread_number) {     #set thread number
            $job_number++;    
            my $input_file = $input[$job_number-1];
            my $output_file_1 = $outfile_1[$job_number-1];
            print "Genbank_extraction: $job_number\n";
            $threads[$job_number-1]=threads->new(\&genbank_protein_sequence_extract, $input_file, $output_file_1); 

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


sub genbank_protein_sequence_extract {
####### Part of script in this subroutine comes from CGView Comparison Tool. Reference: Grant JR, Arantes AS, Stothard P (2012) Comparing thousands of circular genomes using the CGView Comparison Tool. BMC Genomics 13:202.

    my ($input, $output_1)=@_;
    my $file_name = basename($input);
    open(OUTPUT_1, ">$output_1");

        my $in = new Bio::SeqIO( -file => $input, -format => 'genbank' );

        while ( my $seqObject = $in->next_seq() ) {
            my $acc = $seqObject->accession;
            my @features = $seqObject->get_SeqFeatures();
            my $count = 0;       
            
            foreach (@features) {
                my $feat = $_;	        
                unless ($feat->primary_tag =~ /^CDS/) {
                    next;
                }
                    if ( $feat->has_tag('pseudo') ) {

                    }else {$count++;}
                    #print "genome_$count\n";

                    my $locus_tag;
                    if ( $feat->has_tag('locus_tag') ) {
                        ($locus_tag) = $feat->get_tag_values('locus_tag');

                    }else {$locus_tag = "CDS_$acc"."_".$count;}
               
                    my $product;
                    if ( $feat->has_tag('product') ) {
                        ($product) = $feat->get_tag_values('product');
                        $product =~ s/ /_/g;
                    }
 
                    my $protein_seqeunce;
                    if ( $feat->has_tag('translation')) {
                        ($protein_seqeunce) = $feat->get_tag_values('translation');
                      
                            print OUTPUT_1 ">$locus_tag#$product;$file_name\n$protein_seqeunce\n";
                        
                    }
            }
        }

    close OUTPUT_1;

}



sub format_sequence {
    my $makeblastdb = shift;
    my $db_file = shift;
    system ("$makeblastdb -in $db_file -dbtype prot -logfile $db_file.makeblast.log");
}



sub do_blastP {
    my ($input, $db_file, $outfile, $e_value, $blastp) =@_;
    system("$blastp -outfmt '6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore' -evalue $e_value -num_alignments 1 -num_threads 1 -query $input -db $db_file -out $outfile");

}



sub bacth_blast {

    my ($work_number, $thread_number, $subfunction, $sub_function_parameters) = @_;
    my ($input_file, $db_file, $output_file, $e_value, $blastp) = @$sub_function_parameters;
 
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
            $threads[$job_number-1]=threads->new(\&$subfunction, $input_file, $db_file, $output_file, $e_value, $blastp);   
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



sub blast_filter {

    my ($infile, $e_value, $identify, $coverage, $match_length) = @_;
    open (INFILE, $infile);
    my $all_vs_all_cluster = "$blast_fold/all_vs_all.cluster";

    open (OUTFILE, ">$all_vs_all_cluster");

    print OUTFILE "$output_title\n";

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
