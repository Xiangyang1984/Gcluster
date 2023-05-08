#!usr/bin/perl -w

use strict;
use warnings;
use Bio::SeqIO;
use File::Basename qw<basename dirname>;
use threads;

#########################################################################################################################################
# Before using, interested_gene_generation.pl needs to pre-install several Perl Modules.
# Required Perl Modules: threads File::Basename Bio::SeqIO; (part of BioPerl, http://bioperl.org)

#USAGE
#===========
# perl extract_rRNA_dir.pl path_genbank_dir thread_number
# ***Please use the absolute path for script extract_rRNA_dir.pl***
# For example: perl /home/xiangyang/Gcluster_v1.01-master/script/extract_rRNA_dir.pl /home/xiangyang/Gcluster_v1.01-master/test_data/gbk 3
#########################################################################################################################################


my $path_genbank = $ARGV[0];
my $thread_number = $ARGV[1];
my @information;

my $rRNA_gene_dir = (dirname $path_genbank)."/rRNA_gene_dir"; 
mkdir $rRNA_gene_dir;
opendir (RRNA_DIR, $rRNA_gene_dir);
my @rRNA_gene_dir = readdir RRNA_DIR; 
@rRNA_gene_dir = grep ($_!~/^\./ ,@rRNA_gene_dir);
closedir RRNA_DIR;
system ("rm $rRNA_gene_dir/*") if scalar @rRNA_gene_dir > 0;





@information = batch_genbank_16S_rRNA_extract($path_genbank, $rRNA_gene_dir, $thread_number);

chdir $rRNA_gene_dir;
system ("cat *.fasta > all_16S_rRNA_gene.fasta");

my $log = "$rRNA_gene_dir/log.txt";

open (LOG_FILE, ">$log");

print LOG_FILE join ("", @information);
 
sub batch_genbank_16S_rRNA_extract {

    my ($path_genbank, $rRNA_gene_dir, $thread_number)=@_;
    my @information;

    system ("rename 's/ /_/g' $path_genbank/*");  #repalce the space with "_"  for all GenBank files
    system ("rename 's/_=_/_/g' $path_genbank/*");  #repalce "=" with "_"  for all GenBank files
    opendir PATH_GENBANK, $path_genbank or die "could not open $path_genbank";

    my @path_genbank = readdir PATH_GENBANK;
    @path_genbank =grep ($_!~/^\./ ,@path_genbank);  #delete hidden file . .. 
    closedir PATH_GENBANK;

    my $sub_file_number = scalar @path_genbank;

    my @input;
    my @outfile;
 
    foreach my $file_genbank(@path_genbank){ 

        my $input="$path_genbank/$file_genbank"; 
        push (@input,$input);  
        @input = sort @input;  

        my $output="$rRNA_gene_dir/$file_genbank.16S_rRNA.fasta";
        push (@outfile,$output);
        @outfile = sort @outfile;

    }

    my $thread;
    my @threads;
    my $job_number=0;

    while(){ 

        last if ($job_number eq $sub_file_number);  

        while(scalar(threads->list())<$thread_number) {     #set thread number
            $job_number++;    
            my $input_file = $input[$job_number-1];
            my $output_file = $outfile[$job_number-1];
            print "16S rRNA gene is extracted from Genbank $job_number\n";
            $threads[$job_number-1]=threads->new({context=>'array'}, \&genbank_16S_rRNA_sequence_extract, $input_file, $output_file); 

            last if ($job_number eq $sub_file_number);  
        }

        foreach $thread(threads->list(threads::all)){
            if($thread->is_joinable()) {              
              push (@information, join ("-", $thread->join()) );                     
            }
        }
    }

    foreach my $th (threads->list(threads::all)) {
        push (@information, join ("-", $th->join()) );      
    }
    return @information;
}


sub genbank_16S_rRNA_sequence_extract { 
    my ($input, $output)=@_;
    my @information;
    my $file_name = basename($input);
    open(OUTPUT, ">$output");

    my $each_input = new Bio::SeqIO(-file =>$input, -format => 'genbank');

    my %hash_rRNA;
    my $i=0;
    while ( my $seqobject = $each_input->next_seq() ) {
        my $acc = $seqobject->accession;
        my $count = 0; 
        my @features = $seqobject->get_SeqFeatures();
#================================================================ display organism name for 16S rRNA gene		
	my @org;
        foreach my $feat (@features) {

	    if ($feat->primary_tag() eq "source"){
		@org = $feat->get_tag_values('organism');
		next;
	    }	
#================================================================ Xiangyang li addition 2013-01-25
            unless ( $feat->primary_tag eq "rRNA") {
                next;
            }

            $count++;

            my (@fasta_ID, @check_tag);

            if ( $feat->has_tag('locus_tag') ) {
                push(@fasta_ID, $feat->get_tag_values('locus_tag'));
            }else {push(@fasta_ID, "CDS_$acc"."_".$count);}

            if ( $feat->has_tag('product') ) {                    
            	push (@fasta_ID, $feat->get_tag_values('product'));

            	@check_tag = $feat->get_tag_values('product');
            }
            $org[0] =~ s/ = /_/g;
            push (@fasta_ID, $file_name);                       
         
            my $fasta_ID = join(";", @fasta_ID );

            $fasta_ID =~ s/\s+/_/g;
            $fasta_ID =~ s/\n//g;
            $fasta_ID =~ s/\t+/ /g;
            #print " $check_tag[0]\n";
            if ($check_tag[0] eq "16S ribosomal RNA") {
                $i++;
                $hash_rRNA{$fasta_ID} = ($feat->spliced_seq->seq)[0];
            }
        }                                                                                
    }
           
    if ($i > 0) {
        my ($max_ID, $max_seq, $max_length) = maxlength_rRNA_genbank(\%hash_rRNA);                                             
  
        push (@information, "$file_name has $i 16S rRNA gene\(s\), and $max_ID has the max length \($max_length bp\)\n");
        $max_ID =~ s/.*;//g; 
        print OUTPUT ">$max_ID\n$max_seq\n"; 

    } else { 

        push (@information, "$file_name has no 16S rRNA gene\n");

    }
    close OUTPUT;
    return @information;
}                                                                                                 


sub maxlength_rRNA_genbank {

    my $hash_rRNA = shift;
    #my %hash_rRNA = %$hash_rRNA_1;
    my $max_length = 0;
    my $max_ID;
    my $max_seq;
    my %max_hash;

    foreach (keys %$hash_rRNA) {
        chomp;
        
        if ($max_length <= length ($hash_rRNA->{$_})){
            $max_length =length ($hash_rRNA->{$_});
            $max_ID = $_;
            $max_seq = $hash_rRNA->{$_};

        }

    }

    #return $max_ID, $max_seq, $max_length;

}
