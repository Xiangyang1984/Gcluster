#!usr/bin/perl  

use strict;
use warnings;
use FindBin;

###Test-step1: Checks for Gcluster dependencies
print "\n\nTest-step1: Checks for Gcluster dependencies...\n";
print "################################################################\n";
#Checks for perl modules
&check_perl_modules;
print "\n----------------------------------------------------------------\n";
#Chechs for programs
&check_programs;
print "################################################################\n";


my $script_dir = $FindBin::Bin;

###Test-step2: Begin test Gcluster.pl
print "\n\n\n\n\nTest-step2: Begin test Gcluster.pl...\n";
print "################################################################\n";
system ("perl $script_dir/Gcluster.pl -dir $script_dir/test_data/gbk -gene $script_dir/test_data/interested_gene_name.txt -tree $script_dir/test_data/16S_rRNA_tree.nwk -m 9");
print "################################################################\n";
print "Ok, Gcluster.pl works success!\n";



###Test-step3: Begin test interested_gene_generation.pl
print "\n\n\n\n\nTest-step3: Begin test interested_gene_generation.pl...\n";
print "################################################################\n";
system ("perl $script_dir/interested_gene_generation.pl -dir $script_dir/test_data/gbk -db $script_dir/test_data/aioB.fasta -m 9");
print "################################################################\n";
print "Ok, interested_gene_generation.pl works success!\n\n";



#check Perl modules dependencies;
sub check_perl_modules {
 
    my @test_modules = ("GD", "GD::SVG", "SVG", "threads", "File::Basename", "FindBin", "lib", "Getopt::Long", "Math::BigFloat", "Storable", "vars", "File::Spec", "Bio::SeqIO", "Bio::Tree::NodeI", "Bio::TreeIO");

    my $check_number=0;
    foreach (@test_modules){
        my $cmd = "perl -M$_ -e 'print \"***$_ Version\\t \".$_->VERSION.\"\tok.\n\"'";
        #perl -MGD::SVG -e 'print GD::SVG->VERSION. "\n"'
        my $test_out = system ($cmd);

        print  "***Warning:\t$_ is not installed\n" if $test_out ne 0;
        $check_number++ if $test_out ne 0;
    }
    if ($check_number == 0){
        print "!!!Ok, all dependencies Perl modulers are installed*\n";
    }else{
        print "!!!Not ok, it appears that certain dependencies Perl modulers are installed*\n";
    }
}


# checks software dependencies
sub check_programs{

    print "Checking for makeblastdb ... ";
    my $makeblastdb_path = `which makeblastdb`;
    chomp $makeblastdb_path;
    if (not -e $makeblastdb_path){
    	print "error: makeblastdb is not installed\n";
    }
    else{
	    print "OK, makeblastdb is installed at: $makeblastdb_path\n";
    }
    
    print "Checking for blastp ... ";
    my $blastp_path = `which blastp`;
    chomp $blastp_path;
    if (not -e $blastp_path){
    	print "error: blastp is not installed\n";
    }
    else{
    	print "OK, blastp is installed at: $blastp_path\n";
    }
    
    print "Checking for mcl ... ";
    my $mcl_path = `which mcl`;
    chomp $mcl_path;
    if (not -e $mcl_path){
    	print "error: mcl is not installed\n";
    }
    else{
	    print "OK, mcl is installed at: $mcl_path\n";
    }

}
