#!usr/bin/perl -w

use strict;
use warnings;
use FindBin;

my $script_dir = $FindBin::Bin;
my $gbk_dir = "/media/xiangyang/X128/Gcluster-Gcluster-v2.0.6/stutzeri_123";
my $db = "/media/xiangyang/X128/Gcluster-Gcluster-v2.0.6/ectC_FAA.txt";
my $tree = "/media/xiangyang/X128/Gcluster-Gcluster-v2.0.6/127_core_concatenation_542_faa-NJ1000.nwk";

###Test-step3: Begin test interested_gene_generation.pl with the test_data
print "\n\nStep1: run interested_gene_generation.pl...\n";
print "################################################################\n";
my $check_i = system ("perl $script_dir/interested_gene_generation.pl -dir $gbk_dir -db $db -m 3");
print "################################################################\n";
if ($check_i eq 0){
    print "Ok, interested_gene_generation.pl works success!\n";
}else {
    print "Not Ok, interested_gene_generation.pl works with some errors!\n";
}


###Test-step2: Begin test Gcluster.pl with the test_data
print "\n\nStep2: run Gcluster.pl...\n";
print "################################################################\n";
my $check_g = system ("perl $script_dir/Gcluster.pl -dir $gbk_dir -gene $script_dir/interested_gene_workplace/interested_gene_name.txt -tree $tree -m 3 -n 20");
print "################################################################\n";
if ($check_g eq 0){
    print "Ok, Gcluster.pl works success!\n\n";
}else {
    print "Not Ok, Gcluster.pl works with some errors!\n";
}



