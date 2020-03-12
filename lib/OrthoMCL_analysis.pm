package OrthoMCL_analysis;

# This module mainly comes from an old version of ORTHOMCL (Verson 1.4, https://orthomcl.org/common/downloads/software/unsupported/v1.4/) to avoid using mysql, Reference: Li Li, Christian J. Stoeckert, Jr. and David S. Roos. OrthoMCL: Identification of Ortholog Groups for Eukaryotic Genomes. Genome Research.13:2178-2189, 2003. We modified this script to enable Gcluster working with multible threads.

use strict;
use warnings;
use threads;
use Storable;

use vars qw($maximum_weight $e_value $identify $coverage $genome_gene_file $blast_file $bpo_file $bpo_idx_file $bpo_se_file $matrix_file $index_file $mcl_file $mcl_bi_file $orthomcl_log_file $bbh_file @taxa %gindex %gindex2 %blastquery %graph %weight $match_length_single $address_ref $blastquery_ref @mcl_index); # global variables storing data


sub orthoMCL {

my ($directory_sub_blastout, $directory_protein_part, $directory_homologs_cluster, $interested_gene_name, $homologous_gene_cutoff, $thread_number, $blastp, $makeblastdb, $MCL, $inflation) =@_;

my (%connect, %ortho);

($e_value, $identify, $coverage, $match_length_single) = split (",", $homologous_gene_cutoff);

$genome_gene_file        = "$directory_sub_blastout/genome_gene_file.gg.txt";
$blast_file              = "$directory_sub_blastout/all.blast";
$bpo_file                = "$directory_sub_blastout/all.bpo";
$bpo_idx_file            = "$directory_sub_blastout/all_bpo.idx";
$bpo_se_file             = "$directory_sub_blastout/all_bpo.se";
$matrix_file             = "$directory_sub_blastout/all_ortho.mtx";
$index_file              = "$directory_sub_blastout/all_ortho.idx";
$mcl_file                = "$directory_sub_blastout/all_ortho.mcl";

$bbh_file                = "$directory_homologs_cluster/all_blast.bbh";
open (BBH,">$bbh_file") or die "can't create $bbh_file";

$mcl_bi_file             = "$directory_homologs_cluster/all_orthomcl.out";

$orthomcl_log_file       = "$directory_homologs_cluster/orthomcl.log";
open (LOG,">$orthomcl_log_file") or die "can't create $orthomcl_log_file";


blast_homologs_cluster($directory_sub_blastout, $directory_protein_part, $interested_gene_name, $homologous_gene_cutoff, $thread_number, $blastp, $makeblastdb);

read_ggfile($genome_gene_file);


if ($bpo_file =~ m/(\S+)\.(\S+)/) {
    $bpo_idx_file     = $1.'_bpo.idx';
    $bpo_se_file      = $1.'_bpo.se';
} else {
    $bpo_idx_file     = $bpo_file.'_bpo.idx';
    $bpo_se_file      = $bpo_file.'_bpo.se';
}

$maximum_weight = maximum_weight($bpo_file); #caculate the minimum p-value (excluding evulue = 0) using blast software. For example, if 1e-180 is the case you should use -log(1e-181)=181.

&constructIDX_for_bpofile($bpo_file,$bpo_idx_file) unless (-e $bpo_idx_file);

&constructSE_for_bpofile($bpo_file,$bpo_se_file) unless (-e $bpo_se_file);


&open_bpofile($bpo_file);
&retrieve_from_file($bpo_idx_file,$bpo_se_file);


foreach my $taxon (@taxa) {
	write_log("\nIdentifying inparalogs from $taxon\n");
	@{$connect{$taxon.' '.$taxon}}  = &makeInparalog($taxon);                      # identification of inparalogs
}

for(my $i=0;$i<scalar(@taxa)-1;$i++) {
	for(my $j=$i+1;$j<scalar(@taxa);$j++) {
		write_log("\nIdentifying ortholog pairs between $taxa[$i] and $taxa[$j]\n");  

		@{$connect{$taxa[$i].' '.$taxa[$j]}} = &makeOrtholog($taxa[$i],$taxa[$j]); # identification of orthologs
		write_log("Appending co-ortholog pairs between $taxa[$i] and $taxa[$j]: ");  
		my $c_coortholog=0;

		my %e;
		my $edge_ref=$connect{$taxa[$i].' '.$taxa[$j]}->[0];
		foreach my $pi (keys %$edge_ref) {@{$e{$pi}}=@{$edge_ref->{$pi}};}  #make a copy of current edge data structure into %e

		my %w =  %{$connect{$taxa[$i].' '.$taxa[$j]}->[1]};
		my $sumw =  $connect{$taxa[$i].' '.$taxa[$j]}->[3];
		my $c_ortholog = $connect{$taxa[$i].' '.$taxa[$j]}->[4];
		my %p1 = %{$connect{$taxa[$i].' '.$taxa[$i]}->[0]};
		my %p2 = %{$connect{$taxa[$j].' '.$taxa[$j]}->[0]};
		my %para;
		foreach my $p (keys %{$connect{$taxa[$i].' '.$taxa[$i]}->[0]}) {
			push (@{$para{$p}}, @{$p1{$p}}); }
		foreach my $p (keys %{$connect{$taxa[$j].' '.$taxa[$j]}->[0]}) {
			push (@{$para{$p}}, @{$p2{$p}}); }

		
		foreach my $n (keys %e) {
			$ortho{$n} = 1;
			my (@nodes1, @nodes2);

			if (exists($para{$n})) {push (@nodes1, $n, @{$para{$n}});}
			else {push (@nodes1, $n);}

			foreach (@{$e{$n}}) {
				if (exists($para{$_})) {push (@nodes2, $_, @{$para{$_}});}
				else {push (@nodes2, $_);}
			}

			@nodes1=@{nonredundant_list(\@nodes1)}; #can be commented
			@nodes2=@{nonredundant_list(\@nodes2)};

			for(my $k=0;$k<scalar(@nodes1);$k++) {
				for(my $l=0;$l<scalar(@nodes2);$l++) {
					
					next if(exists($w{$nodes1[$k].' '.$nodes2[$l]}));
					my ($pv1, $pv2);

					if (blastqueryab($nodes1[$k],$nodes2[$l])) {
						my ($s,$pm,$pe,$pi)=(blastqueryab($nodes1[$k],$nodes2[$l]))[0,3,4,5];
						next if($pm.'e'.$pe > $e_value || $pi< $identify);
						if($coverage) {
							next if(&simspan($s) < $coverage);
						}
						if($pm==0) { $pv1 = $maximum_weight;} else { $pv1 = -log($pm.'e'.$pe)/log(10); }
					} else {next;}

					if (blastqueryab($nodes2[$l],$nodes1[$k])) {
						my ($s,$pm,$pe,$pi)=(blastqueryab($nodes2[$l],$nodes1[$k]))[0,3,4,5];
						next if($pm.'e'.$pe > $e_value || $pi < $identify);
						if($coverage) {
							next if(&simspan($s) < $coverage);
						}
						if($pm==0) { $pv2 = $maximum_weight;} else { $pv2 = -log($pm.'e'.$pe)/log(10); }
						push (@{$edge_ref->{$nodes1[$k]}}, $nodes2[$l]);
						push (@{$edge_ref->{$nodes2[$l]}}, $nodes1[$k]);
						my $wt = ($pv1+$pv2)/2;
						# use averaged score as edge weight
						$w{$nodes1[$k].' '.$nodes2[$l]} = sprintf("%.3f", $wt);
						$w{$nodes2[$l].' '.$nodes1[$k]} = sprintf("%.3f", $wt);
						$sumw += $wt;
						$c_coortholog++;
					}
				}
			}
		}
		write_log("$c_coortholog pairs\n");
		my $avgw = 'N/A';
		if ($c_ortholog+$c_coortholog) {
			$avgw = $sumw/($c_ortholog+$c_coortholog);
		}
		write_log("$taxa[$i] and $taxa[$j] average weight: $avgw\n");
		foreach my $p (keys %w) {
			$w{$p} = sprintf("%.3f", $w{$p}/$avgw);
		}
		$connect{$taxa[$i].' '.$taxa[$j]}->[1] = \%w;
	}
}

%blastquery=();
%gindex=();

foreach my $taxon (@taxa) {
	write_log("\ncalculate average weight from $taxon\n");
	my %e = %{$connect{$taxon.' '.$taxon}->[0]};
	my %w = %{$connect{$taxon.' '.$taxon}->[1]};

	my $count=0; my $sum=0;
	my $count_all=0; my $sum_all = 0;
	foreach my $pair (keys %w) {
		my ($n,$p) = split(' ',$pair);
		$count_all++; $sum_all += $w{$n.' '.$p};
		if ($ortho{$n} || $ortho{$p}) {
			$count++;
			$sum += $w{$n.' '.$p};
		}
	}
	my $avgw;
	# normalize the in-paralog weights by the average weight of inparalogs which have orthologs in other species
	# common case, for eukaryotes and most prokaryotes
	if ($count) {
		$avgw = $sum/$count;
	}
	# OR normalize the in-paralog weights by the average weight of all inparalogs
	# not common, useful for prokaryotes or pathogens
	elsif ($count_all) {
		$avgw = $sum_all/$count_all;
		write_log("taxon average weight is calculated based on all inparalog pairs\n");
	}
	# OR no normalization since $count_all=0 and there is nothing stored in %weight
	# not common, useful for prokaryotes or pathogens 
	else {
		$avgw = 'N/A';
		write_log("taxon average weight is not calculated because there's no inparalog pairs\n");
	}
	write_log("$taxon average weight: $avgw\n");
	foreach my $p (keys %w) {
		$w{$p} = sprintf("%.3f", $w{$p}/$avgw);
	}
	$connect{$taxon.' '.$taxon}->[1] = \%w; 
}
%ortho=();


foreach my $p (keys %connect) {
	my %e = %{$connect{$p}->[0]};
	my %w =  %{$connect{$p}->[1]};
	
	foreach my $n (keys %e) {
		push(@{$graph{$n}}, @{$e{$n}});
		delete $e{$n};
	}
	%e=();
	foreach my $n (keys %w) {
		$weight{$n} = $w{$n};
		delete $w{$n};
	}
	%w=();
	delete $connect{$p};
}
%connect=();

write_matrix_index($matrix_file,$index_file);

%graph=();
%weight=();

executeMCL($MCL, $matrix_file,$mcl_file,$inflation);
mcl_backindex($mcl_file,$mcl_bi_file);
%gindex2=();
return $mcl_bi_file;
} # orthoMCL subrouting
#######################################SUBROUTINES###########################################
# This subroutine is an important part of OrthoMCL, used to
# look for inparalog (recent paralog) which is defined as 
# reciprocal better hits here.
# Please refer to the OrthoMCL paper for more details.
# One Arguments:
# 1. String Variable: Taxon name
# Last modified: 10/02/06
sub makeInparalog {
	my $taxon = $_[0];
	my (%seqs, %inbest, %pvalue,%sim);
	foreach (@{$gindex{$taxon}}) {$seqs{$_} = 1;}
	foreach my $qid (keys %seqs) {
		my ($sStart,$sEnd);
		if (defined $blastquery{$qid}) {
			($sStart,$sEnd)=split ("\t",$blastquery{$qid});
		} else {next;}
		my @sorted_simid=pvtie_sort($sStart,$sEnd,$taxon);
		LINE:foreach (0..$#sorted_simid) {
			my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($sorted_simid[$_]))[0,3,5,6,7];
			if ($sid ne $qid) {
				last LINE unless ($seqs{$sid});                                  ## better hit not from the same species
				last LINE if($pm.'e'.$pe > $e_value || $pi < $identify);      ## better hit not meet the cutoff
				if($coverage) {
					next LINE if(&simspan($s) < $coverage);
				}
				push(@{$inbest{$qid}}, $sid);
				$pvalue{$qid.' '.$sid} = $pm.'e'.$pe; 
			}
		}
	}
	my @b = keys %inbest;
	write_log(scalar(@b)." sequences have better hits within species\n");
	return &matrix(\%inbest, \%pvalue);
} ##makeInparalog



# This subroutine is an important part of OrthoMCL, used to
# look for ortholog which is defined as the reciprocal best
# hit between two species.
# Please refer to the OrthoMCL paper for more details.
# Two Arguments:
# 1. String Variable: Taxon name
# 2. String Variable: Taxon name
# Last modified: 10/02/06
sub makeOrtholog {
	my ($ta,$tb) = @_;
	my (@seqs,%best,%sim,%pvalue);
	foreach my $qid (@{$gindex{$ta}}) {
		my ($sStart,$sEnd);
		if (defined $blastquery{$qid}) {
			($sStart,$sEnd)=split ("\t",$blastquery{$qid});
		} else {next;}
		my ($lastpm,$lastpe);
		my $hit_id=0;
		LINE:foreach ($sStart..$sEnd) {
			my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($_))[0,3,5,6,7];
			if (defined $gindex2{$sid}) {
				if ($gindex2{$sid} eq $tb) {
					$hit_id++;
					if ($hit_id==1) {
						push(@{$sim{$qid}},"$sid,$pm,$pe,$pi,$s");
						$lastpm=$pm;$lastpe=$pe;
					}
					else {
						if (($lastpm==$pm) && ($lastpe==$pe)) {
							push(@{$sim{$qid}},"$sid,$pm,$pe,$pi,$s");
						}
						else {last LINE;}
					}
				}
			} else {write_log("$sid gindex2 not defined; lineid: $_\n");}
		}
	}
	foreach my $qid (@{$gindex{$tb}}) {
		my ($sStart,$sEnd);
		if (defined $blastquery{$qid}) {
			($sStart,$sEnd)=split ("\t",$blastquery{$qid});
		} else {next;}
		my ($lastpm,$lastpe);
		my $hit_id=0;
		LINE:foreach ($sStart..$sEnd) {
			my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($_))[0,3,5,6,7];
			if (defined $gindex2{$sid}) {
				if ($gindex2{$sid} eq $ta) {
					$hit_id++;
					if ($hit_id==1) {
						push(@{$sim{$qid}},"$sid,$pm,$pe,$pi,$s");
						$lastpm=$pm;$lastpe=$pe;
					}
					else {
						if (($lastpm==$pm) && ($lastpe==$pe)) {
							push(@{$sim{$qid}},"$sid,$pm,$pe,$pi,$s");
						}
						else {last LINE;}
					}
				}
			} else {write_log("$sid gindex2 not defined; lineid: $_\n");}
		}
	}

	foreach my $q (keys %sim) {
		foreach (@{$sim{$q}}) {
			my @bla=split (',',$_);
			next if($bla[1].'e'.$bla[2] > $e_value || $bla[3]< $identify);
			if($coverage) {
				next if(&simspan($bla[4]) < $coverage);
			}
			push(@{$best{$q}}, $bla[0]);
			$pvalue{$q.' '.$bla[0]} = $bla[1].'e'.$bla[2];

		}
	}
	my @b = keys %best;
	write_log(scalar(@b)." sequences have best hits from the other species\n");
	return &matrix(\%best, \%pvalue);
} ## makeOrtholog




# This subroutine is used to choose two-way hits among one-way hits (best
# hits between two species or better hits within one species), 
# calculate the weight between two nodes (minus logrithm of the p-value, 
# or $MAX_WEIGHT_DEFAULT for p-value 0 ), and calculate average
# weight among all inparalogs within one species or all orthologs between
# two species. (Weighting process takes place in the main script)
# Two Arguments:
# 1. Reference Variable: reference to a hash which stores all the possible
#    gene pairs (one-way best hit, or better hit).
# 2. Reference Variable: reference to a hash which stores the pvalue for
#    the gene pairs.
# Last modified: 10/02/06
sub matrix {
	my %best      = %{$_[0]};
	my %pvalue    = %{$_[1]};
	my (%edge, %weight);
	my $count=0;
	my $sumw=0;

	foreach my $query (sort keys %best) {
		foreach my $subject (@{$best{$query}}) {
			next if($weight{$query.' '.$subject});
			my $flag = 0;
			foreach my $q (@{$best{$subject}}) {
				if($q eq $query) { $flag = 1; }
			}
			if($flag == 1) {
				push (@{$edge{$query}}, $subject);
				push (@{$edge{$subject}}, $query);
				#use -logP as weights and treat P=0 as -logP=$maximum_weight (DEFAULT=300)
				my ($pv1, $pv2);
				if($pvalue{$query.' '.$subject} == 0) {
					$pv1 = $maximum_weight;
				}else { 
					$pv1 = -log($pvalue{$query.' '.$subject})/log(10);
				}	    
				if($pvalue{$subject.' '.$query} == 0) {
					$pv2 = $maximum_weight;
				}else {
					$pv2 = -log($pvalue{$subject.' '.$query})/log(10);
				}
				write_bbh("$query	$subject	".$pvalue{$query.' '.$subject}."	".$pvalue{$subject.' '.$query}."\n");
				my $w = ($pv1+$pv2)/2;
				$sumw += $w;
				$count++;
				# use averaged score as edge weight
				$weight{$query.' '.$subject} = sprintf("%.3f", $w);
				$weight{$subject.' '.$query} = sprintf("%.3f", $w);
			}
		}
	}
	my $avgw = 'N/A';
	if ($count) {
		$avgw = $sumw/$count;
	}
	my $no_tmp = scalar(keys %weight)/2;
	write_log("$no_tmp sequence pairs were identified as Reciprocal Better/Best Hit\n");
	return (\%edge, \%weight, $avgw, $sumw, $count);
} ## matrix


# This subroutine is used by the subroutine makeInparalog,
# to solve the pv_tie problem. It rearranges the pv-tied blast hits
# so that the hits from a specific taxon are moved higher than hits from
# other species.
# Three Arguments:
# 1. Number Variable: starting line id (or similarity id) of bpo file (blast
#    parse out file)
# 2. Number Variable: ending line id (or similarity id) of bpo file (blast
#    parse out file)
# 3. String Variable: taxon
# Last modified: 07/20/04
sub pvtie_sort {
	my ($sStart,$sEnd,$taxon)=@_;
	my (@sorted_simid,@tmp);
	my ($lastpm,$lastpe)=(getline_from_bpofile($sStart))[5,6];
	foreach ($sStart..$sEnd) {
		my ($s,$sid,$pm,$pe,$pi)=(&getline_from_bpofile($_))[0,3,5,6,7];
		if (($lastpm==$pm) && ($lastpe==$pe)) {
			if ($gindex2{$sid} eq $taxon) {
				push (@sorted_simid,$s);
			}
			else {
				push (@tmp,$s);
			}
		}
		else {
			if (scalar(@tmp)>0) {
				push (@sorted_simid,@tmp);
				@tmp=();
			} 
			if ($gindex2{$sid} eq $taxon) {
				push (@sorted_simid,$s);
			}
			else {
				push (@tmp,$s);
			}
		}
		$lastpm=$pm;$lastpe=$pe;
	}
	if (scalar(@tmp)>0) {push (@sorted_simid,@tmp);}
	return @sorted_simid;
} ## pvtie_sort



# This subroutine, together with matchlen, are used to calculate
# how much of the query sequences match each other.
# One Argument:
# 1. Number Variable: line id (or similarity id) of bpo file (blast
# parse out file)
# Last modified: 07/21/04
sub simspan {
	my $s = $_[0];
	my (%sub_start, %sub_length, %query_start, %query_length);
	my @hsp=split ('\.',(&getline_from_bpofile($s))[8]);

	foreach (@hsp) {

		if (/(\d+)\:(\d+)\-(\d+)\:(\d+)\-(\d+)/) {
			$sub_start{$1}=$4; 
			$sub_length{$1}=$5-$4+1;
			$query_start{$1}=$2;
			$query_length{$1}=$3-$2+1;
		}
	}
	my $match_lengths = &matchlen(\%sub_start,\%sub_length);
	my $match_lengthq = &matchlen(\%query_start,\%query_length);			
	my ($lengthq,$lengths)=(&getline_from_bpofile($s))[2,4];   # June 3
	if($lengths >= $lengthq) {
		return 100*$match_lengthq/$lengthq;
	}else{
		return 100*$match_lengths/$lengths;
	}
} ##simspan



# This subroutine, together with simspan, are used to calculate
# how much of the query sequences match each other.
# Two Arguments:
# 1. Reference Variable: reference to an hash which stores the starting
#    position of each HSP.
# 2. Reference Variable: reference to an hash which stores the length
#    of each HSP.
# Last modified: 07/19/04
sub matchlen {

	my %start        = %{$_[0]}; 
	my %length       = %{$_[1]};

	my @starts = sort{$start{$a}<=>$start{$b}} (keys %start);
	return $length{$starts[0]} if(scalar(@starts)==1);
	my $i=1; 
	my  $match_length = $length{$starts[0]}; 
	my $pos = $length{$starts[0]} + $start{$starts[0]} ;
	while($i<scalar(@starts)) {

	if($length{$starts[$i]} + $start{$starts[$i]} <= $pos) {
		$i++;
		next;
	}
	if($start{$starts[$i]}> $pos) {
		$match_length += $length{$starts[$i]};
		$pos = $start{$starts[$i]} + $length{$starts[$i]};
	}else {
		$match_length += $length{$starts[$i]} - ($pos - $start{$starts[$i]});
		$pos = $start{$starts[$i]} + $length{$starts[$i]};
	}
	$i++;
	}

	return $match_length;
} ## matchlen



# This subroutine is used to facilitate reading blast information
# from blastparseout file, instead of reading all blast result into
# memory.
# It works by storing every line (of blastparseout file)'s address
# into an array, and then uses Storable module to save this array
# into a file for future reference.
# Two Arguments:
# 1. String Variable: blast parse out file name
# 2. String Variable: index file name, storing index information of
#    blast parse out file.
# Last modified: 07/19/04
sub constructIDX_for_bpofile {
	my ($file,$idxfile)=@_;
	my @address;
	open (FILE,$file) or die "can't open $file file";
	push (@address, tell (FILE));
	while (<FILE>) {push (@address, tell (FILE));}
	store \@address, $idxfile;
	close (FILE);
} ## constructIDX_for_bpofile


# Two arguments:
# 1. String Variable: genome-gene file
# Last modified: 07/20/04
sub read_ggfile {
	my $gg_file=$_[0];
	my $totalgeneno=0;
	open (GG,$gg_file) or die "can't open $gg_file!";
	while (<GG>) {
		$_=~s/\r|\n//g;
		if (/(\S+)\(\d+\):(.*)/) {
			my $taxon=$1;
			my @genes=split (" ",$2);
			push (@taxa,$taxon);
			push (@{$gindex{$taxon}},@genes);
			foreach (@genes) {$gindex2{$_}=$taxon;}
			$totalgeneno+=scalar(@{$gindex{$taxon}});
		}
		elsif (/(\S+):(.*)/) {
			my $taxon=$1;
			my @genes=split (" ",$2);
			push (@taxa,$taxon);
			push (@{$gindex{$taxon}},@genes);
			foreach (@genes) {$gindex2{$_}=$taxon;}
			$totalgeneno+=scalar(@{$gindex{$taxon}});
		}
	}
	write_log("\nThere are ".@taxa." genomes, $totalgeneno sequences in $gg_file !\n\n");
	close (GG);
} ## read_ggfile


# This subroutine is also used to facilitate reading blast information
# from blastparseout file, instead of reading all blast result into
# memory.
# It works by storing Starting and Ending line_ids of specific query_gene_id
# into a hash, and then uses Storable module to save this hash into a file
# for future reference.
# Two Arguments:
# 1. String Variable: blast parse out file name
# 2. String Variable: file name storing the hash
# Last modified: 07/19/04
sub constructSE_for_bpofile {
	my ($bpo_file,$bpo_se_file)=@_;
	my $lastquery='';
        my $lastsimid=0;
	open (BPOUT,$bpo_file) or die "can't open $bpo_file";
	while (<BPOUT>) {
		$_=~s/\r|\n//g;
		next unless ($_ ne '');
		my @bpout=split ("\t",$_);
		if ($lastquery ne $bpout[1]) {
			$blastquery{$bpout[1]}=$bpout[0];
			$blastquery{$lastquery}.="\t".$lastsimid if $lastquery;
		}
		$lastquery=$bpout[1];
		$lastsimid=$bpout[0];
	}
	$blastquery{$lastquery}.="\t".$lastsimid;
	store \%blastquery,$bpo_se_file;
	close(BPOUT);
} ## constructSE_for_bpofile


# Three Arguments:
# 1. String Variable: matrix file name
# 2. String Variable: mcl file name
# 3. Number Variable: inflation parameter for MCL algorithm
# Last modified: 07/19/04
sub executeMCL {

    my ($MCL, $matrix_file, $mcl_file, $inflation) = @_;

    system("$MCL $matrix_file -I $inflation -o $mcl_file -V all");
    if (-e $mcl_file) {
        write_log("\nMCL result $mcl_file generated!\n");
    }else {

        die "$mcl_file failed to be generated!";
    }
} ## executeMCL



# This module is used to open BPO (Blast Parse Out) file, 
# to provide an filehandle for later use in blast query process.
# One Argument:
# 1. String Variable: BPO file name
# Last modified: 07/21/04
sub open_bpofile {
	my $bpo_file = $_[0];
	open (BPOFILE,$bpo_file) or die "can't open $bpo_file file"; # getline_from_bpofile subroutine will use this file handle 
																			
} # open_bpofile


# This subroutine is used to retrieve @address and %blastquery
# from blast_idx_file and blast_SE_file
# Two arguments:
# 1. String Variable: blast_idx_file
# 2. String Variable: blast_SE_file
# Last modified: 07/19/04
sub retrieve_from_file {
	my $blastidxfile  = $_[0];
	my $blastsefile   = $_[1];
	$address_ref=retrieve ($blastidxfile);
	$blastquery_ref=retrieve ($blastsefile);
	%blastquery=%{$blastquery_ref};
} ## retrieve_from_file


# This subroutine is used to retrieve blast information from bpo file
# (blast parse out file), given the line_id.
# One Argument:
# 1. Reference Variable: reference to the address array
# 2. Number Variable: line_id of bpo file
# Last modified: 09/08/04
sub getline_from_bpofile {
	my $lineid         = $_[0];
	seek (BPOFILE,$address_ref->[$lineid-1],0);
	my $line=<BPOFILE>; 
	$line=~s/\r|\n//g;
	chop $line;

	my @bpout=split ("\t",$line);

	my ($pm,$pe);
	if ($bpout[5] eq 0) {
             $pm=0;
             $pe=0;
        }
#	elsif ($bpout[5]=~/(\d+)e(\-\d+)/) {$pm=$1;$pe=$2;}
	elsif ($bpout[5]=~/(\S+)e(\-\S+)/) {$pm=$1;$pe=$2;}  #For WU-BLAST their p-value has the pattern \d\.\de\-\d+  OR 0.
	else {$pm=$bpout[5];$pe=0;}
	return ($bpout[0],$bpout[1],$bpout[2],$bpout[3],$bpout[4],$pm,$pe,$bpout[6],$bpout[7]);
} ## getline_from_bpofile


# This subroutine is used to make a nonredundant list.
# One Argument:
# 1. Reference Variable: reference to an array
# Last modified: 07/19/04
sub nonredundant_list {
	my $list_ref=$_[0];
	my %nr=();
	foreach (@{$list_ref}) {$nr{$_}=1;}
	my @nr=sort (keys %nr);
	return \@nr;
} ## nonredundant_list


# This subroutine is used to generate input file (matrix file)
# for MCL and index file.
# Four arguments:
# 1. String Variable: matrix file name
# 2. String Variable: index file name
# Last modified: 02/24/05
sub write_matrix_index {

	my $matrix_file  = $_[0];
	my $index_file   = $_[1];

	my $size = scalar(keys %graph);
	write_log("\nThere are $size sequences to cluster\n");
	open (MTX,">$matrix_file") or die "cannot write to file $matrix_file";
	print MTX "(mclheader\nmcltype matrix\ndimensions ".$size."x".$size."\n)\n\n(mclmatrix\nbegin\n\n";

	my $i=0;
	my %mcl_index2;
	foreach my $p (sort keys %graph) {
		$mcl_index2{$p}=$i;$mcl_index[$i]=$p;
		$i++;
	}
	foreach my $p (sort keys %graph) {

		print MTX $mcl_index2{$p}."    ";
		foreach my $m (@{$graph{$p}}) {
			print MTX $mcl_index2{$m}.":".$weight{$p.' '.$m}." ";
		}
		print MTX "\$\n";
	}
	print MTX ")\n\n";

	close (MTX);

	write_log("Matrix($size X $size) file $matrix_file generated\n");

	#################################WRITE INDEX FILE######################################
	open(IDX,">$index_file") or die "cannot write to file $index_file";
	foreach my $id (sort { $mcl_index2{$a} <=> $mcl_index2{$b} } keys %mcl_index2) {
		print IDX "$mcl_index2{$id}\t$id\n";
	}
	close(IDX);
	write_log("\nIndex file $index_file generated\n");
}


# This subroutine is used to back index all gene_ids present 
# in MCL output and generate the final result.
# Two arguments:
# 1. String Variable: mcl file name
# 2. String Variable: mcl back_index file name
# Last modified: 11/13/06
sub mcl_backindex {
	
	my $mcl_file      = $_[0];
	my $mcl_bi_file   = $_[1];

	open (MCL,$mcl_file) or die "can't open $mcl_file";
	my $last=0;
	my $lastline='';
	my @mcl;
	while (<MCL>) {
#		chomp;chop;  #bug, reported by Robson Francisco de Souza.    #this may result in wrong gene index by chopping the last digit of line    #for new format of MCL output. Not a bug for older versions, e.g. mcl-02-063    #replaced to the following line: substitute the '$' with nothing, and     #end reading data in the end of data block, preventing loading comment    #rows.
		chomp; s/\$$//; last if (/^\)/ && scalar(@mcl));

		if (/^(\d+)\s+(.*)\$/) {
			$mcl[$last]=$lastline;
			$last=$1;$lastline=$2;
		}
		elsif (/^(\d+)\s+(.*)/) {
			$mcl[$last]=$lastline;
			$last=$1;$lastline=$2;
		}
		elsif (/^\s+/) {$lastline.=$_;}
	}
	$mcl[$last]=$lastline;
	close (MCL);

	open (MCLBI,">$mcl_bi_file") or die "can't write to $mcl_bi_file";
	my $orthomcl_cluster_id=1;
	foreach my $mcl_cluster_id (0..$last) {
		$mcl[$mcl_cluster_id]=~s/\s+/ /g;
		my @g=split (" ",$mcl[$mcl_cluster_id]);
		next unless (scalar(@g)>=2);
		my @taxa=();
		foreach (@g) {
			my $taxon=$gindex2{$mcl_index[$_]};
			my $presence=0;
			foreach (@taxa) {
				if ($_ eq $taxon) {$presence=1;}
			}
			push (@taxa,$taxon) unless ($presence);
		}

		print MCLBI "homologous_gene_cluster_".$mcl_cluster_id."(".scalar(@g)." genes,".scalar(@taxa)." taxa):";
		foreach (@g) {
			print MCLBI " $mcl_index[$_]"; #modify it to product a similary format with group files generated by latest version of OrthoMCL software
			#print MCLBI " $mcl_index[$_]($gindex2{$mcl_index[$_]})";
		}
		print MCLBI "\n";
#No Species Cutoff
	}
        close MCLBI; # addied bu xiangyang 2019-11-25
	write_log("\n\nFinal ORTHOMCL Result: $mcl_bi_file generated!!!\n\n");

    return $mcl_bi_file; # addied bu xiangyang 2019-11-25
} ## mcl_backindex


# This subroutine is used to retrieve blast information from bpo file
# (blast parse out file), given the query gene id and the subject gene
# id. If such information doesn't exist, zero is returned.
# Four Arguments:
# 1. String Variable: query gene id
# 2. String Variable: subject gene id
# Last modified: 07/19/04
sub blastqueryab {
	my $a              = $_[0];
	my $b              = $_[1];
	my ($s,$e);
	if (defined $blastquery_ref->{$a}) {
		($s,$e)=split ("\t",$blastquery_ref->{$a});
	} else {return 0;}
    foreach my $i ($s..$e) {
		my @c=(&getline_from_bpofile($i))[0,1,3,5,6,7];
		if ($c[2] eq $b) {return(@c);}
	}
	return 0;
} ## blastqueryab

sub write_log {
	my $comment = $_[0];
	print LOG $comment;
	#print STDERR $comment;
}


sub write_bbh {
	my $comment = $_[0];
	print BBH $comment;
}


##################################################################################################
###### Subrounting--blast_homologs_cluster
###### Function:
###### Cluster analysis of homologs protein
##################################################################################################
sub blast_homologs_cluster {
    
    my ($directory_sub_blastout, $directory_protein_part, $interested_gene_name, $homologous_gene_cutoff, $thread_number, $blastp, $makeblastdb) = @_;
       

    
    system ("cat $directory_protein_part/*.fasta > $directory_sub_blastout/$interested_gene_name.fasta");
    system ("$makeblastdb -in $directory_sub_blastout/$interested_gene_name.fasta -dbtype prot -logfile $directory_sub_blastout/$interested_gene_name.fasta.makeblast.log");
##################################################################################################

    chdir $directory_protein_part;
    
    my @input;
    my @outfile;
    opendir PATH_PROTEIN_PART, $directory_protein_part or die "could not open $directory_protein_part";
    my @path_protein_part = readdir PATH_PROTEIN_PART; 
    @path_protein_part =grep ($_!~/^\./ ,@path_protein_part);  #delete hidden file . ..
    closedir PATH_PROTEIN_PART;

    open (GG,">$genome_gene_file");

    foreach my $each_protein_file(@path_protein_part){
 
        my $input="$directory_protein_part/$each_protein_file"; 
        push (@input,$input);  
        #print "input: $input\n";
        @input = sort @input;
        my $output="$directory_sub_blastout/$each_protein_file.blastout";
        push (@outfile,$output);
        @outfile = sort @outfile;

######################################################################################
        open (EACHFILE_GG, $input);

	print GG "$each_protein_file:";
	while (<EACHFILE_GG>) {
            chomp;
	    print GG " $1" if $_ =~ />(.*)/;
	} 
	    print GG "\n";

	close EACHFILE_GG;
######################################################################################
    }

    my $sub_file_number = scalar @input;

    my @sub_function_parameters = (\@input, "$directory_sub_blastout/$interested_gene_name.fasta", \@outfile, $e_value, $blastp);

    bacth_blast($sub_file_number, $thread_number, "do_blastP", \@sub_function_parameters);

    system ("cat $directory_sub_blastout/*.blastout > $blast_file");

    blast_parse($directory_sub_blastout, $blast_file, $bpo_file, $e_value, $identify, $coverage, $match_length_single);  #####################################################################

}

##################################################################################################
###### Subrounting--do_blastP
###### Function:
###### do blastp analysis
##################################################################################################
sub do_blastP {
    my ($input, $db_file, $outfile, $e_value, $blastp) =@_;    
    #print "test:----$input\n";
    system("$blastp -outfmt '6 qseqid sseqid pident qlen slen length qstart qend sstart send evalue bitscore' -evalue $e_value -num_threads 1 -query $input -db $db_file -out $outfile -num_alignments 100000");

}


##################################################################################################
###### Subrounting--do_batch_blastP
###### Function:
###### do bacth blastp analysis
##################################################################################################
sub bacth_blast {

    my ($work_number, $thread_number, $subfunction, $sub_function_parameters) = @_;
    my ($input_file, $db_file, $output_file, $e_value, $blastp) = @$sub_function_parameters;
 
    my @input = @$input_file;
    my @outfile = @$output_file;
    my $thread;
    my @threads;
    my $job_number=0;
    print "Blast_perform_percent: ";
    while(){ 
        last if ($job_number eq $work_number);                         
        while(scalar(threads->list())<$thread_number) {     #set threadnumber；
            $job_number++;                                 
            my $input_file = $input[$job_number-1];
            my $output_file = $outfile[$job_number-1];
            my $Blast_progress = int (($job_number/$work_number)*100);
            print "$Blast_progress","...";
            $threads[$job_number-1]=threads->new(\&$subfunction, $input_file, $db_file, $output_file, $e_value, $blastp);  
	    last if ($job_number eq $work_number);                         

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



sub blast_parse {

    my ($directory_sub_blastout, $blastFile, $parseoutfile, $e_value, $identify, $coverage, $match_length) = @_;

    open (PARSEOUT, ">$parseoutfile"); 
    open(F, $blastFile) || die "can't open BLAST file '$blastFile'\n";

    my $prevSubjectId = 'blah';
    my $prevQueryId = 'blah';
    my $subject;  # hash to hold subject info
    my $queryShorter;
    my $similarityid=0;
    my $simspanid;

    while(<F>) {
        chomp;
        my ($queryId, $subjectId, $percentIdentity, $querylen, $subjectlen, $length, $queryStart, $queryEnd, $subjectStart, $subjectEnd, $evalue, $bits) = split;


        if ($queryId ne $prevQueryId || $subjectId ne $prevSubjectId) {
            $simspanid =1;

            my $percentIdent = int($subject->{totalIdentities} / $subject->{totalLength} * 10 + .5)/10 if $subject;
	    print PARSEOUT "$similarityid\t$subject->{queryId}\t$subject->{queryLength}\t$subject->{subjectId}\t$subject->{subjectLength}\t$subject->{evalue}\t$percentIdent\t$subject->{simspan}\n" if $subject;
            $similarityid++;
	    # initialize new one from first HSP
	    $prevSubjectId = $subjectId;
            $prevQueryId = $queryId;

	    $subject = {}; 
	    $subject->{queryId} = $queryId;
	    $subject->{subjectId} = $subjectId;
	    $subject->{queryLength} = $querylen;
	    $subject->{subjectLength} = $subjectlen;
	    $subject->{evalue} = $evalue;
	
        }

        # get additional info from subsequent HSPs

	$subject->{simspan}.="$simspanid:$queryStart-$queryEnd:$subjectStart-$subjectEnd.";
	$simspanid++;

        $subject->{totalIdentities} += $percentIdentity * $length;
        $subject->{totalLength} += $length;
 

    }

        my $percentIdent = int($subject->{totalIdentities} / $subject->{totalLength} * 10 + .5)/10;
	print PARSEOUT "$similarityid\t$subject->{queryId}\t$subject->{queryLength}\t$subject->{subjectId}\t$subject->{subjectLength}\t$subject->{evalue}\t$percentIdent\t$subject->{simspan}\n";
    close PARSEOUT;
}


# this (corrected) version of formatEvalue provided by Robson de Souza
sub formatEvalue {
    my ($evalue) = @_;
    $evalue = '1' . $evalue if ($evalue =~ /^e/);
    $evalue = sprintf("%.3e",$evalue);
    my ($evalue_mant, $evalue_exp) = split(/e/, $evalue);
    $evalue_mant = sprintf("%.2f",$evalue_mant);
    $evalue_mant =~ s/\.0+$//;
    $evalue_exp =~ s/\+//;
    $evalue_exp = 0 if ($evalue_exp eq '00');
    return ($evalue_mant, $evalue_exp);
} #formatEvalue


# Caculate the minimum p-value (excluding evalue = 0) using blast software. For example, if 1e-180 is the case you should use -log(1e-181)=181.
sub maximum_weight {

    my $bpo_file = shift;
    my @min_value;
    open (BPO_FILE_1, $bpo_file);
    while (<BPO_FILE_1>){

        chomp;
        my @blast_parse_array = split ("\t", $_);
        my ($evalue_mant, $evalue_exp) = formatEvalue($blast_parse_array[5]);
        my $value = $blast_parse_array[5];
        push (@min_value, $value) if $evalue_mant ne 0;

    }
    close BPO_FILE_1;

    @min_value =sort{$a<=>$b}@min_value;
    my ($evalue_mant_min, $evalue_exp_min) = formatEvalue($min_value[0]);
    $evalue_exp_min = -($evalue_exp_min-1);
    return $evalue_exp_min;
} #maximum_weight


1;  

__END__   

