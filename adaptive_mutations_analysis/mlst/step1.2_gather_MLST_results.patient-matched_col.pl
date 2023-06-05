#!/usr/local/bin/perl -w

use strict;


my $dir_file = "/project/directory/data/colonising_pairs_patient-matched_gwas_dataset.all_pairs.paths.txt";
my $mlst_dir = "/project/directory/mlst/mlst_results/";
my $mlst_table = "/project/directory/mlst/colonising_pairs_patient-matched_gwas_dataset.all_pairs.mlst.txt";


open A, "<$dir_file";
my @laLines = <A>;
close A;

open O, ">$mlst_table";
my $header = "sample\tST\tNewST\tContamination\tarcc\taroe\tglpf\tgmk\tpta\ttpi\tyqil\n";
print O $header;

foreach my $line (@laLines)
{
	#print $line;
	chomp($line);
	my @a = split("/",$line);
	my $sam = $a[$#a];
	print $sam."\n";
	my $mlst_file = $mlst_dir.$sam."/mlst_results.allele.csv";
	my $aa = $sam."\tND\tND\tND\tND\tND\tND\tND\tND\tND\tND\n";
	if(-e $mlst_file)
	{
		open B, "<$mlst_file";
		my @laLines2 = <B>;
		close B;
		#print $laLines2[1]."\n";
		my @b = split("\t",$laLines2[1]);
		$aa = $sam."\t".join("\t",@b[1..$#b]);
	} else
	{
		print "MLST file not found\n";
	}
	print O $aa;
	print $aa;
}


