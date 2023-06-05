#!/usr/local/bin/perl -w

use strict;
use Bio::SeqIO;
 
# This Perl script creates the CC22 multiple alignment

my $cc_lists = "/project/directory/phylogeny_cc/colonising_collections_list_of_CCs.txt";
my $liExpLength = 2832299; # Expected alleles concatenated length
my $ref = "CC22_EMRSA15.fasta"; # do not change
unless(-e $ref){ print $ref." could not be found\n"; exit(0); }

open A, "<$cc_lists";
my @laCCs = <A>;
close A;

foreach my $cc (@ laCCs)
{
	chomp($cc);
        print $cc."\n";

	my $paths1 = "/project/directory/phylogeny_cc/colonising_collections_".$cc."_samples.txt";
	my $multiAlign = "/project/directory/phylogeny_cc/colonising_collections_".$cc.".EMRSA15.fa";
	my $log = "/project/directory/phylogeny_cc/colonising_collections_".$cc.".EMRSA15.log";
	unless(-e $paths1){ print $paths1." could not be found\n"; exit(0); }

	open A, "<$paths1";
	my @laLines = <A>;
	close A;

	open O, ">$multiAlign";
	open O2, ">$log";

	foreach my $sam (@laLines)
	{
	chomp($sam);
        print $sam."\n";

        my $align_file = "/project/directory/mapping/".$sam."_SMALT/".$sam.".mfa";
	print $align_file."\n";
        if(-e $align_file)
        {
		my $seqio_object = Bio::SeqIO->new(-file => $align_file);
		my $seq_object   = $seqio_object->next_seq;
		my $length = $seq_object->length;
		my $seq = $seq_object->seq;
		print $sam."\t".$length."\n";
		if($length == $liExpLength)
		{
			my $newfasta = ">".$sam."\n".$seq."\n";
			print O $newfasta;
		}else
		{
			print O2 $sam."\tlength_".$length."\n";
		}
	} else
	{
		print O2 $sam."\t.mfa_file_not_found\n";
	}
	}
	# Adding reference genome
	my $seqio_object = Bio::SeqIO->new(-file => $ref);
        my $seq_object   = $seqio_object->next_seq;
        my $length = $seq_object->length;
        my $seq = $seq_object->seq;
        print $ref."\t".$length."\n";
        if($length == $liExpLength)
        {
		my $newfasta = ">EMRSA15\n".$seq."\n";
		print O $newfasta;
        } else
        {
		print O2 "EMRSA15\tlength_".$length."\n";
	}
	close O;
	close O2;
}




