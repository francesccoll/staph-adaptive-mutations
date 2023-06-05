#!/usr/bin/perl -w

# FC notes
# This a modified version of the script by John Lees downloaded from: https://github.com/johnlees/paired-samples/tree/master/
# cite original article: J. A. Lees et al., Large scale genomic analysis shows no evidence for pathogen adaptation between the blood and cerebrospinal fluid niches during bacterial meningitis. Microbial Genomics. 3, 1–12 (2016).
# IMPORTANT NOTE: script tested with bcftools version 1.9 - it may not work with older versions
# The script (function blastn_ref) was modified to make sure that only the best blast hit is saved and to count the number of best hits. Before, the last/worst hit (was saved) if multiple hits found.
# This way only variants that map uniquely to the reference genome are kept.
# Variants whose reference or alternative allele does not match the new reference reference genome are discarded.


use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

# Allows use of perl modules in ./
# use Cwd 'abs_path';
# use File::Basename;
# use lib dirname( abs_path $0 );
# FC: line below commented out to load assembly_common module from same directory
# use lib dirname( abs_path $0 ) . "/../assembly_scripts";

# Perl modules - assembly
# FC note: assembly_common::standardise_contig_names not used
# use assembly_common;
# FC note: the functions defined in compare_variants.pm have been brought to this script to make sure it can run standalone
# use compare_variants;

#
# Help and usage
#
my $usage_message = <<USAGE;
Usage: ./insert_variants.pl --map_ref <assembly.fasta> --new_ref <assembly.fasta> --input_vcf <input_variants.vcf> --output_vcf <output_variants.vcf> --not_mapped_vcf <not_mapped_variants.vcf> --window_size <int> --window_size_per <int>

Changes coordinates of vcf between references

   Options
   --map_ref           Reference vcf is mapped to
   --new_ref           Reference to lift over to
   --input_vcf         Input VCF file to operate on
   --output_vcf        Output VCF file with coordinates on the new reference
   --not_mapped_vcf    Output VCF file with the variants that could not be mapped on the new reference
   --window_size       Window size in bp (e.g. 100)
   --window_size_per   Minimum blast alignment length as a proportion of window size (e.g. 0.5 for 50%)

   --dirty             Don't clean up temporary files

   -h, --help          Shows this help.

USAGE


#****************************************************************************************#
#* Functions                                                                            *#
#****************************************************************************************#


my $blastn_location = "blastn";
my $bcftools_location= "bcftools-1.9";

my %flip = (
   "N" => "N",
   "A" => "T",
   "T" => "A",
   "G" => "C",
   "C" => "G");

sub classify_var($$$$)
{
   my (@variants) = @_;

   # Any site > 1 base must be an indel
   my ($type, $match);

   # Check for indels
   my @del_bases;
   for (my $i = 0; $i<=3; $i+=2)
   {
      if (length($variants[$i]) > 1 || length($variants[$i+1]) > 1)
      {
         $type = "INDEL";
         $del_bases[$i/2] = decompose_indel($variants[$i], $variants[$i+1]);
      }
      else
      {
         $del_bases[$i/2] = "";
      }

   }

   # Classify snps
   if (!defined($type))
   {
      $type = "SNP";

      my ($q_ref, $q_alt, $s_ref, $s_alt) = @variants;
      if ($q_ref eq $s_ref && $q_alt eq $s_alt)
      {
         $match = "EXACT";
      }
      elsif ($q_ref eq $s_alt && $q_alt eq $s_ref)
      {
         $match = "SWITCHED";
      }
      # This distinction can't really be made for A/T and G/C SNPs
      elsif (flip_strand($q_ref) eq $s_ref && flip_strand($q_alt) eq $s_alt)
      {
         $match = "EXACT_FLIP";
      }
      elsif (flip_strand($q_ref) eq $s_alt && flip_strand($q_alt) eq $s_ref)
      {
         $match = "SWITCH_FLIP";
      }
      else
      {
         $match = "MISMATCH";
      }
   }
   # Classify indels
   else
   {
      if ($del_bases[0] eq $del_bases[1])
      {
         $match = "EXACT";
      }
      elsif (flip_strand($del_bases[0]) eq $del_bases[1])
      {
         $match = "EXACT_FLIP";
      }
      elsif (flip_strand($del_bases[0], "reverse") eq $del_bases[1])
      {
         $match = "SWITCH_FLIP";
      }
      else
      {
         $match = "MISMATCH";
      }
   }

   # Work out match status. Assembly strands/contigs may be flipped!
   return($type, $match);
}

# Reverse complement a sequence
sub flip_strand($;$)
{
   my ($sequence, $direction) = @_;

   # Set whether to also reverse sequence
   if (!defined($direction))
   {
      $direction = "forward";
   }

   # Decompose string into bases
   my @bases = split("", $sequence);
   my @flipped;

   # Flip each base, and append to either start or end of array
   foreach my $base (@bases)
   {
      if ($direction eq "reverse")
      {
         unshift(@flipped, $flip{$base});
      }
      else
      {
         push(@flipped, $flip{$base});
      }
   }

   # Join bases back into a string and return
   my $flipped_sequence = join("", @flipped);
   return($flipped_sequence);
}

# Needs left aligned INDELS. Use bcftools norm
sub decompose_indel($$)
{
   my ($ref, $alt) = @_;

   # Treat everything as deletions i.e. move longer alts to be the ref
   if (length($alt) > length($ref))
   {
      my $tmp_ref = $ref;
      $ref = $alt;
      $alt = $tmp_ref;
   }

   my $del_bases;
   # INDEL could go at either end, as either end of ref sequence matches alt
   if (substr($ref, 0, length($alt)) eq substr($ref, length($ref) - length($alt), length($alt)))
   {
      # Homopolymer indels will hit this. They need to be treated
      # as below
      if (length($ref) == 2*length($alt))
      {
         $del_bases = substr($ref, length($alt));
      }
      else
      {
         $del_bases = substr($ref, length($alt), length($ref) - 2*length($alt));
      }
   }
   # For left aligned indels, can just take bases of ref past alt
   else
   {
      $del_bases = substr($ref, length($alt));
   }

   #DEBUG
   #print STDERR "$ref,$alt,$del_bases\n";
   return($del_bases);
}

# list diff, union, intersection
sub list_compare($$$)
{
   my ($mode, $list1, $list2) = @_;

   my %count;
   my (@union, @intersection, @diff);

   foreach my $element (@$list1, @$list2)
   {
      $count{$element}++;
   }
   foreach my $element (keys %count)
   {
      push(@union, $element);

      if ($count{$element} > 1)
      {
         push(@intersection, $element);
      }
      else
      {
         push(@diff, $element);
      }
   }

   my $return_ref;
   if ($mode eq "diff")
   {
      $return_ref = \@diff;
   }
   elsif ($mode eq "intersection")
   {
      $return_ref = \@intersection;
   }
   elsif ($mode eq "union")
   {
      $return_ref = \@union;
   }
   else
   {
      die("Invalid list mode $mode\n");
   }

   return($return_ref);
}


# blastn of all v all. Returns reference to hash of sample by sample scores
sub blastn_pairwise($$)
{
   my ($query_file, $subject_file) = @_;

   my %blast_scores;

   # blastn of all sequences in query with all in subject. Output is tab
   # sepearated query_id, subject_id, e_value and score for each pairwise
   # comparison
   my $blastn_command = "$blastn_location -subject $subject_file -query $query_file -outfmt \"6 qseqid sallseqid evalue score\"";
   my $blast_result = `$blastn_command`;

   # Return data in a hash
   foreach my $pair (split("\n", $blast_result))
   {
      my ($q_id, $s_id, $e_val, $score) = split("\t", $pair);

      $blast_scores{$q_id}{$s_id} = $score;
   }

   return(\%blast_scores);
}

# blastn of all v ref. Returns reference to hash of variant and new position
# FC note: this function has been modified to make sure the top blast hit is selected
sub blastn_ref($$$$)
{
   my ($query_file, $ref_file, $window_size, $window_size_per) = @_;

   my %blast_scores;

   my $min_align_length = $window_size*$window_size_per;

   print "\$min_align_length --> ".$min_align_length."\n";

   # blastn of all sequences in query with all in subject. Output is tab
   # sepearated query_id, subject_id, e_value and score for each pairwise
   # comparison
   my $blastn_command = "$blastn_location -subject $ref_file -query $query_file -outfmt \"6 qseqid sstart send evalue score\"";
   my $blast_result = `$blastn_command`;
   print "\$blast_result:\n".$blast_result."\n";

   # Return data in a hash
   foreach my $pair (split("\n", $blast_result))
   {
      my ($q_id, $start, $end, $evalue, $align_length) = split("\t", $pair);

      # FC: code below modified to make sure only unique blast hit are saved by counting the number of best hits
      # If the blast hit alignment length is greater that the one chosen, the blast hit is counted
      # If the same variant has more than one blast hit counted that means it does not map uniquely to the reference and it will be discarded
      if($align_length > $min_align_length)
      {
        unless(defined $blast_scores{$q_id})
        {
      	  $blast_scores{$q_id}{start} = $start;
      	  $blast_scores{$q_id}{end} = $end;
      	  $blast_scores{$q_id}{count} = 1;
        } else
        {
	   $blast_scores{$q_id}{start} = $start;
      	   $blast_scores{$q_id}{end} = $end;
      	   $blast_scores{$q_id}{count} = $blast_scores{$q_id}{count} + 1;
        }
      }
   }

   return(\%blast_scores);
}

# Extracts windows of a specified size around a specified list of variants,
# creates a multi-fasta file of them
sub variant_windows($$$$;$)
{
   my ($window_size, $variant_list, $sequence_file, $output_file, $before) = @_;

   # Input and output multifasta
   my $sequence_in = Bio::SeqIO->new( -file   => "<$sequence_file",
                                      -format => "fasta" ) || die ($!);
   my $sequence_out = Bio::SeqIO->new( -file   => ">$output_file",
                                       -format => "fasta") || die ($!);

   print "\$sequence_in ".$sequence_file."\n";

   # Each fasta sequence, check all remaining variants
   my @found_variants;
   while (my $sequence = $sequence_in->next_seq())
   {
      foreach my $variant (@$variant_list)
      {
         my ($contig, $position, $ref, $alt, @samples) = split(",", $variant);

         if ($contig eq $sequence->display_id())
         {
            print "\n\n".$variant."\n";
	    print "\$contig ".$contig."\t"."\$sequence->display_id() ".$sequence->display_id()."\n"; # FC line

            # Remove from variant list if found in this contig
	    # FC note: this is to make sure the contig each variant is annotated on is found in the reference genome
            push (@found_variants, $variant);

            # Variant length. 1 for a SNP, > 1 for an insertion, < 1 for
            # a deletion
            my $var_length = length($alt) - length($ref) + 1;
            my $position_start = $position - 1;
            my $position_end;

            # SNPs
            if ($var_length == 1)
            {
               $position_end = $position + $var_length;
            }
            # Deletions
            elsif ($var_length < 1)
            {
               $position_end = $position - $var_length + 2;
            }
            # Insertions
            else
            {
               $position_end = $position + length($ref);
            }

	    print "\$window_size ".$window_size." \$position ".$position." \$position_start ".$position_start." \$position_end ".$position_end." \$sequence->length() ".$sequence->length()."\n";

            # Extract the window, then write it to the output
            my $window;
            if ($before)
            {
               print "Running extract_bp_before function\n";
               $window = extract_bp_before($sequence, $window_size, $position_start);
            }
            else
            {
               print "Running extract_bp_window function\n";
               $window = extract_bp_window($sequence, $window_size, $position_start, $position_end);
            }

	    if(defined $window)
            {
            	$window->display_id($variant);
            	$sequence_out->write_seq($window);
	    }
         }
      }
   }

   # Print a list of missed variants
   my $missed_variants = list_compare("diff", $variant_list, \@found_variants);
   foreach my $variant (@$missed_variants)
   {
      print STDERR "Variant $variant not found in input $sequence_file\n";
   }
}

# Extracts sequence of a window size around a position, returning a Bio:Seq
# object of this
# FC: function not used
sub extract_bp_window($$$;$)
{
   my ($sequence, $window_size, $position_start, $position_end) = @_;

   print "extract_bp_window: \$sequence ".$sequence." \$window_size ".$window_size." \$position_start ".$position_start." \$position_end ".$position_end."\n";

   # Window size must be even, as we will extract a region symmetric around the
   # variant not including it
   if ($window_size % 2 != 0)
   {
      $window_size--;
   }
   my $half_size = ($window_size)/2;

   my ($window, $start, $end);

   if ($sequence->length() < ($window_size + $position_end - $position_start))
   {
      print STDERR "Cannot extract $window_size bp around $position_start as sequence too small\n";
   }
   else
   {
      if (($position_end + $half_size) > $sequence->length())
      {
         # Extract from end backwards
         # FC: if right window longer than contig length, then extract from end backwards
         $start = $sequence->length() - ($window_size + $position_end - $position_start);
         $end = $sequence->length();
      }
      elsif (($position_start - $half_size) <= 0)
      {
         # Extract from start forwards
         # FC: if left window smaller than 0, then extract from start forwards
         $start = 1;
         $end = $window_size + ($position_end - $position_start);
      }
      else
      {
         # FC: else, extract simetric windows at both sides
         $start = $position_start - $half_size + 1;
         $end = $position_end + $half_size - 1;
      }

      my $window_sequence = $sequence->subseq($start, $position_start) . $sequence->subseq($position_end, $end);
      $window = Bio::Seq->new( -seq => $window_sequence);
   }

   return($window);
}

# Extracts sequence of a window size around a position, returning a Bio::Seq
# object of this
# FC: Extracts sequence of a window size BEFORE a position
sub extract_bp_before($$$)
{
   my ($sequence, $seq_size, $position_start) = @_;

   print "extract_bp_before: \$sequence ".$sequence." \$seq_size ".$seq_size." \$position_start ".$position_start."\n";

   my $context;
   if ($sequence->length() < $seq_size || $position_start < ($seq_size + 1) )
   {
      print STDERR "Cannot extract $seq_size bp before $position_start as sequence too small\n";
   }
   else
   {
      my $context_sequence = $sequence->subseq($position_start-$seq_size, $position_start);
      $context = Bio::Seq->new( -seq => $context_sequence);
   }

   return($context);
}

# Takes a list of paired vcfs and reference fastas and creates a multifasta of
# windows around the variants
# FC: function not used
sub create_blastn_input($$$;$$)
{
   my ($vcfs, $refs, $out_prefix, $filter, $type) = @_;

   my $i = 1;
   foreach my $vcf (@$vcfs)
   {
      my $variant_lists = extract_vcf_variants($vcf, $filter, $type);
      variant_windows(300, $variant_lists, shift(@$refs), "$out_prefix.$i.fa");

      $i++;
   }

}

# Get a list of variant sites from a vcf
sub extract_vcf_variants($;$$)
{
   my ($vcf_in, $filter, $type) = @_;

   # Tab delimited contig, position list of all variants in vcf
   my $bcftools_command;
   if (defined($type) && defined($filter))
   {
      $bcftools_command = "$bcftools_location view -f $filter -v $type $vcf_in | $bcftools_location query -f '%CHROM\t%POS\t%REF\t%ALT\n' -";
   }
   else
   {
      $bcftools_command = "$bcftools_location view $vcf_in | $bcftools_location query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n' -";
   }
   my $bcf_return = `$bcftools_command`;

   # Create a hash of variant details
   my @variant_list;
   foreach my $variant (split("\n", $bcf_return))
   {
      $variant =~ s/\t/,/g;
      push (@variant_list, $variant);
   }

   return(\@variant_list);
}


#****************************************************************************************#
#* Main                                                                                 *#
#****************************************************************************************#




#* gets input parameters
my ($map_ref, $new_ref, $vcf_in, $vcf_out, $vcf_nm_out, $window_size, $window_size_per, $dirty, $help);
GetOptions ("map_ref=s"  => \$map_ref,
            "new_ref=s"  => \$new_ref,
            "input_vcf=s"      => \$vcf_in,
	    "output_vcf=s"      => \$vcf_out,
	    "not_mapped_vcf=s"      => \$vcf_nm_out,
            "window_size=s"      => \$window_size,
	    "window_size_per=s"      => \$window_size_per,
            "dirty"      => \$dirty,
            "help|h"     => \$help
		   ) or die($usage_message);


# Globals

my $prefix = $vcf_in; $prefix =~ s/\.vcf$//g;
# my $tmp_ref = $prefix.".reference_renamed.fa";
my $blast_prefix = $prefix.".blast_windows";


my $new_chrom = "NA";

open O, ">$vcf_out";
open O2, ">$vcf_nm_out";

# Parse input
if (defined($help))
{
   print STDERR $usage_message;
}
elsif (!defined($map_ref) || !defined($new_ref) || !defined($vcf_in))
{
   print STDERR $usage_message;
}
else
{
   # FC note: the line below is commented out to make sure contig names matched between map_ref and vcf_in
   # assembly_common::standardise_contig_names($map_ref, $tmp_ref);
   # $tmp_ref = $map_ref; # FC added

   my $blast_output = "$blast_prefix.fa";
   my $variant_lists = extract_vcf_variants($vcf_in);
   print "\$variant_lists\n".join("\n",@{$variant_lists})."\n";
   
   variant_windows($window_size, $variant_lists, $map_ref, $blast_output, 1); # Only function extract_bp_before is used
   print "\n\n\$blast_output: ".$blast_output." \$new_ref ".$new_ref."\n";

   my $blast_scores = blastn_ref($blast_output, $new_ref, $window_size, $window_size_per);

   my $sequence_in = Bio::SeqIO->new( -file   => "<$new_ref",
                                      -format => "fasta" ) || die ($!);
   my $new_ref_sequence = $sequence_in->next_seq();

   $new_chrom = $new_ref_sequence->id();

   #copy header then print:
   #vcf lines, tab separated
   #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample
   #AE007317 start-1 . A C . PASS . GT 1
   #
   #so ends up as its own vcf
   #
   my $vcf_header = `bcftools view -h $vcf_in`;
   print $vcf_header;  print O $vcf_header; print O2 $vcf_header;

   # for each blast line: e.g. contig7,118277,T,C	2819424	2819624	4.34e-104	201
   foreach my $q_id (sort keys %$blast_scores)
   {
      my ($chrom, $pos, $ref, $alt, @samples) = split(",", $q_id);

      my $new_ref = $ref;
      my $new_alt = $alt;
      my $mapped = "no"; # variable to store whether variants could be both uniquely mapped and new_ref matches new reference genome
      my $not_mapped_reason = "not_determined";

      # FC note: make sure variant's window has only one top blast hit, that is, maps uniquely to reference genome
      if($$blast_scores{$q_id}{count} == 1)
      {
	$mapped = "yes";
      } else
      {
        $mapped = "no";
	$not_mapped_reason = "Variant mapped to multiple locations on reference genome ".$q_id;
        print STDERR "WARNING: Variant mapped to multiple locations on reference genome $q_id\n";
      }

      if ($$blast_scores{$q_id}{start} > $$blast_scores{$q_id}{end})
      {
	 # Variants on the negative strand in the original genome need to be passed to the positive strand on the target genome
      	 print "\n\$q_id ".$q_id."\n";

	 # FC: $$blast_scores{$q_id}{end} corresponds to $position_start in extract_bp_before,
	 # which corresponds to ($position - 1) in variant_windows 
        
	 # If SNP 
         if (length($alt) == 1 && length($ref) == 1)
         {
	   print "SNP annotated on the negative strand\n";
           $pos = $$blast_scores{$q_id}{end} - 1;
	   $new_ref = flip_strand($ref, "reverse");
           $new_alt = flip_strand($alt, "reverse");
         }
	 # If insertion
	 if (length($alt) > length($ref))
	 {
	   print "Insertion annotated on the negative strand\n";
           $pos = $$blast_scores{$q_id}{end} - 2;
	   my $inserted_seq = $alt;
	   $inserted_seq =~ s/^$ref//g;
           $new_ref = $new_ref_sequence->subseq($pos, $pos);
	   $new_alt = $new_ref.flip_strand($inserted_seq, "reverse");
	 }
	 # If deletion
	 if (length($alt) < length($ref))
	 {
	   print "Deletion annotated on the negative strand\n";
	   my $deleted_seq = $ref;
	   $deleted_seq =~ s/^$alt//g;
           $pos = $$blast_scores{$q_id}{end} - 1 - length($deleted_seq) - 1;
           $new_alt = $new_ref_sequence->subseq($pos, $pos);
	   $new_ref = $new_alt.flip_strand($deleted_seq, "reverse");
	 }
         # Complex SNP
	 if (length($new_alt) == length($new_ref) && length($new_ref) > 1)
         {
           $pos = $$blast_scores{$q_id}{end} - length($new_alt);
	   $new_ref = flip_strand($ref, "reverse");
           $new_alt = flip_strand($alt, "reverse");
         }
         # indels include an extra base
         # else
         # {
         #   $pos = $$blast_scores{$q_id}{end} - (abs(length($new_alt) - length($new_ref)) + 1) - 1;
         # }

	 print "\$ref ".$ref."\n";
	 print "\$new_ref ".$new_ref."\n";
         print "\$alt ".$alt."\n";
	 print "\$new_alt ".$new_alt."\n";
      }
      else
      {
         print "\n\$q_id ".$q_id."\n";
	 print "\$ref ".$ref."\n";
         print "\$alt ".$alt."\n";
         $pos = $$blast_scores{$q_id}{end} + 1;
	 # FC: $$blast_scores{$q_id}{end} corresponds to $position_start in extract_bp_before,
	 # which corresponds to ($position - 1) in variant_windows 
        
         # if (length($new_alt) == 1 && length($new_ref) == 1)
         # {
         #   $pos = $$blast_scores{$q_id}{start} - 1;
         # }
         # else
         # {
         #   $pos = $$blast_scores{$q_id}{start};
         # }
      }

      # Check which is really the ref
      my $flipped;
      if ($new_ref_sequence->subseq($pos, $pos + length($new_ref) - 1) =~ /^$new_ref$/i)
      {
         $flipped = 0;
      }
      elsif ($new_ref_sequence->subseq($pos, $pos + length($new_alt) - 1) =~ /^$new_alt$/i)
      {
         $flipped = 1;

         my $tmp_store = $new_alt;
         $new_alt = $new_ref;
         $new_ref = $tmp_store;
      }
      else
      {
        $mapped = "no";
	$not_mapped_reason = "Reference or alternative allele does not match sequence in new reference genome ".$q_id;
        print STDERR "WARNING: Could not map $q_id\n";
      }

      my $sample_str = "";
      foreach my $sample_gt (@samples)
      {
         if ((!$flipped && $sample_gt =~ /^$ref(\/|$)/i) || ($flipped && $sample_gt =~ /^$alt(\/|$)/i))
         {
            $sample_str .= "0\t";
         }
         else
         {
            $sample_str .= "1\t";
         }
      }

      if ($mapped eq "yes")
      {
         my $vcf_new_line = join("\t", $new_chrom, $pos, ".", $new_ref, $new_alt, ".", "PASS", ".", "GT", $sample_str . "\n");
         print $vcf_new_line; print O $vcf_new_line;
      }
      else
      {
         my $vcf_new_line2 = join("\t", $new_chrom, $pos, ".", $new_ref, $new_alt, ".", "PASS", ".", "GT", $sample_str, $not_mapped_reason . "\n");
         print $vcf_new_line2; print O2 $vcf_new_line2;
      }

   } # end of foreach my $q_id

   unless($dirty)
   {
      # unlink $tmp_ref, $blast_output;
      unlink $blast_output;
  }
}

exit(0);

