#!/bin/bash

set -e

# This bash script is used to normalise variants produced by Snippy 

vcf_in=$1
vcf_out=$2
ref_fasta=$3
ref_id=$4
ref_length=$5

### Chromosome id

if [ "$ref_id" == "je2" ]
then
	ref_id="USA300_JE2";
fi

### File names

header1=`echo $vcf_in | sed 's/.vcf$/.header1.txt/g'`;

header_line=`echo $vcf_in | sed 's/.vcf$/.header_line.txt/g'`;

vcf_tmp1=`echo $vcf_in | sed 's/.vcf$/.header2.vcf/g'`;

vcf_tmp2=`echo $vcf_in | sed 's/.vcf$/.header2.sorted.vcf/g'`;


### Commands

# First, remove header ##contig lines containing original contig ids

cat $vcf_in | grep "^#" | grep -v "^##contig" > $header1


# Second, creating new header ##contig line (e.g. "##contig=<ID=NCTC8325,length=2821361>")

ref_id=${ref_id^^}; # to upper case

echo "##contig=<ID="$ref_id",length="$ref_length">" > $header_line


# Create new VCF file

num_header_lines=`cat $header1 | wc -l`;

insert_at=`expr $num_header_lines - 2`;

cat $header1 | head -n $insert_at > $vcf_tmp1

cat $header_line >> $vcf_tmp1

cat $header1 | tail -n 2 >> $vcf_tmp1

cat $vcf_in | grep -v "^#" >> $vcf_tmp1


# Sorting and indexing VCF

bcftools-1.9 sort -o $vcf_tmp2 -O v $vcf_tmp1

bgzip-1.9 $vcf_tmp2

tabix-1.9 -p vcf $vcf_tmp2".gz"


# Forth, normalising variants

bcftools-1.9 norm -f $ref_fasta $vcf_tmp2".gz" > $vcf_out

# Finally, removed temporary files

rm $vcf_tmp1 $vcf_tmp2".gz" $vcf_tmp2".gz.tbi" $header1 $header_line
