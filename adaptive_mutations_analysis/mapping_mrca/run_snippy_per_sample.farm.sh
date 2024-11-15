#!/bin/bash

set -e

# This bash script is used to run snippy, copy the final VCF file and delete the whole snippy directory

out_dir=$1;
reference=$2;
fastq1=$3;
fastq2=$4;
vcf_out=$5;

snippy --outdir $out_dir --ref $reference --R1 $fastq1 --R2 $fastq2
cp $out_dir"/snps.filt.vcf" $vcf_out
rm -r $out_dir



