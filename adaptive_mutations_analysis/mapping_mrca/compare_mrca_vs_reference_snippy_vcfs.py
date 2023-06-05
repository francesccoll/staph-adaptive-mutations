#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import vcf
from Bio import SeqIO


# --------------------------------------------------------------------------------------------------------------------
# Notes
# --------------------------------------------------------------------------------------------------------------------

# This Python 3 script is used to remove variants in MRCA Snippy VCF files (resulting from mapping Fastq files to host
# MRCA reconstructed sequence using the identify_host_ancestral_isolate.step[].py pipeline) not present in the reference
# Snippy VCF files (resulting from mapping Fastq files to the complete reference genome), in an attempt to remove
# VCF variants arising from miss-assemblies in the out-group assembly and errors in ancestral reconstruction.


# --------------------------------------------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------------------------------------------


def parse_arguments():
    description = "Script to remove variants at the edge of sequences in VCF files"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-m", "--input_mrca_vcf", action="store", dest="input_mrca_vcf",
        help="input MRCA VCF file",
        required=True, metavar="INPUT_MRCA_VCF")
    group.add_argument(
        "-r1", "--input_ref_vcf1", action="store", dest="input_ref_vcf1",
        help="input reference VCF file",
        required=True, metavar="INPUT_REF_VCF1")
    group.add_argument(
        "-r2", "--input_ref_vcf2", action="store", dest="input_ref_vcf2",
        help="input reference VCF file",
        required=True, metavar="INPUT_REF_VCF2")
    group.add_argument(
        "-g", "--reference_genome", action="store", dest="reference_genome",
        help="FASTA sequence of the reference genome",
        required=True, metavar="REF")
    parser.add_argument(
        "-o", "--output_mrca_rm_vcf", action="store", dest="output_mrca_rm_vcf",
        help="output MRCA VCF file with variants not present in input reference VCF file",
        required=True, metavar="RM_VCF")
    group.add_argument(
        "-k", "--output_mrca_kept_vcf", action="store", dest="output_mrca_kept_vcf",
        help="output MRCA VCF file with variants also present in input reference VCF file",
        required=True, metavar="KEEP_VFC")

    return parser.parse_args()


def check_file_exists(my_file):
    if not os.path.isfile(my_file):
        logging.error(f'File {my_file} not found!')
        sys.exit(-1)


def get_vcf_reader(my_vcf):
    if os.path.splitext(my_vcf)[-1].lower() == '.gz':
        return vcf.Reader(open(my_vcf, 'rb'))
    else:
        return vcf.Reader(open(my_vcf, 'r'))


# ---------------------------------------------------------------------------------------------------------------------
# Main program
# ---------------------------------------------------------------------------------------------------------------------


def _main():
    # Configure logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s',
        level=logging.INFO
    )
    # Get arguments
    args = parse_arguments()

    check_file_exists(args.input_mrca_vcf)
    check_file_exists(args.input_ref_vcf1)
    check_file_exists(args.input_ref_vcf2)
    check_file_exists(args.reference_genome)

    # First, saving all variants in reference VCF file
    vcf_records_ref = dict()
    logging.info(f'Saving all variants in reference VCF file {args.input_ref_vcf1}...')
    vcf_reader1 = get_vcf_reader(args.input_ref_vcf1)
    for record in vcf_reader1:
        vcf_records_ref[record.POS] = record

    logging.info(f'Saving all variants in reference VCF file {args.input_ref_vcf2}...')
    vcf_reader2 = get_vcf_reader(args.input_ref_vcf2)
    for record in vcf_reader2:
        vcf_records_ref[record.POS] = record

    # Second, save reference genome
    logging.info(f'Loading reference genome {args.reference_genome}...')
    reference_records = SeqIO.parse(args.reference_genome, "fasta")
    reference_seq = ""
    chr_id = 'NA'
    for record in reference_records:
        chr_id = record.id
        chr_length = len(record.seq)
        reference_seq = record.seq
    logging.info(f'Extracted chromosome id: {chr_id}')

    #logging.info(f'Reference genome sequence extracted: {reference_seq[0:19]}...')

    # Third, load MRCA VCF variants and save those found in ref VCF or not
    logging.info(f'Loading all variants in MRCA VCF file {args.input_mrca_vcf}...')
    vcf_reader = get_vcf_reader(args.input_mrca_vcf)
    vcf_writer_kept = vcf.Writer(open(args.output_mrca_kept_vcf, 'w'), vcf_reader)
    vcf_writer_rm = vcf.Writer(open(args.output_mrca_rm_vcf, 'w'), vcf_reader)

    for record in vcf_reader:
        if str(record.CHROM) == str(chr_id):
            print("\n\n")
            if vcf_records_ref.get(record.POS) is None:
                print(str(record.POS) + ' not found in reference VCF file')
                # if position not found but sample's allele same as reference allele, keep
                var_start = int(record.POS) - 1
                var_end = var_start + len(record.ALT[0])
                print("var_start ", var_start, " var_end ", var_end)
                var_sequence = str(reference_seq[var_start:var_end])
                print("var_sequence ", var_sequence)
                if record.ALT == var_sequence:
                    record.FILTER = "ALT_MATCH"
                    vcf_writer_kept.write_record(record)
                    print('   but record.ALT ' + str(record.ALT) + ' same as reference allele ' +
                          var_sequence + '. Kept.')
                else:
                    record.FILTER = "NO_POS_MATCH"
                    vcf_writer_rm.write_record(record)
                    print('   and record.ALT ' + str(record.ALT) + ' different from reference allele ' +
                          var_sequence + '. Removed.')
            else:
                print(str(record.POS) + ' found in reference VCF file.')
                # Make sure REF and ALT match
                if vcf_records_ref.get(record.POS).REF == record.REF:
                    if vcf_records_ref.get(record.POS).ALT == record.ALT:
                        record.FILTER = "POS_REF_ALT_MATCH"
                        vcf_writer_kept.write_record(record)
                        print('  record.REF and record.ALT match. Kept.')
        else:
            logging.error(f'record.CHROM {record.CHROM} does not match chromosome id {chr_id} in '
                          f'{args.reference_genome}')
            exit(-1)

    logging.info(f'All done!')


if __name__ == "__main__":
    _main()