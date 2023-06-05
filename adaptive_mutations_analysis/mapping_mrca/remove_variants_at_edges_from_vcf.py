#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import vcf


# --------------------------------------------------------------------------------------------------------------------
# Notes
# --------------------------------------------------------------------------------------------------------------------

# This Python 3 script is used to remove variants at the start and end of sequences.
# This is intended to remove variants called at the edge of contigs that are not considered confident.
# Module PyVCF required. Version PyVCF-0.6.8 used
#   pip3 install PyVCF
# This script assumes the length of contigs are in the VCF header (##contig=<ID=contig4,length=212317>)

# --------------------------------------------------------------------------------------------------------------------
# Global variables
# --------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------------------------------------------


def parse_arguments():
    description = "Script to remove variants at the edge of sequences in VCF files"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-i", "--input_vcf", action="store", dest="input_vcf",
        help="input VCF file",
        required=True, metavar="INPUT_VCF"
    )
    group.add_argument(
        "-o", "--output_vcf", action="store", dest="output_vcf",
        help="output VCF file with kept variants",
        required=True, metavar="OUTPUT_VCF")
    group.add_argument(
        "-w", "--windom_length", action="store", dest="windom_length",
        help="length (in base pairs) from the start and to the end of sequences where "
             "variants should be deleted (default: %(default)s)",
        required=False, default="100", metavar="LENGTH")

    parser.add_argument(
        "-r", "--removed_vcf", action="store", dest="removed_vcf",
        help="output VCF file with removed variants",
        required=False, metavar="RM_VCF"
    )

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

    check_file_exists(args.input_vcf)
    vcf_reader = get_vcf_reader(args.input_vcf)
    contig_lengths = vcf_reader.contigs

    # Output VCFs
    vcf_writer = vcf.Writer(open(args.output_vcf, 'w'), vcf_reader)
    if args.removed_vcf is not None:
        vcf_writer_r = vcf.Writer(open(args.removed_vcf, 'w'), vcf_reader)

    for record in vcf_reader:
        contig_length = contig_lengths[record.CHROM].length  # length of the reference contig
        length_on_ref = len(record.REF) - len(record.ALT)  # variant length of the reference genome
        lower_offset = int(record.POS) - int(args.windom_length)  # allowed offset after start of contig
        upper_offset = int(record.POS) + int(length_on_ref) + int(args.windom_length)  # before end of contig
        if lower_offset <= 0 or upper_offset >= contig_length:
            logging.info(f'{record} removed')
            if args.removed_vcf is not None:
                vcf_writer_r.write_record(record)
        else:
            logging.info(f'{record} kept')
            vcf_writer.write_record(record)


if __name__ == "__main__":
    _main()