#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import Phylo
import dendropy

# ---------------------------------------------------------------------------------------------------------------------
# Notes
# ---------------------------------------------------------------------------------------------------------------------

# This python script is used to identify the most ancestral isolate among multiple isolates of the same host by:
#   1. identifying an appropriate out-group to all host isolates (the closest isolate from the MRCA),
#   2. identifying an appropriate reference genome (a close isolate to both the out-group and host isolates)
#       implemented in identify_host_ancestral_isolate.step1.py
#   3. mapping the out-group and host isolates to this reference genome and calling variants using Snippy,
#       implemented in identify_host_ancestral_isolate.step2.py
#   4. creating a high-resolution phylogenetic tree with these isolates,
#   5. reconstructing the ancestral sequence of all isolates of the same host,
#   6. and identifying the least diverged isolate (with the smallest number of changes) from this ancestral sequence
#       implemented in identify_host_ancestral_isolate.step3.py
#
# NOTES on expected input format and content:
#   - All isolates in INPUT_METADATA are expected in INPUT_TREE and INPUT_DISTANCES files but not the other way round
#   - This script assumes all isolate genomes from the same host are clonal (clonality)
#   - This script expects a rooted core-genome phylogeny in INPUT_TREE
#   - Genetic distances in INPUT_DISTANCES are accepted as a matrix of pairwise SNP distances produced by pairsnp
#   - The script expects the first column in INPUT_ASSEMBLIES and INPUT_FASTQ files to contain isolate Ids
#   - The script expects the second column in INPUT_ASSEMBLIES file to contain full paths to genome assemblies
#   - The script expects the second and third columns in INPUT_FASTQ file to contain full paths to _1.fastq.gz and
#       _2.fastq.gz files
#   - Isolate Ids must match across all input files

# NOTES on modules versions:
#   - Tested on Python 3.6.0 and [GCC 4.6.3] on linux
#   - Phylo Version: 4.2.1
#   - dendropy 4.2.0


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------


def parse_arguments():
    description = "Script to identify the most ancestral isolate among multiple isolates of the same host"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-m", "--input_metadata", action="store", dest="input_metadata",
        help="input csv metadata file with isolate and patient Ids",
        required=True, metavar="INPUT_METADATA"
    )
    group.add_argument(
        "-o", "--output_table", action="store", dest="output_table",
        help="output table produced by identify_host_ancestral_isolate.step1.py",
        required=True, metavar="OUTPUT_TABLE"
    )
    group.add_argument(
        "-a", "--input_assemblies", action="store", dest="input_assemblies",
        help="input csv file with isolate Ids and full path to their assemblies",
        required=True, metavar="INPUT_ASSEMBLIES"
    )
    group.add_argument(
        "-f", "--input_fastq", action="store", dest="input_fastq",
        help="input csv file with isolate Ids and full path to their fastq files",
        required=True, metavar="INPUT_FASTQ"
    )
    group.add_argument(
        "-v", "--vcf_dir", action="store", dest="vcf_dir",
        help="directory where Snippy VCF files are expected/saved",
        required=True, metavar="VCF_DIR"
    )
    group.add_argument(
        "-b", "--output_bash", action="store", dest="output_bash",
        help="Output bash script with snippy farm jobs",
        required=True, metavar="OUTPUT_BASH"
    )

    return parser.parse_args()


# ------------------------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------------------------

def _main():
    # Configure logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s',
        level=logging.INFO
    )
    # Get arguments
    args = parse_arguments()

    # Making sure input files exist
    input_files = [args.input_metadata, args.output_table, args.input_assemblies, args.input_fastq, args.vcf_dir]
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Extracting host and isolate Ids
    host_isolate_ids = dict()  # dictionary to save host and isolate ids as keys
    if args.input_metadata is not None:
        logging.info(f"Opening input file {args.input_metadata}")
        with open(args.input_metadata, 'r') as metadata_file:
            metadata_file.readline()
            for metadata_line in metadata_file:
                    (host_id, isolate_id) = metadata_line.strip().split('\t')
                    if host_isolate_ids.get(host_id) is None:
                        host_isolate_ids[host_id] = {}
                        host_isolate_ids[host_id][isolate_id] = '-'
                    else:
                        host_isolate_ids[host_id][isolate_id] = '-'

    # Extracting fastq directories
    isolate_fastq_dir = dict()  # dictionary to save isolate fastq directories
    if args.input_fastq is not None:
        logging.info(f"Opening input file {args.input_fastq}")
        with open(args.input_fastq, 'r') as input_fastq_file:
            input_fastq_file.readline()
            for input_fastq_line in input_fastq_file:
                    (isolate_id, fastq_dir) = input_fastq_line.strip().split('\t')
                    isolate_fastq_dir[isolate_id] = fastq_dir
                    # print('isolate_fastq_dir[' + isolate_id + '] ' + isolate_fastq_dir[isolate_id])

    # Extracting assemblies locations
    isolate_assemblies = dict()  # dictionary to save isolate fastq directories
    if args.input_assemblies is not None:
        logging.info(f"Opening input file {args.input_assemblies}")
        with open(args.input_assemblies, 'r') as input_assemblies_file:
            input_assemblies_file.readline()
            for input_assemblies_line in input_assemblies_file:
                    (isolate_id, assembly) = input_assemblies_line.strip().split('\t')
                    isolate_assemblies[isolate_id] = assembly
                    # print('isolate_assemblies[' + isolate_id + '] ' + isolate_assemblies[isolate_id])

    output_bash_file = open(args.output_bash, 'w')

    # Opening table with reference and out-group information
    if args.output_table is not None:
        logging.info(f"Opening input file {args.output_table}")
        with open(args.output_table, 'r') as output_file:
            output_file.readline()
            for output_line in output_file:
                line_items = output_line.strip().split('\t')
                host_id = line_items[0]
                host_isolates = line_items[1].split(';')
                outgroup_isolate = line_items[5]
                # reference_isolate = line_items[7]
                print(host_id + ' ' + '-'.join(host_isolates) + ' ' + outgroup_isolate)
                # host isolates need to be mapped against the reference, as wel as outgroup/reference (as control)
                isolates_to_be_mapped = []
                isolates_to_be_mapped.extend(host_isolates)
                isolates_to_be_mapped.append(outgroup_isolate)
                # isolates_to_be_mapped.append(reference_isolate)
                reference_isolate = outgroup_isolate
                for isolate_to_be_mapped in isolates_to_be_mapped:
                    pair_id = isolate_to_be_mapped + '.' + reference_isolate + '.' + host_id
                    snippy_out_dir = args.vcf_dir + pair_id
                    if isolate_assemblies.get(reference_isolate) is None:
                        logging.error(f'Assembly file for isolate {reference_isolate} not found!')
                        sys.exit(-1)
                    reference_assembly = isolate_assemblies[reference_isolate]
                    if isolate_fastq_dir.get(isolate_to_be_mapped) is None:
                        logging.error(f'Fastq dir for isolate {isolate_to_be_mapped} not found!')
                        sys.exit(-1)
                    fastq1_file = isolate_fastq_dir[isolate_to_be_mapped] + '/' + isolate_to_be_mapped + '_1.fastq.gz'
                    fastq2_file = isolate_fastq_dir[isolate_to_be_mapped] + '/' + isolate_to_be_mapped + '_2.fastq.gz'
                    output_snippy_vcf_file = args.vcf_dir + pair_id + '.snippy.vcf'
                    # Creating snippy farm jobs
                    bsub_line = 'bsub -q normal -G team81 -J ' + pair_id + ' -o ' + pair_id + '.out -n8 ' + \
                                '-R\"span[hosts=1] select[mem > 10000] rusage[mem=10000]\" -M 10000 ' + \
                                '\"bash run_snippy_per_sample.farm.sh ' + snippy_out_dir + ' ' + reference_assembly + \
                                ' ' + fastq1_file + ' ' + fastq2_file + ' ' + output_snippy_vcf_file + '\"'
                    print(bsub_line)
                    output_bash_file.write(bsub_line + '\n')


if __name__ == "__main__":
    _main()