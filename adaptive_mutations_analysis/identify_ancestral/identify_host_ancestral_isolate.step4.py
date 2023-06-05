#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
import dendropy
import vcf
from Bio import Phylo


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

# NOTES on software and modules versions. Tested on
#   - Python 3.6.0 and [GCC 4.6.3] on linux
#   - Phylo Version: 4.2.1
#   - dendropy 4.2.0
#   - bcftools 1.9. Change bcftools executable if necessary
#   - RAxML version 8.2.8
#   - snp-sites 2.4.1
#   - pastml 1.9.20

# NOTES on the use of RAxML
#   - The single-threaded (Serial) version of RAxML (raxmlHPC-SSE3) is used as small alignments are expected
#   - The "rapid Bootstrap analysis and search for best-scoring ML tree in one program run" is used ("-f a")


# ------------------------------------------------------------------------------------
# Global variables
# ------------------------------------------------------------------------------------

_DEPENDENCIES = ['pastml']


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to identify the most ancestral isolate among multiple isolates of the same host (step 3)"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
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
        "-v", "--vcf_dir", action="store", dest="vcf_dir",
        help="directory where Snippy VCF files are expected/saved",
        required=True, metavar="VCF_DIR"
    )
    group.add_argument(
        "-r", "--raxml_dir", action="store", dest="raxml_dir",
        help="directory where RAxML files are expected/saved",
        required=True, metavar="RAXML_DIR"
    )

    return parser.parse_args()


# def get_node_labels(node_list):
#     """
#     Given a list of Node object, returns their taxa labels
#     :param node_list: list of Node objects
#     :return: list of taxa labels
#     """
#     node_labels = list()
#     for i, node in enumerate(node_list):
#         node_labels.append(node_list[i].taxon.label)
#     return node_labels


def get_node_labels(node_list):
    """
    Given a list of Node object, returns their taxa labels (Phylo Tree object)
    :param node_list: list of Node objects
    :return: list of taxa labels
    """
    node_labels = list()
    for i, node in enumerate(node_list):
        node_labels.append(node_list[i].name)
    return node_labels


def run_command_shell(command_line):
    """
    This function executes a command line, check for execution errors and but does not return stdout
    This is to be used when the stdout is not needed
    Note: shell=True needs to be set if I/O redirection operators are to be used (e.g. >) in the command line,
    otherwise they will have no special meaning, they are treated as ordinary arguments
    Note: if shell=True is used then the command line must be provided as a string, not a list
    :param command_line: it must be a list not a string
    """
    command_line_string = ' '.join(command_line)
    try:
        subprocess.run(command_line_string,
                       check=True,
                       shell=True,
                       )
    except subprocess.CalledProcessError as err:
        print('ERROR:', err)

#
# def is_tree_valid(input_tree):
#     try:
#         Phylo.read(input_tree, 'newick')
#         tree = dendropy.Tree.get_from_path(input_tree, 'newick')
#     except:
#         print("Error with the input starting tree: Is it a valid Newick file?")
#         return 0
#     return 1


def check_dependency(executable_name):
    """ Returns true if executable exists, else false """
    found = False
    try:
        output = subprocess.check_output(['which', executable_name]).strip()
        if output:
            found = True
    except subprocess.CalledProcessError as err:
        print('ERROR:', err)
    return found


def check_file_exists(my_file):
    if not os.path.isfile(my_file):
        logging.error(f'File {my_file} not found!')
        sys.exit(-1)


def get_vcf_reader(my_vcf):
    if os.path.splitext(my_vcf)[-1].lower() == '.gz':
        return vcf.Reader(open(my_vcf, 'rb'))
    else:
        return vcf.Reader(open(my_vcf, 'r'))


def create_mutation(mutable_seq, mutation_pos, mutation_ref, mutation_alt, mutation_type):
    """
    Function to introduce mutations into a Bio:Seq object
    :param mutable_seq: a MutableSeq object.
    :param mutation_pos: a 1-indexed VCF position
    :param mutation_ref: VCF REF allele
    :param mutation_alt: VCF ALT allele
    :param mutation_type: VCF INFO TYPE (snp, mnp, ins, del or complex)
    :return:
    """

    # make sure mutation_type is one of the five expected
    if mutation_type not in ['snp', 'mnp', 'ins', 'del', 'complex']:
        logging.error(f'Unknown TYPE from VCF file \'{mutation_type}\' at position {mutation_pos}. Exiting...')
        sys.exit(-1)

    # make sure VCF reference allele matches the one extracted from reference genome
    mutation_pos_to = mutation_pos + len(mutation_ref)
    genome_sequence = mutable_seq[mutation_pos-1:mutation_pos_to-1]
    if mutation_ref != str(genome_sequence).upper():
        logging.error(f'Reference allele from VCF file {mutation_ref} does not match that '
                      f'from reference genome {genome_sequence}. Exiting...')

    if mutation_type == 'snp':
        # make sure both REF and ALT alleles have length of 1
        if len(mutation_ref) != 1 & len(mutation_alt) != 1:
            logging.error(f'REF allele from VCF file {mutation_ref} OR ALT allele from VCF file {mutation_alt} '
                          f'have lengths different from 1 at position {mutation_pos}. Exiting...')
            sys.exit(-1)
        print(f"\tReplacing SNP reference allele {mutable_seq[mutation_pos-1]} with {mutation_alt} at position "
              f"{mutation_pos}")
        print('Before: ' + str(mutable_seq[mutation_pos-10:mutation_pos+10]))
        mutable_seq[mutation_pos-1] = mutation_alt
        print('Alter: ' + str(mutable_seq[mutation_pos-10:mutation_pos+10]))

    if mutation_type == 'del':
        # make sure length of REF > length of ALT
        if len(mutation_ref) <= len(mutation_alt):
            logging.error(f'REF allele from VCF file {mutation_ref} must be longer than from ALT {mutation_alt} '
                          f'for the deletion at position {mutation_pos}. Exiting...')
            sys.exit(-1)
        print(f"\tDeleting REF bases {mutable_seq[mutation_pos:mutation_pos_to-1]} at position {mutation_pos}")
        print('Before: ' + str(mutable_seq[mutation_pos-10:mutation_pos+10]))
        mutable_seq = mutable_seq[0:mutation_pos] + mutable_seq[mutation_pos_to-1:]
        print('After: ' + str(mutable_seq[mutation_pos-10:mutation_pos+10]))

    if mutation_type == 'ins':
        # make sure length of ALT > length of REF
        if len(mutation_ref) >= len(mutation_alt):
            logging.error(f'ALT allele from VCF file {mutation_ref} must be longer than from REF {mutation_alt} '
                          f'for the insertion at position {mutation_pos}. Exiting...')
            sys.exit(-1)
        if len(mutation_ref) != 1:
            logging.error(f'REF allele from VCF file {mutation_ref} extected to be of length 1 for insertion at '
                          f'position {mutation_pos}. Exiting...')
            sys.exit(-1)
        print(f"\tInserting ALT bases {mutation_alt} at position {mutation_pos}")
        print('Before: ' + str(mutable_seq[mutation_pos-10:mutation_pos+10]))
        mutable_seq = mutable_seq[0:mutation_pos-1] + mutation_alt + mutable_seq[mutation_pos_to-1:]
        print('After: ' + str(mutable_seq[mutation_pos-10:mutation_pos+10]))

    if mutation_type == 'mnp':
        # make sure both REF and ALT alleles have same lengths
        if len(mutation_ref) != len(mutation_alt):
            logging.error(f'REF allele from VCF file {mutation_ref} AND ALT allele from VCF file {mutation_alt} '
                          f'have different lengths at position {mutation_pos}. Exiting...')
            sys.exit(-1)
        print(f"\tReplacing MNP reference allele {mutable_seq[mutation_pos-1:mutation_pos_to-1]} with {mutation_alt} "
              f"at position {mutation_pos}")
        print('Before: ' + str(mutable_seq[mutation_pos-10:mutation_pos+10]))
        mutable_seq[mutation_pos-1:mutation_pos_to-1] = mutation_alt
        print('After: ' + str(mutable_seq[mutation_pos-10:mutation_pos+10]))

    if mutation_type == 'complex':
        # if REF and ALT alleles have same lengths then same logic as for 'mnp' variants
        if len(mutation_ref) == len(mutation_alt):
            print(f"\tReplacing COMPLEX reference allele {mutable_seq[mutation_pos-1:mutation_pos_to-1]} with "
                  f"{mutation_alt} at position {mutation_pos}")
            print('Before: ' + str(mutable_seq[mutation_pos - 10:mutation_pos + 10]))
            mutable_seq[mutation_pos-1:mutation_pos_to-1] = mutation_alt
            print('After: ' + str(mutable_seq[mutation_pos - 10:mutation_pos + 10]))
        # if length of REF > length of ALT, deletion
        if len(mutation_ref) > len(mutation_alt):
            # e.g. CTGTCTCTTATA	TG
            # First, replace MNP (e.g. CT with TG)
            mutation_pos_to = mutation_pos + len(mutation_alt)
            print(f"\tReplacing COMPLEX reference allele {mutable_seq[mutation_pos-1:mutation_pos_to-1]} with "
                  f"{mutation_alt} at position {mutation_pos}")
            print('Before: ' + str(mutable_seq[mutation_pos - 10:mutation_pos + 20]))
            mutable_seq[mutation_pos - 1:mutation_pos_to - 1] = mutation_alt
            print('After 1: ' + str(mutable_seq[mutation_pos - 10:mutation_pos + 20]))
            # Second, delete deleted sequence (e.g. GTCTCTTATA)
            mutation_pos_to_del = mutation_pos + len(mutation_ref)
            print(f"\tDeleting REF bases {mutable_seq[mutation_pos_to-1:mutation_pos_to_del-1]} at position "
                  f"{mutation_pos_to}")
            mutable_seq = mutable_seq[0:mutation_pos_to - 1] + mutable_seq[mutation_pos_to_del-1:]
            print('After 2: ' + str(mutable_seq[mutation_pos - 10:mutation_pos + 20]))
        if len(mutation_ref) < len(mutation_alt):
            print('len(mutation_ref) < len(mutation_alt)')
            # e.g. AAAGT	GAAAGA
            mutation_pos_to = mutation_pos + len(mutation_ref)
            print(f"\tReplacing COMPLEX reference allele {mutable_seq[mutation_pos-1:mutation_pos_to-1]} with "
                  f"{mutation_alt} at position {mutation_pos}")
            print('Before: ' + str(mutable_seq[mutation_pos - 10:mutation_pos + 20]))
            mutable_seq[mutation_pos - 1:mutation_pos_to - 1] = mutation_alt
            print('After: ' + str(mutable_seq[mutation_pos - 10:mutation_pos + 20]))

    return mutable_seq


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

    # Making sure dependencies exist
    logging.info('Making sure dependencies exist...')
    for dependency in _DEPENDENCIES:
        if check_dependency(dependency):
            logging.info(f'{dependency} is installed!')
        else:
            logging.error(f'{dependency} is NOT installed!')
            sys.exit(-1)

    # Making sure input files exist
    input_files = [args.output_table, args.vcf_dir, args.raxml_dir]
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Extracting assemblies locations
    isolate_assemblies = dict()  # dictionary to save isolate assemblies
    if args.input_assemblies is not None:
        logging.info(f"Opening input file {args.input_assemblies}")
        with open(args.input_assemblies, 'r') as input_assemblies_file:
            input_assemblies_file.readline()
            for input_assemblies_line in input_assemblies_file:
                (isolate_id, assembly) = input_assemblies_line.strip().split('\t')
                isolate_assemblies[isolate_id] = assembly

    # Opening table with reference and out-group information
    if args.output_table is not None:
        logging.info(f"Opening input file {args.output_table}")
        with open(args.output_table, 'r') as output_file:
            output_file.readline()
            for output_line in output_file:
                line_items = output_line.strip().split('\t')
                host_id = line_items[0]
                # if host_id == "mpros-1328":
                # if host_id == "mpros-623":
                # if host_id == "mpros-350":
                host_isolates = line_items[1].split(';')
                outgroup_isolate = line_items[5]
                # reference_isolate = line_items[7]
                reference_isolate = outgroup_isolate
                print(host_id + ' ' + '-'.join(host_isolates) + ' ' + outgroup_isolate)

                # Making sure final host MRCA fasta file does not exist already
                fasta_id = host_id + '.mrca.' + outgroup_isolate
                fasta_file = args.vcf_dir + fasta_id + '.fasta'

                if not os.path.isfile(fasta_file):
                    # Make sure final phylogenetic tree exists (produced by)
                    final_tree_file = args.raxml_dir + "RAxML_bipartitions." + host_id + '.' + reference_isolate + \
                                      '.raxml_tree' + "_rooted"

                    if os.path.isfile(final_tree_file):
                        logging.info(f"Input tree file {final_tree_file} found. Moving on.")
                    else:
                        logging.error(f'Input tree file {final_tree_file} not found. Exiting...')
                        sys.exit(-1)

                    # NOTE: dendropy functions to_outgroup_position and reroot_at_edge were not successful at
                    # rooting the tree properly, so Phylo method root_with_outgroup had to be used instead
                    # NOTE: avoid using both dendropy and Phylo read/write methods in the same script as newick
                    # formats differ slightly leading to errors in read interpretation
                    logging.info(f"Reading input tree file {final_tree_file}")
                    # tree = dendropy.Tree.get_from_path(final_tree_file, 'newick', preserve_underscores=True)
                    tree = Phylo.read(final_tree_file, "newick")

                    # Printing basic information on tree
                    logging.info(f"Printing basic information on input tree")
                    # print('Number of taxa in tree: ' + str(len(tree.taxon_namespace)))
                    # print('Number of leaf nodes in tree: ' + str(len(tree.leaf_nodes())))
                    # print('Number of internal nodes in tree: ' + str(len(tree.internal_nodes())))
                    print('Number of taxa in tree2: ' + str(len(tree.get_terminals())))
                    print('Number of leaf nodes in tree2: ' + str(tree.count_terminals()))
                    print('Number of internal nodes in tree2: ' + str(len(tree.get_nonterminals())))
                    # print(tree.as_ascii_plot(plot_metric='length'))
                    Phylo.draw_ascii(tree)

                    # Removing _dup taxa
                    logging.info(f"Removing duplicated taxa from input tree")
                    # tree_isolates = get_node_labels(tree.leaf_nodes())
                    tree_isolates = get_node_labels(tree.get_terminals())
                    for host_isolate in host_isolates:
                        if host_isolate not in tree_isolates:
                            logging.error(f'Host isolate {host_isolate} not in input phylogeny {args.input_phylogeny}!')
                            sys.exit(-1)
                        else:
                            dup_id = host_isolate + '_dup'
                            # tree.prune_taxa_with_labels([dup_id])
                            tree.prune(dup_id)

                    # Rooting the tree
                    # https://dendropy.org/primer/treemanips.html
                    logging.info(f"Rooting the tree...")
                    tree.root_with_outgroup(outgroup_isolate)
                    Phylo.draw_ascii(tree)
                    # root_node = tree.find_node_with_taxon_label(outgroup_isolate)
                    # tree.to_outgroup_position(root_node, update_bipartitions=False)
                    # tree.reroot_at_edge(root_node.edge, update_bipartitions=False)
                    # print(tree.as_ascii_plot(plot_metric='length'))

                    # Printing basic information on tree after removing duplicated taxa and re-rooting
                    logging.info(f"Printing basic information on input tree after removing duplicated taxa")
                    print('Number of taxa in tree2: ' + str(len(tree.get_terminals())))
                    print('Number of leaf nodes in tree2: ' + str(tree.count_terminals()))
                    print('Number of internal nodes in tree2: ' + str(len(tree.get_nonterminals())))

                    # Writing the tree
                    Phylo.write(tree, final_tree_file+'_rr', 'newick', plain=True)

                    # Assigning unique Id labels to internal nodes
                    tree = Phylo.read(final_tree_file+'_rr', "newick")
                    for idx, node in enumerate(tree.get_nonterminals()):
                        node.name = 'node' + str(idx)
                    # for idx, node in enumerate(tree.internal_nodes()):
                    #     node.label = 'node' + str(idx)

                    # Writing tree with internal nodes annotated
                    Phylo.write(tree, final_tree_file + '_nodes', 'newick', plain=True)
                    # Note: suppress_rooting=True so that the rooting token [&R] is not written
                    # Note: unquoted_underscores=True so that taxa labels are not quoted. If quoted pastml fails.
                    # tree.write(path=final_tree_file + '_nodes', schema="newick", suppress_rooting=True,
                    #            unquoted_underscores=True)

                    # Loading and saving Snippy VCF variants
                    #  VCF variants will be stored in the following variables
                    #    vcf_variants: dict to save what isolates contain what VCF variants/records
                    #         vcf_variants[host_isolate][var_id] = record
                    #    vcf_variants_all: dict to save all VCF variants/records across all isolates
                    #         vcf_variants_all[var_id] = record
                    #    vcf_variants_all_alt: dict to save all ALT to account for multi-allelic variants
                    vcf_variants = dict()
                    vcf_variants_all = dict()
                    vcf_variants_all_alt = dict()
                    for host_isolate in host_isolates:
                        vcf_variants[host_isolate] = {}
                        pair_id = host_isolate + '.' + reference_isolate + '.' + host_id
                        snippy_vcf_file = args.vcf_dir + pair_id + '.snippy.no_ref.vcf'
                        check_file_exists(snippy_vcf_file)
                        vcf_reader = get_vcf_reader(snippy_vcf_file)
                        for record in vcf_reader:
                            var_id = str(record.CHROM) + '.' + str(record.POS)
                            # Saving what isolates contain what variants
                            vcf_variants[host_isolate][var_id] = record
                            vcf_variants_all[var_id] = record
                            if vcf_variants_all_alt.get(var_id) is None:
                                vcf_variants_all_alt[var_id] = {}
                                vcf_variants_all_alt[var_id][str(record.ALT[0])] = record
                            else:
                                vcf_variants_all_alt[var_id][str(record.ALT[0])] = record
                            # print(var_id)

                    # Copying reference genome locally
                    reference_file_local = args.vcf_dir + reference_isolate + '.fasta'
                    if not os.path.exists(reference_file_local):
                        run_command_shell(
                            ["cp",
                             isolate_assemblies[reference_isolate],
                             reference_file_local],
                        )

                    # Creating PastML annotation table, if there are variants in the VCF file
                    # https://pastml.pasteur.fr/help
                    vcf_variants_ids_sorted = sorted(vcf_variants_all.keys())
                    if len(vcf_variants_ids_sorted) > 0:
                        # vcf_variants_ids: dict to save pastml_var_id (simplified id) as keys and var_id as
                        vcf_variants_ids = dict()

                        pastml_table_file = args.vcf_dir + host_id + '.snippy_vcf_variants.pastml_table.csv'
                        pastml_output_table = open(pastml_table_file, 'w')

                        pastml_table_line = 'ID'
                        for var_idx, var_id in enumerate(vcf_variants_ids_sorted):
                            pastml_var_id = 'var' + str(var_idx+1)
                            vcf_variants_ids[pastml_var_id] = var_id
                            pastml_table_line += '\t' + pastml_var_id
                            # print(var_id + ' > ' + pastml_var_id)
                        # print(pastml_table_line)
                        pastml_output_table.write(pastml_table_line + '\n')

                        for host_isolate in host_isolates:
                            pastml_table_line = host_isolate
                            for var_id in vcf_variants_ids_sorted:
                                if vcf_variants.get(host_isolate).get(var_id) is None:
                                    pastml_table_line += '\t' + str(vcf_variants_all[var_id].REF)
                                else:
                                    pastml_table_line += '\t' + str(vcf_variants[host_isolate][var_id].ALT[0])
                            # print(pastml_table_line)
                            pastml_output_table.write(pastml_table_line + '\n')
                        pastml_table_line = reference_isolate
                        for var_id in vcf_variants_ids_sorted:
                            pastml_table_line += '\t' + str(vcf_variants_all[var_id].REF)
                        # print(pastml_table_line)
                        pastml_output_table.write(pastml_table_line + '\n')
                        pastml_output_table.close()

                        # Running PastML (Ancestral reconstruction)
                        pastml_dir = host_id + '.snippy_vcf_variants.pastml'
                        logging.info(f"Running PastML (Ancestral reconstruction)...")
                        run_command_shell(
                            ["pastml",
                             "--tree", final_tree_file+'_nodes',
                             "--data", pastml_table_file,
                             "--work_dir", pastml_dir,
                             "--prediction_method", "ACCTRAN"],
                        )

                        # Identifying the MRCA of all isolates from host
                        logging.info(f"Identifying the MRCA of all host isolates")
                        tree = Phylo.read(final_tree_file+'_nodes', "newick")
                        host_mrca_node = tree.common_ancestor(host_isolates)
                        # tree = dendropy.Tree.get_from_path(final_tree_file+'_nodes', 'newick', preserve_underscores=True)
                        # host_mrca_node = tree.mrca(taxon_labels=host_isolates)
                        # print(host_mrca_node.label)
                        print(host_mrca_node.name)

                        # Reading PastML output file and saving MRCA alleles
                        vcf_variants_all_mrca = vcf_variants_all
                        pastml_output_file = pastml_dir + '/combined_ancestral_states.tab'
                        logging.info(f"Reading PastML output file {pastml_output_file}")
                        pastml_line_tracker = dict()
                        host_mrca_node_pastml_line_found = False
                        with open(pastml_output_file, 'r') as input_file:
                            header = input_file.readline().strip().split('\t')
                            for input_line in input_file:
                                line_items = input_line.strip().split('\t')
                                node_id = line_items[0]
                                # NOTE: PastML may output multiple lines (presumably multiple ancestral character
                                # reconstruction (ACR) solutions) per node. Only the first occurrence/line per node is kept.
                                if pastml_line_tracker.get(node_id) is None:
                                    # Only lines with expected number of items/characters are kept
                                    expected_line_items = len(vcf_variants_ids_sorted) + 1
                                    if len(line_items) == expected_line_items:
                                        # Only the alleles of host_mrca_node are kept
                                        if node_id == host_mrca_node.name:
                                            host_mrca_node_pastml_line_found = True
                                            for idx, allele in enumerate(line_items):
                                                if idx > 0:
                                                    var_id = vcf_variants_ids[header[idx]]
                                                    # printing info
                                                    print(node_id + ' ' + header[idx] + ' ' + allele + ' ' + var_id)
                                                    # replacing VCF ALT allele with MRCA ancestral allele
                                                    vcf_variants_all_mrca[var_id].ALT[0] = str(allele)
                                                    if vcf_variants_all_alt[var_id].get(str(allele)) is not None:
                                                        vcf_variants_all_mrca[var_id] = vcf_variants_all_alt[var_id][str(allele)]
                                pastml_line_tracker[node_id] = ''
                        if not host_mrca_node_pastml_line_found:
                            logging.error(f'Sequence of host MRCA node {host_mrca_node.name} could not be extracted'
                                          f'from PastML output file {pastml_output_file}. Check!')
                            sys.exit(-1)

                        # Re-constructing host isolates MRCA sequence
                        # The outgroup/reference genome used is read and their alleles replaced with those of MRCA

                        mutated_records = []

                        for record in SeqIO.parse(reference_file_local, "fasta"):
                            print(record.id)
                            mutated_record = record
                            mutated_record.seq = mutated_record.seq.tomutable()
                            # For each contig, keep VCF variants from that contig only, as a list of VCF records
                            vcf_variants_all_mrca_contig = []
                            for var_id in vcf_variants_all_mrca.keys():
                                if vcf_variants_all_mrca[var_id].CHROM == record.id:
                                    vcf_variants_all_mrca_contig.append(vcf_variants_all_mrca[var_id])

                            # If VCF variant found in this contig, then modify mutated contig sequence
                            # NOTE: mutations need to be introduced in a decreasing order to make sure original contig
                            # coordinates are not modified when introducing indels
                            if len(vcf_variants_all_mrca_contig) > 0:
                                for a in range(len(vcf_variants_all_mrca_contig) - 1, -1, -1):
                                    vcf_record = vcf_variants_all_mrca_contig[a]
                                    mutation_type = str(vcf_record.INFO["TYPE"][0])
                                    mutation_alt = str(vcf_record.ALT[0])
                                    mutation_pos = vcf_record.POS
                                    mutation_ref = str(vcf_record.REF)
                                    # Make sure reference allele does not contain Ns
                                    # if 'N' not in mutation_ref:
                                    # Make sure reference and alternative alleles are different
                                    if mutation_ref != mutation_alt:
                                        print('chr ' + vcf_record.CHROM + ' ' + str(mutation_pos) + ' ' + mutation_type)
                                        mutated_record.seq = create_mutation(mutated_record.seq, mutation_pos,
                                                                                 mutation_ref, mutation_alt, mutation_type)

                            # Save mutated contig sequence
                            mutated_records.append(mutated_record)

                        logging.info(f"Saving MRCA FASTA file (contigs) into {fasta_file}")
                        SeqIO.write(mutated_records, fasta_file, "fasta")

                        # Remove temporary files
                        logging.info(f"Removing temporary files ")
                        tmp_files = [final_tree_file+'_rr', final_tree_file+'_nodes', pastml_table_file]
                        for tmp_file in tmp_files:
                            logging.info(f"Removing temporary file {tmp_file}")
                            run_command_shell(["rm -f", tmp_file])
                        tmp_dirs = [pastml_dir]
                        for tmp_dir in tmp_dirs:
                            logging.info(f"Removing temporary directory {tmp_dir}")
                            run_command_shell(["rm -r", tmp_dir])
                    else:
                        # if there are not variants in the VCF file, the outgroup/reference genome is not modified
                        # Copying reference file
                        run_command_shell(
                            ["cp",
                             reference_file_local,
                             fasta_file],
                        )
                        # Remove temporary files
                        logging.info(f"Removing temporary files ")
                        tmp_files = [final_tree_file + '_rr', final_tree_file + '_nodes']
                        for tmp_file in tmp_files:
                            logging.info(f"Removing temporary file {tmp_file}")
                            run_command_shell(["rm -f", tmp_file])


if __name__ == "__main__":
    _main()