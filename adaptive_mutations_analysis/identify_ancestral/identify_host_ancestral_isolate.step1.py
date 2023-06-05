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
    description = "Script to identify an appropriate out-group for the multiple clonal isolates of the same host"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-t", "--input_phylogeny", action="store", dest="input_phylogeny",
        help="Input rooted phylogenetic tree",
        required=True, metavar="INPUT_TREE"
    )
    group.add_argument(
        "-d", "--input_distances", action="store", dest="input_distances",
        help="Matrix of pairwise SNP distances produced by pairsnp",
        required=True, metavar="INPUT_DISTANCES"
    )
    group.add_argument(
        "-m", "--input_metadata", action="store", dest="input_metadata",
        help="input csv metadata file with isolate and patient Ids",
        required=True, metavar="INPUT_METADATA"
    )
    group.add_argument(
        "-s", "--min_snp_distance", action="store", dest="min_snp_distance",
        help="Minimum SNP distance to outgroup (e.g. 5 SNPs)",
        required=True, metavar="MIN_DIS"
    )
    group.add_argument(
        "-o", "--output_table", action="store", dest="output_table",
        help="output table with outgroup isolate per host",
        required=True, metavar="OUTPUT_TABLE"
    )

    return parser.parse_args()


def is_integer(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def is_tree_valid(input_tree):
    try:
        Phylo.read(input_tree, 'newick')
        tree = dendropy.Tree.get_from_path(input_tree, 'newick')
    except:
        print("Error with the input starting tree: Is it a valid Newick file?")
        return 0
    return 1


def get_node_labels(node_list):
    """
    Given a list of Node object, returns their taxa labels
    :param node_list: list of Node objects
    :return: list of taxa labels
    """
    node_labels = list()
    for i, node in enumerate(node_list):
        node_labels.append(node_list[i].taxon.label)
    return node_labels


def load_snp_distances(input_distances):
    """
    This function saves pairwise SNP distances in the form of a two dimensional dictionary
    where keys are isolate Ids and values SNP distances
    :param input_distances: pairsnp output file
    :return: two-key dictionary of SNP distances
    """
    distances_dict = dict()
    # First, get isolate Ids from first column
    isolate_ids = []
    with open(input_distances, 'r') as distances_file:
        for distances_line in distances_file:
            distances = distances_line.strip().split(',')
            isolate_ids.append(distances[0])
    # Second, save SNP distances as a two-key dictionary
    with open(input_distances, 'r') as distances_file:
        for isolate_a_idx, distances_line in enumerate(distances_file):
            isolate_a = isolate_ids[isolate_a_idx]
            # print(str(isolate_a_idx) + ' ' + isolate_a + ' ' + distances_line)
            distances = distances_line.strip().split(',')
            distances.pop(0)
            for isolate_b_idx, distance in enumerate(distances):
                isolate_b = isolate_ids[isolate_b_idx]
                if distances_dict.get(isolate_a) is None:
                    distances_dict[isolate_a] = {}
                    distances_dict[isolate_a][isolate_b] = int(distance)
                else:
                    distances_dict[isolate_a][isolate_b] = int(distance)
                # print('distances_dict[' +  isolate_a + '][' + isolate_b + ']' + distances_dict[isolate_a][isolate_b])
                if distances_dict.get(isolate_b) is None:
                    distances_dict[isolate_b] = {}
                    distances_dict[isolate_b][isolate_a] = int(distance)
                else:
                    distances_dict[isolate_b][isolate_a] = int(distance)
    return distances_dict


def get_all_snp_distances(isolates_list, distances_dict):
    """
    This function return a list of all pairwise SNP distances for a list of isolate Ids
    :param isolates_list: list of isolate Ids
    :param distances_dict: two-key dictionary of SNP distances
    :return: list of SNP distances
    """
    isolates_snp_distances = []
    for isolate_a_idx, isolate_a in enumerate(isolates_list):
        for isolate_b_idx, isolate_b in enumerate(isolates_list):
            if isolate_b_idx > isolate_a_idx:
                if distances_dict.get(isolate_a) is None:
                    sys.exit("Error: No SNP distances found for " + isolate_a)
                if distances_dict.get(isolate_b) is None:
                    sys.exit("Error: No SNP distances found for " + isolate_b)
                if distances_dict.get(isolate_a).get(isolate_b) is None:
                    sys.exit("Error: No SNP distance found for pair: " + isolate_a + " and " + isolate_b)
                isolates_snp_distances.append(distances_dict[isolate_a][isolate_b])
    return isolates_snp_distances


def get_min_snp_distances(source_list, target_list, distances_dict):
    """
    This function is used to calculate the maximum SNP distance of isolates in the target list
    against those in the source list
    Given two lists of isolates
    :param source_list: isolate list used for comparison
    :param target_list: isolate list of interest
    :param distances_dict: two-key dictionary of SNP distances
    :return: a list of same length/order as target list
    """
    isolates_snp_distances = []
    for isolate_t_idx, isolate_t in enumerate(target_list):
        tmp_distances = list()
        for isolate_s_idx, isolate_s in enumerate(source_list):
                if distances_dict.get(isolate_t) is None:
                    sys.exit("Error: No SNP distances found for " + isolate_t)
                if distances_dict.get(isolate_s) is None:
                    sys.exit("Error: No SNP distances found for " + isolate_s)
                if distances_dict.get(isolate_t).get(isolate_s) is None:
                    sys.exit("Error: No SNP distance found for pair: " + isolate_t + " and " + isolate_s)
                tmp_distances.append(distances_dict[isolate_t][isolate_s])
        isolates_snp_distances.append(min(tmp_distances))
    return isolates_snp_distances


def substract_a_from_b(list_a, list_b):
    """
    This function is used to substract items in list A from list B and return the items
    remaining in list B
    :param list_a: subset of list B
    :param list_b: longer list containing items of list A
    :return: items remaining in list B
    """
    remaining_items = list(set(list_b).difference(list_a))
    return remaining_items


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
    input_files = [args.input_phylogeny, args.input_distances, args.input_metadata]
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Making sure min_snp_distance is an integer digit
    if not is_integer(args.min_snp_distance):
        logging.error(f'min_snp_distance {args.min_snp_distance} must be an integer!')
        sys.exit(-1)
    else:
        if int(args.min_snp_distance) <= 0:
            logging.error(f'min_snp_distance {args.min_snp_distance} must be greater than 0!')
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

    # Loading pairwise SNP distances
    logging.info(f"Loading SNP distances from input file {args.input_distances}")
    isolate_distances = load_snp_distances(args.input_distances)

    # Loading phylogenetic tree
    if is_tree_valid(args.input_phylogeny) == 0:
        logging.error(f'Input phylogenetic tree {args.input_phylogeny} is invalid!')
        sys.exit(-1)
    tree = dendropy.Tree.get_from_path(args.input_phylogeny, 'newick', preserve_underscores=True)

    # Printing basic information on tree
    print('Number of taxa in tree: ' + str(len(tree.taxon_namespace)))
    print('Number of leaf nodes in tree: ' + str(len(tree.leaf_nodes())))
    print('Number of internal nodes in tree: ' + str(len(tree.internal_nodes())))

    # Assigning unique Id labels to internal nodes
    for idx, node in enumerate(tree):
        node.label = 'node' + str(idx)

    # Write tree with numbered nodes
    tree.write(path="output.tre", schema="newick")

    # Writing into output table file
    output_table_file = open(args.output_table, 'w')
    # items_to_save = ['host_id', 'host_isolates', 'min_host_isolates_dis', 'max_host_isolates_dis', 'host_mrca_node',
    #                  'outgroup_isolate', 'outgroup_snp_dis', 'reference_isolate', 'reference_snp_dis']
    items_to_save = ['host_id', 'host_isolates', 'min_host_isolates_dis', 'max_host_isolates_dis', 'host_mrca_node',
                     'outgroup_isolate', 'outgroup_snp_dis']
    output_table_file.write('\t'.join(items_to_save) + '\n')

    # For each host: extract outgroup isolates
    for host_id in host_isolate_ids.keys():
        # if host_id == 'paterson2015-Staff_A':
        # if host_id == 'young2017-P071':
            print(host_id)
            items_to_save = list()
            items_to_save.append(host_id)
            # Get all isolates for host
            host_isolates = list(host_isolate_ids.get(host_id).keys())
            print('host_isolates ' + str(host_isolates))
            items_to_save.append(';'.join(host_isolates))
            # Get SNP distances for host_isolates
            host_isolates_dis = get_all_snp_distances(host_isolates, isolate_distances)
            items_to_save.append(str(min(host_isolates_dis)))
            items_to_save.append(str(max(host_isolates_dis)))

            # Make sure host isolates exist in tree
            tree_isolates = get_node_labels(tree.leaf_nodes())
            for host_isolate in host_isolates:
                if host_isolate not in tree_isolates:
                    logging.error(f'Host isolate {host_isolate} not in input phylogeny {args.input_phylogeny}!')
                    sys.exit(-1)

            # Get MRCA of all isolates for host
            host_mrca_node = tree.mrca(taxon_labels=host_isolates)
            items_to_save.append(str(host_mrca_node.label))

            # Get outgroup isolate
            # Note: the while look is design to find the ancestor node from host_mrca with an isolate having a
            # genetic distance greater than the specified by the user. This is to avoid choosing outgroup isolates
            # that are genetically identical to host isolates
            snp_distance_to_outgroup = 0
            while int(snp_distance_to_outgroup) <= int(args.min_snp_distance):
                host_outgroup_node = host_mrca_node.parent_node
                host_outgroup_isolates = substract_a_from_b(get_node_labels(host_mrca_node.leaf_nodes()),
                                                            get_node_labels(host_outgroup_node.leaf_nodes()))
                host_outgroup_snp_dis = get_min_snp_distances(host_isolates, host_outgroup_isolates, isolate_distances)
                snp_distance_to_outgroup = min(host_outgroup_snp_dis)
                host_mrca_node = host_outgroup_node

            host_outgroup_idx = host_outgroup_snp_dis.index(min(host_outgroup_snp_dis))
            host_outgroup_isolate = host_outgroup_isolates[host_outgroup_idx]
            items_to_save.append(str(host_outgroup_isolate))
            items_to_save.append(str(min(host_outgroup_snp_dis)))
            print('host_outgroup_isolates ' + str(host_outgroup_isolates))
            print('host_outgroup_snp_dis ' + str(host_outgroup_snp_dis))
            print('host_outgroup_idx ' + str(host_outgroup_idx))
            print('host_outgroup_isolate ' + str(host_outgroup_isolate))

            # Get reference isolate
            # host_reference_node = host_outgroup_node.parent_node
            # if host_reference_node is None:
            #     items_to_save.append(str('not_found'))
            #     items_to_save.append(str('non_applicable'))
            # else:
            #     host_reference_isolates = substract_a_from_b(get_node_labels(host_outgroup_node.leaf_nodes()),
            #                                                  get_node_labels(host_reference_node.leaf_nodes()))
            #     host_reference_snp_dis = get_min_snp_distances(host_isolates, host_reference_isolates, isolate_distances)
            #     host_reference_idx = host_reference_snp_dis.index(min(host_reference_snp_dis))
            #     host_reference_isolate = host_reference_isolates[host_reference_idx]
            #     print('host_reference_isolates ' + str(host_reference_isolates))
            #     print('host_reference_snp_dis ' + str(host_reference_snp_dis))
            #     print('host_reference_idx ' + str(host_reference_idx))
            #     print('host_reference_isolate ' + str(host_reference_isolate))
            #     items_to_save.append(str(host_reference_isolate))
            #     items_to_save.append(str(min(host_reference_snp_dis)))

            # Saving information
            output_table_file.write('\t'.join(items_to_save) + '\n')


if __name__ == "__main__":
    _main()