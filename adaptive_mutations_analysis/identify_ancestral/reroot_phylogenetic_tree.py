#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import Phylo
import dendropy


# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------


def parse_arguments():
    description = "Script to re-root a phylogenetic tree in Newick format"
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "-i", "--input_tree", action="store", dest="input_tree",
        help="Input phylogenetic tree file in Newick format",
        required=True, metavar="INPUT_TREE")
    parser.add_argument(
        "-r", "--out_group", action="store", dest="out_group",
        help="Out-group used to re-root the input tree",
        required=True, metavar="OUTGROUP")
    parser.add_argument(
        "-o", "--output_tree", action="store", dest="output_tree",
        help="Output re-rooted phylogenetic tree file in Newick",
        required=True, metavar="OUTPUT_TREE")

    return parser.parse_args()


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
    input_files = [args.input_tree]
    for input_file in input_files:
        if not os.path.exists(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Loading phylogenetic tree
    logging.info(f"Opening input tree file {args.input_tree}")
    if is_tree_valid(args.input_tree) == 0:
        logging.error(f'Input phylogenetic tree {args.input_tree} is invalid!')
        sys.exit(-1)
    tree = dendropy.Tree.get_from_path(args.input_tree, 'newick', preserve_underscores=True)

    # Printing basic information on tree before removing taxa
    logging.info(f"Input tree information (before re-rooting)")
    print('Number of taxa in tree: ' + str(len(tree.taxon_namespace)))
    print('Number of leaf nodes in tree: ' + str(len(tree.leaf_nodes())))
    print('Number of internal nodes in tree: ' + str(len(tree.internal_nodes())))

    # Extracting all taxa ids in input tree
    tree_taxa = get_node_labels(tree.leaf_nodes())

    # Removing taxa
    if args.out_group in tree_taxa:
        outgroup_node = tree.find_node_with_taxon_label(args.out_group)
        tree.to_outgroup_position(outgroup_node, update_bipartitions=False)
        logging.info(f'Input tree rooted on {args.out_group}')
    else:
        logging.error(f'Taxon {args.out_group} chosen as out-group not found in tree')
        sys.exit(-1)

    # Printing basic information on tree after removing taxa
    logging.info(f"Output tree information (after re-rooting)")
    print('Number of taxa in tree: ' + str(len(tree.taxon_namespace)))
    print('Number of leaf nodes in tree: ' + str(len(tree.leaf_nodes())))
    print('Number of internal nodes in tree: ' + str(len(tree.internal_nodes())))

    # Writing output file
    logging.info(f"Saving output tree file {args.output_tree}")
    tree.write(path=args.output_tree, schema="newick")


if __name__ == "__main__":
    _main()
