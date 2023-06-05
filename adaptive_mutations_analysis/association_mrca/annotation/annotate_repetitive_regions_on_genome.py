#!/usr/bin/env python3

import argparse
import logging
import os
import subprocess
import sys
from pandas import DataFrame
from Bio import SeqIO
from Bio import SeqFeature as sf


# ------------------------------------------------------------------------------------
# Notes
# ------------------------------------------------------------------------------------

# This identified repetitive regions in a reference genome, defined by BLASTing the
# reference genome against itself.

# Notes on software versions used:
# makeblastdb: 2.8.1+
# blastn: 2.8.1+

# ------------------------------------------------------------------------------------
# Global variables
# ------------------------------------------------------------------------------------

_DEPENDENCIES = ['makeblastdb', 'blastn']

# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------


def parse_arguments():
    description = "Script to identify repetitive regions in a bacterial reference genome"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-g", "--input_genome", action="store", dest="input_genome",
        help="input genome in FASTA format",
        required=True, metavar="GENOME"
    )
    group.add_argument(
        "-p", "--prefix", action="store", dest="prefix",
        help="prefix used to name temporary files",
        required=True, metavar="PREFIX")
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s 1.0")
    group = parser.add_argument_group('Blastn arguments (optional)')
    group.add_argument(
        "-l", "--length", action="store", dest="length",
        help="minimum length (in base pairs) of blast hit to be included (default: %(default)s)",
        required=False, default="100", metavar="LENGTH")
    group.add_argument(
        "-i", "--per_identity", action="store", dest="per_identity",
        help="minimum percentage identity of blast hit to be included (default: %(default)s)",
        required=False, default="70", metavar="PERC")

    return parser.parse_args()


def check_dependency(executable_name):
    """ Returns true if executable exists, else false """
    found = False
    output = subprocess.check_output(['which', executable_name]).strip()
    if output:
        found = True
    return found


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


def get_feature_a_overlap_with_feature_b(feature_a, feature_b):
    """
    Given two Bio.SeqFeature objects, return overlap of feature_a with feature_b
    :param feature_a: Bio.SeqFeature object
    :param feature_b: Bio.SeqFeature object
    :return: percentage overlap of feature_a with feature_b
    """
    range_a = set(range(int(feature_a.location.start), int(feature_a.location.end) + 1, 1))
    range_b = set(range(int(feature_b.location.start), int(feature_b.location.end) + 1, 1))
    overlap_a = range_a.intersection(range_b)
    per_ovelap = (len(overlap_a)/len(range_a))*100
    return per_ovelap


def merge_features(feature_a, feature_b, sseqid, type):
    """
    Given two overlapping Bio.SeqFeature objects, returns the union
    The qualifiers of the longest Bio.SeqFeature object are kept
    :param feature_a: Bio.SeqFeature object
    :param feature_b: Bio.SeqFeature object
    :param sseqid: sequence id where Bio.SeqFeature lays
    :param type: type of Bio.SeqFeature
    :return: merged Bio.SeqFeature
    """
    range_a = set(range(int(feature_a.location.start), int(feature_a.location.end) + 1, 1))
    range_b = set(range(int(feature_b.location.start), int(feature_b.location.end) + 1, 1))
    range_union = sorted(list(range_a.union(range_b)))
    new_location = sf.FeatureLocation(
        int(range_union[0]), int(range_union[-1]), ref=sseqid, strand=feature_b.location.strand,
    )
    new_feature = sf.SeqFeature(new_location, type=type)
    new_feature.qualifiers = feature_a.qualifiers if len(range_a) >= len(range_b) else feature_b.qualifiers
    new_feature.qualifiers["locus_tag"] = sseqid + '_' + str(new_location.start) + '_' + str(new_location.end)

    return new_feature


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
    if not os.path.isfile(args.input_genome):
        logging.error(f'Input genome {args.input_genome} not found!')
        sys.exit(-1)

    # Loading input genome and genome length
    logging.info(f"Opening input genome {args.input_genome}")
    input_records = SeqIO.parse(args.input_genome, "fasta")
    chr_length = 0
    chr_id = ''
    for record in input_records:
        chr_length = len(record)
        chr_id = record.id
        print(record.id + ' ' + str(chr_length))

    # Creating input genome blast database
    blast_database = args.input_genome + ".nin"
    if os.path.isfile(blast_database):
        logging.info(f'Index for {args.input_genome} already built...')
    else:
        logging.info(f'Building blast index for {args.input_genome}...')
        run_command_shell(
            ["makeblastdb",
             "-in", args.input_genome,
             "-dbtype", "nucl"],
        )

    # Running blastn
    logging.info(f"Running blastn... this may take a few minutes")
    blastn_output = args.prefix + ".blastn_output.txt"
    if os.path.isfile(blastn_output):
        logging.info(f'File {blastn_output} already found')
    else:
        run_command_shell(
            ["blastn",
             "-query", args.input_genome,
             "-db", args.input_genome,
             "-out", blastn_output,
             "-task", "blastn",
             "-sorthsps", "4",
             "-outfmt", "\"", "6 qseqid sseqid length qlen sstart send pident sstrand qcovs", "\"",
             ],
        )
        if not os.path.exists(blastn_output):
            logging.error(f'{blastn_output} not created after running Blast. Something went wrong!')
            sys.exit(-1)

    # Saving blast hits as a dataframe so that blast hits can be sorted by position and length
    logging.info(f"Saving blast hits as a dataframe so that blast hits can be sorted")
    # Lists used as columns to create dataframe blast_hits_df
    qseqid_list = list()
    length_list = list()
    sstart_list = list()
    send_list = list()
    pident_list = list()
    sstrand_list = list()

    for blast_line in open(blastn_output, "r"):
        if not blast_line.strip() == '':
            # NCTC8325	NCTC8325	5405	2821361	448714	454116	99.593	plus	100
            (qseqid, sseqid, length, qlen, sstart, send, pident, sstrand, qcovs) = blast_line.split("\t")
            # if length < chr_length: is used to remove first line
            if int(length) < int(chr_length):
                # making sure blast hit meets length and percentage identity requirements
                if float(pident) >= float(args.per_identity):
                    if int(length) >= int(args.length):
                        qseqid_list.append(qseqid)
                        length_list.append(int(length))
                        sstart_list.append(int(sstart))
                        send_list.append(int(send))
                        pident_list.append(pident)
                        sstrand_list.append(sstrand)

    blast_hits = {'qseqid': qseqid_list,
                  'length': length_list,
                  'sstart': sstart_list,
                  'send': send_list,
                  'pident': pident_list,
                  'sstrand': sstrand_list,
                  }

    blast_hits_df = DataFrame(blast_hits, columns=['qseqid', 'length', 'sstart', 'send', 'pident', 'sstrand'])

    blast_hits_df.sort_values(by=['sstart','length'], inplace=True, ascending=True)

    # Saving blast hits as repetitive regions
    # Creating repeat features (SeqFeature) from BLAST output
    logging.info(f"Creating repeat features (SeqFeature) from sorted blast hits")
    repeat_features = list()

    for index, row in blast_hits_df.iterrows():
        (qseqid, length, sstart, send, pident, sstrand) = (row['qseqid'], row['length'], row['sstart'], row['send'],
                                                           row['pident'], row['sstrand'])
        strand = +1 if sstrand == "plus" else -1
        start = sstart if sstrand == "plus" else send
        end = send if sstrand == "plus" else sstart
        # create FeatureLocation object, which will be used to create the repeat SeqFeature object
        repeat_location = sf.FeatureLocation(
            int(start), int(end), ref=sseqid, strand=strand,
        )
        repeat_locus_tag = qseqid + '_' + str(sstart) + '_' + str(send)
        repeat_feature = sf.SeqFeature(repeat_location, type="repeat")
        repeat_feature.qualifiers["locus_tag"] = repeat_locus_tag
        repeat_feature.qualifiers["inference"] = "blastn"
        repeat_features.append(repeat_feature)
        print(repeat_feature)

    print("A total of " + str(len(repeat_features)) + " blast hits are saved")

    # Removing repeats contained within other repeats
    logging.info(f"Removing repeats contained within other repeats. This may take a few minuts...")
    # Note: the list repeat_features needs to be iterated backwards for the list.remove to work
    for a in range(len(repeat_features) - 1, -1, -1):
        repeat_feature_a = repeat_features[a]
        for b in range(len(repeat_features) - 1, -1, -1):
            repeat_feature_b = repeat_features[b]
            if repeat_feature_a != repeat_feature_b:
                if get_feature_a_overlap_with_feature_b(repeat_feature_a, repeat_feature_b) >= 100:
                    repeat_features.remove(repeat_feature_a)
                    break

    # Merging overlapping repeats. If A overlaps with B, delete A and B and create a merged A+B
    logging.info(f"Merging overlapping repeats. This may take a few minuts...")
    for a in range(len(repeat_features) - 1, -1, -1):
        repeat_feature_a = repeat_features[a]
        for b in range(len(repeat_features) - 1, -1, -1):
            repeat_feature_b = repeat_features[b]
            if repeat_feature_a != repeat_feature_b:
                if get_feature_a_overlap_with_feature_b(repeat_feature_a, repeat_feature_b) > 1:
                    merged_repeat_feature = merge_features(repeat_feature_a, repeat_feature_b, chr_id, "repeat")
                    repeat_features.append(merged_repeat_feature)
                    repeat_features.remove(repeat_feature_a)
                    repeat_features.remove(repeat_feature_b)
                    break

    print("A total of " + str(len(repeat_features)) + " repeats are left after merging overlapping repeats")

    # Saving repeats
    # Note: the function SeqIO.write consumes the input_records iterator. SeqIO.parse needs to be called again
    input_records = SeqIO.parse(args.input_genome, "fasta")
    output_records = []
    for record in input_records:
        for repeat_feature in repeat_features:
            # repeat_feature.annotations['molecule_type'] = record.annotations['molecule_type']
            record.features.append(repeat_feature)
    record.annotations['molecule_type'] = "DNA"
    output_records.append(record)
    final_output = args.prefix + ".repetitive_regions.gbk"
    SeqIO.write(output_records, final_output, "genbank")


if __name__ == "__main__":
    _main()


