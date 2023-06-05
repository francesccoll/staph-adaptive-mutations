#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
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

# NOTES on software and modules versions. Tested on
#   - Python 3.6.0 and [GCC 4.6.3] on linux
#   - Phylo Version: 4.2.1
#   - dendropy 4.2.0
#   - bcftools 1.9. Change bcftools executable if necessary
#   - RAxML version 8.2.8
#   - snp-sites 2.4.1

# NOTES on the use of RAxML
#   - The single-threaded (Serial) version of RAxML (raxmlHPC-SSE3) is used as small alignments are expected
#   - The "rapid Bootstrap analysis and search for best-scoring ML tree in one program run" is used ("-f a")

# ------------------------------------------------------------------------------------
# Global variables
# ------------------------------------------------------------------------------------

_DEPENDENCIES = ['bcftools-1.9', 'snippy-vcf_extract_subs', 'snp-sites', 'raxmlHPC-SSE3', 'vcfintersect']
# _DEPENDENCIES = ['bcftools', 'snippy-vcf_extract_subs', 'snp-sites']


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


def run_command(command_line):
    """
    This function executes a command line, check for execution errors and returns stdout
    :param command_line: it must be a list not a string
    :return: stdout
    """
    try:
        process_completed = subprocess.run(
            command_line,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
    except subprocess.CalledProcessError as err:
        print('ERROR:', err)
    return process_completed.stdout.decode('utf-8')


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


def concatenate_contigs(ref_file, contigs_file, fasta_id, fasta_file):
    """
    This function is used to read an assembly file and concatenate their contigs
    :param ref_file: reference used for mapping in FASTA format
    :param contigs_file: assembly file with contigs in FASTA format
    :param fasta_id: id of concatenated contigs file in FASTA format
    :param fasta_file: concatenated contigs file in FASTA format
    :return:
    """
    expected_length = 0
    for record in SeqIO.parse(ref_file, "fasta"):
        expected_length = expected_length + len(record)

    concatenated_seq = ""
    contigs = {}
    contig_order = []

    for record in SeqIO.parse(contigs_file, "fasta"):
        contig_order.append(record.id)
        contigs[record.id] = str(record.seq)

    for contig in contig_order:
        concatenated_seq = concatenated_seq + contigs[contig]

    concatenated_seq_record = SeqRecord(Seq(concatenated_seq))
    concatenated_seq_record.id = fasta_id
    print("Length of concatenated sequence: " + str(len(concatenated_seq_record)))
    if len(concatenated_seq_record)!= expected_length:
        logging.error(f'Length of concatenated sequence {str(len(concatenated_seq_record))} different from expected '
                      f'{str(expected_length)}!')
        sys.exit(-1)
    SeqIO.write(concatenated_seq_record, fasta_file, "fasta")


def duplicate_alignment(input_aln_file, taxa_to_duplicate, output_aln_file):
    """
    This function is used duplicate the sequences of selected isolates in a FASTA multiple alignment
    Duplicated sequences are given a new Id that ends with '_dup'
    :param input_aln_file: input multiple alignment
    :param taxa_to_duplicate: labels of taxa to duplicate
    :param output_aln_file: output multiple alignment
    :return: nothing
    """

    output_records = []

    for record in SeqIO.parse(input_aln_file, "fasta"):
        output_records.append(record)
        if record.id in taxa_to_duplicate:
            record_dup = SeqRecord(Seq(str(record.seq)), id=str(record.id) + '_dup', description='')
            output_records.append(record_dup)

    SeqIO.write(output_records, output_aln_file, "fasta")


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
                # if host_id == "harrison2016-CR02":
                host_isolates = line_items[1].split(';')
                outgroup_isolate = line_items[5]
                # reference_isolate = line_items[7]
                reference_isolate = outgroup_isolate
                print(host_id + ' ' + '-'.join(host_isolates) + ' ' + outgroup_isolate)

                # Make sure final phylogenetic tree is not created already, otherwise, continue
                final_tree_file = args.raxml_dir + "RAxML_bipartitions." + host_id + '.' + reference_isolate + \
                                  '.raxml_tree' + "_rooted"
                if os.path.isfile(final_tree_file):
                    logging.info(f"Final tree {final_tree_file} already found.")
                else:
                    # List to save all consensus sequences in fasta format
                    snippy_consensus_files = []

                    # Making sure expected Snippy VCF files exist
                    isolates_in_alignment = []
                    isolates_in_alignment.extend(host_isolates)
                    isolates_in_alignment.append(outgroup_isolate)
                    for isolate_in_alignment in isolates_in_alignment:
                        pair_id = isolate_in_alignment + '.' + reference_isolate + '.' + host_id
                        snippy_vcf_file = args.vcf_dir + pair_id + '.snippy.vcf'
                        if not os.path.exists(snippy_vcf_file):
                            logging.error(f'Expected Snippy VCF file {snippy_vcf_file} not found!')
                            sys.exit(-1)

                        # Copying reference genome locally
                        reference_file_local = args.vcf_dir + reference_isolate + '.fasta'
                        if not os.path.exists(reference_file_local):
                            run_command_shell(
                                ["cp",
                                 isolate_assemblies[reference_isolate],
                                 reference_file_local],
                            )

                        # Removing variants found when mapping reference against its own assembly
                        # This are putative false positive variants arising from missassemblies
                        snippy_vcf_file_no_ref = snippy_vcf_file.replace('.snippy.vcf', '.snippy.no_ref.vcf')
                        if os.path.isfile(snippy_vcf_file_no_ref):
                            logging.info(f"VCF file with no reference variants {snippy_vcf_file_no_ref} already found")
                        else:
                            logging.info(f"Running vcfintersect to remove reference variants...")
                            ref_snippy_vcf_file = args.vcf_dir + reference_isolate + '.' + reference_isolate + '.' +\
                                                  host_id + '.snippy.vcf'
                            run_command_shell(
                                ["vcfintersect",
                                 "-v",
                                 "-i", ref_snippy_vcf_file,
                                 "-r", reference_file_local,
                                 "-w", "1",
                                 snippy_vcf_file,
                                 ">",
                                 snippy_vcf_file_no_ref],
                            )
                            if not os.path.exists(snippy_vcf_file_no_ref):
                                logging.error(f'{snippy_vcf_file_no_ref} not created. Something went wrong!')
                                sys.exit(-1)
                            else:
                                print(snippy_vcf_file_no_ref)

                        # Keeping only substitutions in VCF file, indels are ignored
                        snippy_vcf_file_subs = snippy_vcf_file.replace('.snippy.vcf', '.snippy.subs.vcf')
                        if os.path.isfile(snippy_vcf_file_subs):
                            logging.info(f"snippy-vcf_extract_subs output {snippy_vcf_file_subs} already found")
                        else:
                            logging.info(f"Running snippy-vcf_extract_subs...")
                            run_command_shell(
                                ["snippy-vcf_extract_subs",
                                 snippy_vcf_file_no_ref,
                                 ">",
                                 snippy_vcf_file_subs],
                            )
                            if not os.path.exists(snippy_vcf_file_subs):
                                logging.error(f'{snippy_vcf_file_subs} not created. Something went wrong!')
                                sys.exit(-1)
                            else:
                                print(snippy_vcf_file_subs)

                        # Filtering out SNPs with N REF alleles
                        snippy_vcf_file_noNs = snippy_vcf_file.replace('.snippy.vcf', '.snippy.subs.noNs.vcf')
                        if os.path.isfile(snippy_vcf_file_noNs):
                            logging.info(f"Snippy VCF with no Ns {snippy_vcf_file_noNs} already found")
                        else:
                            logging.info(f"Filtering sites with REF missing alleles (Ns)...")
                            run_command_shell(
                                ["bcftools-1.9",
                                 "filter",
                                 "--exclude",
                                 "\'REF=\"N\"\'",
                                 snippy_vcf_file_subs,
                                 ">",
                                 snippy_vcf_file_noNs],
                            )
                            if not os.path.exists(snippy_vcf_file_noNs):
                                logging.error(f'{snippy_vcf_file_noNs} not created. Something went wrong!')
                                sys.exit(-1)
                            else:
                                print(snippy_vcf_file_noNs)

                        # Creating consensus/pseudo-molecule sequence
                        snippy_consensus_file = snippy_vcf_file.replace('.snippy.vcf', '.snippy.fasta')
                        if os.path.isfile(snippy_consensus_file):
                            logging.info(f"{snippy_consensus_file} already found")
                        else:
                            logging.info(f"Creating Snippy substitutions consensus file...")
                            if isolate_assemblies.get(reference_isolate) is None:
                                logging.error(f'Assembly file for isolate {reference_isolate} not found!')
                                sys.exit(-1)
                            ref = reference_file_local
                            vcf = snippy_vcf_file_noNs
                            consensus = snippy_consensus_file
                            run_command_shell(["bcftools-1.9", "convert", "-Oz", "-o", vcf + '.gz', vcf])
                            run_command_shell(["bcftools-1.9", "index", "-f", vcf + '.gz'])
                            run_command_shell(["bcftools-1.9", "consensus", "-f", ref, "-o", consensus, vcf + '.gz'])
                            run_command_shell(["rm -f", vcf + '.gz', vcf + '.gz.csi', vcf + '.tbi'])
                            if not os.path.exists(snippy_consensus_file):
                                logging.error(f'{snippy_consensus_file} not created. Something went wrong!')
                                sys.exit(-1)
                            else:
                                print(snippy_consensus_file)

                        # Concatenating contigs into a single sequence
                        snippy_consensus_concat_file = snippy_vcf_file.replace('.snippy.vcf', '.snippy.concat.fasta')
                        if os.path.isfile(snippy_consensus_concat_file):
                            logging.info(f"{snippy_consensus_concat_file} already found")
                        else:
                            logging.info(f"Concatenating Snippy consensus FASTA file (contigs)...")
                            ref = reference_file_local
                            concatenate_contigs(ref, snippy_consensus_file, isolate_in_alignment,
                                                snippy_consensus_concat_file)
                            if not os.path.exists(snippy_consensus_concat_file):
                                logging.error(f'{snippy_consensus_concat_file} not created. Something went wrong!')
                                sys.exit(-1)
                            else:
                                print(snippy_consensus_concat_file)
                        # Saving concatenated sequence
                        snippy_consensus_files.append(snippy_consensus_concat_file)

                        # Removing temporary files
                        logging.info(f"Removing temporary files ")
                        tmp_files = [reference_file_local, reference_file_local + '.fai', snippy_vcf_file_subs,
                                     snippy_vcf_file_noNs, snippy_consensus_file]
                        for tmp_file in tmp_files:
                            logging.info(f"Removing temporary file {tmp_file}")
                            run_command_shell(["rm -f", tmp_file])

                    # Creating host multiple alignment of consensus/pseudo-molecule sequence
                    multiple_aln = host_id + '.' + reference_isolate + '.aln'
                    if os.path.isfile(multiple_aln):
                        logging.info(f"Host multiple alignment {multiple_aln} already found")
                    else:
                        logging.info(f"Creating host multiple alignment...")
                        run_command_shell(
                            ["cat",
                             ' '.join(snippy_consensus_files),
                             ">",
                             multiple_aln],
                        )
                        if not os.path.exists(multiple_aln):
                            logging.error(f'{multiple_aln} not created. Something went wrong!')
                            sys.exit(-1)
                        else:
                            print(multiple_aln)

                    # Keeping only SNP sites in multiple alignment
                    multiple_aln_snp = multiple_aln.replace('.aln', '.snp.aln')
                    if os.path.isfile(multiple_aln_snp):
                        logging.info(f"Host SNP multiple alignment {multiple_aln} already found")
                    else:
                        logging.info(f"Finding SNP sites in multiple alignment {multiple_aln}...")
                        run_command_shell(
                            ["snp-sites ",
                             "-c",
                             "-o", multiple_aln_snp,
                             multiple_aln],
                        )
                        if not os.path.exists(multiple_aln_snp):
                            logging.error(f'{multiple_aln_snp} not created. Something went wrong!')
                            sys.exit(-1)
                        else:
                            print(multiple_aln)

                    # Duplicating alignment. As RAxML will fail if only three taxa, host taxa are duplicated
                    multiple_aln_snp_dup = multiple_aln.replace('.aln', '.snp.dup.aln')
                    if os.path.isfile(multiple_aln_snp_dup):
                        logging.info(f"Duplicated host SNP multiple alignment {multiple_aln_snp_dup} already found")
                    else:
                        logging.info(f"Duplicating host SNP multiple alignment {multiple_aln_snp}...")
                        duplicate_alignment(multiple_aln_snp, host_isolates, multiple_aln_snp_dup)

                    # Converting multiple alignment form FASTA to PHYLIP format
                    multiple_aln_phylip = multiple_aln_snp.replace('.snp.aln', '.snp.phy')
                    if os.path.isfile(multiple_aln_phylip):
                        logging.info(f"Host multiple alignment {multiple_aln_phylip} in PHYLIP format already found")
                    else:
                        align = AlignIO.read(multiple_aln_snp_dup, "fasta")
                        SeqIO.write(align, multiple_aln_phylip, "phylip-relaxed")

                    # Running RAxML
                    raxml_tree = host_id + '.' + reference_isolate + '.raxml_tree'
                    raxml_tree_file = "RAxML_bipartitions." + raxml_tree
                    if os.path.isfile(args.raxml_dir + raxml_tree_file):
                        logging.info(f"RAxML tree {raxml_tree_file} already found")
                    else:
                        logging.info(f"Running RAxML for {multiple_aln_phylip}...")
                        run_command_shell(
                            ["raxmlHPC-SSE3",
                             "-f a",
                             "-x 99999",
                             "-p 99999",
                             "-m GTRGAMMA",
                             "-N 100",
                             "-s", multiple_aln_phylip,
                             "-n", raxml_tree,
                             ],
                        )
                        if not os.path.exists(raxml_tree_file):
                            logging.error(f'{raxml_tree_file} not created. Something went wrong!')
                            sys.exit(-1)
                        else:
                            print(raxml_tree_file)
                        # Copying final tree file and removing rest of RAxML files
                        raxml_files_to_remove = ["RAxML_bestTree." + raxml_tree, "RAxML_bootstrap." + raxml_tree,
                                                 "RAxML_bipartitionsBranchLabels." + raxml_tree,
                                                 "RAxML_info." + raxml_tree]
                        for raxml_file_to_remove in raxml_files_to_remove:
                            logging.info(f"Removing temporary file {raxml_file_to_remove}")
                            run_command_shell(["rm -f", raxml_file_to_remove])
                        run_command_shell(["mv", raxml_tree_file, args.raxml_dir + raxml_tree_file])

                    # Rooting the tree on outgroup
                    rooted_raxml_tree_file = args.raxml_dir + raxml_tree_file + "_rooted"
                    if os.path.isfile(rooted_raxml_tree_file):
                        logging.info(f"Rooted RAxML tree {rooted_raxml_tree_file} already found")
                    else:
                        logging.info(f"Rooting RAxML {raxml_tree_file}...")
                        tree = dendropy.Tree.get_from_path(args.raxml_dir + raxml_tree_file, 'newick',
                                                           preserve_underscores=True)
                        outgroup_node = tree.find_node_with_taxon_label(outgroup_isolate)
                        tree.to_outgroup_position(outgroup_node, update_bipartitions=False)
                        tree.write(path=rooted_raxml_tree_file, schema="newick")

                    # Remove temporary files
                    logging.info(f"Removing temporary files ")
                    tmp_files = [multiple_aln, multiple_aln_snp, multiple_aln_snp_dup, multiple_aln_phylip,
                                 multiple_aln_phylip + '.reduced']
                    tmp_files.extend(snippy_consensus_files)
                    for tmp_file in tmp_files:
                        logging.info(f"Removing temporary file {tmp_file}")
                        run_command_shell(["rm -f", tmp_file])


if __name__ == "__main__":
    _main()