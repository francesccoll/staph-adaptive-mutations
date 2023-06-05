#!/usr/bin/env python3

import argparse
import logging
import os
import sys

# ---------------------------------------------------------------------------------------------------------------------
# Notes
# ---------------------------------------------------------------------------------------------------------------------

# This python script is used to create a tab-delimited table of annotated variants by snpEff from annotated snpEff
# VCF files. Both single-sample and multi-sample VCF files are supported
# Only GFF files are supported as genome annotation files

# ---------------------------------------------------------------------------------------------------------------------
# Development Notes
# ---------------------------------------------------------------------------------------------------------------------

# To do:
#   - read vcf.gz files


# ---------------------------------------------------------------------------------------------------------------------
# Global variables
# ---------------------------------------------------------------------------------------------------------------------

GFF_FIELD_COLUMNS = {'seqname': 0, 'source': 1, 'feature': 2, 'start': 3, 'end': 4, 'score': 5, 'strand': 6,
                     'frame': 7, 'attribute': 8, }


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to parse single-sample VCF files annotated by snpEff"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-i", "--input_vcf", action="store", dest="input_vcf",
        help="single-sample snpEff annotated VCF file",
        required=True, metavar="INPUT_VCF"
    )
    group.add_argument(
        "-o", "--output_table", action="store", dest="output_table",
        help="output table with snpEff annotated variants",
        required=True, metavar="OUTPUT_TABLE"
    )

    parser.add_argument(
        "-g", "--gff_file", action="store", dest="gff_file",
        help="GFF file used by snpEff",
        required=False, metavar="GFF_FILE"
    )
    parser.add_argument(
        "-f", "--gff_fields", action="store", dest="gff_fields",
        help="comma-delimited fields from GFF to save.\n"
             "To choose from: seqname, source, feature, start, end, score, strand, frame, attribute"
             " (e.g. feature,start,end)",
        required=False, metavar="GFF_FIELDS"
    )
    parser.add_argument(
        "-a", "--gff_ann_fields", action="store", dest="gff_ann_fields",
        help="comma-delimited fields from GFF attribute field to save (e.g. locus_tag,name)",
        required=False, metavar="GFF_ANN_FIELDS"
    )
    parser.add_argument(
        "-t", "--use_locus_tag", action="store_true", dest="use_locus_tag",
        help="Use GFF locus_tag as gene ID to match VCF SnpEff ANN entries",
        required=False
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------------------------------------------------
# Development notes
# ---------------------------------------------------------------------------------------------------------------------

# Script adapted to process SnpEff annotated VCFs with multiple fields in INFO


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
    if not os.path.isfile(args.input_vcf):
        logging.error(f'Input VCF file {args.input_vcf} not found!')
        sys.exit(-1)

    if args.gff_file is not None:
        if not os.path.isfile(args.gff_file):
            logging.error(f'Input GFF file {args.gff_file} not found!')
            sys.exit(-1)
        if args.gff_ann_fields is None:
            logging.error(f'--gff_ann_fields need to be defined if --gff_file used')
            sys.exit(-1)
        if args.gff_fields is None:
            logging.error(f'--gff_fields need to be defined if --gff_file used')
            sys.exit(-1)

    # Make sure user-selected gff_fields exist and convert to column indices
    if args.gff_file is not None:
        gff_fields = args.gff_fields.split(',')
        gff_fields_indices = list()
        for gff_field in gff_fields:
            if gff_field in GFF_FIELD_COLUMNS:
                gff_fields_indices.append(GFF_FIELD_COLUMNS[gff_field])
            else:
                logging.error(f'GFF field {gff_field} not allowed. Select from this list: '
                              f'seqname, source, feature, start, end, score, strand, frame, attribute')
                sys.exit(-1)

    # Parsing GFF file if chosen
    gff_field_values_for_gene = dict()  # dictionary to save snpeff_gene_id as key and gff fields to save as value
    gff_field_values_not_found = ''
    gff_header_fields = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    gff_field_values_header = ''
    if args.gff_file is not None:
        gff_ann_fields = args.gff_ann_fields.split(',')
        logging.info(f"Opening input GFF file {args.gff_file}")
        with open(args.gff_file, 'r') as gff_file:
            for gff_line in gff_file:
                gff_field_values_header = ''
                if not gff_line.startswith('#'):
                    print(gff_line)
                    gff_field_values_not_found = ''
                    (seqname, source, feature, start, end, score, strand, frame, attribute) = gff_line.strip().split('\t')
                    snpeff_gene_id = feature + '_' + seqname + '_' + start + '_' + end  # CDS_EMRSA15_1671393_1672145
                    gene_gff_fields_to_save = ''
                    for i, s in enumerate(gff_line.strip().split('\t')):
                        if i in gff_fields_indices:
                            gene_gff_fields_to_save += '\t' + s
                            gff_field_values_not_found += '\t' + '-'
                            gff_field_values_header += '\t' + gff_header_fields[i]
                    gff_attributes = attribute.split(';')
                    for gff_field in gff_ann_fields:
                        index = [i for i, s in enumerate(gff_attributes) if gff_field in s]
                        if len(index) == 1:
                            attr = gff_attributes[index[0]].replace(gff_field, '').strip().replace('\"', '').replace('=', '')
                        else:
                            attr = "-"
                        gene_gff_fields_to_save += '\t' + attr
                        gff_field_values_not_found += '\t' + '-'
                        gff_field_values_header += '\t' + gff_field
                        # NOTE: in some SnpEff annotated VCFs, gene names/ids can only be matched via GFF locus_tag
                        if args.use_locus_tag:
                            if gff_field == "locus_tag":
                                snpeff_gene_id = attr
                    gff_field_values_for_gene[snpeff_gene_id] = gene_gff_fields_to_save
                    # print("gff_field_values_for_gene["+snpeff_gene_id+"] --> "+gff_field_values_for_gene[snpeff_gene_id])
    print('GFF fields extracted: ' + gff_field_values_header)

    logging.info(f"Reading input VCF file {args.input_vcf}")
    saved_annotated_variants = dict()  # dictionary to keep track of saved annotated variants
    vcf_tab_file = open(args.output_table, 'w')
    vcf_ann_fields = ["chrm", "pos", "id", "ref", "alt", "qual", "filter", "info_left", "allele", "annotation",
                      "annotation_impact", "gene_name", "gene_id", "feature_type", "feature_id", "transcript_biotype",
                      "rank", "hgvd_c", "hgvd_p", "cdna_pos_length", "cds_pos_length", "aa_pos_length"]
    vcf_tab_file_header = '\t'.join(vcf_ann_fields) + gff_field_values_header + '\n'
    vcf_tab_file.write(vcf_tab_file_header)
    with open(args.input_vcf, 'r') as input_vcf:
        for vcf_line in input_vcf:
            if not vcf_line.startswith('#'):
                (chrm, pos, id, ref, alt, qual, filter, info, *_) = vcf_line.strip().split('\t')
                # Keeping INFO fields that are not ANN
                info_left = ""
                info_ann = ""
                for info_field in info.split(';'):
                    if info_field.startswith("ANN="):
                        info_ann = info_field
                    else:
                        info_left += info_field + ";"
                print('info_left ' + info_left)
                print('info_ann ' + info_ann)
                multiple_annotations = info_ann.split(',')
                for each_annotation in multiple_annotations:
                    if len(each_annotation.split('|')) >= 15:
                        (allele, annotation, annotation_impact, gene_name, gene_id, feature_type, feature_id,
                         transcript_biotype, rank, hgvd_c, hgvd_p, cdna_pos_length, cds_pos_length,
                         aa_pos_length, distance, *_) = each_annotation.split('|')
                        allele = allele.replace("ANN=", '')
                        annotated_variant_id = pos + '_' + allele + '_' + gene_id
                        # Making sure annotated variant is not already saved
                        if annotated_variant_id not in saved_annotated_variants:
                            saved_annotated_variants[annotated_variant_id] = 1
                            gff_field_values = ''
                            if feature_type == "intergenic_region":
                                (transcript_biotype, rank, hgvd_c, hgvd_p, cdna_pos_length, cds_pos_length, aa_pos_length)\
                                    = ("-", "-", "-", "-", "-", "-", "-")
                                (left_gene, right_gene) = gene_name.split('-')
                                (left_gene_gff_field_values, right_gene_gff_field_values) = ([], [])
                                if left_gene in gff_field_values_for_gene:
                                    left_gene_gff_field_values = gff_field_values_for_gene[left_gene].split('\t')
                                else:
                                    left_gene_gff_field_values = gff_field_values_not_found.split('\t')
                                if right_gene in gff_field_values_for_gene:
                                    right_gene_gff_field_values = gff_field_values_for_gene[right_gene].split('\t')
                                else:
                                    right_gene_gff_field_values = gff_field_values_not_found.split('\t')
                                gff_field_values_zipped = [';'.join(map(str, i)) for i in zip(left_gene_gff_field_values,
                                                                                             right_gene_gff_field_values)]
                                gff_field_values = '\t'.join(gff_field_values_zipped)
                                gff_field_values = gff_field_values.replace(';\t', '\t')
                            else:
                                if gene_name in gff_field_values_for_gene:
                                    gff_field_values = gff_field_values_for_gene[gene_name]
                                else:
                                    gff_field_values = gff_field_values_not_found
                            vcf_line_elements = (chrm, pos, id, ref, alt, qual, filter, info_left, allele, annotation,
                                                 annotation_impact, gene_name, gene_id, feature_type, feature_id,
                                                 transcript_biotype, rank, hgvd_c, hgvd_p, cdna_pos_length,
                                                 cds_pos_length, aa_pos_length)
                            vcf_tab_line = '\t'.join(vcf_line_elements) + gff_field_values + '\n'
                            vcf_tab_file.write(vcf_tab_line)
    logging.info(f"Writing output table file {args.output_table}")
    input_vcf.close()
    vcf_tab_file.close()


if __name__ == "__main__":
    _main()
