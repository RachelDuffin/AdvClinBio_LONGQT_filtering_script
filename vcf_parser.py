import pandas as pd
import re
import os
import argparse
from itertools import islice


def arg_parse():
    """
    Creates argument parser, defines command line arguments, then parses
    supplied command line arguments using the created argument parser.
    """
    parser = argparse.ArgumentParser()
    # defines command line arguments
    parser.add_argument('-i', '--input_file', help="The input .txt file " +
                        "within the current directory (the file you " +
                        "downloaded from VEP")
    parser.add_argument('-o', '--output_file', help="Name for the output " +
                        ".txt file")

    return parser.parse_args()


def gnomAD_AF_filter(line, index):
    """
    For variants with a gnomad allele frequency, filters out those with a
    frequency of < 0.05
    """
    gnomad_AF = re.split(r'\t+', line.rstrip('\t'))[index]
    if "-" in gnomad_AF:
        return True
    else:
        if float(gnomad_AF) < 0.05:
            return True
        else:
            return False


def consequence_filter(line, index, consequence_list):
    """
    Filters out all variants that have the consequences supplied in
    consequence_list
    """
    consequence = re.split(r'\t+', line.rstrip('\t'))[index]
    if not any(variant_type in consequence for variant_type in
               consequence_list):
        return True
    else:
        return False


def transcript_filter(line, transcript_list):
    """
    Filters out all variants that do not match the supplied transcript list
    """
    if any(transcript in line for transcript in transcript_list):
        return True
    else:
        return False


def get_index(df, string):
    """
    Gets the column indices for the columns that the variants will be filtered
    on
    """
    return df.columns.to_list().index(string)


def write_columns(infile, outfile, fields):
    """
    Write columns from infile matching supplied fields to outfile
    """
    with open(outfile, 'w') as outfile:
        df = pd.read_csv(infile, delimiter='\t')
        output = df[fields]
        output.to_csv(outfile, sep='\t')
    outfile.close()


def main():
    args = arg_parse()
    df = pd.read_csv(args.input_file, delimiter="\t")
    # gets the required column indexes for applying filters to
    gnomAD_AF_index = get_index(df, string='gnomAD_AF')
    consequence_index = get_index(df, string="Consequence")
    # specify lists of consequences to filter out
    consequence_list_to_ignore = ['downstream_gene_variant', 'intron_variant',
                                  'upstream_gene_variant',
                                  '3_prime_UTR_variant',
                                  'non_coding_transcript_variant',
                                  'non_coding_transcript_exon_variant',
                                  'synonymous_variant', '5_prime_UTR_variant',
                                  'splice_donor_variant']
    # specify list of transcripts to filter on
    supplied_transcripts = ['NM_199460', 'NM_006888', 'NM_001127670',
                            'NM_172201', 'NM_000238', 'NM_000218', 'NM_000335',
                            'NM_001099404', 'NM_198056', 'NM_005751',
                            'NM_000890', 'NM_174934', 'NM_003098', 'NM_001148',
                            'NM_000891']

    # write to intermediate file all variants remaining after applying filters
    with open(args.input_file, 'r') as infile:
        with open("intermediate.txt", 'w') as outfile:
            # write header line
            header_line = infile.readline()
            outfile.write(header_line)
            # start filtering from line 2
            for line in islice(infile, 1, None):
                if (gnomAD_AF_filter(line, gnomAD_AF_index)
                and consequence_filter(line, consequence_index,
                                       consequence_list_to_ignore)
                and transcript_filter(line, supplied_transcripts)):
                    outfile.write(line)
                else:
                    pass
            outfile.close()
            # write to final output all useful columns
            fields_needed = ['Location', 'Allele', 'Consequence', 'SYMBOL',
                             'Gene', 'Feature_type', 'Feature', 'BIOTYPE',
                             'EXON', 'cDNA_position', 'CDS_position',
                             'Protein_position', 'Amino_acids', 'Codons',
                             'Existing_variation', 'HGNC_ID', 'MANE_SELECT',
                             'MANE_PLUS_CLINICAL', 'TSL', 'SIFT', 'PolyPhen',
                             'gnomAD_AF', 'CLIN_SIG', 'PUBMED']
            write_columns('intermediate.txt', args.output_file, fields_needed)
            os.remove("intermediate.txt")

if __name__ == "__main__":
    main()
