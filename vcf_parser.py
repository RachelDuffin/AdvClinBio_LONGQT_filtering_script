import sys
import pandas as pd
import csv
import re
from itertools import islice

def gnomAD_AF_filter(line, index):
    """
    Filter out all variants with gnomAD allele frequency of less than 0.05
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
    Filter out all variants that are not:
    - downstream_gene_variant 
    - intron_variant 
    - upstream_gene_variant 
    - 3_prime_UTR_variant 
    - non_coding_transcript_variant 
    - non_coding_transcript_exon_variant 
    - synonymous_variant
    """
    consequence = re.split(r'\t+', line.rstrip('\t'))[index]
    if not any(variant_type in consequence for variant_type in consequence_list):
        return True
    else:
        return False

def transcript_filter(line, select_index, clinical_index):
    """
    Filter out all variants without a "MANE SELECT" or "MANE PLUS CLINICAL" transcript
    """
    MANE_SELECT_transcript = re.split(r'\t+', line.rstrip('\t'))[select_index]
    MANE_CLINICAL_transcript = re.split(r'\t+', line.rstrip('\t'))[clinical_index]
    if "-" in MANE_SELECT_transcript and "-" in MANE_CLINICAL_transcript:
        return False
    else: 
        return True

def tsl_filter(line, index):
    """
    Filter for minimum transcript support level - filters out transcripts with a value above 2 
    """
    tsl = re.split(r'\t+', line.rstrip('\t'))[index]
    if float(tsl) < 3:
        return True

def get_index(df, string):
    """
    Get the column indices for the columns that the variants will be filtered on
    """
    return df.columns.to_list().index(string)

def write_columns(infile, outfile, fields):
    with open(outfile, 'w') as outfile:
        df = pd.read_csv(infile, delimiter='\t')
        output = df[fields]
        output.to_csv(outfile, sep='\t')


def main():
    df = pd.read_csv(sys.argv[1], delimiter = "\t")
    gnomAD_AF_index = get_index(df, string='gnomAD_AF')   
    consequence_index = get_index(df, string="Consequence")
    MANE_SELECT_index = get_index(df, string="MANE_SELECT")
    MANE_CLINICAL_index = get_index(df, string="MANE_PLUS_CLINICAL")
    tsl_index = get_index(df, string="TSL")
    consequence_list_to_ignore = ['downstream_gene_variant', 'intron_variant', 'upstream_gene_variant', '3_prime_UTR_variant',
            'non_coding_transcript_variant', 'non_coding_transcript_exon_variant', 'synonymous_variant', '5_prime_UTR_variant', 'splice_donor_variant']
    with open(sys.argv[1], 'r') as rfp:
        with open('filtered_variant.txt', 'w') as wfp:
            # write header line
            header_line = rfp.readline().rstrip()
            wfp.write(header_line)
            # start filtering from line 2
            for line in islice(rfp, 1, None):
                if gnomAD_AF_filter(line, gnomAD_AF_index) and consequence_filter(line, consequence_index, consequence_list_to_ignore) and transcript_filter(line, MANE_SELECT_index, MANE_CLINICAL_index) and tsl_filter(line, tsl_index):
                    wfp.write(line)
                else:
                    pass
            wfp.close()
            fields_needed = ['Location', 'Allele', 'Consequence', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation','HGNC_ID', 'MANE_SELECT', 'MANE_PLUS_CLINICAL', 'TSL', 'SIFT', 'PolyPhen', 'gnomAD_AF', 'CLIN_SIG', 'PUBMED']
            write_columns('filtered_variant.txt', 'final_variants.txt', fields_needed)


if __name__ == "__main__":
    main()   
        
