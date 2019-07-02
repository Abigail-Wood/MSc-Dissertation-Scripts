#!usr/bin/env python 3

import argparse

from Bio.Blast import NCBIXML # Biopython version 1.73

def parse_args():
    """ Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input', required = True, help='Input file in BLAST XML format.')
    parser.add_argument(
        '--output', required = True, help='Output file.')
    return parser.parse_args()

def read(file):
    result_handle = open(file)
    return result_handle

def extract_hits(result_handle, output):
    """ Extracts all hits from BLAST XML format output in CSV format.
    Returns:
        Query id
        Hit identifier
        Gene name
        Length
        E value
        Total searches and hits
    """
    blast_records = NCBIXML.parse(result_handle) # returns an iterator.
    blast_records = list(blast_records)
    with open(f"{output}.csv", 'w+') as f:
        f.write("query,hit identifier,gene,length,e value,"
                + f"total searches and hits,{len(blast_records)}")
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                if alignment.title: # if a hit was found for a query
                    for hsp in alignment.hsps:
                        f.write(f"{blast_record.query},{alignment.title},"
                                + f"{alignment.title.split('|')[4]},"
                                + f"{alignment.length},"
                                + f"{hsp.expect}\n")

query_gene_dict = {}
def extract_gene_list(output):
    """ Extracts list of genes and matched query ids from CSV file of hits."""
    with open(f"{output}.csv",'r') as f:
        for line in f.readlines():
            line= line.strip().split(',')
            if line[2] in query_gene_dict: # if gene in dictionary keys
                if line[0] in query_gene_dict[line[2]]: #query is in dictionary
                    pass
                else: # if gene not value for gene in dictionary, add it
                    query_gene_dict[line[2]].append(line[0])
            else: # add gene to dictionary keys with id value
                query_gene_dict[line[2]] = [line[0]]
    with open(f"{output}_unique_genes.csv", 'w+') as f:
        for key, value in query_gene_dict.items():
            f.write(f"{key},{','.join(value)}\n")

def main():
    args = parse_args()
    result_handle = read(args.input)
    extract_hits(result_handle, args.output)
    extract_gene_list(args.output)

if __name__ == "__main__":
    main()
