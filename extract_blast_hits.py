#!usr/bin/env python 3

import argparse

from Bio.Blast import NCBIXML # Biopython version 1.73
import pandas as pd # Pandas version 0.24.2

def parse_args():
    """ Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i','--input', required = True,
        help='Input file in BLAST XML format.')
    parser.add_argument(
        '-o','--output', required = True, help='Output file (CSV format).')
    parser.add_argument(
        '-a','--annotation', help='Annotation file (CSV format).')
    return parser.parse_args()

def read(file):
    result_handle = open(file)
    return result_handle

def extract_hits(result_handle, output):
    """ Extracts all hits from BLAST XML format output in CSV format.
    Output:
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
                + f"total searches and hits,{len(blast_records)}\n")
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
    """ Extracts list of genes and matched query ids from CSV file of hits.

    Returns:
        List of genes, for use in annotation.

    Outputs:
        Output file with genes and matched query ids in CSV format.
    """
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
    list_genes = []
    with open(f"{output}_unique_genes.csv", 'w+') as f:
        for key, value in query_gene_dict.items():
            list_genes.append(key)
            f.write(f"{key},{','.join(value)}\n")
    return list_genes

def annotate(list_genes,annotation_file,output):
    df = pd.read_csv(annotation_file)
    df2 = df[df.name.isin(list_genes)]
    with open(f"{output}.log", 'w+') as f:
        f.write(f"Total unique genes: {len(list_genes)-1}\n")
        f.write("Gene annotations not found:\n")
        for gene in list_genes:
            if gene == 'gene':
                pass
            elif gene not in df['name'].unique():
                f.write(f"{gene}\n")
    return df2

def main():
    args = parse_args()
    result_handle = read(args.input)
    extract_hits(result_handle, args.output)
    list_genes = extract_gene_list(args.output)
    annotate(list_genes, args.annotation, args.output)

if __name__ == "__main__":
    main()
