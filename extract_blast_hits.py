#!usr/bin/env python 3

import argparse

from Bio.Blast import NCBIXML # Biopython version 1.73
import pandas as pd # Pandas version 0.24.2

def parse_args():
    """ Parse command-line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i','--input', required = True, help='Input in BLAST XML format.')
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
        Transcript
        Ensembl ID
        Ensembl Family ID
        Length
        Bits
        Identities
        Percentage Identity (calculated using identities/length * 100)
        E value
        Total searches and hits
    """
    blast_records = NCBIXML.parse(result_handle) # returns an iterator.
    blast_records = list(blast_records)
    with open(f"{output}.csv", 'w+') as f:
        f.write(f"query,hit identifier,gene,transcript,ensembl,length,bits,identities,percent_identity,e value,ensemblfamily,total searches and hits,{len(blast_records)}\n")
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                if alignment.title: # if a hit was found for a query
                    for hsp in alignment.hsps:
                        id = alignment.title.split('|')
                        f.write(f"{blast_record.query},{alignment.title},"
                                + f"{id[3]},"
                                + f"{id[4]},"
                                + f"{id[2].split()[1]},"
                                + f"{alignment.length},"
                                + f"{hsp.bits},"
                                + f"{hsp.identities},"
                                + f"{hsp.identities / alignment.length * 100},"
                                + f"{hsp.expect},")
                        if len(id)>=6:
                            f.write(f"{id[5]}\n")
                        else:
                            f.write("no family\n")

    return 0

query_gene_dict = {}
def extract_gene_list(output):
    """ Extracts list of genes and matched query ids from CSV file of hits.

    Returns:
        List of ensembl IDs, for use in annotation.

    Outputs:
        Output file with genes and matched query ids in CSV format.
    """
    with open(f"{output}.csv",'r') as f:
        for line in f.readlines():
            line= line.strip().split(',')
            if line[4] in query_gene_dict: # if ensembl id in dictionary keys
                if line[0] in query_gene_dict[line[4]]: #query is in dictionary
                    pass
                else: # if ensembl id not in dictionary, add it
                    query_gene_dict[line[4]].append(line[0])
            else: # add ensembl to dictionary keys with id value
                query_gene_dict[line[4]] = [line[0]]
    list_genes = []
    with open(f"{output}_unique_genes.csv", 'w+') as f:
        for key, value in query_gene_dict.items():
            list_genes.append(key)
            f.write(f"{key},{','.join(value)}\n")
    return list_genes

def annotate(list_genes,annotation_file,output):
    """TODO: add docstring"""
    df = pd.read_csv(annotation_file)
    df2 = df[df.ensembl.isin(list_genes)]
    with open(f"{output}.log", 'w+') as f:
        f.write(f"Total unique genes: {len(list_genes)-1}\n")
        f.write("Gene annotations not found:\n")
        for ensemblID in list_genes:
            if ensemblID == 'ensembl':
                pass
            elif ensemblID not in df['ensembl'].unique():
                f.write(f"{ensemblID}\n")
    df2.to_csv(f"{output}_annotated.csv", index=False)
    return df2

def main():
    args = parse_args()
    result_handle = read(args.input)
    extract_hits(result_handle, args.output)
    list_genes = extract_gene_list(args.output)
    annotate(list_genes, args.annotation, args.output)

if __name__ == "__main__":
    main()
