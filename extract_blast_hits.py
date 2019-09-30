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
    Output (unfiltered):
        Query id
        Hit identifier
        Gene name
        Transcript
        Ensembl ID
        Ensembl Family ID (derived from PANTHER14.1 terms)
        Subfamily ID
        Length
        Bits
        Identities
        Percentage Identity (calculated using identities/length * 100)
        E value
    Additional filtered output takes only the top 3 blast hits for each query.
    """
    blast_records = NCBIXML.parse(result_handle) # returns an iterator.
    blast_records = list(blast_records)
    with open(f"{output}_unfiltered.csv", 'w+') as f:
        f.write(f"query,hit identifier,gene,transcript,ensembl,length,bits,identities,percent_identity,e value,ensemblfamily,subfamily\n")
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                if alignment.title: # if a hit was found for a query
                    for hsp in alignment.hsps:
                        if (hsp.identities / alignment.length * 100) >= 10:
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
                                if '_' in id[5]:
                                    f.write(f"{id[5].split('_')[0]},{id[5].split('_')[1]}\n")
                                else:
                                    f.write(f"{id[5].split('_')[0]}\n")
                            else:
                                f.write("no family\n")
    df = pd.read_csv(f"{output}_unfiltered.csv")
    # e-value column take top 3
    df = df.groupby('query', as_index=False).head(3) # works as df descending sorted by e-value
    # sort each group by % identity
    df = df.groupby('query', as_index=False).apply(lambda x: x.sort_values(['percent_identity'], ascending=False))
    # take top 1 by percent identity
    df = df.groupby('query', as_index=False).head(1)
    df.to_csv(f"{output}_filtered.csv", index=False) 
    return 0

query_gene_dict = {}
def extract_gene_list(output):
    """ Extracts list of genes and matched query ids from CSV file of hits.

    Returns:
        List of ensembl IDs, for use in annotation.

    Outputs:
        Output file with genes and matched query ids in CSV format.
    """
    with open(f"{output}_filtered.csv",'r') as f:
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
    with open(f"{output}_matches_per_gene.csv", 'w+') as f:
        for key, value in query_gene_dict.items():
            list_genes.append(key)
            f.write(f"{key},{len(value)},{','.join(value)}\n")
    return list_genes

def annotate(list_genes,annotation_file,output):
    """Annotates a list of genes.

    Writes ensembl IDs of an non-annotated genes to a .log file.
    """
    df = pd.read_csv(annotation_file)
    df2 = df[df.ensembl.isin(list_genes)]
    df3 = df[~df.ensembl.isin(list_genes)]
    for ensemblID in list_genes:
        if ensemblID == 'ensembl':
            pass
        elif ensemblID not in df['ensembl'].unique():
            with open(f"{output}.log", 'w+') as f:
                f.write(f"{ensemblID}\n")
    df2.to_csv(f"{output}_annotated_genes.csv", index=False)
    df3.to_csv(f"{output}_excluded_genes.csv", index=False)
    return df2

def no_kinases(output, df):
    """ Takes an annotated dataframe and a specified output filename."""
    df = df[df['goFunctions'].str.contains('GO:0004672')]
    kinase_list = list(df['ensembl'])
    # kinase_list of ensembl IDs where 'kinase' is in fullname for df2
    df2 = pd.read_csv(f"{output}_filtered.csv")
    df_noKinases = df2[~df2.ensembl.isin(kinase_list)] # tilde implements 'not in'
    df_noKinases.to_csv(f"{output}_no_kinases.csv", index=False)
    # also not kinases output file.
    df_Kinases = df2[df2.ensembl.isin(kinase_list)]
    df_Kinases.to_csv(f"{output}_kinases.csv", index=False)
    return df

def main():
    args = parse_args()
    result_handle = read(args.input)
    extract_hits(result_handle, args.output)
    list_genes = extract_gene_list(args.output)
    df2 = annotate(list_genes, args.annotation, args.output)
    df = no_kinases(args.output, df2)

if __name__ == "__main__":
    main()
