#!usr/bin/env/python 3

import argparse
import sys

import pandas as pd

def parse_args():
    """ Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='''Identify family distribution with unmasked Ensembl ID''')
    parser.add_argument(
        '-i','--input', nargs='?', type=argparse.FileType('r'),
        default=sys.stdin, help='An input file (CSV).')
    parser.add_argument(
        '-o','--output', nargs='?', default=sys.stdout,
        help='Basename for output file (CSV).')
    return parser.parse_args()

def read_csv(file):
    df = pd.read_csv(file, sep=',')
    return df

ids={}
def family_with_genes(df, file):
    #attempt to pull out ensembl id for each family number row and append to list in dictionary.
    for index, row in df.iterrows():
        if ids.get(row['ensemblfamily'],): # if family id is already in dictionary
            if row['gene'] in set(ids.get(row['ensemblfamily'],)):
                pass
            else:
                ids[row['ensemblfamily']].append(row['gene']) # append the new gene id to the list
        else:
            ids[row['ensemblfamily']]=[row['gene']] # assign new list to family id with gene id as identifier/key
    with open(f"{file}.csv", 'w+') as f:
        f.write("family,number_of_members,members_found\n")
        for family, genes in ids.items():
            f.write(f"{family},{len(genes)},")
            for gene in genes:
                f.write(f"{gene},")
            f.write("\n")

ids={}
def family_non_unique(df, file):
    #attempt to pull out ensembl id for each family number row and append to list in dictionary.
    for index, row in df.iterrows():
        if ids.get(row['ensemblfamily'],): # if family id is already in dictionary
            ids[row['ensemblfamily']].append(row['gene']) # append the new gene id to the list
        else:
            ids[row['ensemblfamily']]=[row['gene']] # assign new list to family id with gene id as key
    with open(f"{file}.csv", 'w+') as f:
        f.write("family,number_of_members,members_found\n")
        for family, genes in ids.items():
            f.write(f"{family},{len(genes)},")
            for gene in genes:
                f.write(f"{gene},")
            f.write("\n")

def main():
    args = parse_args()
    df = read_csv(args.input)
    family_non_unique(df, args.output)

if __name__ == "__main__":
    main()
