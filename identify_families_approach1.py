#!usr/bin/env/python 3

import argparse
import sys

import pandas as pd

def parse_args():
    """ Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='''Annotate genes with masked
                                                    family ID''')
    parser.add_argument(
        '-i','--input', nargs='?', type=argparse.FileType('r'),
        default=sys.stdin, help='An input file (CSV).')
    parser.add_argument(
        '-a','--annotation', required=True, type=argparse.FileType('r'),
        help='An annotation file containing gene ensembl ID and family (CSV).')
    parser.add_argument(
        '-o','--output', nargs='?', default=sys.stdout,
        help='Basename for output file (CSV).')
    return parser.parse_args()

def annotate_families(input,annotation,output):
    """TODO: add docstring"""
    df = pd.read_csv(input)
    df2 = pd.read_csv(annotation)
    df3 = df.merge(df2, how='inner', on='ensembl')
    df3.to_csv(f"{output}.csv", index=False)
    list_ids = df.ensembl.unique()
    with open(f"{output}.log", 'w+') as f:
        ids=''
        count=0
        for ensembl in list_ids:
            if ensembl not in df3.ensembl.unique():
                ids = ids + f"{ensembl}\n"
                count+=1
        f.write(f"ensembl IDs not in annotation file = {count}\n")
        f.write(ids)
    return df3

frequencies={}
ids={}
def unique_family_ids(df):
    for family_number in df.Family.unique():
        if frequencies.get(family_number,):
            frequencies[family_number]+=1
        else:
            frequencies[family_number]=1
    for k in frequencies:
        print(k,frequencies[k])
    print(len(df.Family.unique()),'\n')

    #attempt to pull out ensembl id for each family number row and append to list in dictionary.
    for row in df:
        if ids.get(df['Family'],):
            ids[df['Family']].append(df['ensembl'])
        else:
            ids[df['Family']]=[df['ensembl']]
    for k in ids:
        print(k,len(ids[k]))
        for i in ids[k]:
            print(i, end=' ')
    print('\n')

def main():
    args = parse_args()
    ann_df = annotate_families(args.input, args.annotation, args.output)
    unique_family_ids(ann_df)

if __name__ == "__main__":
    main()
