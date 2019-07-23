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

def remove_subfamilyids(df):
    # Apply a function to one column and assign it back to the dataframe
    df['ensemblfamily'] = df['ensemblfamily'].apply(lambda s: s.split('_')[0])
    return df

ids={}
def unique_family_ids(df, file):
    #attempt to pull out ensembl id for each family number row and append to list in dictionary.
    for index, row in df.iterrows():
        if ids.get(row['ensemblfamily'],): # if family id is already in dictionary
            if row['ensembl'] in set(ids.get(row['ensemblfamily'],)):
                pass
            else:
                ids[row['ensemblfamily']].append(row['ensembl']) # append the new ensembl id to the list
        else:
            ids[row['ensemblfamily']]=[row['ensembl']] # assign new list to family id with ensembl id as identifier/key
    with open(file,'w+') as f:
        for family, ensembl_ids in ids.items():
            f.write(f"{family},{len(ids[family])},")
            for id in ensembl_ids:
                f.write(f"{id},")
            f.write("\n")
            

def main():
    args = parse_args()
    df = read_csv(args.input)
    df = remove_subfamilyids(df)
    unique_family_ids(df, args.output)

if __name__ == "__main__":
    main()
