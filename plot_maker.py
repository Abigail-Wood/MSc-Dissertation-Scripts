#!usr/bin/env/python 3

import sys

import pandas as pd
import matplotlib.pyplot as plt

def family_histogram(file,title):
    df = pd.read_csv(file)
    x = df['number_of_members']

    n, bins, patches = plt.hist(x, df.shape[1]-3, color='green')
    plt.title(title)
    plt.xlabel('Family size (number of unique proteins found)')
    plt.ylabel('Number of families')
    plt.xticks(range(1,df['number_of_members'].max()))
    plt.xlim(1,)

    cm = plt.cm.get_cmap('RdYlBu_r')
    col = (n-n.min())/(n.max()-n.min())
    for c, p in zip(col, patches):
        plt.setp(p, 'facecolor', cm(c))
    plt.savefig("test_figure.png")
    plt.show()
    return

def family_overlap_histogram(file1, file2, file3):
    plt.style.use('seaborn-deep')
    
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    df3 = pd.read_csv(file3)
    x = df1['number_of_members']
    y = df2['number_of_members']
    z = df3['number_of_members']
    
    plt.hist([x, y, z], range(y.min(),y.max()+2), label=['Retained (No repeats)', 'Retained (with repeats)', 'Excluded'], align='left')
    plt.legend(loc='upper right')
    plt.xlabel('Family size (number of unique proteins found)')
    plt.ylabel('Number of families')
    plt.xticks(range(y.min(),y.max()+1))
    plt.xlim(0.5,y.max()+0.5)
    
    plt.show()
    return

def matches_per_gene_histogram(file):
    df = pd.read_csv(file)
    x = df['matches']
    # weird what we have to do with the bins here, but without an empty bin on the end the
    # x.max() values get binned into the previous bin.
    n, bins, patches = plt.hist(x, range(1, x.max()+2), align='left')
    plt.xlim(left=0.5, right=x.max()+0.5)
    plt.xticks(range(x.min(),x.max()+1))
    plt.xlabel('Number of query matches per gene')
    plt.ylabel('Genes (with x number of matches)')
    

    cm = plt.cm.get_cmap('RdYlBu_r')
    col = (n-n.min())/(n.max()-n.min())
    for c, p in zip(col, patches):
        plt.setp(p, 'facecolor', cm(c))
    
    plt.show()
    return 

def main():
    # family_histogram("BLAST_results/blastp_trichuris_proteins_combined_family_numbers_filtered.csv", 'All Proteins')
    # family_overlap_histogram("whole_proteome_blast/filtered_innate_non_unique_family_numbers.csv","whole_proteome_blast/no_kinases_non_unique_family_numbers.csv","whole_proteome_blast/kinases_non_unique_family_numbers.csv")
    # matches_per_gene_histogram("whole_proteome_blast/whole_proteome_matches_per_gene_fig.csv")
    family_overlap_histogram("whole_proteome_blast/filtered_innate_family_numbers.csv","whole_proteome_blast/filtered_innate_non_unique_family_numbers.csv","whole_proteome_blast/excluded_family_numbers.csv")

    

if __name__ == '__main__':
    main()
