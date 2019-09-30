import pandas as pd

def read(file):
    df = pd.read_csv(file)
    return df

def excluded(df_inc, df_total):
    list_genes = list(df_inc['ensembl'])
    df_excl = df_total[~df_total.ensembl.isin(list_genes)]
    df_excl['goFunctions'].fillna(value='No terms', inplace=True) # ensures we can index on this column
    df_excl.to_csv("excluded.csv", index=False)
    return df_excl

def separate_kinases(df, output):
    """ Takes an annotated dataframe and a specified output filename."""
    df_Kinases = df[df['goFunctions'].str.contains('GO:0004672')]
    df_Kinases.to_csv(f"{output}_kinases.csv", index=False)
    kinase_list = list(df_Kinases['ensembl'])# ensembl ID where kinase term in goFunctions
    df_noKinases = df[~df.ensembl.isin(kinase_list)] # tilde implements 'not in'
    df_noKinases.to_csv(f"{output}_no_kinases.csv", index=False)

def extract_headers(df_excl, file):
    ensembl_ids = list(df_excl['ensembl'])
    selected_lines = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                line = line.strip('>')
                if (line.split('|')[0] in ensembl_ids):
                    selected_lines.append(line.strip().split('|')) # add line to dataframe
            else:
                pass
    df_out = pd.DataFrame(selected_lines, columns=['ensembl','gene','transcript','ensemblfamily'])
    df_out = df_out.sort_values(by='ensembl')
    df_out.to_csv("excluded_families.csv",index=False)
                                  
def condense_families(file, output):
    gene_id = None
    gene_name = None
    family_ids = None
    with open(file, 'r') as f, open(output,'w+') as o:
        for line in f:
            values = line.split(',')
            if values[0] == gene_id:
                s = values[3].split('_')[0]
                family_id = s.strip()
                if family_id is not '':
                    family_ids.add(family_id)
            else:
                if gene_id is not None: # if there is a gene_id but it is not the current one
                    s = ' '.join(family_ids)
                    o.write(f"{gene_id},{gene_name},{s}\n")
                    family_ids = set()
                gene_id = values[0]
                gene_name = values[1]
                s = values[3].split('_')[0]
                s = s.strip()
                if s is not '':
                    family_ids = set([s])
        s = ' '.join(family_ids)
        o.write(f"{gene_id},{gene_name},{s}\n") # captures last line in file

def families_kinases_condensed(df_condensed, df_kinases):
    kinases = list(df_kinases['ensembl'])
    df_kinases = df_condensed[df_condensed.ensembl.isin(kinases)]
    df_kinases.to_csv("excluded_kinases_families_condensed.csv", index=False)
    df_no_kinases = df_condensed[~df_condensed.ensembl.isin(kinases)]
    df_no_kinases.to_csv("excluded_no_kinases_families_condensed.csv", index=False)

def main():
    df = read("whole_proteome_blast/filtered_innate.csv")
    df2 = read("human_innate_genes/InnateDB_genes_noMiR.csv")
    df_excl = excluded(df, df2)
    separate_kinases(df_excl, "excluded")
    extract_headers(df_excl, "allxms.fa")
    condense_families("excluded_families.csv", "excluded_families_condensed.csv")
    df_condensed = read("excluded_families_condensed.csv")
    df_kinases = read("excluded_kinases_annotated.csv")
    families_kinases_condensed(df_condensed, df_kinases)
        
if __name__ == "__main__":
    main()
