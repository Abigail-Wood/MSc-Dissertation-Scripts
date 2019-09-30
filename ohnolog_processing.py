import pandas as pd

def read(file):
    df = pd.read_csv(file)
    return df

def read_ohnologs(ohnologs_file):
    """ Ohnolog2 format data read in human data (contains tabs)"""
    ohnologs = pd.read_csv(ohnologs_file, sep='\t')
    ohnolog_ids = set(ohnologs['Ohno1']).union(ohnologs['Ohno2'])
    return ohnolog_ids

def read_species(homologs_file, ohnologs_file):
    """ Ohnolog2 format data read in (contains tabs)"""
    homologs = pd.read_csv(homologs_file)
    ohnolog_ids = read_ohnologs(ohnologs_file)
    ohnolog_homologs = set(homologs.loc[homologs[homologs.columns[1]].isin(ohnolog_ids), 'Gene stable ID'])
    return ohnolog_homologs

def huang_haplosufficiency(file, output):
    with open(file, 'r') as f, open(output, 'w+') as o:
        for line in f:
            o.write(line.replace('|',','))
    return

def ohnologs(df, file):
    ohnologs = []
    with open(file, 'r') as f:
        for line in f:
            ohnologs.append(line.strip())
    df.loc[df['Gene stable ID'].isin(ohnologs), 'Paralog status'] = 'Ohnolog'
    df.to_csv("Ohnologs/paralog_classifications_all_ohnologs.csv",index=False)

human_ohnologs = {}
def multi_tier_ohnologs(species_file, human_ohnolog_file,output):
    """ Returns human ohnologs and their frequencies in other species."""
    species = []
    with open(species_file,'r') as f:
        for line in f:
            line = line.strip().split(',')
            species.append(line)
    human_ohnolog_ids = read_ohnologs(human_ohnolog_file)
    for id in human_ohnolog_ids:
        if id not in human_ohnologs:
            human_ohnologs[id] = []
    species_names = []
    for name, homologs_file, ohnologs_file in species:
        species_names.append(name)
        ohnolog_homologs = read_species(homologs_file, ohnologs_file)
        for id in ohnolog_homologs:
            if id in human_ohnologs:
                human_ohnologs[id].append(name)
    species_names.sort()
    with open(output,'w+') as o:
        o.write(f"Gene stable ID,Ohnolog Count,{','.join(species_names)}\n")
        for key in human_ohnologs:
            species = human_ohnologs[key]
            species_values = ','.join('1' if name in species else '0' for name in species_names)
            o.write(f"{key},{len(species)},{species_values}\n")

def classify_dataset(df, df2,output):
    df = df.merge(df2, how='left', on='ensembl')
    df.to_csv(output, index=False)

def count_paralogs(dfs,output):
    counts = pd.Series()
    for df in dfs:
        counts = counts.append(df['Paralog status'].value_counts(normalize=True,dropna=False) * 100)
    counts.to_csv(output,header=False)
    
def main():
    # huang_haplosufficiency("Ohnologs/haplosufficiency_rank_predictions_unseparated.csv","Ohnologs/haplosufficiency_rank_predictions_separated.csv")
    # df = read("Ohnologs/paralog_classifications_all_no_ohnologs.csv")
    # df2 = read("Ohnologs/hsapiens.Pairs.Strict.2R.csv")
    # ohnologs(df, "Ohnologs/ohnolog2_ids.txt")
    # ohnolog_homologs = read_species("Ohnologs/species_data/mmusculus_human_homologs.csv", "Ohnologs/species_data/mmusculus.Pairs.Strict.2R.tsv")
    # multi_tier_ohnologs("Ohnologs/species_files.txt","Ohnologs/species_data/hsapiens.Pairs.Strict.2R.tsv","Ohnologs/output_counting.csv")
    df = read("whole_proteome_blast/excluded_annotated.csv")
    df2 = read("Ohnologs/paralog_classifications_all_ohnologs.csv")
    classify_dataset(df,df2,"Ohnologs/excluded_ohnologs.csv")
    df_files = ["Ohnologs/innate_filtered_ohnologs.csv","Ohnologs/background_innate_ohnologs.csv","Ohnologs/excluded_ohnologs.csv","Ohnologs/filtered_ohnologs.csv","Ohnologs/paralog_classifications_all_ohnologs.csv"]
    dfs = []
    for file in df_files:
        dfs.append(read(file))
    count_paralogs(dfs,"Ohnologs/ohnolog_table.csv")
    
if __name__ == "__main__":
    main()
