# remove the exact matches from the known file
# adds the adjusted ranking based on the FGs

import sys
import argparse
import pandas as pd
from rdkit import Chem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from utils import is_file_match
import multiprocessing as mp
from functools import partial

def get_with_structure(df):
    merged_results_with_gnps = df.copy()
    merged_results_with_gnps = merged_results_with_gnps[merged_results_with_gnps['Smiles'].isna() == False]
    merged_results_with_gnps = merged_results_with_gnps[merged_results_with_gnps['Smiles'] != '']
    merged_results_with_gnps = merged_results_with_gnps[merged_results_with_gnps['Smiles'] != ' ']
    merged_results_with_gnps['Smiles'] = merged_results_with_gnps['Smiles'].apply(lambda x: x.replace('&gt;', ''))
    return merged_results_with_gnps

def handle_group(group):
    to_drop = set()
    if len(group) <= 1:
        return to_drop
    group_mols = {}
    for index, row in group.iterrows():
        group_mols[index] = Chem.MolFromSmiles(row['Smiles'])
    for index1, row1 in group.iterrows():
        for index2, row2 in group.iterrows():
            if index1 == index2:
                continue
            if index1 in to_drop or index2 in to_drop:
                continue
            mol1 = group_mols[index1]
            mol2 = group_mols[index2]
            if mol1 is None or mol2 is None:
                continue
            if mol1.HasSubstructMatch(mol2) and mol2.HasSubstructMatch(mol1):
                if row1['MQScore'] > row2['MQScore']:
                    to_drop.add(index2)
                elif row1['MQScore'] < row2['MQScore']:
                    to_drop.add(index1)
                else:
                    if index1 not in to_drop:
                        to_drop.add(index2)
    return to_drop

def remove_duplicates(with_structure, num_processes=5):
    # with_structure = get_with_structure(df)
    to_drop = set()
    grouped = with_structure.groupby(['#Scan#', 'SpectrumFile'])

    if num_processes == 1:
        for name, group in grouped:
            to_drop = to_drop.union(handle_group(group))
        return with_structure.drop(to_drop)
    

    pool = mp.Pool(num_processes)
    to_drop_sets = pool.map(handle_group, [group for name, group in grouped])
    pool.close()
    pool.join()
    for drop_set in to_drop_sets:
        to_drop = to_drop.union(drop_set)
    return with_structure.drop(to_drop)


def remove_knowns(library_search, knowns, topk=0):
    # for each result of the library search, check if the matched structure is the same as its true structure
    grouped = library_search.groupby(['#Scan#','SpectrumFile'])
    to_drop = []
    for name, group in grouped:
        # find the row in knowns that matches the scan number and spectrum file if it exists
        known = knowns[(knowns['Scan'] == name[0]) & (knowns['mgf_path'].apply(lambda x: is_file_match(x, name[1])))]
        if known.shape[0] == 0:
            continue
        true_struct = known['SMILES'].values[0]
        true_mol = Chem.MolFromSmiles(true_struct)

        for index, row in group.iterrows():
            # check if the structure is the same as the true structure
            try:
                mol = Chem.MolFromSmiles(row['Smiles'])
            except:
                print("doesn't have structure")
                continue
            if mol is None:
                continue
            if mol.HasSubstructMatch(true_mol) and true_mol.HasSubstructMatch(mol):
                to_drop.append(index)
    library_search = library_search.drop(to_drop)
    # sort by mqscore
    library_search_results = library_search.sort_values(by='MQScore', ascending=False)
    # only with structure
    library_search_results = get_with_structure(library_search_results)
    library_search_results.reset_index(drop=True, inplace=True)
    # remove duplicates
    library_search_results = remove_duplicates(library_search_results)
    # filter top k results
    if topk > 0:
        library_search_results = library_search_results.groupby(['#Scan#','SpectrumFile']).head(topk)
    
    return library_search_results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Test write out a file.')
    parser.add_argument('library_search_results')
    parser.add_argument('knowns')
    parser.add_argument('output_filename')
    parser.add_argument('--topk', type=int,  default=0, help='Top k results to consider, default is no filtering')
    args = parser.parse_args()

    library_search = pd.read_csv(args.library_search_results, sep='\t')
    knowns = pd.read_csv(args.knowns)
    result = remove_knowns(library_search, knowns, args.topk)
    result.to_csv(args.output_filename, sep='\t', index=False)