# adds adjusted ranking to the library search results

import sys
import argparse
import pandas as pd
from rdkit import Chem
import numpy as np
import utils as utils
from tqdm import tqdm
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from utils import is_file_match

def handle_group_adjusted_Score(group, possible_reaction_fgs_mols, experiment_reacted_fgs, and_order_strategy):
    group['adjusted_score'] = group['MQScore']
    group['reacted_with'] = ''
    for index, row in group.iterrows():
        try:
            mol = Chem.MolFromSmiles(row['Smiles'])
            mol = Chem.AddHs(mol)
        except:
            print("bad smiles")
            continue
        if mol is None:
            continue
        flag = 0
        for possible_reaction_fgs_mol in possible_reaction_fgs_mols:
            reacted_with = set()
            reaction_name_map = dict()
            reaction_type = set()
            for index2, reaction_mol in enumerate(possible_reaction_fgs_mol):
                if mol.HasSubstructMatch(reaction_mol['mol']):
                    reacted_with.add(reaction_mol['id'])
                    reaction_name_map[reaction_mol['id']] = reaction_mol['name']
                    reaction_type.add(reaction_mol['REACTION'])
            temp_reacted = ','.join([str(x) for x in reacted_with])
            temp_reacted_name = ','.join([reaction_name_map[x] for x in reacted_with])
            temp_reacted_type = ','.join([x for x in reaction_type])
            temp_reacted = temp_reacted + ' (' + temp_reacted_name + ')' + ' [' + temp_reacted_type + ']'
            if group.at[index, 'reacted_with'] == '':
                group.at[index, 'reacted_with'] = temp_reacted 
            else:
                group.at[index, 'reacted_with'] = group.at[index, 'reacted_with'] + " & " + temp_reacted

            if and_order_strategy and reacted_with == experiment_reacted_fgs:
                flag += 1
            
            elif not and_order_strategy and len(reacted_with.intersection(experiment_reacted_fgs)) > 0:
                flag += 1
        
        if flag > 0 and flag == len(possible_reaction_fgs_mols):
            group.at[index, 'adjusted_score'] = 1 + group.at[index, 'MQScore']
    
    # sort by adjusted score
    group = group.sort_values(by='adjusted_score', ascending=False)
    group['rank_after'] = np.arange(1, len(group) + 1)
    return group


def add_adjusted_score(df, reactants, all_fgs, and_order_strategy=True):
    grouped = df.groupby(['#Scan#', 'SpectrumFile'])
    new_df = pd.DataFrame()
    for name, group in tqdm(grouped):
        reacted = reactants[(reactants['Scan'] == name[0]) & (reactants['mgf_path'].apply(lambda x: is_file_match(x, name[1])))]
        if reacted.empty:
            print(name, (reactants['Scan'] == name[0]).sum())
            continue
        reactions = reacted['Reaction'].values[0]
        reacted_fgs = reacted['FGs'].values[0]
        if reacted_fgs == 'nan' or reacted_fgs == 'NONE' or reacted_fgs == '' or reacted_fgs is None:
            continue
        reacted_fgs = set([int(x) for x in reacted_fgs.split(',')])

        temp_df = group.copy()
        temp_df['rank_before'] = np.arange(1, len(temp_df) + 1)

        reaction_array = reactions.split(',')
        possible_reaction_fgs_mols = []
        for reaction in reaction_array:
            possible_reaction_fgs_mols.append([fg for fg in all_fgs if fg['REACTION'] in reaction])

        temp_df = handle_group_adjusted_Score(temp_df, possible_reaction_fgs_mols, reacted_fgs, and_order_strategy)
        if new_df.empty:
            new_df = temp_df
        else:
            new_df = pd.concat([new_df, temp_df])
    return new_df


def main(library_search_results, reactants, fgs, and_order_strategy=True):
    # convert fgs to mols
    all_fgs = []
    for index, row in fgs.iterrows():
        smarts = utils.split_smarts(row['SMARTS'])
        id = row['Substructure ID']
        reaction = row['REACTION']
        f_Name = row['Functional Group']
        mols = [Chem.MolFromSmarts(smart) for smart in smarts]
        for mol in mols:
            all_fgs.append({'id': id, 'REACTION': reaction, 'mol': mol, 'name': f_Name})
    
    # add adjusted score
    result = add_adjusted_score(library_search_results, reactants, all_fgs, and_order_strategy)
    return result



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Test write out a file.')
    parser.add_argument('library_search_results')
    parser.add_argument('reactants')
    parser.add_argument('fgs')
    parser.add_argument('adjusted_ranking_path')
    parser.add_argument('--reactants_only_path', default=None, help='Path to save the reactants only')
    parser.add_argument('--and_order_strategy', action='store_true', help='Use the AND order strategy, default is OR')
    args = parser.parse_args()

    library_search = pd.read_csv(args.library_search_results, sep='\t')
    reactants = pd.read_csv(args.reactants)
    fgs = pd.read_csv(args.fgs)
    and_order_strategy = args.and_order_strategy
    result = main(library_search, reactants, fgs, and_order_strategy)
    result.to_csv(args.adjusted_ranking_path, sep='\t', index=False)

    # select only the reactants
    if args.reactants_only_path is not None:
        grouped = result.groupby(['#Scan#', 'SpectrumFile'])
        reactants_only = pd.DataFrame()
        for name, group in grouped:
            reacted = reactants[(reactants['Scan'] == name[0]) & (reactants['mgf_path'].apply(lambda x: is_file_match(x, name[1])))]
            if reacted.empty:
                continue
            else:
                if reactants_only.empty:
                    reactants_only = group
                else:
                    reactants_only = pd.concat([reactants_only, group])
        reactants_only.to_csv(args.reactants_only_path, sep='\t', index=False)