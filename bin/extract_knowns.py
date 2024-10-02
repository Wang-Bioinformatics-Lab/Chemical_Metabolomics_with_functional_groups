import pandas as pd
import os
from utils import split_smarts, _read_mgf
from pyteomics import mzml, auxiliary, mgf
from rdkit import Chem
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
from rdkit.Chem import AllChem
from rdkit import DataStructs
import sys
import argparse
import utils as utils
import json

def extract_knowns(df, fgs, mgf_path):
    knowns = pd.DataFrame(columns=['mgf_path','Name','SMILES','FGs', 'true_fgs', 'Scan','Adduct','Reaction', 'FEATURE_ID'])
    known_scans = df[df['name'].notnull() & df['smiles'].notnull()]
    for index, row in known_scans.iterrows():
        try:
            reactivity = json.loads(row['online_reactivity'])
        except:
            print('Could not parse reactivity: ', row['online_reactivity'])
            continue
        scan_fgs = []
        reactions = []
        for r in reactivity:
            if r['educt_smarts'] in fgs['SMARTS'].values:
                # add the id of the functional group
                scan_fgs.append(fgs[fgs['SMARTS'] == r['educt_smarts']]['Substructure ID'].values[0])
                reactions.append(fgs[fgs['SMARTS'] == r['educt_smarts']]['REACTION'].values[0])
        
        try:
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol is None:
                print('Could not parse smiles: ', row['smiles'])
                continue
        except:
            print('Could not parse smiles: ', row['smiles'])
            continue
        try:
            mol = Chem.AddHs(mol)
        except:
            pass
        true_fgs = []
        for i, fg in fgs[fgs['REACTION'].isin(reactions)].iterrows():
            smarts = fg['SMARTS']
            smarts = split_smarts(smarts)
            mols = [Chem.MolFromSmarts(smart) for smart in smarts]
            for m in mols:
                if mol.HasSubstructMatch(m):
                    true_fgs.append(fg['Substructure ID'])
                    break
        
        true_fgs = list(set(true_fgs))
        scan_fgs = list(set(scan_fgs))

        true_fgs = ','.join([str(fg) for fg in true_fgs])
        scan_fgs = ','.join([str(fg) for fg in scan_fgs])

        reactions = list(set(reactions))
        reactions = ','.join(reactions)

        if len(knowns) == 0:
            knowns = pd.DataFrame({'mgf_path': mgf_path,
                                   'Name': row.get('name', None),
                                   'SMILES': row['smiles'],
                                   'FGs': scan_fgs,
                                   'true_fgs': true_fgs,
                                   'Scan': row['scans'],
                                   'Adduct': row.get('adduct', None),
                                   'Reaction': reactions,
                                   'FEATURE_ID': row.get('feature_id', None)
                                   }, index=[0])
        else:
            knowns = pd.concat([knowns, pd.DataFrame({'mgf_path': mgf_path,
                                   'Name': row.get('name', None),
                                   'SMILES': row['smiles'],
                                   'FGs': scan_fgs,
                                   'true_fgs': true_fgs,
                                   'Scan': row['scans'],
                                   'Adduct': row.get('adduct', None),
                                   'Reaction': reactions,
                                   'FEATURE_ID': row.get('feature_id', None)
                                   }, index=[len(knowns)]),])
    
    return knowns

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts knowns from a dataset')
    parser.add_argument('input', help='Path to the input file')
    parser.add_argument('fgs', help='Path to the fgs file')
    parser.add_argument('output', help='Path to the output file')
    args = parser.parse_args()

    df = utils._read_mgf(args.input)
    fgs = pd.read_csv(args.fgs)
    mgf_path = args.input
    knowns = extract_knowns(df, fgs, mgf_path)
    file_name = os.path.basename(args.input)
    file_name = file_name.split('.')[0]
    knowns.to_csv(os.path.join(args.output, file_name + '_knowns.csv'), index=False)