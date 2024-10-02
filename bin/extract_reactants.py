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

def extract_reactants(df, fgs, mgf_path):
    reactants = pd.DataFrame(columns=['mgf_path','FGs','Scan','Reaction', 'FEATURE_ID'])
    reactant_scans = df[df['online_reactivity'].notnull()]
    for index, row in reactant_scans.iterrows():
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
        
        scan_fgs = list(set(scan_fgs))
        scan_fgs = ','.join([str(fg) for fg in scan_fgs])
        reactions = list(set(reactions))
        if len(reactions) == 0:
            continue
        reactions = ','.join(reactions)

        if len(reactants) == 0:
            reactants = pd.DataFrame({'mgf_path': mgf_path,
                                   'FGs': scan_fgs,
                                   'Scan': row['scans'],
                                   'Reaction': reactions,
                                   'FEATURE_ID': row.get('feature_id', None)}, index=[0])
        else:
            reactants = pd.concat([reactants, pd.DataFrame({'mgf_path': mgf_path,
                                   'FGs': scan_fgs,
                                   'Scan': row['scans'],
                                   'Reaction': reactions,
                                   'FEATURE_ID': row.get('feature_id', None)}, index=[len(reactants)]),])
    
    return reactants

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts reactants from a dataset')
    parser.add_argument('input', help='Path to the input file')
    parser.add_argument('fgs', help='Path to the fgs file')
    parser.add_argument('output', help='Path to the output file')
    args = parser.parse_args()

    df = utils._read_mgf(args.input)
    fgs = pd.read_csv(args.fgs)
    mgf_path = args.input
    reactants = extract_reactants(df, fgs, mgf_path)
    file_name = os.path.basename(args.input)
    file_name = file_name.split('.')[0]
    reactants.to_csv(os.path.join(args.output, file_name + '_reactants.csv'), index=False)