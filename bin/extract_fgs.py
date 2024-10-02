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

def extract_fgs(df):
    reactants = df[df['online_reactivity'].notnull()]

    fgs = pd.DataFrame(columns=['REACTION','Functional Group','SMARTS','delta_mz'])
    for index, row in reactants.iterrows():
        if row['online_reactivity'] is not None:
            try:
                reactivity = json.loads(row['online_reactivity'])
            except:
                # if reactivity is a list then it is already in the correct format
                if isinstance(row['online_reactivity'], list):
                    reactivity = row['online_reactivity']
                else:
                    print('Error in row', index)
                    continue
            for r in reactivity:
                # if reactivity not in fgs then add it
                if r['educt_smarts'] not in fgs['SMARTS'].values:
                    if len(fgs) == 0:
                        fgs = pd.DataFrame({'REACTION': r['reaction_name'], 
                                            'Functional Group': r['functional_group_name'],
                                            'SMARTS': r['educt_smarts'],
                                            'delta_mz': r['delta_mz']}, index=[0])
                    else:
                        fgs = pd.concat([fgs, pd.DataFrame({'REACTION': r['reaction_name'], 
                                                            'Functional Group': r['functional_group_name'],
                                                            'SMARTS': r['educt_smarts'],
                                                            'delta_mz': r['delta_mz']}, index=[len(fgs)]),])
    return fgs


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extracts functional groups and knowns from a dataset')
    parser.add_argument('input', help='Path to the input file')
    parser.add_argument('output', help='Path to the output file')
    args = parser.parse_args()

    df = utils._read_mgf(args.input)
    fgs = extract_fgs(df)
    file_name = os.path.basename(args.input)
    file_name = file_name.split('.')[0]
    fgs.to_csv(os.path.join(args.output, file_name + '_fgs.csv'), index=False)
    