import pandas as pd
import numpy as np
import sys
import os
from tqdm import tqdm
import re
from rdkit import Chem
import json
import argparse
import multiprocessing
from functools import partial
import utils as utils

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Check the format of the mgf file to ensure it is in the correct format')
    parser.add_argument('mgf_file', help='Path to the mgf file')
    args = parser.parse_args()

    to_lower = lambda x: x.lower()
    bad_formats = ['na', 'n/a', 'none', 'n.a', 'n.a.', 'none', '', 'null', 'nan']
    to_contain = ["compound_type", "reaction_name", "filename_contains", "functional_group_name", "educt_smarts", "reaction_smarts", "delta_mz"]
    to_contain = set(to_contain)
    # read the mgf file
    mgf = utils._read_mgf(args.mgf_file)

    # check if the MGF has the SCANS field and they are correct format
    if 'scans' not in mgf.columns:
        print('SCANS field is missing')
        sys.exit(1)
    
    try:
        SCANS = mgf['scans'].values.astype(int) 
    except:
        print('SCANS should be integer')
        sys.exit(1)
    # sort the SCANS
    SCANS = np.sort(SCANS)
    # check if SCANS are unique and starting from 1 and incrementing by 1
    if not len(SCANS) == len(np.unique(SCANS)):
        print('SCANS are not unique')
        sys.exit(1)
    if not SCANS[0] == 1:
        print('SCANS should start from 1')
        sys.exit(1)
    if not all(np.diff(SCANS) == 1):
        print('SCANS should increment by 1')
        sys.exit(1)
    
    # check PEPMASS field is present and it floats
    if 'pepmass' not in mgf.columns:
        print('PEPMASS field is missing')
        sys.exit(1)
    
    PEPMASS = mgf['pepmass']
    try:
        PEPMASS = PEPMASS.apply(lambda x: (float(x[0]), x[1]))
    except:
        print('PEPMASS should be float')
        # print which row is not float
        for i in range(len(PEPMASS)):
            try:
                float(PEPMASS[i])
            except:
                print(f'PEPMASS at row {i} is not float', PEPMASS[i])
                sys.exit(1)
        sys.exit(1)
    
    if 'online_reactivity' not in mgf.columns:
        print('WARNING!! ONLINE_REACTIVITY field is missing')
    
    else:
        ONLINE_REACTIVITY = mgf['online_reactivity']
        # only get non-empty ONLINE_REACTIVITY
        ONLINE_REACTIVITY = ONLINE_REACTIVITY[ONLINE_REACTIVITY.notnull()]
        ONLINE_REACTIVITY = ONLINE_REACTIVITY[ONLINE_REACTIVITY.apply(to_lower).isin(bad_formats) == False]
        # check if ONLINE_REACTIVITY is in json format
        try:
            ONLINE_REACTIVITY = ONLINE_REACTIVITY.apply(json.loads)
        except:
            # print which row is not json
            print('ONLINE_REACTIVITY should be in json format')
            for i in range(len(ONLINE_REACTIVITY)):
                try:
                    json.loads(ONLINE_REACTIVITY[i])
                except:
                    print(f'ONLINE_REACTIVITY at row {i} is not in json format', ONLINE_REACTIVITY[i])
                    sys.exit(1)
            sys.exit(1)
        
        # check if the keys are in the correct format
        for i in range(len(ONLINE_REACTIVITY)):
            try:
                row = ONLINE_REACTIVITY.iloc[i]
                if not isinstance(row, list):
                    print('ONLINE_REACTIVITY should be a list')
                    sys.exit(1)
                for reactivity in row:
                    keys = set(reactivity.keys())
                    # if not all keys are present
                    if not keys.issuperset(to_contain):
                        print(f'ONLINE_REACTIVITY should contain {to_contain} but found {keys.difference(to_contain)} missing')
                        sys.exit(1)
            except:
                print(f'ONLINE_REACTIVITY at row {i} is not in correct format', ONLINE_REACTIVITY)
                sys.exit(1)
    
    if 'smiles' not in mgf.columns:
        print('WARNING!! SMILES field is missing, no known compounds is identified')
    else:
        SMILES = mgf['smiles']
        # if SMILES is present, check if it has any value
        count_good_smiles = 0
        # check if SMILES is in correct format
        for i in range(len(SMILES)):
            if SMILES[i].lower() in bad_formats:
                continue
            if not Chem.MolFromSmiles(SMILES[i]):
                print('SMILES is not in correct format', SMILES[i])
                sys.exit(1)
            count_good_smiles += 1
        if count_good_smiles == 0:
            print('WARNING!! No valid SMILES found, no known compounds is identified')

    
    print('The mgf file is in correct format')
    
