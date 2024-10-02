import numpy as np
import os
import pickle
import pandas as pd
from rdkit.Chem import AllChem

from tqdm import tqdm 
import pyteomics
from pyteomics import mzml, auxiliary, mgf
import matplotlib.pyplot as plt
import argparse
import json
from rdkit import Chem
from rdkit import DataStructs
import re

def split_smarts(smarts):
    """
    This function is to fix an error in smart parsing of rdkit
    """
    
    # Regular expression to find patterns like [X,Y,Z,...]
    pattern = r'\[([^]]+)\]'
    
    # Function to process the content within brackets
    def process_options(options):
        elements = options.split(',')
        if 'H' in elements and len(elements) > 1:
            without_h = ','.join([elem for elem in elements if elem != 'H'])
            return [f'{without_h}', 'H']
        return [f'{options}']
    
    # Initialize output lists
    res = []
    
    last_end = 0
    # get the first match
    first_match = re.search(pattern, smarts)
    if first_match:
        start, end = first_match.span()
        starting_string = smarts[:start]
        solved_after = split_smarts(smarts[end:])
        inbetween = process_options(smarts[start+1:end-1])
        for item in inbetween:
            if solved_after is None or len(solved_after) == 0:
                solved_after = ['']
            for after in solved_after:
                res.append(starting_string + '[' + item +']' + after)
    else:
        res.append(smarts)
    
    return res

def _read_mgf(mgf_path: str) -> pd.DataFrame:
    msms_df = []
    with mgf.MGF(mgf_path) as reader:
        for spectrum in reader:
            try:
                d = spectrum['params']
                d['spectrum'] = np.array([spectrum['m/z array'],
                                        spectrum['intensity array']])
                if 'precursor_mz' not in d:
                    d['precursor_mz'] = d['pepmass'][0]
                else:
                    d['precursor_mz'] = float(d['precursor_mz'])
                msms_df.append(d)
            except Exception as e:
                print(e)

    msms_df = pd.DataFrame(msms_df)
    if 'precursor_mz' in msms_df.columns and 'scans' in msms_df.columns:
        msms_df['precursor_mz'] = msms_df['precursor_mz'].astype(float)
        msms_df['scans'] = msms_df['scans'].astype(int)
    return msms_df

def _write_mgf(msms_df: pd.DataFrame, mgf_path: str):
    specs = []
    for i, row in msms_df.iterrows():
        spectrum = {
            'params': row.drop('spectrum').to_dict(),
            'm/z array': row['spectrum'][0],
            'intensity array': row['spectrum'][1]
        }
        specs.append(spectrum)
    with open(mgf_path, 'w') as out:
        mgf.write(specs, out)

def calculate_tanimoto(smiles1, smiles2, method = "morgan", radius = 2, maxPath=3, fpSize=2048):
    """
    Calculate the tanimoto similarity between two smiles
    params:
    smiles1: str
    smiles2: str
    method: str, "morgan" or "rdkit"
    radius: int, radius for morgan fingerprints
    maxPath: int, maxPath for rdkit fingerprints
    fpSize: int, fpSize for rdkit fingerprints
    """
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 is None or mol2 is None:
        return 0
    if method == "morgan":
        fpgen = AllChem.GetMorganGenerator(radius=radius, fpSize=fpSize)
    else:
        fpgen = AllChem.GetRDKitFPGenerator(maxPath=maxPath, fpSize=fpSize)
    fps = [fpgen.GetFingerprint(x) for x in [mol1, mol2]]
    tanimoto = DataStructs.TanimotoSimilarity(fps[0],fps[1])
    return tanimoto

def is_file_match(name1, name2):
    if name1 in name2 or name2 in name1:
        return True
    return False