import argparse
import utils as utils
import os
import pandas as pd
from rdkit import Chem
from multiprocessing import Pool
import numpy as np
from functools import partial

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='convert_knowns_to_correct_format')
    parser.add_argument('--mgfs', type=str, help='Path to the root of mgfs', default='data/mgfs')
    parser.add_argument('--knowns', type=str, help='Path to the known compounds', default='data/knowns')
    parser.add_argument('--output', type=str, help='Path to the output folder', default='data/cleaned')
    args = parser.parse_args()

    # Load the libraries
    mgf_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(args.mgfs) for f in filenames if f.endswith('.mgf')]
    known_files = [f for f in os.listdir(args.knowns) if f.endswith('.csv')]

    reaction_types = ['AQC', 'cysteine', 'hydroxylamine']

    columns = ['mgf_path','Name','SMILES','FGs','Scan','Adduct', 'Reaction']
    result_df = pd.DataFrame(columns=columns)

    skip_count = 0
    for file in known_files:
        file_df = pd.DataFrame()
        knowns = pd.read_csv(os.path.join(args.knowns, file))
        reaction_type = None
        for rt in reaction_types:
            rt_lower = rt.lower()
            file_lower = file.lower()
            if rt_lower in file_lower:
                reaction_type = rt
                break
        if reaction_type is None:
            raise ValueError(f"Could not find reaction type in the name of the file {file}")
        
        for index, row in knowns.iterrows():
            folder_name = row['Sample name']
            sub_folder_name = row['Sample']
            reaction_name = folder_name.split('_')[1]
            if reaction_name != reaction_type and reaction_name != None:
                raise ValueError(f"Reaction name {reaction_name} does not match reaction type {reaction_type}")
            folder_name = folder_name.split('_')[0]
            if sub_folder_name == folder_name:
                full_mfg_path = os.path.join(args.mgfs, folder_name, reaction_name, f'{sub_folder_name}_{reaction_name}.mgf')
            else:
                full_mfg_path = os.path.join(args.mgfs, folder_name, reaction_name, sub_folder_name, f'{sub_folder_name}_{reaction_name}.mgf')
            if mgf_files.count(full_mfg_path) == 0:
                raise ValueError(f"Could not find mgf file {full_mfg_path}, error in naming")
            
            scan_ids = []
            if 'Educt ID H+' in row:
                scan = row['Educt ID H+'] 
                if type(scan) == float or type(scan) == int:
                    if np.isnan(scan):
                        scan = None
                    else:
                        scan = str(int(scan))

                if scan is None or scan == '':
                    skip_count += 1
                else:
                    scan = scan.replace(' ', '')
                    for s in scan.split(';'):
                        scan_ids.append((s, 'H+'))
            
            if 'Educt ID Na+' in row:
                scan = row['Educt ID Na+'] 
                if type(scan) == float or type(scan) == int:
                    if np.isnan(scan):
                        scan = None
                    else:
                        scan = str(int(scan))
                
                if scan is None or scan == '':
                    skip_count += 1
                else:
                    scan = scan.replace(' ', '')
                    for s in scan.split(';'):
                        scan_ids.append((s, 'Na+'))
            
            for scan_id, adduct in scan_ids:
                adjusted_row = {
                    'mgf_path': full_mfg_path,
                    'Name': row['Name'],
                    'SMILES': row['SMILE'],
                    'FGs': row['FGs'],
                    'Scan': scan_id,
                    'Adduct': adduct,
                    'Reaction': reaction_name
                }
                if len(file_df) == 0:
                    file_df = pd.DataFrame(adjusted_row, index=[0])
                else:
                    file_df = pd.concat([file_df, pd.DataFrame(adjusted_row, index=[0])], ignore_index=True)
            
        
        # remove duplicates from the file
        print(f"Before removing duplicates from {file}, it has {len(file_df)} scans and {len(knowns)} features")
        file_df = file_df.drop_duplicates(subset=['mgf_path', 'Scan', 'Adduct'], keep='first')
        print(f"After removing duplicates from {file}, it has {len(file_df)} scans and {len(file_df.groupby('Name'))} features")
        result_df = pd.concat([result_df, file_df], ignore_index=True)
    
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    result_df.to_csv(os.path.join(args.output, 'knowns.csv'), index=False)
