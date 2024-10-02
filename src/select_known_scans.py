import argparse
import utils as utils
import os
import pandas as pd
from rdkit import Chem
from multiprocessing import Pool
import numpy as np
from functools import partial
import json

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create libraries from MGF files')
    parser.add_argument('--mgf', type=str, help='Path to the MGF file', default='data/mgfs')
    parser.add_argument('--knowns', type=str, help='Path to the known compounds csv_file', default='data/cleaned/knowns.csv')
    parser.add_argument('--output', type=str, help='Path to the output folder', default='data/cleaned')
    args = parser.parse_args()

    knowns = pd.read_csv(args.knowns)
    mgf_files = knowns['mgf_path'].unique().tolist()
    old_mgfs = dict()
    for file in mgf_files:
        mgf = utils._read_mgf(file)
        old_mgfs[file] = mgf
        old_mgfs[file]['scans'] = old_mgfs[file]['scans'].astype(int)

    new_mgfs = {file: pd.DataFrame() for file in mgf_files}
    scan_counter = 1
    reindexed_knowns = pd.DataFrame()
    merged_mgfs = pd.DataFrame()

    fgs = pd.read_csv('data/cleaned/all_reactions_smarts.csv')
    all_fgs = []
    for index, row in fgs.iterrows():
        smarts = utils.split_smarts(row['SMARTS'])
        id = row['Substructure ID']
        reaction = row['REACTION']
        Functional_Group = row['Functional Group']
        smart = row['SMARTS']
        mols = [Chem.MolFromSmarts(smart) for smart in smarts]
        for mol in mols:
            all_fgs.append({'id': id, 'REACTION': reaction, 'mol': mol,
                            'Functional_Group': Functional_Group, 'SMARTS': smart})

    for index, row in knowns.iterrows():
        mgf_path = row['mgf_path']
        scan_number = row['Scan']
        scan_number = int(float(scan_number))
        k_Name = row['Name']
        k_SMILES = row['SMILES']
        k_Adduct = row['Adduct']
        scan = old_mgfs[mgf_path][old_mgfs[mgf_path]['scans'] == scan_number]
        new_scan = scan.copy()
        new_scan['Name'] = np.nan
        new_scan['SMILES'] = np.nan
        new_scan['Adduct'] = np.nan
        # change the dtype of the columns
        new_scan['Name'] = new_scan['Name'].astype('object')
        new_scan['SMILES'] = new_scan['SMILES'].astype('object')
        new_scan['Adduct'] = new_scan['Adduct'].astype('object')

        if new_scan.shape[0] == 0:
            name = row['Name']
            reaction = row['Reaction']
            adduct = row['Adduct']
            print("ERROR!", f"Could not find scan {scan_number} with name {name} and reaction {reaction} and adduct {adduct} in {mgf_path}")
            continue
        else:
            # if new row is known, update the online_reactivity
            for scan_index, scan_row in new_scan.iterrows():
                if 'online_reactivity' not in scan_row or pd.isna(scan_row['online_reactivity']):
                    continue
                online_reactivity = json.loads(scan_row['online_reactivity'])
                for educt_index, educt in enumerate(online_reactivity):
                    temp_educt = educt.copy()
                    educt_mol = Chem.MolFromSmarts(educt['educt_smarts'])
                    educt_smart = educt['educt_smarts']
                    # print(educt_mol)
                    flag = False
                    for fg in all_fgs:
                        if educt_smart == fg['SMARTS']:
                            temp_educt['functional_group_name'] = fg['Functional_Group']
                            # temp_educt['extra_info'] = temp_educt['reaction_name']
                            temp_educt['reaction_name'] = fg['REACTION']
                            flag = True
                            break
                    if not flag:
                        print("no fg found")
                        raise ValueError("no fg found", index, scan_index, educt_index)
                    online_reactivity[educt_index] = temp_educt
                scan_row_copy = scan_row.copy()
                scan_row_copy['online_reactivity'] = json.dumps(online_reactivity)
                scan_row_copy['Name'] = k_Name
                scan_row_copy['SMILES'] = k_SMILES
                scan_row_copy['Adduct'] = k_Adduct
                new_scan.loc[scan_index] = scan_row_copy
                # print new_scan columns and their dtypes
                # print(new_scan.dtypes)
                # print(scan_row_copy.dtypes)
                # raise ValueError("stop here")


        new_mgfs[mgf_path] = pd.concat([new_mgfs[mgf_path], new_scan])
        
        new_row = row.copy()
        new_row['old_scan'] = scan_number
        new_row['Scan'] = scan_counter
        
        reindexed_knowns = pd.concat([reindexed_knowns, pd.DataFrame([new_row])])

        if len(new_scan) > 1:
            print("ERROR!", f"Found multiple scans for {scan_number} in {mgf_path}")
            raise ValueError(f"Found multiple scans for {scan_number} in {mgf_path}")
        else:
            scan_dict = new_scan.iloc[0].to_dict()
            scan_dict['old_scan'] = scan_number
            scan_dict['scans'] = scan_counter
            merged_mgfs = pd.concat([merged_mgfs, pd.DataFrame([scan_dict])])
        scan_counter += 1

        

    df = pd.DataFrame()
    for file in mgf_files:
        item = {
            'file': os.path.basename(file),
            'old_size': old_mgfs[file].shape[0],
            'new_size': new_mgfs[file].shape[0]
        }
        df = pd.concat([df, pd.DataFrame([item])])
    df.to_csv(os.path.join(args.output, 'scan_comparison.csv'), index=False)
    
    for file in mgf_files:
        new_file_path = file.replace(args.mgf, '')
        if new_file_path[0] == '/':
            new_file_path = new_file_path[1:]
        new_file_path_old_format = os.path.join(args.output, 'old_format', new_file_path)
        if not os.path.exists(os.path.dirname(new_file_path_old_format)):
            os.makedirs(os.path.dirname(new_file_path_old_format))
        utils._write_mgf(new_mgfs[file], new_file_path_old_format)
    
    # updated_knowns = knowns.copy()
    # updated_knowns['old_scan'] = updated_knowns['Scan']
    # merged_mgfs = pd.DataFrame()
    # scan_number = 1
    # for file in mgf_files:
    #     mgf = new_mgfs[file]
    #     for index, row in mgf.iterrows():
    #         updated_knowns.loc[(updated_knowns['mgf_path'] == file) & (updated_knowns['Scan'] == row['scans']), 'Scan'] = scan_number
    #         row['old_scan'] = row['scans']
    #         row['scans'] = scan_number
    #         scan_number += 1
    #         merged_mgfs = pd.concat([merged_mgfs, pd.DataFrame([row])])
    reindexed_knowns['mgf_path'] = os.path.join(args.output, 'merged.mgf')
    
    reindexed_knowns.to_csv(os.path.join(args.output, 'reindexed_knowns.csv'), index=False)
    utils._write_mgf(merged_mgfs, os.path.join(args.output, 'merged.mgf'))


        