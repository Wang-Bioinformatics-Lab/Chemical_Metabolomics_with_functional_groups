"""
This code recieves a directory with ms files and the fg groups and returns a directory with mgf files.
"""

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

def fix_fgs(fgs):
    """
    Fix the functional groups in the dataframe as the SMARTS are not correctly formatted.
    return: dict with the functional groups
    """
    fgs_dict = dict()
    for i, row in fgs.iterrows():
        temp_dict = dict()
        temp_dict['REACTION'] = row['REACTION']
        temp_dict['Functional Group'] = row['Functional Group']
        temp_dict['delta_mz'] = row['delta_mz']
        temp_dict['SMARTS'] = row['SMARTS']
        smarts = utils.split_smarts(row['SMARTS'])
        temp_dict['mols'] = [Chem.MolFromSmarts(smart) for smart in smarts]
        fgs_dict[row['Substructure ID']] = temp_dict
    return fgs_dict

def parse_file_to_sections(path, file):
    with open(os.path.join(path, file), "r") as f:
        lines = f.readlines()
        
        sections = []
        # section based on empty line
        section = []
        for line in lines:
            if line == "\n":
                sections.append(section)
                section = []
            else:
                # trim the line
                # remove the newline character
                line = line.strip()
                # remove the leading and trailing quote characters
                if line[0] == "\"" or line[0] == "\'":
                    line = line[1:]
                if line[-1] == "\"" or line[-1] == "\'":
                    line = line[:-1]
                section.append(line)
        if len(section) > 0:
            sections.append(section)
        
        return sections

def parse_key(key):
    if key.startswith(">"):
        key = key.split(" ")[0][1:]
    return key

def get_peaks(section):
    mzs = []
    intensities = []
    for line in section:
        mz, intensity = line.split(" ")
        mzs.append(float(mz))
        intensities.append(float(intensity))
    
    # sort the peaks by mz
    indices = np.argsort(mzs)
    mzs = np.array(mzs)[indices]
    intensities = np.array(intensities)[indices]
    peaks = [mzs, intensities]
    peaks = np.array(peaks)
    return peaks

def merge_peaks(array_of_peaks):
    """
    Merge the peaks from all the sections (coolision energies)
    """
    bin_size = 6
    merged_intensities = dict()

    for peaks in array_of_peaks:
        mzs1 = peaks[0]
        intensities1 = peaks[1]
        for i in range(len(mzs1)):
            mz = mzs1[i]
            intensity = intensities1[i]
            mz = round(mz, bin_size)
            if mz in merged_intensities:
                merged_intensities[mz].append(intensity)
            else:
                merged_intensities[mz] = [intensity]
    
    mzs = []
    intensities = []
    for mz, intensity in merged_intensities.items():
        mzs.append(mz)
        intensities.append(np.mean(intensity))
    
    indices = np.argsort(mzs)
    mzs = np.array(mzs)[indices]
    intensities = np.array(intensities)[indices]
    merged_peaks = [mzs, intensities]
    merged_peaks = np.array(merged_peaks)
    return merged_peaks

def get_online_ractivity(smiles, fgs):
    """
    Get the online reactivity of the molecule
    input:
    smiles: str, smiles of the molecule
    fgs: dict of functional group id to its information
    """
    online_reactivity = []
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    reaction_set = set()
    for fg in fgs:
        fg_data = fgs[fg]
        for fg_mol in fg_data['mols']:
            if mol.HasSubstructMatch(fg_mol):
                reaction_set.add(fg_data['REACTION'])
                break
                # # in case we can seprate them, we can just add the reaction to the list here
                # temp_reaction = {
                #     "reaction_name": fg_data['REACTION'],
                #     "filename_contains": fg_data['REACTION'],
                #     "functional_group_name": fg_data['Functional Group'],
                #     "compound_type": "Educt",
                #     "educt_smarts": fg_data['SMARTS'],
                #     "reaction_smarts": "N/A",
                #     "delta_mz": "N/A",
                #     "linked_ids": [],
                # }
                # online_reactivity.append(temp_reaction)
                # break
    
    for reaction in reaction_set:
        for fg in fgs:
            fg_data = fgs[fg]
            if fg_data['REACTION'] == reaction:
                temp_reaction = {
                    "reaction_name": fg_data['REACTION'],
                    "filename_contains": fg_data['REACTION'],
                    "functional_group_name": fg_data['Functional Group'],
                    "compound_type": "Educt",
                    "educt_smarts": fg_data['SMARTS'],
                    "reaction_smarts": "N/A",
                    "delta_mz": "N/A",
                    "linked_ids": [],
                }
                online_reactivity.append(temp_reaction)
    return online_reactivity

def handle_file(file, path, fgs_dict):
    temp_data = dict()
    sections = parse_file_to_sections(path, file)
    for item in sections[0]:
        if not item.startswith(">"):
            continue
        vals = item.split(" ")
        key = vals[0]
        value = " ".join(vals[1:])
        key = parse_key(key)
        temp_data[key] = value

    array_of_peaks = []
    for section in sections[1:]:
        key = section[0].split(" ")[0]
        key = parse_key(key)
        peaks = get_peaks(section[1:])
        array_of_peaks.append(peaks)

    # make one merged spectrum for all peaks
    merged_peaks = merge_peaks(array_of_peaks)
    temp_data['spectrum'] = merged_peaks
    temp_data['mslevel'] = 2
    temp_data['filename'] = file
    temp_data['FEATURE_ID'] = file
    if 'precursor_mz' not in temp_data:
        temp_data['precursor_mz'] = temp_data['parentmass']
    if 'PEPMASS' not in temp_data:
        temp_data['PEPMASS'] = temp_data['precursor_mz']
    if 'Adduct' not in temp_data:
        temp_data['Adduct'] = temp_data['ionization']
    match = re.search(r"^([a-zA-Z]+)", file)
    if match:
        temp_data['group'] = match.group(1)
    else:
        temp_data['group'] = file[0:4]
    
    if 'NAME' not in temp_data and 'compound' in temp_data:
        temp_data['NAME'] = temp_data['compound']

    online_reactivity = get_online_ractivity(temp_data['smiles'], fgs_dict)
    temp_data['ONLINE_REACTIVITY'] = json.dumps(online_reactivity)
    return temp_data

def handle_files(files, path, fgs_dict):
    data = []
    for file in tqdm(files):
        temp_data = handle_file(file, path, fgs_dict)
        data.append(temp_data)
    return data

def ms_to_mgf(path, files, fg_groups):
    df = pd.DataFrame()
    peaks_keys = ['ms2peaks', 'collision']

    num_cores = 10
    chunks = np.array_split(files, num_cores)
    pool = multiprocessing.Pool(num_cores)
    results = pool.map(partial(handle_files, path=path, fgs_dict=fg_groups), chunks)
    pool.close()
    pool.join()

    for result in results:
        if len(df) == 0:
            df = pd.DataFrame(result)
        else:
            df = pd.concat([df, pd.DataFrame(result)], ignore_index=True)
    return df
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert ms files to mgf files')
    parser.add_argument('ms_dir', help='Directory with ms files')
    parser.add_argument('output_mgf_dir', help='Directory to save mgf files')
    parser.add_argument('fg_groups', help='File with fg groups')
    parser.add_argument('--number_of_output_files', help='Number of output files', default=1, type=int)
    args = parser.parse_args()

    fgs = pd.read_csv(args.fg_groups)
    fgs = fix_fgs(fgs)

    files = os.listdir(args.ms_dir)

    result = ms_to_mgf(args.ms_dir, files, fgs)

    # check if the output directory exists
    if not os.path.exists(args.output_mgf_dir):
        os.makedirs(args.output_mgf_dir)

    num_sections = int(args.number_of_output_files)
    section_size = len(result) // num_sections
    section_start = 0
    for i in tqdm(range(num_sections)):
        section_end = section_start + section_size
        if i == num_sections - 1:
            section_end = len(result)
        section = result.iloc[section_start:section_end]
        # make section a new dataframe
        section = pd.DataFrame(section)
        # reset the index
        section.reset_index(drop=True, inplace=True)
        section['SCANS'] = np.arange(1, len(section) + 1)
        section_start = section_end
        utils._write_mgf(section, os.path.join(args.output_mgf_dir, f"output_{i}.mgf"))
