import os
import pandas as pd
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='merge_all_fgs_to_one_file')
    parser.add_argument('input', help='Path to the input file')
    parser.add_argument('output', help='Path to the output file')
    parser.add_argument('--images_dir', help='Path to the directory to save images', default='images')
    args = parser.parse_args()
    
    df = pd.read_csv(args.input)
    # unique based on SMARTS and keep the first one
    df = df.drop_duplicates(subset='SMARTS', keep='first')

    # reindex
    df = df.reset_index(drop=True)

    # add 'Substructure ID' column
    df['Substructure ID'] = np.arange(1, len(df)+1)

    # if the images directory does not exist, create it
    if not os.path.exists(args.images_dir):
        os.makedirs(args.images_dir)

    # draw the structure
    for i in range(len(df)):
        smart = df.loc[i, 'SMARTS']
        mol = Chem.MolFromSmarts(smart)
        img = Draw.MolToImage(mol)
        img.save(os.path.join(args.images_dir, f'FG_{i+1}.png'))
    
    df['image_dir'] = [os.path.join(args.images_dir, f'FG_{i+1}.png') for i in range(len(df))]
    df.to_csv(args.output, index=False)