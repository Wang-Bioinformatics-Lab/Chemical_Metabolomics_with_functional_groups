import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import argparse
import os
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
from utils import is_file_match
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial


def enrich_group_with_tanimoto(grouped_scan, knowns):
    name, group = grouped_scan
    known = knowns[(knowns['Scan'] == name[0]) & (knowns['mgf_path'].apply(lambda x: is_file_match(x, name[1])))]
    if known.shape[0] == 0:
        return pd.DataFrame()
    true_struct = known['SMILES'].values[0]
    true_mol = Chem.MolFromSmiles(true_struct)

    temp_group = group.copy()
    # rest index
    temp_group.reset_index(drop=True, inplace=True)
    temp_group['tanimoto'] = np.nan
    temp_group['FEATURE_ID'] = known['FEATURE_ID'].values[0]

    for index, row in temp_group.iterrows():
        # check if the structure is the same as the true structure
        mol = Chem.MolFromSmiles(row['Smiles'])
        if mol is None:
            continue
            # calculate the tanimoto similarity
        fpgen = AllChem.GetRDKitFPGenerator()
        fps = [fpgen.GetFingerprint(x) for x in [true_mol, mol]]
        tanimoto = DataStructs.TanimotoSimilarity(fps[0],fps[1])
        temp_group.loc[index, 'tanimoto'] = tanimoto
    
    temp_group = temp_group[temp_group['tanimoto'].isna() == False]
    return temp_group

def enrich_with_tanimoto(library_search, knowns):
    grouped = library_search.groupby(['#Scan#','SpectrumFile'])
    num_processes = min(5, os.cpu_count())

    with Pool(num_processes) as p:
        func = partial(enrich_group_with_tanimoto, knowns=knowns)
        results = list(tqdm(p.imap(func, grouped), total=len(grouped)))
        p.close()
        p.join()
    
    enriched_library_search = pd.concat(results, ignore_index=True)
    return enriched_library_search

def get_analysis_result(enriched_library_search, col_name):
    analysis_results = pd.DataFrame()
    not_found_scans = []
    groups = enriched_library_search.groupby(['#Scan#', 'SpectrumFile'])
    for name, group in groups:
        row = {'scan': name[0], 'file': name[1], 'FEATURE_ID': group['FEATURE_ID'].values[0]}
        highest_tanimoto = group[col_name].max()
        row['highest_tanimoto'] = highest_tanimoto
        for k in [1, 5, 10]:
            group_before = group[group['rank_before'] <= k]
            group_after = group[group['rank_after'] <= k]
            average_tanimoto_before = group_before[col_name].mean()
            average_tanimoto_after = group_after[col_name].mean()
            row[f'top_{k}_average_tanimoto_before'] = average_tanimoto_before
            row[f'top_{k}_average_tanimoto_after'] = average_tanimoto_after
            row[f'top_{k}_best_tanimoto_before'] = group_before[col_name].max()
            row[f'top_{k}_best_tanimoto_after'] = group_after[col_name].max()
        try:
            row[f'rank_highest_tanimoto_before'] = group[group[col_name] == highest_tanimoto].iloc[0]['rank_before']
        except:
            print(f'Error in {name}, col_name: {col_name}, highest_tanimoto: {highest_tanimoto}')
            print(group)
            raise Exception(f'Error in {name}')
        row[f'rank_highest_tanimoto_after'] = group[group[col_name] == highest_tanimoto].iloc[0]['rank_after']
        analysis_results = pd.concat([analysis_results, pd.DataFrame([row])], ignore_index=True)
    return analysis_results


def calculate_improved_result(filtered_analysis_results):
    improvement_result = pd.DataFrame(columns=['group', 'improved', 'worsened', 'same'])

    for k in [1, 5, 10]:
        tanimoto_improved = filtered_analysis_results[f'top_{k}_average_tanimoto_after'] > filtered_analysis_results[f'top_{k}_average_tanimoto_before']
        tanimoto_worsened = filtered_analysis_results[f'top_{k}_average_tanimoto_after'] < filtered_analysis_results[f'top_{k}_average_tanimoto_before']
        same = filtered_analysis_results[f'top_{k}_average_tanimoto_after'] == filtered_analysis_results[f'top_{k}_average_tanimoto_before']

        changed = filtered_analysis_results[filtered_analysis_results[f'top_{k}_average_tanimoto_after'] != filtered_analysis_results[f'top_{k}_average_tanimoto_before']]
        improvement_result = pd.concat([improvement_result, pd.DataFrame([{'group': f'top {k}',
                                                                           'average_before': changed[f'top_{k}_average_tanimoto_before'].mean(),
                                                                            'average_after': changed[f'top_{k}_average_tanimoto_after'].mean(),
                                                                            'improved': sum(tanimoto_improved),
                                                                            'worsened': sum(tanimoto_worsened), 
                                                                            'same': sum(same)}])], ignore_index=True)        

    highest_rank_improved = filtered_analysis_results[f'rank_highest_tanimoto_after'] < filtered_analysis_results[f'rank_highest_tanimoto_before']
    highest_rank_worsened = filtered_analysis_results[f'rank_highest_tanimoto_after'] > filtered_analysis_results[f'rank_highest_tanimoto_before']
    highest_rank_same = filtered_analysis_results[f'rank_highest_tanimoto_after'] == filtered_analysis_results[f'rank_highest_tanimoto_before']
    changed = filtered_analysis_results[filtered_analysis_results[f'rank_highest_tanimoto_after'] != filtered_analysis_results[f'rank_highest_tanimoto_before']]
    improvement_result = pd.concat([improvement_result,
                                    pd.DataFrame([{'group': 'Highest tanimoto rank',
                                                    'average_before': changed['rank_highest_tanimoto_before'].mean(),
                                                    'average_after': changed['rank_highest_tanimoto_after'].mean(), 
                                                    'improved': sum(highest_rank_improved), 
                                                    'worsened': sum(highest_rank_worsened),
                                                        'same': sum(highest_rank_same)}])], ignore_index=True)

    improvement_result['total'] = improvement_result['improved'] + improvement_result['worsened'] + improvement_result['same']

    # add percentage
    for col in ['improved', 'worsened', 'same']:
        improvement_result[f'{col}_percentage'] = (improvement_result[col] / improvement_result['total'] * 100)
        improvement_result[f'{col}_percentage'] = improvement_result[f'{col}_percentage'].apply(lambda x: round(x, 2))
    return improvement_result

def calculate_improved_result_best(filtered_analysis_results):
    improvement_result = pd.DataFrame(columns=['group', 'improved', 'worsened', 'same'])
    for k in [1, 5, 10]:
        tanimoto_improved = filtered_analysis_results[f'top_{k}_best_tanimoto_after'] > filtered_analysis_results[f'top_{k}_best_tanimoto_before']
        tanimoto_worsened = filtered_analysis_results[f'top_{k}_best_tanimoto_after'] < filtered_analysis_results[f'top_{k}_best_tanimoto_before']
        same = filtered_analysis_results[f'top_{k}_best_tanimoto_after'] == filtered_analysis_results[f'top_{k}_best_tanimoto_before']

        changed = filtered_analysis_results[filtered_analysis_results[f'top_{k}_best_tanimoto_after'] != filtered_analysis_results[f'top_{k}_best_tanimoto_before']]
        improvement_result = pd.concat([improvement_result, pd.DataFrame([{'group': f'top {k} (best)',
                                                                           'best_before': changed[f'top_{k}_best_tanimoto_before'].mean(),
                                                                            'best_after': changed[f'top_{k}_best_tanimoto_after'].mean(),
                                                                            'improved': sum(tanimoto_improved),
                                                                            'worsened': sum(tanimoto_worsened), 
                                                                            'same': sum(same)}])], ignore_index=True)
    

    improvement_result['total'] = improvement_result['improved'] + improvement_result['worsened'] + improvement_result['same']

    # add percentage
    for col in ['improved', 'worsened', 'same']:
        improvement_result[f'{col}_percentage'] = (improvement_result[col] / improvement_result['total'] * 100)
        improvement_result[f'{col}_percentage'] = improvement_result[f'{col}_percentage'].apply(lambda x: round(x, 2))
    return improvement_result

def draw_scatter(analysis_results, name = 1, output_path='.'):
    color_worsened = '#f08511'
    color_improved = '#3b7e23'
    font_size = 28
    font_size2 = 22
    tanimoto_improved = analysis_results[analysis_results[f'{name}_after'] > analysis_results[f'{name}_before']]
    tanimoto_worsened = analysis_results[analysis_results[f'{name}_after'] < analysis_results[f'{name}_before']]
    tanimoto_same = analysis_results[analysis_results[f'{name}_after'] == analysis_results[f'{name}_before']]
    # merged = pd.concat([tanimoto_worsened, tanimoto_improved, tanimoto_same], ignore_index=True)

    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    ax.scatter(tanimoto_improved[f'{name}_before'], tanimoto_improved[f'{name}_after'], color=color_improved, label='Improved')
    ax.scatter(tanimoto_worsened[f'{name}_before'], tanimoto_worsened[f'{name}_after'], color=color_worsened, label='Worsened')
    ax.scatter(tanimoto_same[f'{name}_before'], tanimoto_same[f'{name}_after'], color='gray', label='Same')
    ax.set_xlabel('Tanimoto Score w/o FG', fontsize=font_size2)
    ax.set_ylabel('Tanimoto Score w/ FG', fontsize=font_size2)
    # fix ticks font size
    ax.tick_params(axis='both', which='major', labelsize=font_size2*0.8)

    name_count = name.split('_')[1]
    ax.set_title(f'Top {name_count} Tanimoto Score', fontsize=font_size, fontweight='bold')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xticks(np.arange(0, 1.1, 0.25))
    ax.set_yticks(np.arange(0, 1.1, 0.25))
    # save the figure
    fig.savefig(os.path.join(output_path, f'{name}.png'), bbox_inches='tight')
    fig.savefig(os.path.join(output_path, f'{name}.svg'), bbox_inches='tight')
    plt.close(fig)


    fig2, ax2 = plt.subplots(figsize=(5, 5), dpi=300)
    for index_j, j in enumerate(['before', 'after']):
        changed = analysis_results[analysis_results[f'{name}_after'] != analysis_results[f'{name}_before']]
        data = [changed[f'{name}_{j}'].values]
        if data[0].shape[0] == 0:
            continue
        ax2.violinplot(data, positions=[index_j*0.7], showmeans=True, showmedians=False)
        ax2.text(index_j*0.7, changed[f'{name}_{j}'].mean(), f'{changed[f"{name}_{j}"].mean():.2f}', ha='center', va='bottom', fontsize=font_size)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.set_xticks([0, 1*0.7])
    ax2.set_xticklabels(['Before', 'After'], fontsize=font_size)
    # remove y ticks
    ax2.set_yticks([])
    fig2.savefig(os.path.join(output_path, f'{name}_violin.png'), bbox_inches='tight')
    fig2.savefig(os.path.join(output_path, f'{name}_violin.svg'), bbox_inches='tight')



def draw_scatter_rank(analysis_results, name = 1, output_path='.'):
    color_worsened = '#f08511'
    color_improved = '#3b7e23'
    font_size = 28
    font_size2 = 22
    tanimoto_improved = analysis_results[analysis_results[f'{name}_after'] > analysis_results[f'{name}_before']]
    tanimoto_worsened = analysis_results[analysis_results[f'{name}_after'] < analysis_results[f'{name}_before']]
    tanimoto_same = analysis_results[analysis_results[f'{name}_after'] == analysis_results[f'{name}_before']]
    # merged = pd.concat([tanimoto_worsened, tanimoto_improved, tanimoto_same], ignore_index=True)

    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    ax.scatter(tanimoto_improved[f'{name}_before'], tanimoto_improved[f'{name}_after'], color=color_worsened, label='Improved')
    ax.scatter(tanimoto_worsened[f'{name}_before'], tanimoto_worsened[f'{name}_after'], color=color_improved, label='Worsened')
    ax.scatter(tanimoto_same[f'{name}_before'], tanimoto_same[f'{name}_after'], color='gray', label='Same')
    ax.set_xlabel('Rank of Most Similar w/o FG', fontsize=font_size2)
    ax.set_ylabel('Rank of Most Similar w/ FG', fontsize=font_size2)

    ax.tick_params(axis='both', which='major', labelsize=font_size2*0.8)

    name_count = name.split('_')[1]
    ax.set_title(f'Rank of Most Similar', fontsize=font_size, fontweight='bold')
    # remove right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    max_rank = max(analysis_results[f'{name}_before'].max(), analysis_results[f'{name}_after'].max())

    ax.set_xticks(np.arange(0, max_rank + 0.8, 5))
    ax.set_yticks(np.arange(0, max_rank + 0.8, 5))
    # save the figure
    fig.savefig(os.path.join(output_path, 'rank_highest_tanimoto.png'), bbox_inches='tight')
    fig.savefig(os.path.join(output_path, 'rank_highest_tanimoto.svg'), bbox_inches='tight')
    plt.close(fig)
    
    fig2, ax2 = plt.subplots(figsize=(5, 5), dpi=300)
    for index_j, j in enumerate(['before', 'after']):
        changed = analysis_results[analysis_results[f'{name}_after'] != analysis_results[f'{name}_before']]
        data = [changed[f'{name}_{j}'].values]
        if data[0].shape[0] == 0:
            continue
        ax2.violinplot(data, positions=[index_j*0.7], showmeans=True, showmedians=False)
        ax2.text(index_j*0.7, changed[f'{name}_{j}'].mean(), f'{changed[f"{name}_{j}"].mean():.2f}', ha='center', va='bottom', fontsize=font_size)
    ax2.set_xticks([0, 1*0.7])
    ax2.set_xticklabels(['Before', 'After'], fontsize=font_size)
    # remove y ticks
    ax2.set_yticks([])
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    fig2.savefig(os.path.join(output_path, 'rank_tanimoto_violin.png'), bbox_inches='tight')
    fig2.savefig(os.path.join(output_path, 'rank_tanimoto_violin.svg'), bbox_inches='tight')

def main(search_results, knowns, output_path, images_path, filter_threshold=0.5):
    enriched_search_results = enrich_with_tanimoto(search_results, knowns)
    enriched_search_results.to_csv(os.path.join(output_path, 'enriched_search_results.csv'), index=False)
    analysis_results = get_analysis_result(enriched_search_results, 'tanimoto')
    for col in analysis_results.columns:
        if "tanimoto" in col and "rank" not in col:
            # convert to float
            analysis_results[col] = analysis_results[col].astype(float)
            analysis_results[col] = analysis_results[col].apply(lambda x: round(x, 4))
    
    analysis_results.to_csv(os.path.join(output_path, 'analysis_results.csv'), index=False)

    analysis_results = analysis_results[analysis_results['highest_tanimoto'] >= filter_threshold]
    improved_results = calculate_improved_result(analysis_results)
    draw_scatter(analysis_results, 'top_1_average_tanimoto', images_path)
    draw_scatter(analysis_results, 'top_5_average_tanimoto', images_path)
    draw_scatter(analysis_results, 'top_10_average_tanimoto', images_path)
    draw_scatter_rank(analysis_results, 'rank_highest_tanimoto', images_path)
    # add image path to the improved results
    improved_results['image_path'] = ['images/top_1_average_tanimoto.png', 'images/top_5_average_tanimoto.png', 'images/top_10_average_tanimoto.png', 'images/rank_highest_tanimoto.png']
    improved_results['violin_path'] = ['images/top_1_average_tanimoto_violin.png', 'images/top_5_average_tanimoto_violin.png', 'images/top_10_average_tanimoto_violin.png', 'images/rank_tanimoto_violin.png']
    improved_results.to_csv(os.path.join(output_path, 'improved_results.csv'), index=False)

    improved_results_best = calculate_improved_result_best(analysis_results)
    draw_scatter(analysis_results, 'top_1_best_tanimoto', images_path)
    draw_scatter(analysis_results, 'top_5_best_tanimoto', images_path)
    draw_scatter(analysis_results, 'top_10_best_tanimoto', images_path)
    improved_results_best['image_path'] = ['images/top_1_best_tanimoto.png', 'images/top_5_best_tanimoto.png', 'images/top_10_best_tanimoto.png']
    improved_results_best['violin_path'] = ['images/top_1_best_tanimoto_violin.png', 'images/top_5_best_tanimoto_violin.png', 'images/top_10_best_tanimoto_violin.png']
    # copy last row of improved results to the best results
    improved_results_best = pd.concat([improved_results_best, improved_results.iloc[-1:].copy()], ignore_index=True)
    improved_results_best.to_csv(os.path.join(output_path, 'improved_results_best.csv'), index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get analysis result')
    parser.add_argument('search_results', help='Path to the search results')
    parser.add_argument('knowns', help='Path to the knowns')
    parser.add_argument('output', help='Path to the output file')
    parser.add_argument('images', help='Path to the output file')
    parser.add_argument('--filter_threshold', help='Threshold to filter the search results', default=0.5, type=float)

    args = parser.parse_args()
    search_results = pd.read_csv(args.search_results, sep='\t')
    knowns = pd.read_csv(args.knowns)

    # make sure the output path exists
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    if not os.path.exists(args.images):
        os.makedirs(args.images)

    main(search_results, knowns, args.output, args.images, args.filter_threshold)
    
#     df = enrich_with_tanimoto(df, knowns)
