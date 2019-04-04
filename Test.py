# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 10:55:15 2019

@author: hcji
"""

import os
import pandas as pd
import numpy as np
from PyCFMID.PyCFMID import cfm_id_database, search_pubchem, parser_cfm_id
from joblib import Parallel, delayed
from tqdm import tqdm
import multiprocessing

num_cores = min(multiprocessing.cpu_count(), 10)
result = pd.read_excel(os.path.join('Output', 'result.xlsx'))

'''
os.mkdir('Candidate')
for i in tqdm(range(len(result))):
    formula = result['formula'][i]
    kegg = result['kegg'][i]
    candidate = search_pubchem(formula)
    candidate['ID'] = candidate.index
    candidate.to_csv(os.path.join(os.getcwd(), 'Candidate', str(kegg)+'.csv'), index=False)
''' 

def process_one_sample(i, database='biodb'):
    kegg = result['kegg'][i]
    formula = result['formula'][i]
    chebi = result['chebi'][i]
    pubchem = str(result['pubchem'][i])
    spectrum_file = os.path.join(os.getcwd(), 'Spectra', kegg + '.csv')
    spectrum_dataframe = pd.read_csv(spectrum_file) 
    candidate_file = os.path.join(os.getcwd(), 'Candidate', kegg + '.csv')
    try:
        os.mkdir(os.path.join('Input', str(kegg)))
    except:
        pass
    input_dir = os.path.join(os.getcwd(), 'Input', kegg)
    output_file = os.path.join(os.getcwd(), 'Output', str(kegg) + '.txt')
    if database == 'biodb':
        result_biodb = cfm_id_database(spectrum_dataframe, formula, database='biodb', input_dir=input_dir, output_file=output_file)
    else:
        if str(kegg) + '.txt' in os.listdir(os.path.join(os.getcwd(), 'Output')):
            result_biodb = {}
            result_biodb['result'] = parser_cfm_id(os.path.join(os.getcwd(), 'Output', str(kegg) + '.txt'))
            result_biodb['candidates'] = pd.read_csv(os.path.join(os.getcwd(), 'Candidate', str(kegg)+'.csv'))
        else:    
            result_biodb = cfm_id_database(spectrum_dataframe, formula, database=candidate_file, input_dir=input_dir, output_file=output_file)
    for j in result_biodb['candidates'].index:
        if database == 'biodb':
            x = result_biodb['candidates']['ChEBI'][j]
            if chebi in str(x):
                whtrue = j
        else:
            x = result_biodb['candidates']['PubChem'][j]
            if pubchem in str(x):
                whtrue = j            
    score_of_true = max(result_biodb['result']['Score'][ result_biodb['result']['ID']==whtrue])
    rank = len(np.where(result_biodb['result']['Score'] > score_of_true)[0]) + 1
    return rank

rank = Parallel(n_jobs=num_cores, verbose=5)(delayed(process_one_sample)(i) for i in range(len(result)))
result['rank'] = rank
result.to_excel(os.path.join('Output', 'result.xlsx'), index=False)

'''
for i in tqdm(range(len(result))):
    kegg = result['kegg'][i]
    formula = result['formula'][i]
    chebi = result['chebi'][i]
    pubchem = result['pubchem'][i]
    spectrum_file = os.path.join(os.getcwd(), 'Spectra', kegg + '.csv')
    spectrum_dataframe = pd.read_csv(spectrum_file)
    result_biodb = cfm_id_database(spectrum_dataframe, formula, database='biodb')
    for j in result_biodb['candidates'].index:
        x = result_biodb['candidates']['ChEBI'][j]
        if chebi in str(x):
            whtrue = j
    score_of_true = max(result_biodb['result']['Score'][ result_biodb['result']['ID']==whtrue])
    rank = len(np.where(result_biodb['result']['Score'] > score_of_true)[0]) + 1
    result['rank'][i] = rank
    if i % 20 == 19:
        result.to_excel(os.path.join('Output', 'result.xlsx'), index=False)
'''