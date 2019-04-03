# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 10:55:15 2019

@author: hcji
"""

import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from PyCFMID.PyCFMID import cfm_id_database

result = pd.read_excel(os.path.join('Output', 'result.xlsx'))
for i in tqdm(range(len(result))):
    kegg = result['kegg'][i]
    formula = result['formula'][i]
    chebi = result['chebi'][i]
    pubchem = result['pubchem'][i]
    spectrum_file = os.path.join(os.getcwd(), 'Spectra', kegg + '.csv')
    spectrum_dataframe = pd.read_csv(spectrum_file)
    result_biodb = cfm_id_database(spectrum_dataframe, formula, database='biodb')
    for i in result_biodb['candidates'].index:
        x = result_biodb['candidates']['ChEBI'][i]
        if chebi in str(x):
            whtrue = i
    score_of_true = max(result_biodb['result']['Score'][ result_biodb['result']['ID']==whtrue])
    rank = len(np.where(result_biodb['result']['Score'] > score_of_true)[0]) + 1
    result['rank'][i] = rank
    