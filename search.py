#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 14:17:54 2020

@author: warren

Search ChEMBL DB
"""

# Import library to connect Postgresql
import psycopg2
import pandas as pd
from rdkit import Chem

# Connect to chembl_26 DB
try:
    connect_str = "chembl_26"
    
    # use our connection values to establish a connection
    conn = psycopg2.connect(database = connect_str)
    
    # create a psycopg2 cursor that can execute queries
    cursor = conn.cursor()
    
    # query ChemBL database
    print('PostgreSQL database version:')
    cursor.execute('SELECT version()')
    
    # display the PostgreSQL database server version
    db_version = cursor.fetchone()
    print(db_version)
except Exception as e:
    print("Uh oh, can't connect. Invalid dbname, user or password?")
    print(e)


def sql_query(chembl_id):
    # attempt SQL query with help from ChemBL blog
    qtext = ("""
    SELECT m.chembl_id AS compound_chembl_id,   
    s.canonical_smiles,
    d.year,
    a.description                   AS assay_description,   act.standard_type,
    act.standard_relation,   
    act.standard_value,
    act.standard_units
    FROM compound_structures s,   
    molecule_dictionary m,   
    compound_records r,
    docs d,
    activities act,   
    assays a,
    target_dictionary t
    WHERE s.molregno     = m.molregno 
    AND m.molregno       = r.molregno 
    AND r.record_id      = act.record_id
    AND r.doc_id         = d.doc_id
    AND act.assay_id     = a.assay_id
    AND a.tid            = t.tid """
    "AND t.chembl_id     = '{}'".format(chembl_id))
    
    # cursor execute and fetch
    cursor.execute(qtext)
    query_result = cursor.fetchall()
    
    # create dataframe 
    return pd.DataFrame(query_result, columns=['Chembl id', 'SMILES',
                                               'Pub_year', 'Assay_des', 
                                               'Activity_type',
                                               'Relation', 'Activity_value',
                                               'Activity_units'])

def standardise_smiles(df):
    std_mols = [Chem.MolFromSmiles(SMILES) for SMILES in df.SMILES]
    std_smiles = [Chem.MolToSmiles(mol) for mol in std_mols]
    df.SMILES = std_smiles    
    
def find_matches(chembl_df, compare_df):
    # Merge dfs to see what is available in Chembl db
    matches_df = pd.merge(chembl_df, compare_df, 
                          how='inner', on=['SMILES'])
    
    return matches_df
    
def clean_matches(matches_df):   
    matches_df = matches_df[(matches_df.Activity_type == 'Ki') 
                            | (matches_df.Activity_type == 'IC50')]
    
    # Remove row if no activity value given
    matches_df = matches_df[matches_df.Activity_value.notna()]
    
    return matches_df

def filter_matches(matches_df, target_name):
    groups = matches_df.groupby(by=['Mol'])
    data_to_add = []
    for group_name, df_group in groups:
        if len(df_group) == 1:
            data_to_add.append(df_group)
        if len(df_group) >= 2:
            if len(df_group[df_group.Activity_type== 'Ki']) == 0:
                df_group = df_group[df_group.Activity_type == 'IC50']
                if len(df_group[df_group.Pub_year == df_group.Pub_year. max()]) == 0:
                    df_group = df_group.drop_duplicates(subset=['Activity_type'], keep='first')
                    data_to_add.append(df_group)
                if len(df_group[df_group.Pub_year == df_group.Pub_year. max()]) == 1:
                    df_group = df_group[df_group.Pub_year == df_group.Pub_year. max()]
                    data_to_add.append(df_group)
                if len(df_group[df_group.Pub_year == df_group.Pub_year. max()]) >= 2:
                    df_group = df_group[df_group.Pub_year == df_group.Pub_year. max()]
                    df_group = df_group.drop_duplicates(subset=['Activity_type'], keep='first')
                    data_to_add.append(df_group)     
            if len(df_group[df_group.Activity_type== 'Ki']) == 1:
                df_group = df_group[df_group.Activity_type == 'Ki']
                data_to_add.append(df_group)    
            if len(df_group[df_group.Activity_type== 'Ki']) >= 2:
                df_group = df_group[df_group.Activity_type == 'Ki']
                if len(df_group[df_group.Pub_year == df_group.Pub_year.max()]) == 0:
                    df_group = df_group.drop_duplicates(subset=['Activity_type'], keep='first')
                    data_to_add.append(df_group)
                if len(df_group[df_group.Pub_year == df_group.Pub_year.max()]) == 1:
                    df_group = df_group[df_group.Pub_year == df_group.Pub_year. max()]
                    data_to_add.append(df_group)
                if len(df_group[df_group.Pub_year == df_group.Pub_year. max()]) >= 2:
                    df_group = df_group[df_group.Pub_year == df_group.Pub_year. max()]
                    df_group = df_group.drop_duplicates(subset=['Activity_type'], keep='first')
                    data_to_add.append(df_group)
    
    # Create DF from list of dfs
    filtered_df = pd.concat(data_to_add)
    
    # Add df target name
    filtered_df.insert(loc=0, 
                       column='Target_name', 
                       value=target_name) 
    
    # Save to .csv
    filtered_df.to_csv("Filtered_lists/{}_matches.csv".format(target_name), 
                      index=False)
    
    return filtered_df

def save_one_csv(list_all_dfs):
    # Concat list dfs onto a single df    
    one_df = pd.concat(list_all_dfs)
    
    # Make more modelling friendly
    one_df = one_df[['Target_name','Chembl id',
                    'SMILES', 'Activity_type',
                    'Relation', 'Activity_value',
                    'Activity_units']]
    
    # Save to .csv
    one_df.to_csv("Filtered_lists/All_target_matches.csv", 
                      index=False)
    

chembl_smi_list = ['chembl_28', 'chembl_219', 'chembl_276',
                    'chembl_10378', 'chembl_10498', 'chembl_10752',
                    'chembl_11279', 'chembl_11359', 'chembl_11534',
                    'chembl_12670', 'chembl_12968',
                    'chembl_18061', 'chembl_20014']
chembl_targets = ['CHEMBL1952', 'CHEMBL245', 'CHEMBL254',
                  'CHEMBL4072','CHEMBL3837','CHEMBL1991',
                  'CHEMBL3772','CHEMBL288','CHEMBL2954',
                  'CHEMBL1974','CHEMBL4792',
                  'CHEMBL4296','CHEMBL4722']

list_all_dfs = []

for chembl_smiles, chembl_target in zip(chembl_smi_list,chembl_targets):
    chembl_id = chembl_target         
    chembl_df = sql_query(chembl_id) 
    compare_df = pd.read_csv('SMILES_data/{}_actives.smi'.format(chembl_smiles), 
                             delim_whitespace=True, header=None,
                             names=['SMILES', 'Mol'])

    # Rename Mol in compare_df to only have numbers
    compare_df.Mol = [int(name.strip('mol')) for name in compare_df.Mol] 
    
    # Standardise smiles
    standardise_smiles(chembl_df)
    standardise_smiles(compare_df)
    
    # Find matches
    matches_df = find_matches(chembl_df, compare_df)
    
    # Filter matches
    filt_matches = clean_matches(matches_df)
    
    # Clean up and get Ki or IC50 with most recent chosen
    filt_matches_final = filter_matches(filt_matches,chembl_id)
    
    # Append to list for saving as one csv
    list_all_dfs.append(filt_matches_final)

save_one_csv(list_all_dfs)

# close database connection and psycopg2 cursor
cursor.close()
conn.close()

