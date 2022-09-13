"""
Functions needed to prepare the input fasta for Alphafold2 modeling.
Author: Xiaozhe Ding
Email: xding@caltech.edu, dingxiaozhe@gmail.com
"""


import os
import pandas as pd

def save_fasta(jobname, sequence, folder_path='./input_fasta/'):
    """
    Save a fasta file given the job name and sequence.
    """
    if folder_path[-1] != '/':
        folder_path += '/'

    if not os.path.exists(folder_path):
        # Create a new directory because it does not exist
        os.makedirs(folder_path)
        print("> Created folder {}".format(folder_path))

    with open(folder_path + jobname + '.fasta', 'w') as f:
        f.write('> {}\n{}'.format(jobname, sequence))
    print('> Generated ' + folder_path + jobname + '.fasta')

def get_complex_fastas(receptor_name, receptor_seq, list_peptide1_names, \
    list_peptide1_seqs, mode='pairwise', square_matrix=True, \
    list_peptide2_names=[], list_peptide2_seqs=[], pool_size=4, folder_path='./input_fasta/'):
    """
    Create and save input fastas for AlphaFold2.
    """

    list_query_sequence = []
    list_jobname = []

    if mode == 'pairwise':
        if square_matrix:
            list_peptide2_names = list_peptide1_names
            list_peptide2_seqs = list_peptide1_seqs

        for j, peptide1_name in enumerate(list_peptide1_names):
            if len(list_peptide2_names) == 0:
                query_sequence =  list_peptide1_seqs[j] + ":" + receptor_seq
                jobname = receptor_name + "_and_" + peptide1_name
                list_query_sequence += [query_sequence]
                list_jobname += [jobname]
                save_fasta(jobname, query_sequence, folder_path)
            else:
                for k, peptide2_name in enumerate(list_peptide2_names):
                    query_sequence =  list_peptide1_seqs[j] + ":" + list_peptide2_seqs[k] + ":"+ receptor_seq
                    jobname = receptor_name + "_and_" + peptide1_name + "_vs_" + peptide2_name
                    list_query_sequence += [query_sequence]
                    list_jobname += [jobname]
                    save_fasta(jobname, query_sequence, folder_path)
    elif mode == 'single':
        for j, peptide1_name in enumerate(list_peptide1_names):
            query_sequence =  list_peptide1_seqs[j] + ":" + receptor_seq
            jobname = receptor_name + "_and_" + peptide1_name
            list_query_sequence += [query_sequence]
            list_jobname += [jobname]
            save_fasta(jobname, query_sequence, folder_path)

    elif mode == 'pooled':
        df_peptides_to_model = pd.DataFrame({'peptide_name':list_peptide1_names, 'peptide_seq':list_peptide1_seqs})
        for i in range(int(len(df_peptides_to_model) / pool_size)):
            # Randomly choose a pool and remove it from the original library
            df_pool = df_peptides_to_model.sample(n=pool_size)
            df_peptides_to_model = df_peptides_to_model.drop(df_pool.index)
            df_pool = df_pool.reset_index()

            # Generate query sequence and standard names with the pool
            query_sequence = ''
            jobname = receptor_name + "_and_"
            for j in range(len(df_pool)):
                jobname += df_pool.loc[j]['peptide_name'] + "_vs_"
                query_sequence += df_pool.loc[j]['peptide_seq'] + ":"
            jobname = jobname[0:-4]
            query_sequence += receptor_seq
            list_query_sequence += [query_sequence]
            list_jobname += [jobname]
            save_fasta(jobname, query_sequence, folder_path)

        if len(df_peptides_to_model) > 0:
            # Use the rest variants that were not enough to form a standard-sized pool to form a pool
            df_pool = df_peptides_to_model.reset_index()

            # Generate query sequence and standard names with the pool
            query_sequence = ''
            jobname = receptor_name + "_and_"
            for j in range(len(df_pool)):
                jobname += df_pool.loc[j]['peptide_name'] + "_vs_"
                query_sequence += df_pool.loc[j]['peptide_seq'] + ":"
            jobname = jobname[0:-4]
            query_sequence += receptor_seq
            list_query_sequence += [query_sequence]
            list_jobname += [jobname]
            save_fasta(jobname, query_sequence, folder_path)

    return list_query_sequence, list_jobname


def load_peptides(csv_file_path):
    """
    Load peptide names and sequences from a .csv file.
    """
    df_peptide_lists = pd.read_csv(csv_file_path)
    list_peptide_names = df_peptide_lists['peptide_name'].to_list()
    list_peptide_seqs = df_peptide_lists['peptide_seq'].to_list()
    return list_peptide_names, list_peptide_seqs
