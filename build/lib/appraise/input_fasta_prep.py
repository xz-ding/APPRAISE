"""
Functions needed to prepare the input fasta for Alphafold2 modeling.
Author: Xiaozhe Ding
Email: xding@caltech.edu, dingxiaozhe@gmail.com
"""

import os
import pandas as pd
import numpy as np


def save_fasta(jobname, sequence, folder_path='./input_fasta/'):
    """
    Save a fasta file given the job name and sequence. The function creates a
    new folder at the specified folder_path, if it does not already exist, and
    then saves a fasta file with the given jobname and sequence in the folder.

    Args:
        jobname (str): The name of the job for which the fasta file will be saved.
        sequence (str): The sequence to be saved in the fasta file.
        folder_path (str, optional): The path to the folder where the fasta file will be saved. Defaults to './input_fasta/'.
    """
    if folder_path[-1] != '/':
        folder_path += '/'

    if not os.path.exists(folder_path):
        # Create a new directory because it does not exist
        os.makedirs(folder_path)
        print("$ Created folder {}".format(folder_path))

    with open(folder_path + jobname + '.fasta', 'w') as f:
        f.write('> {}\n{}'.format(jobname, sequence))
    print('$ Generated ' + folder_path + jobname + '.fasta')


def get_complex_fastas(receptor_name, receptor_seq, list_peptide1_names,
                       list_peptide1_seqs, mode='pairwise', square_matrix=True,
                       list_peptide2_names=[], list_peptide2_seqs=[], pool_size=4,
                       folder_path='./input_fasta/', use_glycine_linker=False,
                       glycine_linker_length=30, random_seed=42,
                       prepare_receptor_model=True):
    """
    Create and save input fastas for structural modeling in ColabFold (with
    either AlphaFold-multimer or ESMfold). More specifically, the function generates fasta files for the specified receptor and peptides according to the specified mode and saves them in the specified folder_path. The mode argument can be 'pairwise', in which case fasta files will be generated

    Args:
        receptor_name (str): The name of the receptor for the complex structure.
        receptor_seq (str): The sequence of the receptor.
        list_peptide1_names (list of str): A list of names for the first set of peptides.
        list_peptide1_seqs (list of str): A list of sequences for the first set of peptides.
        mode (str, optional): The mode in which the complex-modeling fasta files will be generated. Can be 'pairwise', 'single', or 'pooled'. Defaults to 'pairwise'.
        square_matrix (bool, optional): If True, the list_peptide1_names and list_peptide1_seqs arguments will be used for both the first and second sets of peptides. If False, the list_peptide2_names and list_peptide2_seqs arguments must be provided. Defaults to True.
        list_peptide2_names (list of str, optional): A list of names for the second set of peptides. Defaults to [].
        list_peptide2_seqs (list of str, optional): A list of sequences for the second set of peptides. Defaults to [].
        pool_size (int, optional): The number of peptides to include in each fasta file when the mode argument is set to 'pooled'. Defaults to 4.
        folder_path (str, optional): The path to the folder where the fasta files will be saved. Defaults to './input_fasta/'.
        use_glycine_linker (bool, optional): If True, the peptides and receptor will be linked with a glycine linker of a specified length. If False, the peptides and receptor will be linked with a colon (:) character. Defaults to False.
        glycine_linker_length (int, optional): The length of the glycine linker to be used when the use_glycine_linker argument is set to True. Defaults to 30.
        prepare_receptor_model (bool, optional): If True, an additional fasta file with receptor sequence only will be prepared. Defaults to True.
    """
    list_query_sequence = []
    list_jobname = []

    if use_glycine_linker is False:
        split_linker = ':'
    else:
        split_linker = 'G' * glycine_linker_length

    if mode == 'pairwise':
        if square_matrix:
            list_peptide2_names = list_peptide1_names
            list_peptide2_seqs = list_peptide1_seqs

        for j, peptide1_name in enumerate(list_peptide1_names):
            if len(list_peptide2_names) == 0:
                query_sequence = list_peptide1_seqs[j] + split_linker + receptor_seq
                jobname = receptor_name + "_and_" + peptide1_name
                list_query_sequence += [query_sequence]
                list_jobname += [jobname]
                save_fasta(jobname, query_sequence, folder_path)
            else:
                for k, peptide2_name in enumerate(list_peptide2_names):
                    query_sequence = list_peptide1_seqs[j] + split_linker + list_peptide2_seqs[k] + split_linker + receptor_seq
                    jobname = receptor_name + "_and_" + peptide1_name + "_vs_" + peptide2_name
                    list_query_sequence += [query_sequence]
                    list_jobname += [jobname]
                    save_fasta(jobname, query_sequence, folder_path)

        # Prepare a fasta for receptor-only model if requested
        if prepare_receptor_model:
            query_sequence = receptor_seq
            jobname = receptor_name + "_receptor_model"
            list_query_sequence += [query_sequence]
            list_jobname += [jobname]
            save_fasta(jobname, query_sequence, folder_path)

    elif mode == 'single':
        for j, peptide1_name in enumerate(list_peptide1_names):
            query_sequence = list_peptide1_seqs[j] + split_linker + receptor_seq
            jobname = receptor_name + "_and_" + peptide1_name
            list_query_sequence += [query_sequence]
            list_jobname += [jobname]
            save_fasta(jobname, query_sequence, folder_path)

        # Prepare a fasta for receptor-only model if requested
        if prepare_receptor_model:
            query_sequence = receptor_seq
            jobname = receptor_name + "_receptor_model"
            list_query_sequence += [query_sequence]
            list_jobname += [jobname]
            save_fasta(jobname, query_sequence, folder_path)

    elif mode == 'pooled':
        df_peptides_to_model = pd.DataFrame({'peptide_name': list_peptide1_names, 'peptide_seq': list_peptide1_seqs})
        #set random state
        np.random.seed(seed=random_seed)

        for i in range(int(len(df_peptides_to_model) / pool_size)):
            # Randomly choose a pool and remove it from the original library
            df_pool = df_peptides_to_model.sample(pool_size)
            df_peptides_to_model = df_peptides_to_model.drop(df_pool.index)
            df_pool = df_pool.reset_index()

            # Generate query sequence and standard names with the pool
            query_sequence = ''
            jobname = receptor_name + "_and_"
            for j in range(len(df_pool)):
                jobname += df_pool.loc[j]['peptide_name'] + "_vs_"
                query_sequence += df_pool.loc[j]['peptide_seq'] + split_linker
            jobname = jobname[0:-4]
            query_sequence += receptor_seq

            #record the generated results
            list_query_sequence += [query_sequence]
            list_jobname += [jobname]
            save_fasta(jobname, query_sequence, folder_path)

        # Deal with the remaining leftover sequences
        if len(df_peptides_to_model) > 0:
            # Use the rest variants that were not enough to form a standard-sized pool to form a pool
            df_pool = df_peptides_to_model.reset_index()

            # Generate query sequence and standard names with the pool
            query_sequence = ''
            jobname = receptor_name + "_and_"
            for j in range(len(df_pool)):
                jobname += df_pool.loc[j]['peptide_name'] + "_vs_"
                query_sequence += df_pool.loc[j]['peptide_seq'] + split_linker
            jobname = jobname[0:-4]
            query_sequence += receptor_seq

            #record the generated results
            list_query_sequence += [query_sequence]
            list_jobname += [jobname]
            save_fasta(jobname, query_sequence, folder_path)

        # Prepare a fasta for receptor-only model if requested
        if prepare_receptor_model:
            query_sequence = receptor_seq
            jobname = receptor_name + "_receptor_model"
            list_query_sequence += [query_sequence]
            list_jobname += [jobname]
            save_fasta(jobname, query_sequence, folder_path)

    elif mode == 'single_chains' or mode == 'single_chain':
        # Prepare single chain fastas for the petides and receptors
        for j, peptide1_name in enumerate(list_peptide1_names):
            query_sequence = list_peptide1_seqs[j]
            jobname = receptor_name + "_targeting_peptide_" + peptide1_name
            list_query_sequence += [query_sequence]
            list_jobname += [jobname]
            save_fasta(jobname, query_sequence, folder_path)

        # Prepare a fasta for receptor-only model if requested
        if prepare_receptor_model:
            query_sequence = receptor_seq
            jobname = receptor_name + "_receptor_model"
            list_query_sequence += [query_sequence]
            list_jobname += [jobname]
            save_fasta(jobname, query_sequence, folder_path)

    return list_query_sequence, list_jobname


def load_peptides(csv_file_path):
    """
    Load peptide names and sequences from a .csv file.

    csv_file_path (str): The path to the .csv file containing the peptide names and sequences.
    """
    df_peptide_lists = pd.read_csv(csv_file_path)
    list_peptide_names = df_peptide_lists['peptide_name'].to_list()
    list_peptide_seqs = df_peptide_lists['peptide_seq'].to_list()
    return list_peptide_names, list_peptide_seqs
