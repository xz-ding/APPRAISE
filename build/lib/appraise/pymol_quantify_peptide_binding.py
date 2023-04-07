'''
Pymol script for quantifying peptide-receptor models.
Author: Xiaozhe Ding
Email: xding@caltech.edu, dingxiaozhe@gmail.com
'''

import platform
import sys
import argparse
import glob
import os

import csv
from datetime import datetime
import numpy as np
import numpy.linalg as LA


from itertools import groupby
from operator import itemgetter



# Analyze a folder using the default setting (activated when running the script from the command line)
use_relaxed_global = False
mod_start_resi_global = 0
mod_end_resi_global = 0
pLDDT_threshold_global = 0
receptor_chain_global = 'last'
anchor_site_global = 'C-term'
pdb_path = ' '
database_path = ' '




def count_clash(selection='(all)', name='bump_check', quiet=1, clash_distance_threshold=1):
    '''
    A function to count the number of clashing atoms.

    Adapted from show_bumps function by Thomas Holder.

    License info for original show_bumps function:

    http://pymolwiki.org/index.php/count_clash

    (c) 2011 Thomas Holder, MPI for Developmental Biology

    License: BSD-2-Clause

    DESCRIPTION

    Visualize VDW clashes

    ARGUMENTS

    selection = string: atom selection {default: all}

    name = string: name of CGO object to create {default: bump_check}
    '''
    if cmd.count_atoms(selection) > 0:
        cmd.delete(name)
        cmd.create(name, selection, zoom=0)
        cmd.set('sculpt_vdw_vis_mode', 1, name)
        cmd.set('sculpt_field_mask', 0x020)  # cSculptVDW
        for state in range(1, 1 + cmd.count_states('%' + name)):
            cmd.sculpt_activate(name, state)
            vdw_strain = cmd.sculpt_iterate(name, state, cycles=0)
            if not int(quiet):
                print('VDW Strain in state %d: %f' % (state, vdw_strain))
        return vdw_strain
    else:
        return 0


def average_b(selection):
    """
    Function for extracting the b value (pLDDT score in AlphaFold-predicted
    modelss) from residues of selection.
    Source: https://pymolwiki.org/index.php/Average_b
    Author: Gregor Hagelueken
    """
    stored.tempfactor = 0
    stored.atomnumber = 0
    cmd.iterate(selection, "stored.tempfactor = stored.tempfactor + b")
    cmd.iterate(selection, "stored.atomnumber = stored.atomnumber + 1")
    if stored.atomnumber > 0:
        averagetempfactor = stored.tempfactor / stored.atomnumber
    else:
        averagetempfactor = 0
    return averagetempfactor

def get_pLDDT_weighted_coordinates(selection):
    """
    Calculate pLDDT_weighted coordinates from residues of selection.
    """
    stored.tempfactor = 0
    stored.atomnumber = 0
    stored.x_product = 0
    stored.y_product = 0
    stored.z_product = 0
    cmd.iterate_state(1, selection, "stored.tempfactor = stored.tempfactor + b")
    cmd.iterate_state(1, selection, "stored.x_product = stored.x_product + x * b")
    cmd.iterate_state(1, selection, "stored.y_product = stored.y_product + y * b")
    cmd.iterate_state(1, selection, "stored.z_product = stored.z_product + z * b")
    cmd.iterate_state(1, selection, "stored.atomnumber = stored.atomnumber + 1")
    weighted_x = stored.x_product / stored.tempfactor
    weighted_y = stored.y_product / stored.tempfactor
    weighted_z = stored.z_product / stored.tempfactor

    return [weighted_x, weighted_y, weighted_z]

def get_pLDDT_weighted_linear_center(selection):
    """
    Find the coordinates of the center of the primary sequences
    of residues of selection.
    """
    stored.tempfactor = 0
    stored.rank_product = 0

    cmd.iterate_state(1, selection, "stored.tempfactor = stored.tempfactor + b")
    cmd.iterate_state(1, selection, "stored.rank_product = stored.rank_product + rank * b")

    weighted_rank = round(stored.rank_product / stored.tempfactor)
    return cmd.get_coords("rank {} and ({})".format(str(weighted_rank), selection))[0]

def generate_pdb_path_list(AF2_results_path, use_relaxed=use_relaxed_global):
    """
    get a list of pdb paths from colabfold-alphafold output results.
    """
    global use_relaxed_global
    # automatically determine if the models were amber-relaxed and whether relaxed models should be used for analysis.
    if use_relaxed_global == 'auto':
        list_all_pdb = glob.glob(AF2_results_path + '*.pdb')
        use_relaxed_global = False
        for pdb_path_to_check in list_all_pdb:
            if '_relaxed_' in pdb_path_to_check:
                use_relaxed_global = True

    # find pdb files with matching names
    use_relaxed = use_relaxed_global
    if use_relaxed:
        list_pdb_path = glob.glob(AF2_results_path + '*_relaxed_*.pdb')
    else:
        list_pdb_path = glob.glob(AF2_results_path + '*_unrelaxed_*.pdb')
    return list_pdb_path

#to do next
def parse_pdb_file_name(pdb_path):
    """
    Parse out the peptide name and receptor name from pdb file names.
    """
    if use_relaxed_global:
        pdb_string = pdb_path.split('/')[-1][0:-4].split('_relaxed_')[0]
    else:
        pdb_string = pdb_path.split('/')[-1][0:-4].split('_unrelaxed_')[0]
    receptor_name = pdb_string.split('_and_')[0]
    list_peptide_name = pdb_string.split('_and_')[1].split('_vs_')
    return receptor_name, list_peptide_name


def find_chain_IDs(model_name, receptor_chain=receptor_chain_global,
    glycine_linkers='auto', min_linker_length=10):
    """
    Find out the receptor chain ID and a list of peptide chain ID.

    model_name(string): name of the model.

    receptor_chain(string): can be 'last' or 'A', 'B', 'C' etc.. If 'last', then
    the last chain in the model will be automatically assigned to be the receptor.

    glycine_linkers: (boolean or string) If True, the chains will be split by
    tandem glycine linkers longer than 10aa before processing. If 'auto', then
    a model containing only one single chain will be split, while a model
    containing multiple chains will be processed as is.
    """
    # find out the receptor chain ID and generate a list of chain IDs for peptides
    chain_list = cmd.get_chains(model_name)

    if glycine_linkers == True or (glycine_linkers == 'auto' and len(chain_list) == 1):
        _ = split_by_glycine_linkers(chain_list)
        chain_list = cmd.get_chains(model_name)

    if receptor_chain == 'last':
        receptor_chain = chain_list[-1]
        chain_list.remove(receptor_chain)
        list_peptide_chain = chain_list
    else:
        list_peptide_chain = []
        for chain in chain_list:
            if chain != receptor_chain:
                list_peptide_chain += [chain]

    return receptor_chain, list_peptide_chain


def find_residue_indices(selection_string='sele'):
    """
    A utility function to find out the residue IDs of a selection.
    """
    stored.residues = set()
    cmd.iterate(selector.process(selection_string), 'stored.residues.add(resv)')

    return list(sorted(stored.residues))

def find_glycine_linkers(selection_string='all', min_linker_length=25):
    """
    Return a list of start and end indices of glycine linkers that are longer
    than the min_linker_length.
    """
    selection_string = 'pepseq {}'.format('G' * min_linker_length)
    glycine_linker_residue_indices = find_residue_indices(selection_string)

    list_glycine_linker_ranges = []
    for k, g in groupby(enumerate(glycine_linker_residue_indices), lambda ix: ix[0] - ix[1]):
        individual_glycine_linker_indices = list(map(itemgetter(1), g))
        list_glycine_linker_ranges += [[individual_glycine_linker_indices[0], individual_glycine_linker_indices[-1]]]
    return list_glycine_linker_ranges


def split_by_glycine_linkers(list_old_chain_ids, min_linker_length=10):
    """
    Split all chains that are linked by glycines.
    """
    print('APPRAISE> Splitting the model using glycine linkers.')
    list_new_chain_indices = []
    for old_chain_id in list_old_chain_ids:
        list_glycine_linker_ranges = find_glycine_linkers(old_chain_id, \
            min_linker_length=min_linker_length)
        if len(list_glycine_linker_ranges) > 0:
            # Identify the ranges for the first sub chain
            if list_glycine_linker_ranges[0][0] != 1:
                list_new_chain_indices += [[old_chain_id, 1, list_glycine_linker_ranges[0][0] - 1]]
            #Identify the ranges for the other sub chains
            for i, glycine_linkers in enumerate(list_glycine_linker_ranges):
                cmd.do('remove chain {} and resi {}-{}'.format(old_chain_id, \
                    glycine_linkers[0], glycine_linkers[1]))
                if i != len(list_glycine_linker_ranges) - 1:
                    list_new_chain_indices += [[old_chain_id, \
                        list_glycine_linker_ranges[i][1] + 1, \
                        list_glycine_linker_ranges[i + 1][0] - 1]]
                else:
                    list_new_chain_indices += [[old_chain_id, \
                        list_glycine_linker_ranges[i][1] + 1, \
                        100000]]

    if len(list_new_chain_indices) > 0:
        # get a list of new chain IDs
        list_new_chain_ids = []
        for i in range(len(list_new_chain_indices)):
            candidate_char = 'A'
            while (candidate_char in list_old_chain_ids) or (candidate_char in list_new_chain_ids):
                 candidate_char = chr(ord(candidate_char) + 1)
            list_new_chain_ids += [candidate_char]

        # alter the IDs of new chains
        for i, new_chain_indces in enumerate(list_new_chain_indices):
            new_chain_id = list_new_chain_ids[i]
            old_chain_id = new_chain_indces[0]
            start_index = new_chain_indces[1]
            end_index = new_chain_indces[2]
            cmd.do('alter chain {} and resi {}-{}, chain=\"{}\"'.format(old_chain_id, \
                start_index, end_index, new_chain_id))
            cmd.do('alter chain {} and resi {}-{}, resi=str(int(resi) - {})'.format(new_chain_id, \
                start_index, end_index, (start_index - 1)))

        #delete the empty old chains
        for old_chain_id in list_old_chain_ids:
            cmd.do('remove chain {}'.format(old_chain_id))

        return list_new_chain_ids
    else:
        print('APPRAISE> No glycine linker detected.')
        return list_old_chain_ids



def quantify_contact_atom(peptide_chain, receptor_chain=receptor_chain_global, \
    pep_mod_start_resi=0, pep_mod_end_resi=0, b_threshold=0, b_weighted=False):
    """
    A script to quantify the number of contact atoms beween a peptide and the receptor.

    peptide_chain (string): chain ID of the peptide.

    receptor_chain (string): chain ID of the receptor.

    pep_mod_start_resi (int): the start residue index of the range of peptide to
     be considered. If 0, the whole peptide will be considered.

    pep_mod_end_resi (int): the end residue index of the range of peptide to be
    considered. If 0, the whole peptide will be considered.

    b_threshold (float): the pLDDT threshold for an residue to be considered.

    b_weighted (boolean): if True, the a pLDDT-weighted number of contact atoms
    will be counted.

    """

    if pep_mod_start_resi == 0 or pep_mod_end_resi == 0:
        # Default: count contacting atoms in the  whole peptide
        contact_atom_in_peptide = cmd.count_atoms('(chain {} and b > {}) within 5 of chain {}'.format(peptide_chain, str(b_threshold), receptor_chain))
        contact_atom_in_receptor = cmd.count_atoms('chain {} within 5 of (chain {} and b > {})'.format(receptor_chain, peptide_chain, str(b_threshold)))
        if b_weighted:
            contact_atom_in_peptide = contact_atom_in_peptide * average_b('(chain {} and b > {}) within 5 of chain {}'.format(peptide_chain, str(b_threshold), receptor_chain))
            contact_atom_in_receptor = contact_atom_in_receptor * average_b('chain {} within 5 of (chain {} and b > {})'.format(receptor_chain, peptide_chain, str(b_threshold)))

        total_contact_atom_in_interface = contact_atom_in_peptide + contact_atom_in_receptor
        return total_contact_atom_in_interface
    else:
        # count contacting atoms with the insertion only
        contact_atom_in_peptide_ins_only = cmd.count_atoms('(chain {} and resi {}-{} and b > {}) within 5 of chain {}'.format(peptide_chain, str(pep_mod_start_resi), str(pep_mod_end_resi), str(b_threshold), receptor_chain))
        contact_atom_in_receptor_ins_only = cmd.count_atoms('chain {} within 5 of (chain {} and resi {}-{} and b > {})'.format(receptor_chain, peptide_chain, str(pep_mod_start_resi), str(pep_mod_end_resi), str(b_threshold)))
        if b_weighted:
            contact_atom_in_peptide_ins_only = contact_atom_in_peptide_ins_only * average_b('(chain {} and resi {}-{} and b > {}) within 5 of chain {}'.format(peptide_chain, str(pep_mod_start_resi), str(pep_mod_end_resi), str(b_threshold), receptor_chain))
            contact_atom_in_receptor_ins_only = contact_atom_in_receptor_ins_only * average_b('chain {} within 5 of (chain {} and resi {}-{} and b > {})'.format(receptor_chain, peptide_chain, str(pep_mod_start_resi), str(pep_mod_end_resi), str(b_threshold)))


        total_contact_atom_in_interface_ins_only = contact_atom_in_peptide_ins_only + contact_atom_in_receptor_ins_only
        return total_contact_atom_in_interface_ins_only




def quantify_peptide_binding_in_pdb(pairwise_mode=True, \
    receptor_chain=receptor_chain_global, glycine_linkers='auto'):
    """
    The main function to measure different metrics for quantifying
    peptide-receptor binding within a single pdb file.

    pairwise_mode: (boolean) if true, the first and last peptide will be treated
     as competitors and analyzed in pairs; the difference between the two peptides
     will be caulculated as well.

    receptor_chain: (string) can be 'last' or 'A', 'B', 'C' etc.. If 'last', then
    the last chain in the model will be automatically assigned to be the receptor.

    glycine_linkers: (boolean or string) If True, the chains will be split by
    tandem glycine linkers longer than 10aa before processing. If 'auto', then
    a model containing only one single chain will be split, while a model
    containing multiple chains will be processed as is.

    """
    #Clean up and load the new pdb file
    cmd.do('delete all')
    cmd.load(pdb_path)
    object_list = cmd.get_object_list('all')
    model_name = object_list[0]
    print('\n\n > Processing model {}'.format(model_name))



    #Split chains
    receptor_chain, list_peptide_chain = find_chain_IDs(model_name, \
        receptor_chain=receptor_chain, glycine_linkers=glycine_linkers)
    # if it is a single-chain model (without glycien linkers), skip.
    if len(list_peptide_chain) == 0:
        print('[Warning] Skipped! Only one chain was present in the model. APPRAISE analysis requires models containing at least two chains or a single chain separated by glycine linkers. ')
        return

    # find out the receptor chain ID and generate a list of chain IDs for peptides


    # get metainfo from the pdb file name
    receptor_name, list_peptide_name = parse_pdb_file_name(pdb_path)

    # ### debugging ###
    # print(receptor_chain)
    # print(list_peptide_chain)
    # ### debugging ###

    # find the receptor center and membrane_anchor (last 30 atoms)
    receptor_coordinates = np.array(cmd.get_coords('{} and chain {}'.format(model_name, receptor_chain )))
    receptor_center = np.mean(receptor_coordinates, axis=0)
    receptor_Rg = np.sqrt(np.sum((receptor_coordinates - receptor_center)**2) / len(receptor_coordinates))

    if anchor_site_global == 'C-term' or anchor_site_global == 'C':
        anchor_site_coordinates = np.mean(np.array(cmd.get_coords('{} and chain {}'.format(model_name, receptor_chain ))[-30:]), axis=0)
    elif anchor_site_global == 'N-term' or anchor_site_global == 'N':
        anchor_site_coordinates = np.mean(np.array(cmd.get_coords('{} and chain {}'.format(model_name, receptor_chain ))[0:30]), axis=0)
    elif type(anchor_site_global) == int or type(anchor_site_global) == float:
        anchor_site_coordinates = np.mean(np.array(cmd.get_coords('{} and chain {} and resi {}'.format(model_name, receptor_chain, anchor_site_global))), axis=0)
    else:
        print("APPRAISE> Unrecognized membrane anchor site. Using C terminus of the receptor by default.")
        anchor_site_coordinates = np.mean(np.array(cmd.get_coords('{} and chain {}'.format(model_name, receptor_chain ))[-30:]), axis=0)

    weighted_receptor_center = np.mean(np.array(get_pLDDT_weighted_coordinates('{} and chain {}'.format(model_name, receptor_chain ))), axis=0)

    # Loop through each chain in the model
    for j, peptide_chain in enumerate(list_peptide_chain):

        # Clean up and load the same model to avoid some memory issue
        cmd.do('delete all')
        cmd.load(pdb_path)
        receptor_chain, list_peptide_chain = find_chain_IDs(model_name, \
            receptor_chain=receptor_chain, glycine_linkers=glycine_linkers)

        # Measure peptide liength
        peptide_length = len(cmd.get_model('{} and chain {}'.format(model_name, peptide_chain)).get_residues())

        # find the name of the peptide
        if len(list_peptide_name) == len(list_peptide_chain):
            peptide_name = list_peptide_name[j]
        else:
            peptide_name = "NA"

        # An AAV9-specific modification to account for shorter peptide length (because of lack of insertion)
        if peptide_name == 'AAV9' and mod_start_resi_global > 2:
            pep_mod_start_resi_global = max(mod_start_resi_global - 1, 0)
            pep_mod_end_resi_global = max(mod_end_resi_global - 6, 0)
        else:
            pep_mod_start_resi_global = mod_start_resi_global
            pep_mod_end_resi_global = mod_end_resi_global

        # get the competitors
        list_competitor_name = list_peptide_name[0:j] + list_peptide_name[j+1:]
        list_competitor_chains = list_peptide_chain[0:j] + list_peptide_chain[j+1:]

        # # If in pairwise mode, find the name of the competitor peptide
        # if pairwise_mode and len(list_peptide_chain) == 2:
        #     competitor_name = list_peptide_name[1-j]
        # else:
        #     competitor_name = 'NA'

        # find the peptide center or weighted peptide center

        # ### debugging ###
        # print(model_name)
        # print(peptide_chain)
        # print(peptide_length)
        # print(list_competitor_chains)
        # ### debugging ###


        ar_coordinates = np.array(cmd.get_coords('{} and chain {}'.format(model_name, peptide_chain)))
        peptide_center = np.mean(ar_coordinates, axis=0)
        peptide_n_receptor_ar_coordinates = np.array(cmd.get_coords('{} and (chain {} or chain {})'.format(model_name, peptide_chain, receptor_chain)))
        peptide_n_receptor_center = np.mean(peptide_n_receptor_ar_coordinates, axis=0)
        ar_mod_coordinates = np.array(cmd.get_coords('{} and chain {} and resi {}-{}'.format(model_name, peptide_chain, str(pep_mod_start_resi_global), str(pep_mod_end_resi_global))))
        #peptide_mod_center = np.mean(ar_mod_coordinates, axis=0)

        ar_contacting_coordinates = np.array(cmd.get_coords('({} and chain {}) within 5 of chain {}'.format(model_name, peptide_chain, receptor_chain)))
        if ar_contacting_coordinates.size > 1:
            contacting_center = np.mean(ar_contacting_coordinates, axis=0)
            pLDDT_weighted_linear_contacting_center = get_pLDDT_weighted_coordinates('({} and chain {}) within 5 of chain {}'.format(model_name, peptide_chain, receptor_chain))

        #peptide_N_end_residue_center = np.mean(np.array(cmd.get_coords('{} and chain {} and resi {}'.format(model_name, peptide_chain, str(1)))), axis=0)
        #peptide_C_end_residue_center = np.mean(np.array(cmd.get_coords('{} and chain {} and resi {}'.format(model_name, peptide_chain, str(peptide_length)))), axis=0)

        weighted_peptide_center = np.mean(np.array(get_pLDDT_weighted_coordinates('{} and chain {}'.format(model_name, peptide_chain, str(pep_mod_start_resi_global), str(pep_mod_end_resi_global)))), axis=0)

        pLDDT_weighted_linear_center = get_pLDDT_weighted_linear_center('{} and chain {}'.format(model_name, peptide_chain))


        # Calculate the average pLDDT score
        average_pLDDT = average_b('{} and chain {} and resi {}-{}'.format(model_name, peptide_chain, str(pep_mod_start_resi_global), str(pep_mod_end_resi_global)))
        interface_pLDDT = average_b('({} and chain {}) within 5 of chain {}'.format(model_name, peptide_chain, receptor_chain))
        n_atom_above_threshold = cmd.count_atoms('{} and chain {} and resi {}-{} and b > {}'.format(model_name, peptide_chain, str(pep_mod_start_resi_global), str(pep_mod_end_resi_global), str(pLDDT_threshold_global)))


        # calculate the angle between receptor membrane_anchor (GPI-anchor) and the peptide center (releative to receptor center)
        if ar_contacting_coordinates.size > 1:
            unit_vector_membrane_anchor = (anchor_site_coordinates - receptor_center)/LA.norm(anchor_site_coordinates - receptor_center)
            unit_vector_peptide = (contacting_center - receptor_center)/LA.norm(contacting_center - receptor_center)
            angle_between_membrane_anchor_and_peptide = np.arccos(np.dot(unit_vector_membrane_anchor, unit_vector_peptide))
            contact_point_distance_to_membrane = LA.norm(anchor_site_coordinates - receptor_center) - np.cos(angle_between_membrane_anchor_and_peptide) * LA.norm(contacting_center - receptor_center)
            insert_contact_distance = LA.norm(peptide_center-contacting_center)
            weighted_peptide_contact_distance = LA.norm(weighted_peptide_center-contacting_center)
            pLDDT_weighted_linear_center_contact_distance = LA.norm(pLDDT_weighted_linear_center - pLDDT_weighted_linear_contacting_center)

        else:
            angle_between_membrane_anchor_and_peptide = 0
            contact_point_distance_to_membrane = 100
            insert_contact_distance = 100
            weighted_peptide_contact_distance = 100
            pLDDT_weighted_linear_center_contact_distance = 100

        # measure distances
        peptide_receptor_distance = LA.norm(receptor_center - peptide_center)
        chang_et_al_distance = LA.norm(receptor_center - peptide_n_receptor_center)
        weighted_peptide_receptor_distance = LA.norm(weighted_receptor_center - weighted_peptide_center)
        peptide_tip_receptor_distance = np.amin(LA.norm(ar_coordinates-receptor_center, axis=1))
        #end_to_end_distance = LA.norm(peptide_N_end_residue_center - peptide_C_end_residue_center)

        # Get peptide direction
        #peptide_direction = LA.norm(receptor_center - peptide_C_end_residue_center) - LA.norm(receptor_center - peptide_N_end_residue_center)

        # get peptide sequence
        peptide_seq = cmd.get_fastastr('{} and chain {}'.format(model_name, peptide_chain ), -1, 1).split('\n')[1]

        # count contact atom number
        total_contact_atom_in_interface_thresholded = quantify_contact_atom(peptide_chain, receptor_chain, b_threshold=pLDDT_threshold_global)
        total_contact_atom_in_interface_ins_only = quantify_contact_atom(peptide_chain, receptor_chain, pep_mod_start_resi_global, pep_mod_end_resi_global, b_threshold=0)
        total_contact_atom_in_interface_weighted = quantify_contact_atom(peptide_chain, receptor_chain, b_threshold=0, b_weighted=True)
        total_contact_atom_in_interface = quantify_contact_atom(peptide_chain, receptor_chain, 0, 0, 0, False)

        #calculate angle-factored contact atom number using a logistic function
        binding_angle_factor = 1 / (1 + np.exp(np.pi - 6 * np.absolute(angle_between_membrane_anchor_and_peptide)))
        total_contact_atom_in_interface_angle_factored = binding_angle_factor * total_contact_atom_in_interface
        #folded_factor = 1 / (1 + np.exp((end_to_end_distance - 20)/2))
        distance_to_membrane_factor = 1 / (1 + np.exp((3 - contact_point_distance_to_membrane)))

        # measure clashes
        #vdw_strain = count_clash('{} and (chain {} or chain {})'.format(model_name, peptide_chain, receptor_chain))
        vdw_strain = count_clash('(({} and chain {}) within 5 of ({} and chain {}) or ({} and chain {}) within 5 of ({} and chain {}))'.format(model_name, peptide_chain, model_name, receptor_chain, model_name, receptor_chain, model_name, peptide_chain))
        clash_number = cmd.count_atoms('(chain {} ) within {} of chain {}'.format(peptide_chain, 1, receptor_chain))


        if pairwise_mode:
            # Clean up and load the same model to avoid some internal pymol bug
            cmd.do('delete all')
            cmd.load(pdb_path)
            receptor_chain, list_peptide_chain = find_chain_IDs(model_name, \
                receptor_chain=receptor_chain, glycine_linkers=glycine_linkers)

            # Measure competitor peptide liength
            competitor_peptide_length = len(cmd.get_model('{} and chain {}'.format(model_name, list_competitor_chains[0])).get_residues())

            #Find out the modification start sites of the peptide
            competitor_peptide_name = list_peptide_name[1-j]

            if competitor_peptide_name == 'AAV9':
                pep_mod_start_resi_global_competitor = mod_start_resi_global - 1
                pep_mod_end_resi_global_competitor = mod_end_resi_global - 6
            else:
                pep_mod_start_resi_global_competitor = mod_start_resi_global
                pep_mod_end_resi_global_competitor = mod_end_resi_global

            # find the peptide center or weighted peptide center in the competitor
            ar_coordinates_competitor = np.array(cmd.get_coords('{} and chain {}'.format(model_name, list_competitor_chains[0])))
            peptide_center_competitor = np.mean(ar_coordinates_competitor, axis=0)
            peptide_n_receptor_ar_coordinates_competitor = np.array(cmd.get_coords('{} and (chain {} or chain {})'.format(model_name, list_competitor_chains[0], receptor_chain)))
            peptide_n_receptor_center_competitor = np.mean(peptide_n_receptor_ar_coordinates_competitor, axis=0)
            ar_mod_coordinates_competitor = np.array(cmd.get_coords('{} and chain {} and resi {}-{}'.format(model_name, list_competitor_chains[0], str(pep_mod_start_resi_global_competitor), str(pep_mod_end_resi_global_competitor))))
            #peptide_mod_center_competitor = np.mean(ar_mod_coordinates_competitor, axis=0)
            ar_contacting_coordinates_competitor = np.array(cmd.get_coords('({} and chain {}) within 5 of chain {}'.format(model_name, list_competitor_chains[0], receptor_chain)))
            if ar_contacting_coordinates_competitor.size > 1 :
                contacting_center_competitor = np.mean(ar_contacting_coordinates_competitor, axis=0)
                pLDDT_weighted_linear_contacting_center_competitor = get_pLDDT_weighted_coordinates('({} and chain {}) within 5 of chain {}'.format(model_name, list_competitor_chains[0], receptor_chain))

            weighted_peptide_center_competitor = np.mean(np.array(get_pLDDT_weighted_coordinates('{} and chain {}'.format(model_name, list_competitor_chains[0], str(pep_mod_start_resi_global_competitor), str(pep_mod_end_resi_global_competitor)))), axis=0)
            # competitor_peptide_N_end_residue_center = np.mean(np.array(cmd.get_coords('{} and chain {} and resi {}'.format(model_name, list_competitor_chains[0], str(1)))), axis=0)
            # competitor_peptide_C_end_residue_center = np.mean(np.array(cmd.get_coords('{} and chain {} and resi {}'.format(model_name, list_competitor_chains[0], str(competitor_peptide_length)))), axis=0)

            pLDDT_weighted_linear_center_competitor = get_pLDDT_weighted_linear_center('{} and chain {}'.format(model_name, list_competitor_chains[0]))


            interface_pLDDT_competitor = average_b('({} and chain {}) within 5 of chain {}'.format(model_name, list_competitor_chains[0], receptor_chain))


            # calculate the angle between receptor membrane_anchor (GPI-anchor) and the peptide center (releative to receptor center)
            if ar_contacting_coordinates_competitor.size > 1:
                unit_vector_membrane_anchor = (anchor_site_coordinates - receptor_center)/LA.norm(anchor_site_coordinates - receptor_center)
                unit_vector_competitor_peptide = (contacting_center_competitor - receptor_center)/LA.norm(contacting_center_competitor - receptor_center)
                angle_between_membrane_anchor_and_competitor_peptide = np.arccos(np.dot(unit_vector_membrane_anchor, unit_vector_competitor_peptide))
                contact_point_distance_to_membrane_competitor = LA.norm(anchor_site_coordinates - receptor_center) - np.cos(angle_between_membrane_anchor_and_competitor_peptide) * LA.norm(contacting_center_competitor - receptor_center)
                insert_contact_distance_competitor = LA.norm(peptide_center_competitor-contacting_center_competitor)
                weighted_peptide_contact_distance_competitor = LA.norm(weighted_peptide_center_competitor-contacting_center_competitor)
                pLDDT_weighted_linear_center_contact_distance_competitor = LA.norm(pLDDT_weighted_linear_center_competitor - pLDDT_weighted_linear_contacting_center_competitor)
            else:
                angle_between_membrane_anchor_and_competitor_peptide = 0
                contact_point_distance_to_membrane_competitor = 100
                insert_contact_distance_competitor = 100
                weighted_peptide_contact_distance_competitor = 100
                pLDDT_weighted_linear_center_contact_distance_competitor = 100

            # measure distances
            peptide_receptor_distance_competitor = LA.norm(receptor_center - peptide_center_competitor)
            chang_et_al_distance_competitor = LA.norm(receptor_center - peptide_n_receptor_center_competitor)
            weighted_peptide_receptor_distance_competitor = LA.norm(weighted_receptor_center - weighted_peptide_center_competitor)
            peptide_tip_receptor_distance_competitor = np.amin(LA.norm(ar_coordinates_competitor-receptor_center, axis=1))
            #end_to_end_distance_competitor = LA.norm(competitor_peptide_N_end_residue_center - competitor_peptide_C_end_residue_center)

            # Get peptide direction
            #peptide_direction_competitor = LA.norm(receptor_center - competitor_peptide_C_end_residue_center) - LA.norm(receptor_center - competitor_peptide_N_end_residue_center)


            # count contacting atoms in the competitors
            total_contact_atom_in_interface_competitor_thresholded = quantify_contact_atom(list_competitor_chains[0], receptor_chain, b_threshold=pLDDT_threshold_global)
            total_contact_atom_in_interface_competitor_ins_only = quantify_contact_atom(list_competitor_chains[0], receptor_chain, pep_mod_start_resi_global_competitor, pep_mod_end_resi_global_competitor, b_threshold=0)
            total_contact_atom_in_interface_weighted_competitor = quantify_contact_atom(list_competitor_chains[0], receptor_chain, b_threshold=0, b_weighted=True)
            total_contact_atom_in_interface_competitor = quantify_contact_atom(list_competitor_chains[0], receptor_chain, 0, 0, 0, False)

            #calculate angle-factored contact atom number using a logistic function
            binding_angle_factor_competitor = 1 / (1 + np.exp(np.pi - 6 * np.absolute(angle_between_membrane_anchor_and_competitor_peptide)))
            total_contact_atom_in_interface_competitor_angle_factored = binding_angle_factor_competitor * total_contact_atom_in_interface_competitor
            #folded_factor_competitor = 1 / (1 + np.exp((end_to_end_distance_competitor - 20)/2))
            distance_to_membrane_factor_competitor = 1 / (1 + np.exp((3 - contact_point_distance_to_membrane_competitor)))

            # measure clashes
            #vdw_strain_competitor = count_clash('{} and (chain {} or chain {})'.format(model_name, list_competitor_chains[0], receptor_chain))
            vdw_strain_competitor = count_clash('(({} and chain {}) within 5 of ({} and chain {}) or ({} and chain {}) within 5 of ({} and chain {}))'.format(model_name, list_competitor_chains[0], model_name, receptor_chain, model_name, receptor_chain, model_name, list_competitor_chains[0]))

            clash_number_competitor = cmd.count_atoms('(chain {} ) within {} of chain {}'.format(list_competitor_chains[0], 1, receptor_chain))


            # calculate differences
            peptide_receptor_distance_difference = peptide_receptor_distance_competitor - peptide_receptor_distance
            weighted_peptide_receptor_distance_difference = weighted_peptide_receptor_distance_competitor - weighted_peptide_receptor_distance
            peptide_tip_receptor_distance_difference = peptide_tip_receptor_distance_competitor - peptide_tip_receptor_distance

            total_contact_atom_in_interface_thresholded_difference = total_contact_atom_in_interface_thresholded - total_contact_atom_in_interface_competitor_thresholded
            total_contact_atom_in_interface_ins_only_difference = total_contact_atom_in_interface_ins_only - total_contact_atom_in_interface_competitor_ins_only
            total_contact_atom_in_interface_difference = total_contact_atom_in_interface - total_contact_atom_in_interface_competitor
            total_contact_atom_in_interface_difference_angle_factored =  total_contact_atom_in_interface_angle_factored - total_contact_atom_in_interface_competitor_angle_factored



            #calculate pLDDT-thresholded version of distances
            pLDDT_threshold_globaled_peptide_coordinates = np.array(cmd.get_coords('{} and chain {} and resi {}-{} and b > {}'.format(model_name, peptide_chain, str(pep_mod_start_resi_global), str(pep_mod_end_resi_global), str(pLDDT_threshold_global))))
            pLDDT_threshold_globaled_peptide_coordinates_competitor = np.array(cmd.get_coords('{} and chain {} and resi {}-{} and b > {}'.format(model_name, list_competitor_chains[0], str(pep_mod_start_resi_global_competitor), str(pep_mod_end_resi_global_competitor), str(pLDDT_threshold_global))))

            if pLDDT_threshold_globaled_peptide_coordinates.size > 1:
                pLDDT_threshold_globaled_peptide_center = np.mean(pLDDT_threshold_globaled_peptide_coordinates, axis=0)
                pLDDT_threshold_globaled_peptide_receptor_distance = LA.norm(receptor_center - pLDDT_threshold_globaled_peptide_center)
                if pLDDT_threshold_globaled_peptide_coordinates_competitor.size > 1:
                    pLDDT_threshold_globaled_peptide_center_competitor = np.mean(pLDDT_threshold_globaled_peptide_coordinates_competitor, axis=0)
                    pLDDT_threshold_globaled_peptide_receptor_distance_competitor = LA.norm(receptor_center - pLDDT_threshold_globaled_peptide_center_competitor)
                    pLDDT_threshold_globaled_peptide_receptor_distance_difference = pLDDT_threshold_globaled_peptide_receptor_distance_competitor - pLDDT_threshold_globaled_peptide_receptor_distance
                else:
                    pLDDT_threshold_globaled_peptide_receptor_distance_difference = 0
            else:
                pLDDT_threshold_globaled_peptide_receptor_distance = 100
                pLDDT_threshold_globaled_peptide_receptor_distance_difference = 0
        else:
            total_contact_atom_in_interface_difference = None

        list_to_append = [model_name, receptor_name, peptide_chain, peptide_name, \
                list_competitor_name, peptide_seq, peptide_length, receptor_Rg, anchor_site_global,\
                str(pep_mod_start_resi_global), str(pep_mod_end_resi_global), peptide_receptor_distance, chang_et_al_distance, chang_et_al_distance_competitor,\
                weighted_peptide_receptor_distance, peptide_tip_receptor_distance, peptide_tip_receptor_distance_competitor, pLDDT_threshold_globaled_peptide_receptor_distance, \
                peptide_receptor_distance_difference, weighted_peptide_receptor_distance_difference, peptide_tip_receptor_distance_difference, pLDDT_threshold_globaled_peptide_receptor_distance_difference,\
                total_contact_atom_in_interface_thresholded, total_contact_atom_in_interface_ins_only,\
                total_contact_atom_in_interface_thresholded_difference, total_contact_atom_in_interface_ins_only_difference, \
                total_contact_atom_in_interface, total_contact_atom_in_interface_competitor,\
                total_contact_atom_in_interface_weighted, total_contact_atom_in_interface_weighted_competitor,\
                angle_between_membrane_anchor_and_peptide, angle_between_membrane_anchor_and_competitor_peptide,\
                binding_angle_factor, binding_angle_factor_competitor,\
                contact_point_distance_to_membrane, contact_point_distance_to_membrane_competitor, \
                distance_to_membrane_factor, distance_to_membrane_factor_competitor, \
                insert_contact_distance, insert_contact_distance_competitor,\
                weighted_peptide_contact_distance, weighted_peptide_contact_distance_competitor,\
                pLDDT_weighted_linear_center_contact_distance, pLDDT_weighted_linear_center_contact_distance_competitor,\
                #folded_factor, folded_factor_competitor,\
                total_contact_atom_in_interface_difference, total_contact_atom_in_interface_angle_factored, \
                total_contact_atom_in_interface_difference_angle_factored, \
                angle_between_membrane_anchor_and_peptide, \
                average_pLDDT, pLDDT_threshold_global, n_atom_above_threshold,\
                interface_pLDDT, interface_pLDDT_competitor,\
                vdw_strain, clash_number, vdw_strain_competitor, clash_number_competitor]
                # peptide_direction, peptide_direction_competitor] #end_to_end_distance, \
        print("APPRAISE> New measurements added for {} (peptide: {})".format(model_name, peptide_name))

        # Open file in append mode
        with open(database_path, 'a') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow(list_to_append)
    #Clean up all pdbs
    cmd.do('delete all')
    return

def quantify_results_folder(AF2_results_path='./*result*/', \
    receptor_chain='last', anchor_site='C-term', use_relaxed='auto', time_stamp=True,
    mod_start_resi=3, mod_end_resi=9, pLDDT_threshold=0, output_path='auto',
    glycine_linkers='auto', dt_string='now'):
    """
    AF2_results_path: (str) path to the folder(s) that contain the AF2 modeling
    results. Wildcard is allowed. Default: all folders containing "results" in
    their name under the current working directory.
        # Format requirement of the pdb file (pdb files generated by the pipeline should already meet the requirements):
        # Filname should look like ReceptorName_and_Peptide1Name_vs_Peptide2Name(_vs_...PeptideXName)_relaxed_rank_X_model_X
        # The pdb should have more than 1 chains (('A'), 'B', 'C'...),
        # The last chain should be the receptor, with all the chains ahead being peptide(s)
        # Example 1: chain A = peptide 1, chain B = peptide 2, chain C = receptor
        # Example 2: chain B = peptide 2, chain C = receptor

    use_relaxed: (string or boolean) whether to use Amber-relaxed models for the
    quantification. If 'auto', APPRAISE will use relaxed models if relaxed models
     are available and use unrelaxed models otherwise.

    time_stamp: (boolean) whether to add timestamp to the output file name
    (to avoid complications between multiple measurements).

    receptor_chain: (string) chain ID of the receptor. Can be 'last' or 'A', 'B'
    , 'C' etc.. If 'last', then the last chain will be asigned as the receptor
    chain.

    anchor_site: (string or integer) anchor site used for calculating binding
    angle theta. Can be 'N-term', 'N', 'C-term', 'C', or an integer. If it's an
    integer, the reisdue with that index will be assigned as the anchor site.

    mod_start_resi: (int) Define the start of the modified residues (insert or
    substitution) within the peptide (for new features being developped).

    mod_end_resi: (int) Define the end of the modified residues (insert or
    substitution) within the peptide (for new features being developped).

    pLDDT_threshold_global: (int) pLDDT threshold for the thresholed
    measurements (for new features being developped).

    output_path: (str) the folder that will store the database file. If auto,
    it will be the same folder that stores the AF2_results folder.

    glycine_linkers: (boolean or string) If True, the chains will be split by
    tandem glycine linkers longer than 10aa before processing. If 'auto', then
    a model containing only one single chain will be split, while a model
    containing multiple chains will be processed as is.

    dt_string: (str) the date-time string to be attached to the database file as
    time stamp.
    """

    global use_relaxed_global
    global mod_start_resi_global
    global mod_end_resi_global
    global pLDDT_threshold_global
    global receptor_chain_global
    global anchor_site_global
    global pdb_path
    global database_path

    #Update the global variables according to function argument
    use_relaxed_global = use_relaxed
    mod_start_resi_global = mod_start_resi
    mod_end_resi_global = mod_end_resi
    pLDDT_threshold_global = pLDDT_threshold
    anchor_site_global = anchor_site
    receptor_chain_global = receptor_chain


    #clean up format of the paths
    if AF2_results_path[-1] != '/':
        AF2_results_path += '/'
    if AF2_results_path[0] != '/' and AF2_results_path[0] != '~' and AF2_results_path[0] != '.':
        AF2_results_path = './' + AF2_results_path

    #generate output path
    if output_path == 'auto':
        output_path = '/'.join(AF2_results_path.split('/')[0:-2])

    #generate time-stampped output file name
    if time_stamp:
        if dt_string == 'now':
            dt_string = datetime.now().strftime("%m%d%H%M")
        database_path = output_path + '/database_APPRAISE_measurements_{}.csv'.format(dt_string)
    else:
        database_path = output_path + '/database_APPRAISE_measurements.csv'


    # list out all pdb files in the folder
    list_pdb_path = generate_pdb_path_list(AF2_results_path, use_relaxed=use_relaxed_global)

    # Create a new row
    list_to_append = ["model_name", "receptor_name", "peptide_chain", "peptide_name", \
            'competitors', "peptide_seq", "peptide_length", "receptor_Rg", "anchor_site",\
            "pep_mod_start_resi", "pep_mod_end_resi", "peptide_receptor_distance", "chang_et_al_distance", "chang_et_al_distance_competitor",\
            "pLDDT_weighted_peptide_receptor_distance", "peptide_tip_receptor_distance", "peptide_tip_receptor_distance_competitor", "pLDDT_threshold_globaled_peptide_receptor_distance", \
            "peptide_receptor_distance_difference", "weighted_peptide_receptor_distance_difference", "peptide_tip_receptor_distance_difference", "pLDDT_threshold_globaled_peptide_receptor_distance_difference",\
            "total_contact_atom_in_interface_thresholded", "total_contact_atom_in_interface_ins_only", \
            "total_contact_atom_in_interface_thresholded_difference", "total_contact_atom_in_interface_ins_only_difference", \
            "total_contact_atom_in_interface", "total_contact_atom_in_interface_competitor", \
            "total_contact_atom_in_interface_weighted", "total_contact_atom_in_interface_weighted_competitor",\
            "angle_between_membrane_anchor_and_peptide", "angle_between_membrane_anchor_and_competitor_peptide",\
            "binding_angle_factor", "binding_angle_factor_competitor",\
            "contact_point_distance_to_membrane", "contact_point_distance_to_membrane_competitor",
            "distance_to_membrane_factor", "distance_to_membrane_factor_competitor", \
            "insert_contact_distance", "insert_contact_distance_competitor",\
            "weighted_peptide_contact_distance", "weighted_peptide_contact_distance_competitor",\
            "pLDDT_weighted_linear_center_contact_distance", "pLDDT_weighted_linear_center_contact_distance_competitor",\
            #"folded_factor", "folded_factor_competitor",\
            "total_contact_atom_in_interface_difference", "total_contact_atom_in_interface_angle_factored", \
            "total_contact_atom_in_interface_difference_angle_factored", \
            "angle_between_membrane_anchor_and_peptide",
            "average_pLDDT", "pLDDT_threshold_global", "n_atom_above_threshold", \
            "interface_pLDDT", "interface_pLDDT_competitor",\
            "vdw_strain", "clash_number", "vdw_strain_competitor", "clash_number_competitor"]
            #"peptide_direction", "peptide_direction_competitor"] #"end_to_end_distance", \


    # Check if the database file exists
    if not os.path.exists(database_path):
        # The file does not exist, create a new one with headers
        with open(database_path, 'a') as write_obj:
            # Create a writer object from csv module
            csv_writer = csv.writer(write_obj)
            # Add contents of list as last row in the csv file
            csv_writer.writerow(list_to_append)


    #  (redundant block) automatically determine if the models were amber-relaxed and whether relaxed models should be used for analysis.
    if use_relaxed_global == 'auto':
        use_relaxed_global = False
        for pdb_path_loaded in list_pdb_path:
            if '_relaxed_' in pdb_path_loaded:
                use_relaxed_global = True


    # measure the pdb files one by one
    for pdb_path_loaded in list_pdb_path:
        #cmd.load(pdb_path)
        pdb_path = pdb_path_loaded
        quantify_peptide_binding_in_pdb(glycine_linkers=glycine_linkers)
        #cmd.do('delete all')

    print("APPRAISE> Finished! Results are saved in {}".format(database_path))

    return





def main():
    # Identify the location of the additional arguments
    for i, arg in enumerate(sys.argv):
        if '.py' in arg:
            starting_index = i

    # Read the additional argments as input folders
    if len(sys.argv) > starting_index + 1:
        folder_path_list = sys.argv[starting_index + 1:]

    # get a time stamp to name the output file
    dt_string = datetime.now().strftime("%m%d%H%M")

    # Loop through the folders
    for folder_path in folder_path_list:
        print("APPRAISE> Processing folder {} using default settings.".format(folder_path))
        quantify_results_folder(folder_path, dt_string=dt_string)

    # Report success
    if len(folder_path_list) > 1:
        print("APPRAISE> Finished processing all folders in {}.".format(folder_path_list))
    print("APPRAISE> Note: These PDBs were analyzed using default settings. If you"
        + " wish to use non-standard settings, please run the script in PyMOL GUI "
        + "and call function quantify_results_folder() with appropriate arguments.")
            # # Read output path
            # if len(sys.argv) > starting_index + 2:
            #     output_path = sys.argv[starting_index + 2]
            # else:
            #     output_path = 'auto'

            # # Read output path
            # if len(sys.argv) > starting_index + 2:
            #     output_path = sys.argv[starting_index + 2]
            # else:
            #     output_path = 'auto'



if __name__ == '__main__':
    print("""This is a PyMOL script that cannot be directly launched with Python.
If you have installed PyMOL, below are the recommended ways to runn the script:
*******************************
Option 1: Run the quantification script in pymol prompt
    In PyMOL, use the following lines (the can be adjusted):

    # Load the script (replace "/path/to/APPRAISE" with the actual path)
    run /path/to/APPRAISE/appraise/pymol_quantify_peptide_binding.py

    # Call the quantification function (change the parameters as needed)
    quantify_binding('path_to_results_folder/', use_relaxed=False, time_stamp=True, mod_start_resi=3, mod_end_resi=9, pLDDT_threshold=0, membrane_anchor_site='C-term')
*******************************

*******************************
Option 2: Run the quantification script in terminal or JuPyter notebook.
    The pymol script can also be launched from the terminal (https://www.pymolwiki.org/index.php/Launching_From_a_Script).
    General syntax:

    pymol -cq /path/to/APPRAISE/appraise/pymol_quantify_peptide_binding.py /path/to/results/folder/

    You will need to replace "/path/to/APPRAISE" with the actual path. You might also need to change "pymol" to the actual location of the executable, depending on your operation system and PyMOL release.
*******************************

""")

if __name__ == 'pymol':
    from pymol import cmd
    from pymol import stored

    #If the script is run from command line with "pymol -cq"
    if len(sys.argv) > 1:
        main()
        cmd.quit(0)
    #If the script is run from PyMOL GUI
    else:
        # Allow the main functions to be run in PyMol syntax
        cmd.extend("quantify_results_folder", quantify_results_folder)
        cmd.extend('count_clash', count_clash)



# def main():
#     # Create an ArgumentParser object
#     parser = argparse.ArgumentParser()
#
#     # Add arguments to the parser
#     parser.add_argument('-i', '--input', type=str, help='input directory(ies) to be processed')
#     parser.add_log_argument('-o', '--output', type=str, help='output file', required=False)
#
#     # Parse the arguments
#     args = parser.parse_args()
#
#     # Expand the wildcard and get a list of all matching directories
#     directories = glob.glob(args.input)
#
#     # Print the list of files
#     print(f'> directories to be processed: {directories}')
#     if args.output:
#         print(f'> output: {args.output}')
#
#     # Call the quantification function
#     for directory in directories:
#         quantify_results_folder(directory, output_path = args.output)
#
#
