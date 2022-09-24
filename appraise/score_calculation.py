"""
Score calculation functions for appraise package and the demo notebook.
Author: Xiaozhe Ding
Email: xding@caltech.edu, dingxiaozhe@gmail.com
"""
import pandas as pd
import numpy as np


def calculate_B_energetic(N_contact, N_clash=0, penalize_clash=True):
    """
    calcualte B_energetic
    """
    return (N_contact - 1000 * N_clash * penalize_clash) * (N_contact - 1000 * N_clash * penalize_clash> 0)

def calculate_B_angle(theta):
    """
    calcualte B_angle
    """
    return -1 * (theta < np.pi/2) * (1 - theta / (np.pi/2))**10 * 1000

def calculate_relative_depth(distance, R_minor):
    return (R_minor-distance)/R_minor

def calculate_B_depth(d):
    """
    calcualte B_depth
    """
    return d**3 * 100




def calculate_scores(df_measurements, version=1.2, penalize_clash=True, angle_constraint=True, \
    direction_constraint=True, depth_constraint=True):
    """
    Function for getting interactive input in the notebook.
    """

    # For backward compatibility with measuredments from previous versions
    if 'angle_between_C_terminus_and_peptide' in df_measurements.columns:
        df_measurements['angle_between_membrane_anchor_and_peptide'] = df_measurements['angle_between_C_terminus_and_peptide']
        df_measurements['angle_between_membrane_anchor_and_competitor_peptide'] = df_measurements['angle_between_C_terminus_and_competitor_peptide']


    if version == 1:

        ## Calculate interface energy score (APPRAISE-1.0)
        df_measurements['interface_energy_score'] = calculate_B_energetic(df_measurements['total_contact_atom_in_interface'], df_measurements['clash_number'], penalize_clash)
        df_measurements['interface_energy_score_competitor'] = calculate_B_energetic(df_measurements['total_contact_atom_in_interface_competitor'], df_measurements['clash_number_competitor'], penalize_clash)
        df_measurements['interface_energy_score_difference'] = df_measurements['interface_energy_score'] - df_measurements['interface_energy_score_competitor']

        #Summarize the results
        df_measurements['B_POI'] = df_measurements['interface_energy_score']
        df_measurements['B_competitor'] = df_measurements['interface_energy_score_competitor']
        df_measurements['Delta_B'] = df_measurements['interface_energy_score_difference']

    if version == 1.1:
        #Development note: v1.1 ranks 9P36 to #2, Ly6a ROC AUC=0.92, but it doesn't give much signal in TTD rankings.

        ## Calculate interface energy score (APPRAISE-1.0)
        df_measurements['interface_energy_score'] = calculate_B_energetic(df_measurements['total_contact_atom_in_interface'], df_measurements['clash_number'], penalize_clash)
        df_measurements['interface_energy_score_competitor'] = calculate_B_energetic(df_measurements['total_contact_atom_in_interface_competitor'], df_measurements['clash_number_competitor'], penalize_clash)
        df_measurements['interface_energy_score_difference'] = df_measurements['interface_energy_score'] - df_measurements['interface_energy_score_competitor']

        ## Calculate constrained interface energy score (APPRAISE-1.1)
        #calculate constrained interface energy score for the main peptide
        df_measurements['constrained_interface_energy_score'] = df_measurements['interface_energy_score']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score'] += calculate_B_angle(df_measurements['angle_between_membrane_anchor_and_peptide'])
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] += calculate_B_angle(df_measurements['angle_between_membrane_anchor_and_competitor_peptide'])
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']

        #Summarize the results
        df_measurements['B_POI'] = df_measurements['constrained_interface_energy_score']
        df_measurements['B_competitor'] = df_measurements['constrained_interface_energy_score_competitor']
        df_measurements['Delta_B'] = df_measurements['constrained_interface_energy_score_difference']

    if version == 1.2:

        # Works pretty well across receptors -- except PHP.C1 ranking #2
        ## Calculate interface energy score (APPRAISE-1.0)
        df_measurements['interface_energy_score'] = calculate_B_energetic(df_measurements['total_contact_atom_in_interface'], df_measurements['clash_number'], penalize_clash)
        df_measurements['interface_energy_score_competitor'] = calculate_B_energetic(df_measurements['total_contact_atom_in_interface_competitor'], df_measurements['clash_number_competitor'], penalize_clash)
        df_measurements['interface_energy_score_difference'] = df_measurements['interface_energy_score'] - df_measurements['interface_energy_score_competitor']

        ## Calculate constrained interface energy score (APPRAISE-1.1 and 1.2)
        #calculate constrained interface energy score for the main peptide
        df_measurements['constrained_interface_energy_score'] = df_measurements['interface_energy_score']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score'] += calculate_B_angle(df_measurements['angle_between_membrane_anchor_and_peptide'])
        if depth_constraint == True:
            df_measurements['pocket_relative_depth'] = calculate_relative_depth(df_measurements['peptide_tip_receptor_distance'], df_measurements['receptor_Rminor'])
            df_measurements['constrained_interface_energy_score'] += calculate_B_depth(df_measurements['pocket_relative_depth'])
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] += calculate_B_angle(df_measurements['angle_between_membrane_anchor_and_competitor_peptide'])
        if depth_constraint == True:
            df_measurements['peptide_tip_receptor_distance_competitor'] = df_measurements['peptide_tip_receptor_distance'] + df_measurements['peptide_tip_receptor_distance_difference']
            df_measurements['pocket_relative_depth_competitor'] = calculate_relative_depth(df_measurements['peptide_tip_receptor_distance_competitor'], df_measurements['receptor_Rminor'])
            df_measurements['constrained_interface_energy_score_competitor'] += calculate_B_depth(df_measurements['pocket_relative_depth_competitor'])
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']

        #Summarize the results
        df_measurements['B_POI'] = df_measurements['constrained_interface_energy_score']
        df_measurements['B_competitor'] = df_measurements['constrained_interface_energy_score_competitor']
        df_measurements['Delta_B'] = df_measurements['constrained_interface_energy_score_difference']

    return df_measurements.copy()
