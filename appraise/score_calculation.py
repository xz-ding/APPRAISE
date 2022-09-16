"""
Score calculation functions for appraise package and the demo notebook.
Author: Xiaozhe Ding
Email: xding@caltech.edu, dingxiaozhe@gmail.com
"""
import pandas as pd
import numpy as np

def calculate_scores(df_measurements, version=1.2, angle_constraint=True, \
    insert_constraint=True, direction_constraint=True, depth_constraint=True):
    """
    Function for getting interactive input in the notebook.
    """

    # For backward compatibility with measuredments from previous versions
    if 'angle_between_C_terminus_and_peptide' in df_measurements.columns:
        df_measurements['angle_between_membrane_anchor_and_peptide'] = df_measurements['angle_between_C_terminus_and_peptide']
        df_measurements['angle_between_membrane_anchor_and_competitor_peptide'] = df_measurements['angle_between_C_terminus_and_competitor_peptide']


    if version == 1:

        ## Calculate interface energy score (APPRAISE-1.0)
        df_measurements['interface_energy_score'] = df_measurements['total_contact_atom_in_interface'] - 1000 * df_measurements['clash_number']
        df_measurements['interface_energy_score'] = df_measurements['interface_energy_score'] * (df_measurements['interface_energy_score'] > 0)
        df_measurements['interface_energy_score_competitor'] = df_measurements['total_contact_atom_in_interface_competitor'] - 1000 * df_measurements['clash_number_competitor']
        df_measurements['interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor'] * (df_measurements['interface_energy_score_competitor'] > 0)
        df_measurements['interface_energy_score_difference'] = df_measurements['interface_energy_score'] - df_measurements['interface_energy_score_competitor']

        #Summarize the results
        df_measurements['B_POI'] = df_measurements['interface_energy_score']
        df_measurements['B_competitor'] = df_measurements['interface_energy_score_competitor']
        df_measurements['Delat_B'] = df_measurements['interface_energy_score_difference']

    if version == 1.1:
        #Development note: v1.1 ranks 9P36 to #2, Ly6a ROC AUC=0.92, but it doesn't give much signal in TTD rankings.

        ## Calculate interface energy score (APPRAISE-1.0)
        df_measurements['interface_energy_score'] = df_measurements['total_contact_atom_in_interface'] - 1000 * df_measurements['clash_number']
        df_measurements['interface_energy_score'] = df_measurements['interface_energy_score'] * (df_measurements['interface_energy_score'] > 0)
        df_measurements['interface_energy_score_competitor'] = df_measurements['total_contact_atom_in_interface_competitor'] - 1000 * df_measurements['clash_number_competitor']
        df_measurements['interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor'] * (df_measurements['interface_energy_score_competitor'] > 0)
        df_measurements['interface_energy_score_difference'] = df_measurements['interface_energy_score'] - df_measurements['interface_energy_score_competitor']

        ## Calculate constrained interface energy score (APPRAISE-1.1)
        #calculate constrained interface energy score for the main peptide
        df_measurements['constrained_interface_energy_score'] = df_measurements['interface_energy_score']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score'] -= (np.cos(df_measurements['angle_between_membrane_anchor_and_peptide']) + 1)**6
        if insert_constraint == True:
            df_measurements['constrained_interface_energy_score'] -= (df_measurements['peptide_tip_receptor_distance'] / 5)**6
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['peptide_tip_receptor_distance_competitor'] = df_measurements['peptide_tip_receptor_distance'] + df_measurements['peptide_tip_receptor_distance_difference']
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (np.cos(df_measurements['angle_between_membrane_anchor_and_competitor_peptide']) + 1)**6
        if insert_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (df_measurements['peptide_tip_receptor_distance_competitor'] / 5)**6
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']

        #Summarize the results
        df_measurements['B_POI'] = df_measurements['constrained_interface_energy_score']
        df_measurements['B_competitor'] = df_measurements['constrained_interface_energy_score_competitor']
        df_measurements['Delat_B'] = df_measurements['constrained_interface_energy_score_difference']

    if version == 1.2:

        # Works pretty well across receptors -- except PHP.C1 ranking #2
        ## Calculate interface energy score (APPRAISE-1.0)
        df_measurements['interface_energy_score'] = df_measurements['total_contact_atom_in_interface'] - 1000 * df_measurements['clash_number']
        df_measurements['interface_energy_score'] = df_measurements['interface_energy_score'] * (df_measurements['interface_energy_score'] > 0)
        df_measurements['interface_energy_score_competitor'] = df_measurements['total_contact_atom_in_interface_competitor'] - 1000 * df_measurements['clash_number_competitor']
        df_measurements['interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor'] * (df_measurements['interface_energy_score_competitor'] > 0)
        df_measurements['interface_energy_score_difference'] = df_measurements['interface_energy_score'] - df_measurements['interface_energy_score_competitor']

        ## Calculate constrained interface energy score (APPRAISE-1.1 and 1.2)
        #calculate constrained interface energy score for the main peptide
        df_measurements['constrained_interface_energy_score'] = df_measurements['interface_energy_score']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score'] -= (df_measurements['angle_between_membrane_anchor_and_peptide'] < np.pi/3) * (1 - df_measurements['angle_between_membrane_anchor_and_peptide'] / 3 * np.pi)**6 * 1000
        if depth_constraint == True:
            df_measurements['pocket_relative_depth'] = df_measurements['receptor_Rminor'] - df_measurements['peptide_tip_receptor_distance']
            df_measurements['constrained_interface_energy_score'] += df_measurements['pocket_relative_depth']**3 / 10
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (df_measurements['angle_between_membrane_anchor_and_competitor_peptide'] < np.pi/3) * (1 - df_measurements['angle_between_membrane_anchor_and_competitor_peptide'] / 3 * np.pi)**6 * 1000
        if depth_constraint == True:
            df_measurements['peptide_tip_receptor_distance_competitor'] = df_measurements['peptide_tip_receptor_distance'] + df_measurements['peptide_tip_receptor_distance_difference']
            df_measurements['pocket_relative_depth_competitor'] = df_measurements['receptor_Rminor'] - df_measurements['peptide_tip_receptor_distance_competitor']
            df_measurements['constrained_interface_energy_score_competitor'] += df_measurements['pocket_relative_depth_competitor']**3 / 10
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']

        #Summarize the results
        df_measurements['B_POI'] = df_measurements['constrained_interface_energy_score']
        df_measurements['B_competitor'] = df_measurements['constrained_interface_energy_score_competitor']
        df_measurements['Delat_B'] = df_measurements['constrained_interface_energy_score_difference']

    return df_measurements.copy()
