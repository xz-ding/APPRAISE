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


    if version == 0:

        ## Use peptide_receptor distance as surrogate control
        df_measurements['constrained_interface_energy_score_difference'] = 1000 / df_measurements['peptide_receptor_distance']

    if version == 1:

        ## Calculate interface energy score (APPRAISE-1.0)
        df_measurements['interface_energy_score'] = df_measurements['total_contact_atom_in_interface'] - 1000 * df_measurements['clash_number']
        df_measurements['interface_energy_score'] = df_measurements['interface_energy_score'] * (df_measurements['interface_energy_score'] > 0)
        df_measurements['interface_energy_score_competitor'] = df_measurements['total_contact_atom_in_interface_competitor'] - 1000 * df_measurements['clash_number_competitor']
        df_measurements['interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor'] * (df_measurements['interface_energy_score_competitor'] > 0)
        df_measurements['interface_energy_score_difference'] = df_measurements['interface_energy_score'] - df_measurements['interface_energy_score_competitor']


    if version == 1.05:
        #Development note: v1.05: uses angle factor only

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
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (np.cos(df_measurements['angle_between_membrane_anchor_and_competitor_peptide']) + 1)**6
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']



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

    if version == 1.12:
        # Use hydrodynamic radius (based on 1.1)
        # May be overfitted
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
            df_measurements['constrained_interface_energy_score'] -= (df_measurements['angle_between_membrane_anchor_and_peptide'] < np.pi * 3 / 4) * (np.pi * 3 / 4 - df_measurements['angle_between_membrane_anchor_and_peptide'])**6
        if insert_constraint == True:
            df_measurements['pocket_relative_depth'] = df_measurements['receptor_Rminor'] - df_measurements['peptide_tip_receptor_distance']
            df_measurements['constrained_interface_energy_score'] -= (4 - df_measurements['pocket_relative_depth']/5)**6
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (df_measurements['angle_between_membrane_anchor_and_competitor_peptide'] < np.pi * 3 / 4) * (np.pi * 3 / 4 - df_measurements['angle_between_membrane_anchor_and_competitor_peptide'])**6
        if insert_constraint == True:
            df_measurements['peptide_tip_receptor_distance_competitor'] = df_measurements['peptide_tip_receptor_distance'] + df_measurements['peptide_tip_receptor_distance_difference']
            df_measurements['pocket_relative_depth_competitor'] = df_measurements['receptor_Rminor'] - df_measurements['peptide_tip_receptor_distance_competitor']
            df_measurements['constrained_interface_energy_score_competitor'] -= (4 - df_measurements['pocket_relative_depth_competitor']/5)**6
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']


    if version == 1.15:
        #Development note: v1.15 works well for Ly6a, giving 0.95 ROC AUC (same as angle factor only), but it doesn't rank 9P36 high, possibly due to folded conformation

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
            df_measurements['constrained_interface_energy_score'] -= ((df_measurements['weighted_peptide_contact_distance'])/10)**6
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (np.cos(df_measurements['angle_between_membrane_anchor_and_competitor_peptide']) + 1)**6
        if insert_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= ((df_measurements['weighted_peptide_contact_distance_competitor'])/10)**6
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']

    if version == 1.16:
        # Peptide direction factor: useful: improved ROC AUC from 0.92 (v1.05) to 0.95 without sacraficing anything
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
        if direction_constraint == True:
            df_measurements['constrained_interface_energy_score'] -= ((df_measurements['peptide_direction'] < 0) * df_measurements['peptide_direction'])**6
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (np.cos(df_measurements['angle_between_membrane_anchor_and_competitor_peptide']) + 1)**6
        if direction_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= ((df_measurements['peptide_direction_competitor'] < 0) * df_measurements['peptide_direction_competitor'])**6

        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']

    if version == 1.17:
        # Use absolute depth
        # May be overfitted
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
            df_measurements['pocket_relative_depth'] = df_measurements['receptor_Rg'] - df_measurements['peptide_tip_receptor_distance']
            df_measurements['constrained_interface_energy_score'] -= ((df_measurements['pocket_relative_depth'] < 10) * (10 - df_measurements['pocket_relative_depth']))**2
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (np.cos(df_measurements['angle_between_membrane_anchor_and_competitor_peptide']) + 1)**6
        if insert_constraint == True:
            df_measurements['peptide_tip_receptor_distance_competitor'] = df_measurements['peptide_tip_receptor_distance'] + df_measurements['peptide_tip_receptor_distance_difference']
            df_measurements['pocket_relative_depth_competitor'] = df_measurements['receptor_Rg'] - df_measurements['peptide_tip_receptor_distance_competitor']
            df_measurements['constrained_interface_energy_score_competitor'] -= ((df_measurements['pocket_relative_depth_competitor'] < 10) * (10 - df_measurements['pocket_relative_depth_competitor']))**2

        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']

    if version == 1.175:
        # Use hydrodynamic radius
        # May be overfitted
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
            df_measurements['constrained_interface_energy_score'] -= (df_measurements['angle_between_membrane_anchor_and_peptide'] < np.pi * 3 / 4) * (np.pi * 3 / 4 - df_measurements['angle_between_membrane_anchor_and_peptide'])**6
        if insert_constraint == True:
            df_measurements['pocket_relative_depth'] = df_measurements['receptor_Rminor'] - df_measurements['peptide_tip_receptor_distance']
            df_measurements['constrained_interface_energy_score'] -= (1/np.maximum(0.1 * np.ones(len(df_measurements)), df_measurements['pocket_relative_depth']) + 1)**6
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (df_measurements['angle_between_membrane_anchor_and_competitor_peptide'] < np.pi * 3 / 4) * (np.pi * 3 / 4 - df_measurements['angle_between_membrane_anchor_and_competitor_peptide'])**6
        if insert_constraint == True:
            df_measurements['peptide_tip_receptor_distance_competitor'] = df_measurements['peptide_tip_receptor_distance'] + df_measurements['peptide_tip_receptor_distance_difference']
            df_measurements['pocket_relative_depth_competitor'] = df_measurements['receptor_Rminor'] - df_measurements['peptide_tip_receptor_distance_competitor']
            df_measurements['constrained_interface_energy_score_competitor'] -= (1/np.maximum(0.1 * np.ones(len(df_measurements)), df_measurements['pocket_relative_depth_competitor']))**6

        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']


    if version == 1.18:
        # Gives 1.00 ROC AUC for Ly6a, ranks 9P36 to #2, but a bit arbitray -- but bad at Ly6c1,  also gives low-confidence ranking in VYGR variants
        # May be overfitted
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
            df_measurements['pocket_relative_depth'] = 1 - df_measurements['peptide_tip_receptor_distance'] / df_measurements['receptor_Rg']
            df_measurements['constrained_interface_energy_score'] -= ((1 - df_measurements['pocket_relative_depth'])*4)**6
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (np.cos(df_measurements['angle_between_membrane_anchor_and_competitor_peptide']) + 1)**6
        if insert_constraint == True:
            df_measurements['peptide_tip_receptor_distance_competitor'] = df_measurements['peptide_tip_receptor_distance'] + df_measurements['peptide_tip_receptor_distance_difference']
            df_measurements['pocket_relative_depth_competitor'] = 1 - df_measurements['peptide_tip_receptor_distance_competitor'] / df_measurements['receptor_Rg']
            df_measurements['constrained_interface_energy_score_competitor'] -= ((1 - df_measurements['pocket_relative_depth_competitor'])*4)**6

        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']

    if version == 1.19:
        # Combine v1.16 and 1.18
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
        if direction_constraint == True:
            df_measurements['constrained_interface_energy_score'] -= ((df_measurements['peptide_direction'] < 0) * df_measurements['peptide_direction'])**6
        if insert_constraint == True:
            df_measurements['pocket_relative_depth'] = 1 - df_measurements['peptide_tip_receptor_distance'] / df_measurements['receptor_Rg']
            df_measurements['constrained_interface_energy_score'] -= ((1 - df_measurements['pocket_relative_depth'])*4)**6

        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (np.cos(df_measurements['angle_between_membrane_anchor_and_competitor_peptide']) + 1)**6
        if direction_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= ((df_measurements['peptide_direction_competitor'] < 0) * df_measurements['peptide_direction_competitor'])**6
        if insert_constraint == True:
            df_measurements['peptide_tip_receptor_distance_competitor'] = df_measurements['peptide_tip_receptor_distance'] + df_measurements['peptide_tip_receptor_distance_difference']
            df_measurements['pocket_relative_depth_competitor'] = 1 - df_measurements['peptide_tip_receptor_distance_competitor'] / df_measurements['receptor_Rg']
            df_measurements['constrained_interface_energy_score_competitor'] -= ((1 - df_measurements['pocket_relative_depth_competitor'])*4)**6


        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']

    if version == 1.2 or version == 1.22:

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

    if version == 1.21:
        # Works pretty well across receptors -- except PHP.C1 ranking #2
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
            df_measurements['constrained_interface_energy_score'] -= (df_measurements['angle_between_membrane_anchor_and_peptide'] < np.pi/2) * (np.pi/2 - df_measurements['angle_between_membrane_anchor_and_peptide'])**6
        if insert_constraint == True:
            df_measurements['pocket_relative_depth'] = df_measurements['receptor_Rminor'] - df_measurements['peptide_tip_receptor_distance']
            df_measurements['constrained_interface_energy_score'] += df_measurements['pocket_relative_depth']**3 / 10
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (df_measurements['angle_between_membrane_anchor_and_competitor_peptide'] < np.pi/2) * (np.pi/2 - df_measurements['angle_between_membrane_anchor_and_competitor_peptide'])**6
        if insert_constraint == True:
            df_measurements['peptide_tip_receptor_distance_competitor'] = df_measurements['peptide_tip_receptor_distance'] + df_measurements['peptide_tip_receptor_distance_difference']
            df_measurements['pocket_relative_depth_competitor'] = df_measurements['receptor_Rminor'] - df_measurements['peptide_tip_receptor_distance_competitor']
            df_measurements['constrained_interface_energy_score_competitor'] += df_measurements['pocket_relative_depth_competitor']**3 / 10
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']


    if version == 1.23:
        # Works pretty well across receptors -- except PHP.C1 ranking #2
        ## Calculate interface energy score (APPRAISE-1.0)
        df_measurements['interface_energy_score'] = df_measurements['total_contact_atom_in_interface'] - 1000 * df_measurements['clash_number']
        df_measurements['interface_energy_score'] = df_measurements['interface_energy_score'] * (df_measurements['interface_energy_score'] > 0)
        df_measurements['interface_energy_score_competitor'] = df_measurements['total_contact_atom_in_interface_competitor'] - 1000 * df_measurements['clash_number_competitor']
        df_measurements['interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor'] * (df_measurements['interface_energy_score_competitor'] > 0)
        df_measurements['interface_energy_score_difference'] = df_measurements['interface_energy_score'] - df_measurements['interface_energy_score_competitor']

        df_measurements['Rmajor'] = df_measurements['Dmax'] / 2
        df_measurements['Rminor'] = df_measurements['Dmax'] / df_measurements['AxialRatio'] / 2

        ## Calculate constrained interface energy score (APPRAISE-1.1)
        #calculate constrained interface energy score for the main peptide
        df_measurements['constrained_interface_energy_score'] = df_measurements['interface_energy_score']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score'] -= (df_measurements['angle_between_membrane_anchor_and_peptide'] < np.pi/3) * (1 - df_measurements['angle_between_membrane_anchor_and_peptide'] / 3 * np.pi)**6 * 1000
        if depth_constraint == True:
            df_measurements['Rtheta'] = df_measurements['Rminor'] * df_measurements['Rmajor'] / np.sqrt(np.sin(df_measurements['angle_between_membrane_anchor_and_peptide'])**2 * df_measurements['Rminor']**2 + np.cos(df_measurements['angle_between_membrane_anchor_and_peptide'])**2 * df_measurements['Rmajor']**2)
            df_measurements['pocket_relative_depth'] = df_measurements['Rtheta'] - df_measurements['peptide_tip_receptor_distance']
            df_measurements['constrained_interface_energy_score'] += df_measurements['pocket_relative_depth']**3 / 10
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] -= (df_measurements['angle_between_membrane_anchor_and_competitor_peptide'] < np.pi/3) * (1 - df_measurements['angle_between_membrane_anchor_and_competitor_peptide'] / 3 * np.pi)**6 * 1000
        if depth_constraint == True:
            df_measurements['Rtheta_competitor'] = df_measurements['Rminor'] * df_measurements['Rmajor'] / np.sqrt(np.sin(df_measurements['angle_between_membrane_anchor_and_competitor_peptide'])**2 * df_measurements['Rminor']**2 + np.cos(df_measurements['angle_between_membrane_anchor_and_competitor_peptide'])**2 * df_measurements['Rmajor']**2)
            df_measurements['peptide_tip_receptor_distance_competitor'] = df_measurements['peptide_tip_receptor_distance'] + df_measurements['peptide_tip_receptor_distance_difference']
            df_measurements['pocket_relative_depth_competitor'] = df_measurements['receptor_Rminor'] - df_measurements['peptide_tip_receptor_distance_competitor']
            df_measurements['constrained_interface_energy_score_competitor'] += df_measurements['pocket_relative_depth_competitor']**3 / 10
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']

    return df_measurements.copy()
