"""
Score calculation functions for appraise package and the demo notebook.
Author: Xiaozhe Ding
Email: xding@caltech.edu, dingxiaozhe@gmail.com
"""
import pandas as pd
import numpy as np
from appraise.utilities import *

def calculate_B_energetic(N_contact, N_clash=0, clash_penalty_weight=1000):
    """
    calcualte B_energetic score.

    N_contact: (float) number of contacting atoms.

    N_clash: (float) number of clashing atoms.

    clash_penalty_weight: (float) the relative weight of the penalty against clashing
    atoms vs the reward for contacting atoms.

    """
    return (N_contact - clash_penalty_weight * N_clash) * (N_contact - clash_penalty_weight * N_clash > 0)

def calculate_B_angle(theta, critical_angle = np.pi/2, angle_exponent=10, angle_penalty_weight=1000):
    """
    calcualte B_angle score.

    critical_angle: (float) the angle below which penalty starts to apply.

    angle_exponent: (float) the exponent (steepness) of the B_angle function.

    angle_penalty_weight: (float) the relative weight of the penalty against the
    impossible angles.

    """
    return -1 * (theta < critical_angle) * (1 - theta / critical_angle)**angle_exponent * angle_penalty_weight

def calculate_relative_depth(distance, R_minor):
    """
    Calcualte the relative depth used for calculating B_depth score.

    distance: (float) the distance between the atom on the peptide that is the
    closest to the receptor and the center of mass of the receptor.

    R_minor: (float) the minor radius of the ellipsoid hull of the recptor.

    """
    return (R_minor-distance)/R_minor

def calculate_B_depth(d, depth_exponent=3, depth_weight=100):
    """
    Calcualte B_depth score.

    depth_exponent: (float) the exponent  (steepness) of the B_depth function.

    angle_penalty_weight: (float) the relative weight of the depth score.

    """
    return d**depth_exponent * depth_weight




def calculate_scores(df_measurements, version=1.2, angle_constraint=True, \
    depth_constraint=True, clash_penalty_weight=1000, critical_angle = np.pi/2,
    angle_exponent=10, angle_penalty_weight=1000, depth_exponent=3,
    depth_weight=100):
    """
    Calculate the APPRAISE scores based on a dataframe of measurements and
    returns a new dataframe with new columns containing the scores.

    df_measurements: (Pandas dataframe) the input dataframe containing the
    necessary measurements.

    version: (float) the version of APPRAISE to be used for calculation.

    angle_constraint: (boolean) if True and APPRAISE version >= 1.1, then
    B_angle will be considered.

    depth_constraint: (boolean) if True and APPRAISE version >=1.2, then
    B_depths will be considered.

    clash_penalty_weight: (float) the relative weight of the penalty against
    clashing atoms vs the reward for contacting atoms.

    angle_exponent: (float) the exponent (steepness) of the B_angle function.

    angle_penalty_weight: (float) the relative weight of the penalty against the
    impossible angles.

    depth_exponent: (float) the exponent  (steepness) of the B_depth function.

    angle_penalty_weight: (float) the relative weight of the depth score.

    """

    # For backward compatibility with measuredments from previous versions
    if 'angle_between_C_terminus_and_peptide' in df_measurements.columns:
        df_measurements['angle_between_membrane_anchor_and_peptide'] = df_measurements['angle_between_C_terminus_and_peptide']
        df_measurements['angle_between_membrane_anchor_and_competitor_peptide'] = df_measurements['angle_between_C_terminus_and_competitor_peptide']


    if version == 1:

        ## Calculate interface energy score (APPRAISE-1.0)
        df_measurements['interface_energy_score'] = calculate_B_energetic(df_measurements['total_contact_atom_in_interface'], df_measurements['clash_number'], clash_penalty_weight=clash_penalty_weight)
        df_measurements['interface_energy_score_competitor'] = calculate_B_energetic(df_measurements['total_contact_atom_in_interface_competitor'], df_measurements['clash_number_competitor'], clash_penalty_weight=clash_penalty_weight)
        df_measurements['interface_energy_score_difference'] = df_measurements['interface_energy_score'] - df_measurements['interface_energy_score_competitor']

        #Summarize the results
        df_measurements['B_POI'] = df_measurements['interface_energy_score']
        df_measurements['B_competitor'] = df_measurements['interface_energy_score_competitor']
        df_measurements['Delta_B'] = df_measurements['interface_energy_score_difference']

    if version == 1.1:
        #Development note: v1.1 ranks 9P36 to #2, Ly6a ROC AUC=0.92, but it doesn't give much signal in TTD rankings.

        ## Calculate interface energy score (APPRAISE-1.0)
        df_measurements['interface_energy_score'] = calculate_B_energetic(df_measurements['total_contact_atom_in_interface'], df_measurements['clash_number'], clash_penalty_weight=clash_penalty_weight)
        df_measurements['interface_energy_score_competitor'] = calculate_B_energetic(df_measurements['total_contact_atom_in_interface_competitor'], df_measurements['clash_number_competitor'], clash_penalty_weight=clash_penalty_weight)
        df_measurements['interface_energy_score_difference'] = df_measurements['interface_energy_score'] - df_measurements['interface_energy_score_competitor']

        ## Calculate constrained interface energy score (APPRAISE-1.1)
        #calculate constrained interface energy score for the main peptide
        df_measurements['constrained_interface_energy_score'] = df_measurements['interface_energy_score']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score'] += calculate_B_angle(df_measurements['angle_between_membrane_anchor_and_peptide'], critical_angle=critical_angle, angle_exponent=angle_exponent, angle_penalty_weight=angle_penalty_weight)
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['constrained_interface_energy_score_competitor'] += calculate_B_angle(df_measurements['angle_between_membrane_anchor_and_competitor_peptide'], critical_angle=critical_angle, angle_exponent=angle_exponent, angle_penalty_weight=angle_penalty_weight)
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
        df_measurements['interface_energy_score'] = calculate_B_energetic(df_measurements['total_contact_atom_in_interface'], df_measurements['clash_number'], clash_penalty_weight=clash_penalty_weight)
        df_measurements['interface_energy_score_competitor'] = calculate_B_energetic(df_measurements['total_contact_atom_in_interface_competitor'], df_measurements['clash_number_competitor'], clash_penalty_weight=clash_penalty_weight)
        df_measurements['interface_energy_score_difference'] = df_measurements['interface_energy_score'] - df_measurements['interface_energy_score_competitor']

        ## Calculate constrained interface energy score (APPRAISE-1.1 and 1.2)
        #calculate constrained interface energy score for the main peptide
        df_measurements['constrained_interface_energy_score'] = df_measurements['interface_energy_score']
        if angle_constraint == True:
            df_measurements['B_angle'] = calculate_B_angle(df_measurements['angle_between_membrane_anchor_and_peptide'], critical_angle=critical_angle, angle_exponent=angle_exponent, angle_penalty_weight=angle_penalty_weight)
            df_measurements['constrained_interface_energy_score'] += df_measurements['B_angle']
        if depth_constraint == True:
            df_measurements['pocket_relative_depth'] = calculate_relative_depth(df_measurements['peptide_tip_receptor_distance'], df_measurements['receptor_Rminor'])
            df_measurements['B_depth'] = calculate_B_depth(df_measurements['pocket_relative_depth'], depth_exponent=depth_exponent, depth_weight=depth_weight)
            df_measurements['constrained_interface_energy_score'] += df_measurements['B_depth']
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score'] = df_measurements['constrained_interface_energy_score'] * (df_measurements['constrained_interface_energy_score'] > 0)

        #calculate constrained interface energy score for the competitor peptide
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['interface_energy_score_competitor']
        if angle_constraint == True:
            df_measurements['B_angle_competitor'] = calculate_B_angle(df_measurements['angle_between_membrane_anchor_and_competitor_peptide'], critical_angle=critical_angle, angle_exponent=angle_exponent, angle_penalty_weight=angle_penalty_weight)
            df_measurements['constrained_interface_energy_score_competitor'] += df_measurements['B_angle_competitor']
        if depth_constraint == True:
            df_measurements['peptide_tip_receptor_distance_competitor'] = df_measurements['peptide_tip_receptor_distance'] + df_measurements['peptide_tip_receptor_distance_difference']
            df_measurements['pocket_relative_depth_competitor'] = calculate_relative_depth(df_measurements['peptide_tip_receptor_distance_competitor'], df_measurements['receptor_Rminor'])
            df_measurements['B_depth_competitor'] = calculate_B_depth(df_measurements['pocket_relative_depth_competitor'], depth_exponent=depth_exponent, depth_weight=depth_weight)
            df_measurements['constrained_interface_energy_score_competitor'] += df_measurements['B_depth_competitor']
        #Set minimal score to be zero
        df_measurements['constrained_interface_energy_score_competitor'] = df_measurements['constrained_interface_energy_score_competitor'] * (df_measurements['constrained_interface_energy_score_competitor'] > 0)

        #calculate the difference
        df_measurements['constrained_interface_energy_score_difference'] = df_measurements['constrained_interface_energy_score'] - df_measurements['constrained_interface_energy_score_competitor']

        #Summarize the results
        df_measurements['B_POI'] = df_measurements['constrained_interface_energy_score']
        df_measurements['B_competitor'] = df_measurements['constrained_interface_energy_score_competitor']
        df_measurements['Delta_B'] = df_measurements['constrained_interface_energy_score_difference']

    return df_measurements.copy()

def find_rank(string, search_list):
  try:
    return search_list.index(string) + 1
  except ValueError:
    return 0

def get_APPRAISE_ranking(df_average, feature_of_interest='Delta_B', \
    receptor_of_interest='receptor', feature_to_rank_with='auto', \
    tie_threshold='auto', p_value_threshold=0.05, number_of_repeats=10):

    # Determine the feature to rank with
    if feature_to_rank_with == 'auto':
        feature_to_rank_with = feature_of_interest

    #Get a ranked square matrix
    list_peptide_order, tie_threshold, list_match_points = rank_tournament_results(df_average, feature_to_rank_with, by_match_points=True, tie_threshold=tie_threshold, p_value_threshold=p_value_threshold, number_of_repeats=number_of_repeats)
    # df_average = sort_df_by_peptides_and_cleanup(df_average, list_peptide_order)

    # # Add the match points to the dataframe
    # df_average['match_points'] = 0
    # for i, peptide in enumerate(list_peptide_order):
    #     df_average.loc[df_average['peptide_name'] == peptide]['match_points'] = list_match_points[0][i]

    # Calculate winning percentage
    df_average['winning_percentage_by_{}'.format(feature_of_interest)] = df_average['match_points'] / len(list_peptide_order) / 2 + 0.5
    df_average['APPRAISE_ranking_by_{}'.format(feature_of_interest)] = df_average['peptide_name'].apply(find_rank, args=(list_peptide_order,)).astype(int)

    return df_average
