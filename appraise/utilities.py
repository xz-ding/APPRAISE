"""
Utility functions for appraise package and the demo notebook.
Author: Xiaozhe Ding
Email: xding@caltech.edu, dingxiaozhe@gmail.com
"""
import pandas as pd
import numpy as np
import scipy
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def interactive_input(var_name='', default_value=''):
    """
    Function for getting interactive input in the notebook.
    """
    if isinstance(default_value, str):
        print('$ Variable <{}> defaults to <\'{}\'>. Change it? Type new value or hit Enter to skip.'.format(var_name, default_value))
    else:
        print('$ Variable <{}> defaults to <{}>. Change it? Type new value or hit Enter to skip.'.format(var_name, default_value))

    # Determine if the frontend allows interactive input
    if '__file__' in locals():
        # code is being run by nbconvert, use the default value
        var_value = default_value
        print('$ Frontend does not allow interactive input. Used default value.\n')
    else:
        # If the frontend allows, use interactive input
        var_value =  input("$ ")

        if var_value == '' or var_value.lower() == 'n' or var_value.lower() == 'no':
            var_value = default_value
            print('$ Used default value.\n')

    return var_value

def get_peptide_list_from_model_names(df, sorting_metric_name = 'peptide_name'):
    """
    Get a list of peptide from a database dataframe.

    """
    df['competition_name'] = df['model_name'].str.slice_replace(-15, None, '')

    list_peptide_names = []
    list_peptide_seqs = []

    for competition_name in set(df['competition_name']):
        df_competition = df[df['competition_name'] == competition_name].copy()
        df_competition_mean = df_competition.groupby(by=['peptide_name', 'peptide_seq']).mean(numeric_only=True).reset_index()
        df_competition_mean = df_competition_mean.sort_values(by=sorting_metric_name, ascending=False).reset_index()
        for i in range(len(df_competition_mean)):
            peptide_name = df_competition_mean.loc[i]['peptide_name']
            peptide_seq = df_competition_mean.loc[i]['peptide_seq']
            if peptide_name not in list_peptide_names:
                list_peptide_names.append(peptide_name)
                list_peptide_seqs.append(peptide_seq)

    return list_peptide_names, list_peptide_seqs


def sort_df_by_peptides_and_cleanup(df, list_peptide_order, consider_competitor_name=True):
    """
    Sort the a dataframe (averaged or not) by a given order of peptide names.
    """
    df['peptide_name'] = pd.Categorical(df['peptide_name'], list_peptide_order)
    if consider_competitor_name:
        if 'competitor' in df.columns or 'competitors' in df.columns:
            if 'competitors' in df.columns:
                df['competitor'] = [l[2:-2] for l in df['competitors']]
            df['competitor'] = pd.Categorical(df['competitor'], list_peptide_order)
            df_sorted = df.sort_values(by=['peptide_name', 'competitor'])
            df_sorted = df_sorted.dropna(subset=['peptide_name', 'competitor'])
        else:
            df_sorted =  df.sort_values(by=['peptide_name'])
            df_sorted = df_sorted.dropna(subset=['peptide_name'])
    else:
        df_sorted =  df.sort_values(by=['peptide_name'])
        df_sorted = df_sorted.dropna(subset=['peptide_name'])

    return df_sorted

def rank_tournament_results(df_average, metric_name='interface_energy_score_difference',\
    by_match_points=True, tie_threshold='auto', p_value_threshold=0.05, points=[1, 0, -1],\
    number_of_repeats=10):
    """
    Rank a dataframe (averaged between replicates) by counting match results and
    calculating total points.
    """

    df_average[metric_name+'_tstat'] = (df_average[metric_name] - np.mean(df_average[metric_name]))/(np.sqrt(np.var(df_average[metric_name])/number_of_repeats))

    if tie_threshold == 'auto':
        #Get the critical threshold corresponding to the t value
        degree_of_freedom = number_of_repeats - 1
        tie_threshold = scipy.stats.t.ppf(1 - p_value_threshold/2, degree_of_freedom)

        print('Used p-value threshold of {:.3f}'.format(p_value_threshold))
    print('Tie threshold to be {:.2f} of standard deviation: {:.2f}'.format(tie_threshold, np.mean(df_average[metric_name]) + tie_threshold*np.sqrt(np.var(df_average[metric_name])/number_of_repeats)))

    if by_match_points:
        if len(points) != 3:
            print("Error: Format of points has to be a list of three numbers [point_for_winner, point_for_tie, point_for_loser](e.g. [3, 1, 0])")
            return None


        list_peptide = list(set(df_average['peptide_name']))

        for peptide_name in list_peptide:
            df_peptide = df_average[df_average['peptide_name'] == peptide_name]

            #Count number of winning, losing or tie matches with threshold
            n_contact_win = np.sum(df_peptide[metric_name+'_tstat'] > tie_threshold)
            n_contact_tie = np.sum(df_peptide[metric_name+'_tstat']**2 <= tie_threshold**2)
            n_contact_lose = np.sum(df_peptide[metric_name+'_tstat'] < -tie_threshold)


            #Count number of winning, losing or tie matches with higher threshold for breaking the ties
            n_contact_win_tie_breaker = np.sum(df_peptide[metric_name+'_tstat'] > 2 * tie_threshold)
            n_contact_tie_tie_breaker = np.sum(df_peptide[metric_name+'_tstat'] == tie_threshold**2 * 4)
            n_contact_lose_tie_breaker = np.sum(df_peptide[metric_name+'_tstat'] < -2 * tie_threshold)


            #calculate points
            match_points =  points[0] * n_contact_win + points[1] * n_contact_tie + points[2] * n_contact_lose
            match_points_tie_breaker = points[0] * n_contact_win_tie_breaker + points[1] * n_contact_tie_tie_breaker + points[2] * n_contact_lose_tie_breaker

            #record the points
            df_average.loc[df_average['peptide_name'] == peptide_name, 'match_points'] = match_points
            df_average.loc[df_average['peptide_name'] == peptide_name, 'match_points_tie_breaker'] = match_points_tie_breaker


        list_peptide_order = df_average.groupby(by=['peptide_name']).mean(numeric_only=True).sort_values(by=['match_points', 'match_points_tie_breaker', metric_name], ascending=False).reset_index()['peptide_name'].to_list()
        list_match_points = [df_average.groupby(by=['peptide_name']).mean(numeric_only=True).sort_values(by=['match_points', 'match_points_tie_breaker', metric_name], ascending=False).reset_index()['match_points'].to_list(), \
                            df_average.groupby(by=['peptide_name']).mean(numeric_only=True).sort_values(by=['match_points', 'match_points_tie_breaker', metric_name], ascending=False).reset_index()['match_points_tie_breaker'].to_list()]

        return list_peptide_order, tie_threshold, list_match_points
    else:
        return df_average.groupby(by=['peptide_name']).mean(numeric_only=True).sort_values(by=[metric_name], ascending=False).reset_index()['peptide_name'].to_list(), tie_threshold, None


def plot_heatmap(df_average, feature_of_interest='Delta_B', receptor_of_interest='receptor', \
    feature_to_plot_with='auto', feature_to_rank_with='auto', fig_size='auto', \
    tie_threshold='auto', p_value_threshold=0.05, number_of_repeats=10, \
    vmin='auto', vmax='auto', title='auto', palette = "vlag_r",\
    save_figure=True, rank_by_tournament=True, print_label=False, label_dictionary=None,\
    xlabel='Competitor peptide', ylabel='Peptide of interest (POI)',\
    xticklabels='auto', yticklabels='auto'):
    """
    A function to rank and generate result heatpmaps for APPRAISE competitions.

    df_average: (pandas dataframe) a dataframe containing the mean scores.
    The data frame has to contain columns 'peptide_name', 'competitor', as well
    as the column(s) storing the scores needed for ranking and plotting. The
    score columns are specified by feature_of_interest or specified by
    feature_to_plot_with and feature_to_rank_with separately.


    feature_of_interest: (string) The feature / score that is going to be used
    to rank the candidate peptides and plot the heatmap in case
    feature_to_plot_with or feature_to_rank_with is "auto".

    receptor_of_interest: (string) The name of the receptor. Currently this is
    only used for naming the output figure.

    feature_to_plot_with: (string) The feature / score that is going to be used
    to plot the candidate peptides. If 'auto', feature_of_interest will be used.

    feature_to_rank_with: (string) The feature / score that is going to be used
    to rank the candidate peptides. If 'auto', feature_of_interest will be used.

    fig_size: (float or 'auto') figure size of square matrix. If 'auto', it will
    be determined based on the number of peptides.

    tie_threshold: (float or 'auto') the numerical threshold to call a tie for a
    competition. If 'auto', the threshold will be determined using a p value.
    If using the 'auto' mode, make sure the p_value_threshold and
    number_of_repeats are accurate.

    p_value_threshold: (float) The p value threshold (two-sided, single-sample
    Student's t-test) used to determine whether a competition is a tie.

    number_of_repeats: (float) Number of replicate models used to calculate the
    mean scores. This number will be used in the Student's t-test for tie
    threshold determination.

    vmin: (float or 'auto') lower limit for the scale bar to pass to seaborn.

    vmax: (float or 'auto') upper limit for the scale bar to pass to seaborn.

    title: (float or 'auto') title of the plot to pass to seaborn.

    palette = "vlag_r": (string) color palette of the scale bar to pass to
    seaborn.

    save_figure: (boolean) if True, the resulting heatmap will be saved as file.

    rank_by_tournament: (boolean) if True, the competition results will be
    ranked as a point-based round-robin tournament. If False, the mean value
    on each row will be used for ranking.

    print_label: (boolean) If True, the peptide names will be followed by a
    bracketed labels. The labels need to be specified by a label_dictionary.

    label_dictionary: (dictionary) A dictionary that speicifies the labels for
    each peptide.

    xlabel: (string) the label of the x axis to pass to seaborn.

    ylabel: (string) the label of the y axis to pass to seaborn.

    xticklabels: (string) the label of the x axis ticks to pass to seaborn. if
    'auto', the peptide names (optionally with bracketed labels if print_label
    is True) will be used as tick labels.

    yticklabels: (string) the label of the y axis ticks to pass to seaborn. if
    'auto', the peptide names will be used as tick labels.

    """
    sns.set_context("poster")

    #Get the feature to be used to rank and plot
    if feature_to_plot_with == 'auto':
        feature_to_plot_with = feature_of_interest
    if feature_to_rank_with == 'auto':
        feature_to_rank_with = feature_of_interest

    #Get a ranked square matrix
    list_peptide_order, tie_threshold, list_match_points = rank_tournament_results(df_average, feature_to_rank_with, by_match_points=rank_by_tournament, tie_threshold=tie_threshold, p_value_threshold=p_value_threshold, number_of_repeats=number_of_repeats)
    df_average = sort_df_by_peptides_and_cleanup(df_average, list_peptide_order)
    square_matrix = df_average.pivot(index='peptide_name', columns='competitor', values=feature_to_plot_with)

    if xticklabels == 'auto':
        xticklabels = list_peptide_order
    if yticklabels == 'auto':
        if not print_label:
            yticklabels = list_peptide_order
        else:
            if label_dictionary == None:
                print('Warning: print_label option was set to be True but no label was detected.')
                yticklabels = list_peptide_order
            else:
                sr_peptide_name = pd.Series(list_peptide_order)
                sr_labels = sr_peptide_name.map(label_dictionary)
                sr_peptide_name_with_label = sr_peptide_name + ' [' + sr_labels.astype(str) + ']'
                yticklabels = list(sr_peptide_name_with_label)

    #Plot the heat map
    if fig_size == 'auto':
        fig_size = len(square_matrix)/2.5#/3

    if vmin == 'auto':
        vmin = -2 * np.sqrt(np.var(df_average[feature_to_plot_with]))
    if vmax == 'auto':
        vmax = 2 * np.sqrt(np.var(df_average[feature_to_plot_with]))


    sns.set(rc={'figure.figsize':(fig_size,fig_size)})
    sns.heatmap(square_matrix, cmap=palette, square=True, vmin=vmin, vmax=vmax,\
                yticklabels=yticklabels, xticklabels=xticklabels)
    if title=='auto':
        title = feature_to_plot_with

    plt.title(title, size=fig_size * 2)
    plt.xlabel(xlabel, size=fig_size * 1.8)
    plt.ylabel(ylabel, size=fig_size * 1.8)
    #plt.tight_layout()
    if save_figure:
        plt.savefig('{}_ranked_by_{}.png'.format(receptor_of_interest, feature_to_rank_with), bbox_inches = 'tight', dpi=300)

    return xticklabels, yticklabels, list_match_points

def database_quality_check(df):
    # Quality check
    print('\nQuality check: \n \
    The following plot shows the number of models for each peptide variant in the database. \n \
    You should expect to see an equal number of models in all variants. \n \
    If the numbers are not equal, there might some models missing during structure prediction or quantification. You can still proceed with the risk of getting biased rankings.\n \
    Double click on the plot to Zoom in. ')

    plt.rcParams["figure.figsize"] = (len(set(df['peptide_name']))*1.5, 3)
    sns.histplot(df, x='peptide_name')
    plt.tight_layout()
