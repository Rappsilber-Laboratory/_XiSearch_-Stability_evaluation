"""
author: Henning Schiebenhoefer

This Script takes the internal TT and in between TT from each xifdr summary file and plots them in a stacked plot.
For values from the same sample set, averages and standard deviation of internal and between are calculated.

"""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def get_list_of_files(location, file_regex=r"FDR_.*_false_summary_xiFDR(\d+\.)*csv"):
    """generates a list of files, that satisfy specific conditions, such as filename and location
    INPUT: constraints
    RETURNS a list with all the experiment files as values"""
    list_of_files = []
    regex = re.compile(file_regex)
    for rel_dir, sub_dirs, files in os.walk(location):
        for f in files:
            if regex.match(f):
                list_of_files.append(os.path.join(rel_dir, f))
    return list_of_files

# # Test cases
# files = get_list_of_files("../../170317_OCCM_random_fasta_analysis/", r"FDR_.*_false_summary_xiFDR(\d+\.)*csv")
# for f in files:
#     print f
# print os.path.commonprefix(files)
# raw_path = r"../../occm-proteins_14_e-coli-proteins_10/xifdr_output/3/FDR_1.000000_0.000500_1.000000_1.000000_1.000000_10000.000000_false_summary_xiFDR1.0.12.csv"
# result_path = raw_path
# for i in range(3):
#     result_path = os.path.split(result_path)[0]
# print result_path


def convert_list_to_dict_of_files(list_of_files):
    exp_name_last_run = ""
    dict_of_files = {}
    for f in list_of_files:
        # read out folder of experiment
        result_path = f
        for i in range(3):
            result_path = os.path.split(result_path)[0]
        exp_name_this_run = os.path.split(result_path)[1]
        if exp_name_this_run == exp_name_last_run:
            dict_of_files[exp_name_this_run].append(f)
        else:
            dict_of_files[exp_name_this_run] = [f]
        exp_name_last_run = exp_name_this_run
    return dict_of_files


# # test cases
# file_list = get_list_of_files("../../170317_OCCM_random_fasta_analysis/", r"FDR_.*_false_summary_xiFDR(\d+\.)*csv")
# file_dict = convert_list_to_dict_of_files(file_list)
# print file_dict
# for key in file_dict:
#     print "key: %s , value: %s" % (key, file_dict[key])


def read_values_out_of_files(dict_of_files):
    """reads each file into a pandas object, from this object desired values are read and put into an object
    INPUT: dict of files, specifier for values to read
    OUTPUT: an object that contains: experiment, run number, internal TT and between TT"""
    # find out number of runs
    number_of_runs = 0
    for key in dict_of_files:
        if len(dict_of_files[key]) > number_of_runs:
            number_of_runs = len(dict_of_files[key])
    # construct name of columns
    column_names = []
    regex_exp_specifier = re.compile(r"(.+)_(\d+)_(.+)_(\d+)")
    exp_specifier = regex_exp_specifier.match(dict_of_files.keys()[0])
    column_names.extend(exp_specifier.group(1, 3))
    internal_tt_raw_string = "run {}: internal TT"
    between_tt_raw_string = "run {}: between TT"
    for i in range(1, number_of_runs+1):
        column_names.extend([internal_tt_raw_string.format(i), between_tt_raw_string.format(i)])

    # construction of data frame to calculate on
    value_data_frame = pd.DataFrame(index=dict_of_files, columns=column_names)
    # value_data_frame.sort_index(inplace=True)
    # print value_data_frame
    for exper_key in dict_of_files:
        # read out numbers of proteins and write to their specific column
        exp_specifier = regex_exp_specifier.match(exper_key)
        value_data_frame.loc[exper_key, exp_specifier.group(1)] = exp_specifier.group(2)
        value_data_frame.loc[exper_key, exp_specifier.group(3)] = exp_specifier.group(4)
        for i, f in enumerate(dict_of_files[exper_key]):
            csv_file = pd.read_csv(f, names=range(10))
            value_data_frame.loc[exper_key, internal_tt_raw_string.format(i + 1)] = csv_file.iloc[18, 2]
            value_data_frame.loc[exper_key, between_tt_raw_string.format(i + 1)] = csv_file.iloc[18, 5]
            # print i
            # TODO
    # conversion of columns to numeric
    value_data_frame = value_data_frame.apply(pd.to_numeric)
    value_data_frame.sort_values([exp_specifier.group(1), exp_specifier.group(3)], inplace=True)
    # for i in value_data_frame.columns:
    #     print value_data_frame[i].dtype
    return value_data_frame


# # test cases
# file_list = get_list_of_files("../../170317_OCCM_random_fasta_analysis/")
# file_dict = convert_list_to_dict_of_files(file_list)
# # print file_dict.keys()[0]
# values_data_frame = read_values_out_of_files(file_dict)


def calculate_values_of_interest(value_data_frame):
    """takes the object generated by 'read_values_out_of_files' and calculates mean and sdev for internal and between TT
    for each experiment.
    """
    # for each row calculate sd and mean of internal and between TT
    internal_tt_df = value_data_frame.filter(regex=re.compile("internal tt", re.I))
    between_tt_df = value_data_frame.filter(regex=re.compile("between tt", re.I))
    value_data_frame['mean internal TT'] = internal_tt_df.mean(1)
    value_data_frame['std internal TT'] = internal_tt_df.std(1)
    value_data_frame['mean between TT'] = between_tt_df.mean(1)
    value_data_frame['std between TT'] = between_tt_df.std(1)
    return value_data_frame


def plot_stacked_barplot(df_calculated_values):
    """
    y-axis: number of internal and between TT crosslinks stacked on each other
    x-axis: sample ID
    :return:
    """
    # N = 5
    n = len(df_calculated_values.index)
    # menMeans = (20, 35, 30, 35, 27)
    # print df_calculated_values.columns
    internal_tt_means = df_calculated_values['mean internal TT']
    # womenMeans = (25, 32, 34, 20, 25)
    between_tt_means = df_calculated_values['mean between TT']
    # menStd = (2, 3, 4, 1, 2)
    internal_tt_std = df_calculated_values['std internal TT']
    # womenStd = (3, 5, 2, 3, 3)
    between_tt_std = df_calculated_values['std between TT']
    # ind = np.arange(n)  # the x locations for the groups
    ind = np.arange(n)
    # print range(n)
    width = 0.35  # the width of the bars: can also be len(x) sequence
    #
    errorbar_capsize = 2
    # p1 = plt.bar(ind, menMeans, width, color='#d62728', yerr=menStd)
    p1 = plt.bar(ind, internal_tt_means, width, color='#d62728', yerr=internal_tt_std, capsize=errorbar_capsize)
    # p2 = plt.bar(ind, womenMeans, width,
    #              bottom=menMeans, yerr=womenStd)
    p2 = plt.bar(ind, between_tt_means, width, bottom=internal_tt_means, yerr=between_tt_std, capsize=errorbar_capsize)
    #
    plt.xlabel('OCCM DB size:\nE.coli DB size:')
    # plt.ylabel('Scores')
    plt.ylabel('no. of peptide links')
    # plt.title('Scores by group and gender')
    plt.title('xiFDR peptide links as a function of fasta DB size (5% pepfdr, n=3)')
    # plt.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
    xticks = []
    # print df_calculated_values.columns
    for index, row in df_calculated_values.iterrows():
        # xticks.append("{} OCCM, {} E.coli".format(int(row['occm-proteins']), int(row['e-coli-proteins'])))
        xticks.append("{}\n{}".format(int(row['occm-proteins']), int(row['e-coli-proteins'])))
    plt.xticks(ind, xticks, rotation=0)
    # plt.yticks(np.arange(0, 81, 10))
    # plt.legend((p1[0], p2[0]), ('Men', 'Women'))
    plt.legend((p1[0], p2[0]), ('internal TT', 'between TT'))
    plt.subplots_adjust(bottom=0.2)
    #
    # plt.show()
    plt.show()


def runner(file_location, file_regex=r"FDR_.*_false_summary_xiFDR(\d+\.)*csv"):
    # get list of files
    file_list = get_list_of_files(file_location, file_regex)
    # convert list to dict
    file_dict = convert_list_to_dict_of_files(file_list)
    # read out values from files
    values_data_frame = read_values_out_of_files(file_dict)
    # calculate mean and sd for internal and between TT for each experiment
    df_with_calculated_values = calculate_values_of_interest(values_data_frame)
    plot_stacked_barplot(df_with_calculated_values)


# # Test cases
# runner("../../170317_OCCM_random_fasta_analysis/")
runner("/home/henning/mnt/xitu/170323_OCCM_xi_and_xifdr_pipeline/results_with_optional_runs_and_larger_decoyDB/")

# todo
