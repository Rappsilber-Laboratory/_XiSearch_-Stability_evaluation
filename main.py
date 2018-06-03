"""
author: 'Henning Schiebenhoefer'

Input:
Peak Lists
fasta files
Xi cfg
Output:
XiFDR results
"""

from lib import pipeline
from lib.XiWrapper import XiSearchOutOfMemoryException, XiSearchDaemoniseFailureException
import time
import os
import fasta_randomizer
import subprocess
import pandas as pd
import sys  # Anzahl der Parameter beim Aufruf
import argparse     # parsen von Befehlen aus der Kommandozeile
import logging


class Experiment:
    """
    one object per MS analysis

    unique for each object:
    peak files
    iBAQ list
    optimization iBAQ quantiles
    fasta quantile files (these are different due to different iBAQ lists)

    shared between objects:
    fasta ressource file

    """
    def __init__(
            self,
            exp_name,
            list_of_peak_files,
            fasta_target="",
            fasta_decoy="",
            df_rnd_fasta_number_pairs=[]
    ):
        self.name = exp_name
        self.fasta_target = fasta_target
        self.fasta_decoy = fasta_decoy
        self.peak_files = list_of_peak_files
        self.df_rnd_fasta_number_pairs = df_rnd_fasta_number_pairs

    def __str__(self):
        return self.name


def execute_xi_xifdr_pipeline(experiment, lst_fastas, str_run, path_run, xi_xifdr_settings_dict):
    """
    helper function for pipeline execution
    """

    single_pipeline_run_starttime = time.time()
    logging.info("starting xi/xifdr pipeline for {} of experiment '{}'"
                 .format(str_run, experiment.name))
    global error_free_execution
    try:
        pipeline.execute_pipeline(
            # optional general settings
            output_basedir=path_run,
            # xi settings
            list_of_fasta_dbs=lst_fastas,
            xi_config=xi_xifdr_settings_dict['xi_config'],
            peak_files=experiment.peak_files,
            # optional xi settings
            xi_memory=xi_xifdr_settings_dict['xi_memory'],
            additional_xi_parameters=xi_xifdr_settings_dict['additional_xi_parameters'],
            # optional xifdr settings
            pepfdr=xi_xifdr_settings_dict['xifdr_settings']['pepfdr'],
            reportfactor=xi_xifdr_settings_dict['xifdr_settings']['reportfactor'],
            additional_xifdr_arguments=xi_xifdr_settings_dict['xifdr_settings']['additional_xifdr_arguments'])
    except XiSearchDaemoniseFailureException as e:
        logging.error("XiSearch threw daemonising exception for run '{}'"
                      .format(str_run, experiment.name))
        logging.error("XiSearch command: {}".format(" ".join(e.cmd)))
        logging.error("Keeping result.")
        error_free_execution = False
    except XiSearchOutOfMemoryException as e:
        logging.error("XiSearch did not complete '{}' of Experiment '{}'. Please check log."
                      .format(str_run, experiment.name))
        logging.error("XiSearch command: {}".format(" ".join(e.cmd)))
        logging.error("removing file stub {}".format(e.out_file))
        os.remove(e.out_file)
        error_free_execution = False
    logging.info("pipeline execution for {} of experiment {} took {}"
                 .format(str_run, experiment.name, pipeline.calculate_elapsed_time(single_pipeline_run_starttime)))


def generate_random_fasta_files(input_file, output_dir, desired_no):
    """
    takes a fasta file and generates a random sample file from it
    returns name of random sample file
    :param input_file:
    :param output_dir:
    :param desired_no: int
    :return: string: name of random sample file
    """
    # generate output file name based on input file name and output dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    input_file_basename = os.path.basename(input_file)
    output_file_basename = "random_fasta_" + input_file_basename
    output_file = os.path.join(output_dir, output_file_basename)
    # call fasta file generator
    generated_fasta_file = fasta_randomizer.generate_random_protein_list(input_file, output_file, desired_no)
    return generated_fasta_file

# # Test cases
# print generate_random_fasta_files(r"../fasta_randomizer/uniprot-proteome%25253AUP000000625.bin", "random_db", 100, 10)
# print os.path.splitext("OCCM_Scerevisiae.fasta")[0]
# print generate_random_fasta_files("Data/OCCM_Scerevisiae-11_17_09-09_Feb_2016/OCCM_Scerevisiae.fasta", "testing", 5)


def execute_pipeline(
        experiment,
        str_run,
        xi_xifdr_settings_dict,
        amount_of_repeats,  # general settings
        random_fasta_db_size, random_decoy_fasta_db_size,  # random fasta settings
        output_basedir=".",  # optional general settings
):
    """
    This pipeline executes:
    fasta_randomizer
    xiSearch
    and Xifdr
    for amount_of_repeats times
    :param fasta_db_target:
    :param random_fasta_db_size:
    :param fasta_db_decoy:
    :param random_decoy_fasta_db_size:
    :param amount_of_repeats:
    :param xi_config:
    :param peak_files:
    :param pepfdr:
    :param memory:
    :param reportfactor:
    :param additional_xifdr_arguments:
    :param additional_xi_parameters:
    :return:
    """
    # dir structure: exper_name/repeat/[xi_out, xifdr_out, random_fasta_dbs]
    # base dirs for all the files
    # pre_list_of_dirs = ["random_fasta_dbs", "xi_output", "xifdr_output"]
    # loop through generation of all files amount_of_repeats times
    # for loop adjusted to make i represent run number
    for i in range(1, amount_of_repeats+1):
        # generate necessary dirs for this run
        i_out = os.path.join(output_basedir, str(i))
        fasta_out = os.path.join(i_out, "random_fasta_dbs")

        lst_fastas = []
        # list_of_dirs = fun_generate_list_of_dirs(output_basedir, pre_list_of_dirs, i)
        # fun_makedirs(list_of_dirs)

        # generate necessary fastas for this run: one call for true_db, one call for random_db
        if random_fasta_db_size > 0:
            random_fasta_db_target = generate_random_fasta_files(
                input_file=experiment.fasta_target,
                output_dir=fasta_out,
                desired_no=random_fasta_db_size
            )
            lst_fastas.append(random_fasta_db_target)

        if random_decoy_fasta_db_size > 0:
            random_fasta_db_decoy = generate_random_fasta_files(
                input_file=experiment.fasta_decoy,
                output_dir=fasta_out,
                desired_no=random_decoy_fasta_db_size
            )
            lst_fastas.append(random_fasta_db_decoy)

        # call xi/xifdr pipeline
        execute_xi_xifdr_pipeline(
            experiment=experiment,
            lst_fastas=lst_fastas,
            str_run=str_run,
            path_run=i_out,
            xi_xifdr_settings_dict=xi_xifdr_settings_dict
        )

        # # call xisearch
        # xi_result = xi_execution(
        #     xi_config=xi_config,
        #     peak_files=peak_files,
        #     fasta_files=[random_fasta_db_target, random_fasta_db_decoy],
        #     output_folder=list_of_dirs[1],
        #     additional_parameters=additional_xi_parameters
        # )
        #
        # # xi_result is a string but xifdr needs a list as input
        # xifdr_input = []
        # xifdr_input.append(xi_result)
        # # call xifdr
        # xifdr_execution(
        #     xifdr_input_csv=xifdr_input,
        #     xifdr_output_dir=list_of_dirs[2],
        #     pepfdr=pepfdr,
        #     memory=memory,
        #     reportfactor=reportfactor,
        #     additional_xifdr_arguments=additional_xifdr_arguments
        # )
    return


# # Test cases
# xifdr_execution(xifdr_input_csv=['results/occm-proteins_1_e-coli-proteins_0/xi_output/1/xi_results.csv'],
#                 xifdr_output_dir='results/occm-proteins_1_e-coli-proteins_0/xifdr_output/1',
#                 pepfdr='0.05',
#                 memory='1G',
#                 reportfactor='10000',
#                 additional_xifdr_arguments=[])
#
# for i in list('results/occm-proteins_1_e-coli-proteins_0/xi_output/1/xi_results.csv'):
#     print i


def read_experiment_file(input_txt):
    """
    reads the input file
    first line is stored as description
    from second line on, first line element is stored as fasta_db size, second element as decoy_db size.
     if third number exists this is set to "override_repeats"
    reading of file with pandas
    :param input_txt:
    :return: pandas object
    """
    exp_instr = pd.read_csv(input_txt)
    return exp_instr


# test cases
# print read_experiment_file("OCCM-pipeline_settings.txt")
# exp_instr = read_experiment_file("OCCM-pipeline_settings.txt")
# print exp_instr.columns.values[0]


def fun_call_pipeline_for_each_list_element(
        experiment,
        xi_xifdr_settings_dict,
        amount_of_repeats,  # general settings
        output_basedir="."  # optional general settings
):
    """
    call pipeline once for each line (except for the first) of the input file
    call pipeline with number of repeats

    :param input_txt: txt file with first line containing descriptor of DBs, separated by commas
    each consecutive line contains
    :return:
    """
    # read the experiment describers from the input.txt, exp_instr is a pandas object
    exp_instr = experiment.df_rnd_fasta_number_pairs
    # read out column names for later construction of experiment folders
    col_name_T = exp_instr.columns.values[0]
    col_name_D = exp_instr.columns.values[1]

    # call pipeline for each row=experiment_index of input file
    for experiment_index, row in exp_instr.iterrows():
        random_fasta_db_size = row.iloc[0]
        random_decoy_fasta_db_size = row.iloc[1]
        override_repeats = row.iloc[2]
        # check whether amount_of_repeats shall be overridden for this experiment
        if override_repeats == 0:
            actual_repeats = amount_of_repeats
        else:
            actual_repeats = override_repeats
        # construct output basedir consisting of:
        # basedir / column name 1 + db size 1 + column name 2 + db size 2
        str_run = str(col_name_T) + "_" + str(random_fasta_db_size) + "_" + \
                  str(col_name_D) + "_" + str(random_decoy_fasta_db_size)
        exp_folder = os.path.join(output_basedir, str_run)
        if not os.path.exists(exp_folder):
            os.makedirs(exp_folder)

        # calling of the pipeline
        execute_pipeline(
            experiment=experiment,
            str_run=str_run,
            xi_xifdr_settings_dict=xi_xifdr_settings_dict,
            amount_of_repeats=actual_repeats,  # general settings
            random_fasta_db_size=random_fasta_db_size,
            random_decoy_fasta_db_size=random_decoy_fasta_db_size,  # random fasta settings
            output_basedir=exp_folder,  # optional general settings
        )


if __name__ == "__main__":

    error_free_execution = True

    # print help message if script is called without argument
    if len(sys.argv) != 2:
        print \
"""
Script has to be called with config file as argument.
The directory of the config file will be the output dir.
"""
        sys.exit(1)

    rel_config_file = sys.argv[1]

    # # testing
    # print "testing func active"
    # rel_config_file = r"testing/myconfig.py"

    config_file = os.path.abspath(rel_config_file)
    execfile(config_file)
    # set output dir
    out_dir = os.path.split(config_file)[0]

    log_file = os.path.join(out_dir, __name__ + '.log')

    logging.basicConfig(filename=log_file, level=logging.DEBUG,
                        format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')

    repeats = repeats

    list_of_experiment_dicts = list_of_experiments

    xi_xifdr_settings_dict = xi_xifdr_settings_dict

    list_of_experiments = []
    for experiment_dict in list_of_experiment_dicts:
        list_of_experiments.append(
            Experiment(
                list_of_peak_files=experiment_dict['list_of_peak_files'],
                exp_name=experiment_dict['exp_name'],
                fasta_target=experiment_dict['fasta_target'],
                fasta_decoy=experiment_dict['fasta_decoy'],
                df_rnd_fasta_number_pairs=experiment_dict['df_rnd_fasta_number_pairs']
            )
        )

    for exp in list_of_experiments:
        fun_call_pipeline_for_each_list_element(
            experiment=exp,
            xi_xifdr_settings_dict=xi_xifdr_settings_dict,
            amount_of_repeats=repeats,
            output_basedir=out_dir
        )
    if not error_free_execution:
        print "Some experiments failed with errors. Please check log for further details."
    logging.shutdown()


"""running of the entire script"""
# xifdr_filename = "xiFDRDB-1.0.14.34-jar-with-dependencies.jar"
#
# fun_call_pipeline_for_each_list_element(
#         input_txt="OCCM-pipeline_settings_larger_decoyDB.txt",
#         amount_of_repeats=3,  # general settings
#         fasta_db=r"Data/fastas/OCCM_Scerevisiae-11_17_09-09_Feb_2016/OCCM_Scerevisiae.fasta",
#         decoy_fasta_db=r"Data/fastas/Decoy_170316_UPID_UP000000625/ecoli-proteome.bin",        # random fasta settings
#         xi_config="Data/xi_config.cfg",
#         peak_files=[r"Data/peak_files/B161022_OCCM_BS3_Tryp_SECFr16-15_14_05-31_Oct_2016/B161022_05_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr16.HCD.FTMS.peak.apl",
#                     r"Data/peak_files/B161022_OCCM_BS3_Tryp_SECFr16-15_14_05-31_Oct_2016/B161022_05_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr16.HCD.FTMS.sil0.apl",
#                     r"Data/peak_files/B161022_OCCM_BS3_Tryp_SECFr17-15_15_18-31_Oct_2016/B161022_06_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr17.HCD.FTMS.peak.apl",
#                     r"Data/peak_files/B161022_OCCM_BS3_Tryp_SECFr17-15_15_18-31_Oct_2016/B161022_06_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr17.HCD.FTMS.sil0.apl",
#                     r"Data/peak_files/B161022_OCCM_BS3_Tryp_SECFr18-15_16_42-31_Oct_2016/B161022_07_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr18.HCD.FTMS.peak.apl",
#                     r"Data/peak_files/B161022_OCCM_BS3_Tryp_SECFr18-15_16_42-31_Oct_2016/B161022_07_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr18.HCD.FTMS.sil0.apl"],  # xi settings
#         out_dir="results_larger_decoyDB",  # optional general settings
#         pepfdr="5", memory="1G", reportfactor="10000",
#         additional_xi_parameters=["--xiconf=TOPMATCHESONLY:true"],  # optional xi settings
#         additional_xifdr_arguments=list())  # optional xifdr settings


# # Test Cases
# generate_random_fasta_files(input_file="Data/OCCM_Scerevisiae-11_17_09-09_Feb_2016/OCCM_Scerevisiae.fasta",
#                             output_file="trial01/testset.fasta",
#                             size=9,
#                             number_of_random_sample_files=5)


# todo
# Annahme der Argumente ueber die commandline
# das muss dringend uebersichtlicher gestaltet werden...
#   create object class "experiment" that has all the necessary values stored
