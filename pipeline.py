"""
author: 'Henning Schiebenhoefer'

Input:
Peak Lists
fasta files
Xi cfg
Output:
XiFDR results
"""

import os
import fasta_randomizer
import subprocess
import pandas as pd
import sys  # Anzahl der Parameter beim Aufruf
import argparse     # parsen von Befehlen aus der Kommandozeile
import logging


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
    input_file_basename = os.path.basename(input_file)
    output_file_basename = "random_fasta_" + input_file_basename
    output_file = os.path.join(output_dir, output_file_basename)
    # call fasta file generator
    fasta_randomizer.generate_random_protein_list(input_file, output_file, desired_no)
    generated_fasta_file = output_file
    return generated_fasta_file
    # todo

# # Test cases
# print generate_random_fasta_files(r"../fasta_randomizer/uniprot-proteome%25253AUP000000625.bin", "random_db", 100, 10)
# print os.path.splitext("OCCM_Scerevisiae.fasta")[0]
# print generate_random_fasta_files("Data/OCCM_Scerevisiae-11_17_09-09_Feb_2016/OCCM_Scerevisiae.fasta", "testing", 5)


def build_xi_arguments(xi_config, peak_files, fasta_files, output, additional_parameters=()):
    cmd = []
    cmd.extend(["java", "-cp", "XiSearch.jar", "rappsilber.applications.Xi"])
    # config
    cmd.append("--config=" + xi_config)
    # additional parameters
    for par in additional_parameters:
        cmd.append(par)
    # peak_files
    for peak_file in peak_files:
        cmd.append("--peaks=" + peak_file)
    # fasta_files
    for fasta_file in fasta_files:
        cmd.append("--fasta=" + fasta_file)
    # output
    cmd.append("--output=" + output)
    return cmd

# print build_xi_arguments("das_ist_die_cfg", ["erste_peaks", "zweite_peaks"], ["fasta1", "fasta2"], "hier_gehts_hin")
# print build_xi_arguments("das_ist_die_cfg", ["zweite_peaks"], ["fasta1", "fasta2"], "hier_gehts_hin",
#                          additional_parameters=["--xiconfig=TOPMATCHESONLY:true"])


def xi_execution(xi_config, peak_files, fasta_files, output_folder, additional_parameters=list()):
    """
    Calls Xi and gives back Xi result
    Xi gives back a single csv file

    Example cmd line:
    # java -cp XiSearch.jar rappsilber.applications.Xi --conf=[path to config]
    # --xiconfig=TOPMATCHESONLY:true --peaks=[path to peaklist1] --peaks=[path to peaklist2] --peaks=[path to peaklist3]
    # --fasta=[path to fasta file1] --fasta=[path to fasta file2] --fasta=[path to fasta file3]
    # --output=[path to result file]

    :param xi_config: string
    :param peak_files: list of strings
    :param fasta_files: list of strings
    :param output_folder: string
    :param additional_parameters: list of strings
    :return: string: output file name
    """
    assert type(peak_files) == type(fasta_files) == type(additional_parameters) == type([]), \
        """type of the following files needs to be list. It is actually:
        peak_files: {}
        fasta_files: {}
        additional_parameters: {}""".format(type(peak_files), type(fasta_files), type(additional_parameters))
    # generate output_file name
    output_file = os.path.join(output_folder, "xi_results.csv")  # filename example: "xi_results.csv"
    # generate xi commands
    xi_cmd = build_xi_arguments(xi_config, peak_files, fasta_files, output_file, additional_parameters)
    # call xi
    try:
        print subprocess.check_output(xi_cmd)
    except subprocess.CalledProcessError:
        raise AttributeError('XiSearch exited with error message!')
    return output_file


# # Test Cases
# xi_execution("Data/xi_config.cfg", ["Data/B161022_OCCM_BS3_Tryp_SECFr18-15_16_42-31_Oct_2016/B161022_07_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr18.HCD.FTMS.peak.apl",
#              "Data/B161022_OCCM_BS3_Tryp_SECFr18-15_16_42-31_Oct_2016/B161022_07_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr18.HCD.FTMS.sil0.apl"],
#              ["Data/OCCM_Scerevisiae-11_17_09-09_Feb_2016/OCCM_Scerevisiae.fasta"], "test_Xi_results.csv",
#              additional_parameters=["--xiconf=TOPMATCHESONLY:true"])

# xi_cmd = ["java", "-cp", "XiSearch.jar", "-Xmx1G", "rappsilber.applications.Xi", "--help"]
# print subprocess.check_output(xi_cmd)


def build_xifdr_arguments(fdr_input_csv, fdr_output_dir, pepfdr, memory="1G", reportfactor="10000",
                          additional_xifdr_arguments=list(), xifdr_filename="xiFDRDB-1.0.14.34-jar-with-dependencies.jar"):
    assert type(fdr_input_csv) == type([]), """type of fdr_input_csv needs to be list but is: {}""".format(type(fdr_input_csv))
    # Example cmd line:
    # java -Xmx1g -cp xiFDRDB-1.0.13.32-jar-with-dependencies.jar org.rappsilber.fdr.CSVinFDR --psmfdr=X --pepfdr=X
    # --proteinfdr=X --reportfactor=X --linkfdr=X --ppifdr=X --csvOutDir=X --csvBaseName=X csv-file1 csv-file2
    cmd = []
    # memory
    cmd.extend(["java", "-Xmx" + memory, "-cp", xifdr_filename,
                "org.rappsilber.fdr.CSVinFDR", '--reportfactor=10000'])
    # pepfdr
    cmd.append("--pepfdr=" + pepfdr)
    # additional arguments
    for par in additional_xifdr_arguments:
        cmd.append(par)
    # reportfactor
    # taken into default config as it will probably never be changed
    # cmd.append("--reportfactor=" + reportfactor)
    # csvOutDir
    cmd.append("--csvOutDir=" + fdr_output_dir)
    # input csv
    for i in fdr_input_csv:
        cmd.append(i)
    return cmd


# print build_xifdr_arguments(["Xi_results.csv"], "xifdr_test", "5", "1G", "1000")

def xifdr_execution(xifdr_input_csv, xifdr_output_dir, pepfdr="5", memory="1G", reportfactor="10000",
                    additional_xifdr_arguments=list()):
    """
    takes XiSearch output and gives back fdr csv
    :param xifdr_input_csv:
    :param xifdr_output_dir:
    :param additional_xifdr_arguments:
    :return:
    """
    xifdr_cmd = build_xifdr_arguments(xifdr_input_csv, xifdr_output_dir, pepfdr, memory, reportfactor,
                                      additional_xifdr_arguments)
    try:
        print subprocess.check_output(xifdr_cmd)
    except subprocess.CalledProcessError:
        raise AttributeError('XiFDR exited with error message!')
    return


def fun_generate_list_of_dirs(output_basedir, pre_list_of_dirs, n):
    """ creates list of directory names
    for each element in pre_list_of dirs:
    joins output_basedir / element / n
    to a path name"""
    list_of_dirs = []
    for directory in pre_list_of_dirs:
        list_of_dirs.append(os.path.join(output_basedir, directory, str(n)))
    return list_of_dirs

# # Test Cases
# print fun_generate_list_of_dirs("OCCM=1_Ecoli=5",
#                                 ["random_fasta_dbs", "xi_output", "xifdr_output"],
#                                 3)


def fun_makedirs(list_of_dirs):
    """ checks whether directory exists for each directory in list.
    Creates each that does not exist."""
    for directory in list_of_dirs:
        if not os.path.exists(directory):
            os.makedirs(directory)


def execute_pipeline(amount_of_repeats,  # general settings
                     fasta_db, random_fasta_db_size, decoy_fasta_db, random_decoy_fasta_db_size,  # random fasta settings
                     xi_config, peak_files,  # xi settings
                     output_basedir=".",  # optional general settings
                     pepfdr="5", memory="1G", reportfactor="10000", additional_xi_parameters=list(),  # optional xi settings
                     additional_xifdr_arguments=list()):  # optional xifdr settings
    """
    This pipeline executes:
    fasta_randomizer
    xiSeqrch
    and Xifdr
    for amount_of_repeats times
    :param fasta_db:
    :param random_fasta_db_size:
    :param decoy_fasta_db:
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
    # base dirs for all the files
    pre_list_of_dirs = ["random_fasta_dbs", "xi_output", "xifdr_output"]
    # loop through generation of all files amount_of_repeats times
    # for loop adjusted to make i represent run number
    for i in range(1, amount_of_repeats+1):
        # generate necessary dirs for this run
        list_of_dirs = fun_generate_list_of_dirs(output_basedir, pre_list_of_dirs, i)
        fun_makedirs(list_of_dirs)

        # generate necessary fastas for this run: one call for true_db, one call for random_db
        random_fasta_db = generate_random_fasta_files(input_file=fasta_db,
                                                      output_dir=list_of_dirs[0],
                                                      desired_no=random_fasta_db_size)
        random_decoy_fasta_db = generate_random_fasta_files(input_file=decoy_fasta_db,
                                                            output_dir=list_of_dirs[0],
                                                            desired_no=random_decoy_fasta_db_size)

        # call xisearch
        xi_result = xi_execution(xi_config=xi_config,
                                 peak_files=peak_files,
                                 fasta_files=[random_fasta_db, random_decoy_fasta_db],
                                 output_folder=list_of_dirs[1],
                                 additional_parameters=additional_xi_parameters)

        # xi_result is a string but xifdr needs a list as input
        xifdr_input = []
        xifdr_input.append(xi_result)
        # call xifdr
        xifdr_execution(xifdr_input_csv=xifdr_input,
                        xifdr_output_dir=list_of_dirs[2],
                        pepfdr=pepfdr,
                        memory=memory,
                        reportfactor=reportfactor,
                        additional_xifdr_arguments=additional_xifdr_arguments)
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


def fun_create_dir_for_experiment(basedir, column_name_1, random_fasta_db_size, column_name_2, random_decoy_fasta_db_size):
    """constructs a folder name and creates the dir based on the input
    dirname: 'column_name_1'_'random_fasta_db_size'_'column_name_2'_'random_decoy_fasta_db_size'
    returns created folder name"""
    exp_base_folder = str(column_name_1) + "_" + str(random_fasta_db_size) + "_" + \
                      str(column_name_2) + "_" + str(random_decoy_fasta_db_size)
    exp_folder = os.path.join(basedir, exp_base_folder)
    if not os.path.exists(exp_folder):
        os.makedirs(exp_folder)
    return exp_folder


# # test cases
# print fun_create_dir_for_experiment("testing", "occm", 5, "Ecoli", 7)

def fun_call_pipeline_for_each_list_element(
        input_txt,
        amount_of_repeats,  # general settings
        fasta_db, decoy_fasta_db,        # random fasta settings
        xi_config, peak_files,  # xi settings
        output_basedir=".",  # optional general settings
        pepfdr="5", memory="1G", reportfactor="10000",
        additional_xi_parameters=list(),  # optional xi settings
        additional_xifdr_arguments=list()):  # optional xifdr settings
    """
    call pipeline once for each line (except for the first) of the input file
    call pipeline with number of repeats

    :param input_txt: txt file with first line containing descriptor of DBs, separated by commas
    each consecutive line contains
    :return:
    """
    # read the experiment describers from the input.txt, exp_instr is a pandas object
    exp_instr = read_experiment_file(input_txt)
    # read out column names for later construction of experiment folders
    column_name_1 = exp_instr.columns.values[0]
    column_name_2 = exp_instr.columns.values[1]

    # call pipeline for each row=experiment_index of input file
    for experiment_index, row in exp_instr.iterrows():
        random_fasta_db_size = exp_instr.iloc[experiment_index, 0]
        random_decoy_fasta_db_size = exp_instr.iloc[experiment_index, 1]
        override_repeats = exp_instr.iloc[experiment_index, 2]
        # check whether amount_of_repeats shall be overridden for this experiment
        if override_repeats == 0:
            actual_repeats = amount_of_repeats
        else:
            actual_repeats = override_repeats
        # construct output basedir consisting of:
        # basedir / column name 1 + db size 1 + column name 2 + db size 2
        exp_folder = fun_create_dir_for_experiment(basedir=output_basedir,
                                                   column_name_1=column_name_1,
                                                   random_fasta_db_size=random_fasta_db_size,
                                                   column_name_2=column_name_2,
                                                   random_decoy_fasta_db_size=random_decoy_fasta_db_size)
        # calling of the pipeline
        execute_pipeline(amount_of_repeats=actual_repeats,  # general settings
                         fasta_db=fasta_db, random_fasta_db_size=random_fasta_db_size,
                         decoy_fasta_db=decoy_fasta_db, random_decoy_fasta_db_size=random_decoy_fasta_db_size,
                         # random fasta settings
                         xi_config=xi_config, peak_files=peak_files,  # xi settings
                         output_basedir=exp_folder,  # optional general settings
                         pepfdr=pepfdr, memory=memory, reportfactor=reportfactor, additional_xi_parameters=additional_xi_parameters,
                         # optional xi settings
                         additional_xifdr_arguments=additional_xifdr_arguments)  # optional xifdr settings
        # todo


class Runner(object):
    """
    class that initializes the parser and executes the search
    """
    def __init__(self):
        self._args = self.init_parser()

    def print_the_help(self):
        parser.print_help()

    def init_parser(self):
        parser = argparse.ArgumentParser(
            description="""Script to execute Xi Searches and consecutive XiFDR analyses on randomly generated fasta
            files.""")
        parser.add_argument(peak_files, type=list,
                            help="peak files from MaxQuant")
        parser.add_argument(input_cfg,
                            help="""csv file that contains the sizes of the random DBs to generate.
            An additional row states whether standard repeats should be overridden""")
        parser.add_argument(fasta_db,
                            help="fasta protein database to search on")
        parser.add_argument(-r, --repeats, default=1,
                            help="how often should the search be repeated with newly generated databases [default=1]")
        parser.add_argument(-v, --verbose, default=1, type=int, help='set level of logging')
        xisearch_argument_group = parser.add_argument_group('xiSearch', 'parameters for the Xi crosslinking search')
        xisearch_argument_group.add_argument(xi_config,
                                             help="config file for the Xi search")
        xisearch_argument_group.add_argument(--xisearchmem, type=string, default=None,
                                             help="""how much memory to allocate to xiSearch. [default=let java decide]""")
        xisearch_argument_group.add_argument(--add_xisearch_cmd, type=list, nargs='*',
                                             help="""additional parameters to hand to xi [default=--xiconf=TOPMATCHESONLY:true]""")
        xifdr_argument_group = parser.add_argument_group('xiFDR', 'parameters for the analysis with xiFDR')
        xifdr_argument_group.add_argument(--pepfdr, default="5", type=float,
                                          help="""FDR in percent in terms of peptide links [default: 5]""")
        xifdr_argument_group.add_argument(--xifdrmem, default="1G",
                                          help="""how much memory to allocate to xiFDR. [default=1G]""")
        random_fasta_argument_group = parser.add_argument_group('random fasta generator',
                                                                """parameters for the execution of fasta DB randomization""")
        random_fasta_argument_group.add_argument(--decoydb, type=str,
                                                 help="decoy fasta DB to generate false targets from")


        return parser.parse_args()

    def run(self):

if __name__ == "__main__":
    runner = Runner()
    if len(sys.argv) <= 1:  # spuckt die hilfe aus, wenn ein argument oder weniger mitgegeben wurde
        runner.print_the_help()
        sys.exit(1)

    runner.run()


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
#         output_basedir="results_larger_decoyDB",  # optional general settings
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
