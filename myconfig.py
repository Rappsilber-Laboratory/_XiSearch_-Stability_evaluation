import pandas as pd

repeats = 3

# how to construct the fastas. override repeat number? if yes, with what?
labels = ['occm-proteins', 'e-coli-proteins', 'override_repeats']
data = [
    #(1, 0, 10),
    #(2, 0, 10),
    (3, 0, 10)
]
df_rnd_fasta_number_pairs = pd.DataFrame.from_records(data, columns=labels)


# initiate object for each experiment
list_of_experiments = []

fasta_target = r'../../Data/Input/170322_OCCM_random_DB_test/fastas/OCCM_Scerevisiae-11_17_09-09_Feb_2016/OCCM_Scerevisiae.fasta'
fasta_decoy = r'../../Data/Input/170322_OCCM_random_DB_test/fastas/Decoy_170316_UPID_UP000000625/ecoli-proteome.bin'

exp = {
    'exp_name': "OCCM 1710",
    'list_of_peak_files': [
        r"../../Data/Input/170322_OCCM_random_DB_test/peak_files/180403-preprocessed/recalibrated_files.zip"
    ],
    'lst_rnd_fasta_number_pairs': df_rnd_fasta_number_pairs,
    'fasta_target': fasta_target,
    'fasta_decoy': fasta_decoy,
    'df_rnd_fasta_number_pairs': df_rnd_fasta_number_pairs
}
list_of_experiments += [exp]


xi_xifdr_settings_dict = {
    # same config file for all searches
    'xi_config': r"../../Data/Input/170322_OCCM_random_DB_test/xi_config.cfg",
    'xi_memory': '100G',
    'xi_path': r"/data/rappstore/users/hschiebenhoefer/scripts/xiSearch-versions/XiSearch_1.6.731.jar",
    'additional_xi_parameters': [
        "--xiconf=TOPMATCHESONLY:true",
        "--xiconf=tolerance:fragment:10ppm",
        "--xiconf=missedcleavages:2",
        "--xiconf=missing_isotope_peaks:2"
    ],
    'xifdr_settings': {
        'pepfdr': "10",
        'additional_xifdr_arguments': [],
        'reportfactor': "10000",
        'xifdr_path': r'/data/rappstore/users/hschiebenhoefer/scripts/xiFDR-versions/xiFDRDB-1.1.25.55-jar-with-dependencies.jar'
    }
}
