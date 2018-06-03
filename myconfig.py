import pandas as pd

repeats = 3

# how to construct the fastas. override repeat number? if yes, with what?
labels = ['occm-proteins', 'e-coli-proteins', 'override_repeats']
data = [
    (1, 0, 0),
    (2, 0, 0),
    (3, 0, 0),
    (7, 0, 0),
    (13, 0, 0),
    (14, 0, 1),
    (14, 1, 0),
    (14, 2, 0),
    (14, 3, 0),
    (14, 7, 0),
    (14, 10, 0),
    (14, 14, 0),
    (14, 28, 0)
]
df_rnd_fasta_number_pairs = pd.DataFrame.from_records(data, columns=labels)


# initiate object for each experiment
list_of_experiments = []

fasta_target = r'../../Data/Input/170322_OCCM_random_DB_test/fastas/OCCM_Scerevisiae-11_17_09-09_Feb_2016/OCCM_Scerevisiae.fasta'
fasta_decoy = r'../../Data/Input/170322_OCCM_random_DB_test/fastas/Decoy_170316_UPID_UP000000625/ecoli-proteome.bin'

exp = {
    'exp_name': "OCCM 1710",
    'list_of_peak_files': [
        r"../../Data/Input/170322_OCCM_random_DB_test/peak_files/B161022_OCCM_BS3_Tryp_SECFr16-15_14_05-31_Oct_2016/B161022_05_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr16.HCD.FTMS.peak.apl",
        r"../../Data/Input/170322_OCCM_random_DB_test/peak_files/B161022_OCCM_BS3_Tryp_SECFr16-15_14_05-31_Oct_2016/B161022_05_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr16.HCD.FTMS.sil0.apl",
        r"../../Data/Input/170322_OCCM_random_DB_test/peak_files/B161022_OCCM_BS3_Tryp_SECFr17-15_15_18-31_Oct_2016/B161022_06_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr17.HCD.FTMS.peak.apl",
        r"../../Data/Input/170322_OCCM_random_DB_test/peak_files/B161022_OCCM_BS3_Tryp_SECFr17-15_15_18-31_Oct_2016/B161022_06_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr17.HCD.FTMS.sil0.apl",
        r"../../Data/Input/170322_OCCM_random_DB_test/peak_files/B161022_OCCM_BS3_Tryp_SECFr18-15_16_42-31_Oct_2016/B161022_07_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr18.HCD.FTMS.peak.apl",
        r"../../Data/Input/170322_OCCM_random_DB_test/peak_files/B161022_OCCM_BS3_Tryp_SECFr18-15_16_42-31_Oct_2016/B161022_07_Lumos_ML_IN_205_OCCM_BS3_Tryp_SECFr18.HCD.FTMS.sil0.apl"
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
    'additional_xi_parameters': [
        "--xiconf=TOPMATCHESONLY:true",
        "--xiconf=tolerance:fragment:10ppm",
        "--xiconf=missedcleavages:2"
    ],
    'xifdr_settings': {
        'pepfdr': "10",
        'additional_xifdr_arguments': [],
        'reportfactor': "10000"
    }
}
