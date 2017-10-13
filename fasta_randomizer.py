"""
author: Henning Schiebenhoefer

"""

import random
import re
import argparse
# import os

# regular expression to catch protein start line
# example line that needs to be caught by regular expression:
# >sp|Q47274|REQ1_ECOLI Prophage antitermination protein Q homolog QuuD OS=Escherichia coli (strain K12)
# GN=quuD PE=3 SV=1
regex_protein_start = re.compile(r'^(>.*)')


def get_protein_lines(input_file):
    """Returns list of line numbers, which contain protein starts"""
    candidate_list = []
    with open(input_file) as input_fasta:
        for line_i, line in enumerate(input_fasta):
            if regex_protein_start.search(line):
                candidate_list.append(line_i)
    return candidate_list


def generate_output_file(input_file, output_file, chosen_candidates):
    """
    writes proteins belonging to chosen candidate lines in input file to output file.
    :param input_file:
    :param output_file:
    :param chosen_candidates:
    :return:
    """
    with open(output_file, 'w+') as o:
        with open(input_file) as i:
            lines = i.readlines()
            for line in chosen_candidates:
                o.write(lines[line])
                line_belonging_to_this_prot = line
                # write next line, but only:
                # if it does not contain beginning of next protein and
                # if line is not out of line range of the input file
                while True:
                    line_belonging_to_this_prot += 1
                    # break, if line is out of range
                    if line_belonging_to_this_prot == len(lines):
                        break
                    # break if line is beginning of next protein
                    elif regex_protein_start.match(lines[line_belonging_to_this_prot]):
                        break
                    # write another line of chosen candidate protein
                    else:
                        o.write(lines[line_belonging_to_this_prot])
    return


def generate_random_protein_list(input_file, output_file="chosen_random_candidates.fasta", desired_no=1000):
    list_of_candidate_line_numbers = get_protein_lines(input_file)
    if len(list_of_candidate_line_numbers) < desired_no:
        raise ValueError('Candidate list does not contain enough proteins! ' +
                         '\nNumber of proteins in input file: \t' + str(len(list_of_candidate_line_numbers)) +
                         '\nDesired number of output proteins: \t' + str(desired_no))
    # for testing: set deterministic behaviour
    # random.seed(1234)
    chosen_candidates = random.sample(list_of_candidate_line_numbers, desired_no)
    # make sure that identical chosen candidates yield identical output (same order of proteins)
    chosen_candidates.sort()
    generate_output_file(input_file, output_file, chosen_candidates)


# Command line control arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='script to generate random sample of proteins from input fasta file')
    parser.add_argument("input", help="fasta file from which the random sample is generated")
    parser.add_argument("-o", "--output", help="output file. [chosen_random_candidates.fasta]",
                        default="chosen_random_candidates.fasta")
    parser.add_argument("-s", "--size", help="Number of proteins to be in output file. [1000]",
                        type=int, default=1000)
    args = parser.parse_args()

    generate_random_protein_list(args.input, args.output, args.size)
    print("random sample written to {}".format(args.output))

# # Test cases
# generate_random_protein_list(r"test-set.fasta", desired_no=5)
# generate_random_protein_list(r"Data/OCCM_Scerevisiae-11_17_09-09_Feb_2016/OCCM_Scerevisiae.fasta",
#                                output_file="output_uniprot_testset.fasta", desired_no=9)
# # mit parser:
# parser.print_help()
# parser.parse_args(['uniprot-proteome%25253AUP000000625.bin', '-o', 'output_uniprot_testset.fasta', '-s', '9'])

# todo
