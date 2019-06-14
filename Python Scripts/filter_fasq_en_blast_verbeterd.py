from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from datetime import datetime
import time
from urllib import error
from collections import deque


def main():
    """ Calls all other functions."""
    reads_list = get_inhoud()
    seqs_dict = filter_score(reads_list)
    seqs_left_dict = max_amount(seqs_dict)
    blast_file(seqs_left_dict)
    blast_sequences()


def get_inhoud():
    """ Opens a csv file containing Fastq sequences and returns them."""
    file_name = 'data_groep8tsv.csv'
    file = open(file_name, 'r')
    reads_list = []

    for line in file:
        line = line.split('\t')
        reads_list.append(line)

    return reads_list


def filter_score(reads_list):
    """ Filters fastq sequences based on their phred score."""
    max_score_read_1 = 0                # minimal score per nucleotide read 1
    max_score_read_2 = 0                # minimal score per nucleotide read 2
    counter = 0
    seq_list = []
    seqs_dict = {}

    # Loops over every sequence (read 1 + read 2 each loop).
    for read in reads_list:
        for read_nr in 1, 2:
            seq_list.append(">Seq " + str(counter + 1) + " Read " + str(read_nr))
            # Determines if the max score needed is for read 1 or read 2.
            if read_nr == 1:
                max_score = max_score_read_1
                avg_score = determine_fastq_score(read[2])
            else:
                max_score = max_score_read_2
                avg_score = determine_fastq_score(read[5])

            if avg_score >= max_score:
                seq_list.append(reads_list[counter][read_nr * read_nr])
                seq_list.append(avg_score)
            else:
                seq_list.append("-")
                seq_list.append(0)
            # The information about each sequence is stored.
            seqs_dict[counter + 1 + (read_nr / 10)] = seq_list
            seq_list = []

        counter += 1

    print("minimal score per nucleotide read 1: " + str(max_score_read_1))
    print("minimal score per nucleotide read 2: " + str(max_score_read_2))
    print("-"*70)

    return seqs_dict


def determine_fastq_score(fastq_score):
    """ Determines the average phred score per nucleotide of a sequence."""
    dic_ascii = {
        '!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6,
        '(': 7, ')': 8, '*': 9, '+': 10, ',': 11, '-': 12, '.': 13,
        '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
        '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27,
        '=': 28, '>': 29, '?': 30, '@': 31, 'A': 32, 'B': 33, 'C': 34,
        'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40, 'J': 41,
        'K': 42, 'L': 43, 'M': 44
    }
    avg_score = 0

    for score_character in fastq_score.strip():
        avg_score += dic_ascii[score_character]

    # each fastq score contains 301 characters
    avg_score = round(avg_score/301, 3)

    return avg_score


def max_amount(seqs_dict):
    """ Filters out sequences based on the highest phred scores found."""
    seqs_left_dict = {}
    highest_scores = []
    # max amount of sequences to BLAST
    max_seqs = 120
    counter = 0

    for value in seqs_dict.values():
        # If the seqs left contains less than the indicated max amount of sequences it just adds them.
        if len(highest_scores) < max_seqs:
            highest_scores.append(value[2])
        # It looks through all sequences to determine if this sequence has a higher score then the lowest
        # in the list contains the sequences with the highest scores.
        else:
            lowest_score = min(highest_scores)
            if value[2] > lowest_score:
                # Replaces the new score with the lowest score in the list containing the highest scores
                for score in highest_scores:
                    if score == lowest_score:
                        highest_scores[counter] = value[2]
                        break
                    counter += 1
            counter = 0

    print("max amount of sequences: " + str(max_seqs))

    # Makes a dictionary containing the sequences with the highest scores
    for key, value in seqs_dict.items():
        if value[2] != 0:
            if value[2] in highest_scores:
                seqs_left_dict[key] = [value[0], value[1]]

    return seqs_left_dict


def blast_file(seqs_left_dict):
    """ Puts sequences in a fasta file."""
    total = 0
    file_to_blast = open('seqs_to_blast.fasta', 'w')
    for value in seqs_left_dict.values():
        # Every sequence where the sequence isn't like '-'. This was given in the function filter_score
        if value[1] != "-":
            total += 1
            file_to_blast.write(value[0] + '\n')
            file_to_blast.write(value[1] + '\n')

    print("Total amount of sequences left: " + str(total))


def blast_sequences():
    """ Calls on all functions used for running a BLAST."""
    file_to_blast = open('seqs_to_blast.fasta', "r")
    file_to_blast = iter(file_to_blast)
    save_file = open("BLAST_results_goed.txt", "a")
    header = ''
    sequence = ''
    header_list = deque([])
    sequence_list = deque([])

    no_seq_left = False
    seq_nr = 1
    print("-"*70)
    time_start = datetime.now()

    # Loops as long as there are sequences left
    while no_seq_left is False:
        blasting_seqs, no_seq_left, file_to_blast, sequences_list, header_list = sequences_to_blast(file_to_blast,
                                                                                                    header_list,
                                                                                                    sequence_list,
                                                                                                    sequence,
                                                                                                    header,
                                                                                                    no_seq_left)
        if blasting_seqs != "":
            print("Starting next BLASTs:")
            for headers in header_list:
                print(headers)
            blast_results = blast(blasting_seqs, time_start)
            save_blast_results(blast_results, header_list, sequence_list, save_file, seq_nr)

            print("The BLASTs are done!")
            if no_seq_left is False:
                print("+" * 70)

    time_stop = datetime.now()
    save_file.close()
    print("-" * 70)
    print("All BLASTs are done")
    print("BLAST starting time: ", str(time_start))
    print("BLAST done time: ", str(time_stop))
    print("Total BLAST time: ", str(time_stop - time_start))


def sequences_to_blast(file_to_blast, header_list, sequences_list, sequence, header, no_seq_left):
    """ Takes the first sequences from a fasta file."""
    blasting_seqs = ""
    seqs_per_blast = 40             # number of sequences per BLAST
    for _ in range(seqs_per_blast):
        if sequence != "last sequence" or header != "last sequence":
            # Adds the first indicated amount of sequences used for the BLAST to a list.
            header = next(file_to_blast, "last sequence")
            sequence = next(file_to_blast, "last sequence")
            header_list.append(header.strip())
            sequences_list.append(sequence.strip())
            blasting_seqs = blasting_seqs + header + sequence
        else:
            break
    blasting_seqs = blasting_seqs.replace("last sequence", "")
    if header == 'last sequence' or sequence == 'last sequence':
        no_seq_left = True

    return blasting_seqs, no_seq_left, file_to_blast, sequences_list, header_list


def blast(blasting_seqs, time_start):
    """ Runs a BLAST."""
    blast_results_handle = ''
    try:
        blast_results_handle = NCBIWWW.qblast("blastx", "nr", blasting_seqs, matrix_name="BLOSUM62", hitlist_size=75,
                                              expect=10, gapcosts="11 2", word_size=6, filter=True)

    except error.URLError:
        time_stop = datetime.now()
        print("-" * 70)
        print("BLAST failed")
        print("Reason: No internet connection")
        print("BLAST starting time: ", str(time_start))
        print("Failed at: " + str(datetime.now()))
        print("Total BLAST time: ", str(time_stop - time_start))
        exit()
    time.sleep(10)
    blast_results = NCBIXML.parse(blast_results_handle)

    return blast_results


def save_blast_results(blast_results, header_list, sequence_list, save_file, seq_nr):
    """ Saves all results from a BLAST."""
    for blast_result in blast_results:
        header_result = header_list.popleft()
        sequence_result = sequence_list.popleft()
        # Saving 'input' data.
        save_file.write("~" + header_result + "\t" + sequence_result + "\t" + str(seq_nr) + '\n')
        hit_nr = 1

        result_id = header_result.split(" ")
        for alignment in blast_result.alignments:
            # Saving result id
            save_file.write("$" + result_id[1] + "." + result_id[3] + "*" + str(hit_nr) + '\n')
            save_file.write(blast_result.descriptions[hit_nr - 1].title + '\n')
            save_file.write(str(alignment.hsps[0].score) + '\n')
            save_file.write(str(alignment.hsps[0].expect) + '\n')
            # Saving percent identity
            save_file.write(str(round((alignment.hsps[0].identities * 100) /
                                      (alignment.hsps[0].align_length - alignment.hsps[0].gaps), 2)) + '\n')
            # Saving query coverage
            save_file.write(str(round(((alignment.hsps[0].query_end - alignment.hsps[0].query_start)
                                       / alignment.hsps[0].query_end), 3)) + '\n')
            save_file.write(alignment.accession + '\n')
            hit_nr += 1
        seq_nr += 1


main()
