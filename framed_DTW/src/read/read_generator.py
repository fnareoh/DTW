from experiment import compute_matrices, print_details, reverse, complement
from sys import argv
from random import randint, random, choice
import csv


class Error_rate:
    def __init__(self, ed, hom):
        self.inser_del_sub = ed
        self.homopoly = hom


def load_genome():
    genome_file = argv[2]
    file_G = open(genome_file, "r")
    G = "".join(["".join(l.split()) for l in file_G.readlines()[1:]])
    file_G.close()
    return G


def evaluate_dtw_ed(R, G, pos, qual):
    ldtw, ldtw_rev_comp, led, led_rev_comp = compute_matrices(R, G)
    score_dtw, pos_dtw = ldtw.min_last_row_val_index()
    start_dtw, _, _ = ldtw.compute_origin_min_position()
    score_dtw_rev_comp, pos_dtw_rev_comp = ldtw_rev_comp.min_last_row_val_index()
    if score_dtw_rev_comp < score_dtw:
        score_dtw, pos_dtw = score_dtw_rev_comp, pos_dtw_rev_comp
        start_dtw, _, _ = ldtw_rev_comp.compute_origin_min_position()

    score_ed, pos_ed = led.min_last_row_val_index()
    start_ed, _, _ = led.compute_origin_min_position()
    score_ed_rev_comp, pos_ed_rev_comp = led_rev_comp.min_last_row_val_index()
    if score_ed_rev_comp < score_ed:
        score_ed, pos_ed = score_ed_rev_comp, pos_ed_rev_comp
        start_ed, _, _ = led_rev_comp.compute_origin_min_position()

    if abs(pos - start_dtw) > 10 or abs(pos - start_ed) > 10:
        print("************* Read **************")
        print(f"starting position: {pos}")
        print(R)
        print(qual)
        print_details(R, G, ldtw, ldtw_rev_comp, led, led_rev_comp)

    return score_dtw, score_ed


def add_err(S, err):
    nucleotides = ["A", "C", "G", "T"]
    list_read = []
    list_qual = []
    last_c = ""
    nb_IDS = 0
    nb_homopoly = 0
    for c in S:
        if random() < err.inser_del_sub:
            err_type = choice(["substitution", "insertion", "deletion"])
            nb_IDS += 1
            if err_type == "substitution":
                c = choice([x for x in nucleotides if x != c])
                list_qual.append("S")
            if err_type == "insertion":
                tmp_c = choice(nucleotides)
                list_read.append(tmp_c)
                list_qual.append("I")
            if err_type == "deletion":
                c = ""
                list_qual.append("D")
        else:
            list_qual.append("_")
        list_read.append(c)
        if c == last_c:
            if random() < err.homopoly:  # homopoly extension
                list_read.append(c)
                list_qual.append("H")
                nb_homopoly += 1
        last_c = c
    return "".join(list_read), "".join(list_qual), nb_IDS, nb_homopoly


def generate_read(G, N, err):
    read_length = 500
    prob_rev_complement = 0.5
    genome_name = argv[2].split("/")[-1].split(".")[0]
    output_file_name = (
        f"results/{genome_name}_N_{N}_IDS_{err.inser_del_sub}_H_{err.homopoly}"
    )
    read_file = open(output_file_name + ".fastq", "w")
    output_file = open(output_file_name + ".csv", "w")
    csv_output = csv.writer(output_file)
    csv_output.writerow(
        [
            "file",
            output_file_name,
            "Nb subseq",
            N,
            "proba insertion Deletion Substitution",
            err.inser_del_sub,
            "proba Homopolymer extension",
            err.homopoly,
        ]
    )
    csv_output.writerow(["read_id", "nb_IDS", "nb_H", "min DTW", "min ED"])
    for i in range(N):
        r = randint(0, len(G) - read_length - 1)
        subseq = G[r : r + read_length]
        if random() < prob_rev_complement:
            subseq = reverse(complement(subseq))
        read, qual, nb_ids, nb_homopoly = add_err(subseq, err)
        id = (
            f"read_nb_{i}_pos_{r}_length_{read_length}_"
            + f"homopoly_err_{err.homopoly}_err_inser_del_sub_{err.inser_del_sub}_"
            + f"nb_H_{nb_homopoly}_nb_IDS_{nb_ids}"
        )
        read_file.write(f"@{id}\n{read}\n+{id}\n{qual}\n")
        print(id)
        score_dtw, score_ed = evaluate_dtw_ed(read, G, r, qual)
        csv_output.writerow([id, nb_ids, nb_homopoly, score_dtw, score_ed])
    read_file.close()
    output_file.close()


def main():
    err_inser_del_sub_list = [round(0.06 + 0.02 * i, 2) for i in range(6)]
    err_homopoly_list = [round(0.05 + 0.05 * i, 2) for i in range(6)]
    N = int(argv[1])
    G = load_genome()
    for err_ed in err_inser_del_sub_list:
        for err_hom in err_homopoly_list:
            err = Error_rate(err_ed, err_hom)
            generate_read(G, N, err)


if __name__ == "__main__":
    main()
