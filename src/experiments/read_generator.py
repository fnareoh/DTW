from dynamic_prog.local_alignement import LocalED, LocalDTW
from sys import argv
import argparse
from random import randint, random, choice, choices
import csv


class Error_rate:
    def __init__(self, hom):
        self.snp = 1  # percentage of substitution
        self.indel = 0.05  # percentage of indel
        self.max_len_ID = 10
        self.seq_S = 0.001
        self.homopoly = hom

    def IDS_distrib(self):
        return f"I:{round(self.inser,2)},D:{round(self.dele,2)},S:{round(self.sub,2)}"


def evaluate_dtw_ed(R, G, pos, qual, biological_var, bio_qual):
    """ Alligns the sequence R to G w.r. to ED and DTW """
    ldtw = LocalDTW(R, G)
    led = LocalED(R, G)
    score_dtw, pos_dtw = ldtw.min_last_row_val_index()
    start_dtw, al_P, al_T = ldtw.compute_origin_min_position()

    score_ed, pos_ed = led.min_last_row_val_index()
    start_ed, _, _ = led.compute_origin_min_position()

    if False:  # ugly debuging
        print("************* Read **************")
        print(f"starting position: {pos}")
        print(f"biological_var: {biological_var}")
        print(R)
        print(bio_qual)
        print(qual)
        print("******** Dynamic Time Warp **********")
        print(
            f"Alignement (cost {score_dtw}) of Q ending at {pos_dtw} in T starts at {start_dtw} in T."
        )
        print(al_P)
        print(al_T)
        print("******** Edit distance **********")
        print(
            f"Alignement (cost {score_ed}) of Q ending at {pos_ed} in T starts at {start_ed} in T."
        )

    return score_dtw, score_ed


def add_IDS_err(S, err):
    """ Adds Insertion Deletion and Substitution meant to represent a biological variant """
    nucleotides = ["A", "C", "G", "T"]
    list_S = []
    list_qual = []
    nb_IDS = 0
    i = 0
    while i < len(S):
        c = S[i]
        if random() < err.snp / 100:
            nb_IDS += 1
            c = choice([x for x in nucleotides if x != c])
            list_qual.append("S")
        elif random() < err.indel / 100:
            err_type = choice(["insertion", "deletion"])
            if err_type == "insertion":
                len_i = randint(1, err.max_len_ID)
                nb_IDS += len_i
                tmp_c = choices(nucleotides, k=len_i)
                list_S += tmp_c
                list_qual.append("I" * len_i)
            if err_type == "deletion":
                len_d = randint(1, err.max_len_ID)
                nb_IDS += len_d
                c = ""
                i += len_d - 1
                list_qual.append("D" * len_d)
        else:
            list_qual.append("_")
        list_S.append(c)
        i += 1
    return "".join(list_S), "".join(list_qual), nb_IDS


def add_seq_err(S, err):
    """ Adds Substitution and Homopolymer extension meant to represent sequencing errors """
    nucleotides = ["A", "C", "G", "T"]
    list_read = []
    list_qual = []
    last_c = ""
    nb_S = 0
    nb_homopoly = 0
    for c in S:
        if random() < err.seq_S:
            nb_S += 1
            c = choice([x for x in nucleotides if x != c])
            list_qual.append("S")
        else:
            list_qual.append("_")
        list_read.append(c)
        if c == last_c:
            if random() < err.homopoly:  # homopoly extension
                list_read.append(c)
                list_qual.append("H")
                nb_homopoly += 1
        last_c = c
    return "".join(list_read), "".join(list_qual), nb_S, nb_homopoly


def deletion_before(qual):
    """ Compute for each position i the number of deletion in qual[:i] """
    dele_prec = [0] * len(qual)
    nb_del = 0
    for i, c in enumerate(qual):
        dele_prec[i] = nb_del
        if c == "D":
            nb_del += 1
    return dele_prec


def theoretical_distance(r, l, qual, dele_prec):
    nb_match = qual[
        r + dele_prec[r] : min(r + l + dele_prec[r + l - 1], len(qual))
    ].count("_")
    nb_mismatch = (
        min(r + l + dele_prec[r + l], len(qual)) - (r + dele_prec[r]) - nb_match
    )
    if nb_mismatch > 30:
        print(
            "start:",
            r + dele_prec[r + 1],
            "end:",
            r + l + dele_prec[r + l],
            "len(qual):",
            len(qual),
        )
        print("nb match:", nb_match, "nb nb_mismatch:", nb_mismatch)
    return nb_mismatch


def generate_read(Go, Gi, qual_Gi, N, read_length, err):
    dele_prec = deletion_before(qual_Gi)
    genome_name = argv[1].split("/")[-1].split(".")[0]
    output_file_name = f"results/{genome_name}_N_{N}_ID_{err.indel}_SNP_{err.snp}_H_{err.homopoly}_seqS_{err.seq_S}_max_len_indel_{err.max_len_ID}"
    read_file = open(output_file_name + ".fastq", "w")
    output_file = open(output_file_name + ".csv", "w")
    csv_output = csv.writer(output_file)
    csv_output.writerow(
        [
            "file",
            output_file_name,
            "Nb subseq",
            N,
            "percentage of Indel",
            err.indel,
            "percentage of SNP",
            err.snp,
            "proba Homopolymer extension",
            err.homopoly,
            "proba substitution sequencing",
            err.seq_S,
        ]
    )
    csv_output.writerow(
        ["read_id", "biological_var", "seq_S", "seq_H", "min DTW", "min ED"]
    )
    for i in range(N):
        r = randint(0, len(Gi) - read_length - 1)
        biological_var = theoretical_distance(r, read_length, qual_Gi, dele_prec)
        subseq = Gi[r : r + read_length]
        read, qual, nb_s, nb_homopoly = add_seq_err(subseq, err)
        id = (
            f"read_nb_{i}_pos_{r}_length_{read_length}_"
            + f"homopoly_err_{err.homopoly}_sub_err_{err.seq_S}_"
            + f"Indel_{err.indel}_SNP_{err.snp}"
        )
        bio_qual = qual_Gi[
            r
            + dele_prec[r + 1] : min(
                r + read_length + dele_prec[r + read_length - 1], len(qual_Gi)
            )
        ]
        read_file.write(f"@{id}\n{read}\n+{id}\n{bio_qual}\n")
        print(id)
        # print(r, dele_prec[r])
        score_dtw, score_ed = evaluate_dtw_ed(
            read, Go, r, qual, biological_var, bio_qual
        )
        csv_output.writerow(
            [id, biological_var, nb_s, nb_homopoly, score_dtw, score_ed]
        )
    read_file.close()
    output_file.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("genome_file", help="Genome file", type=str)

    # Optional arguments
    parser.add_argument(
        "-N", "--N", help="Number of read per set of parameters", type=int, default=10
    )
    parser.add_argument(
        "-RL", "--read_length", help="Read length", type=int, default=200
    )
    parser.add_argument(
        "-IDS",
        nargs=3,
        metavar=("I", "D", "S"),
        help="Insertion Deletion Substitution Probability",
        type=float,
        default=(1 / 3, 1 / 3, 1 / 3),
    )

    args = parser.parse_args()

    return args.genome_file, args.N, args.read_length, args.IDS


def load_genome(genome_file):
    file_G = open(genome_file, "r")
    G = "".join(["".join(l.split()) for l in file_G.readlines()[1:]])
    file_G.close()
    return G


def main():
    err_inser_del_sub_list = [round(0.01 * i, 2) for i in range(11)]
    err_homopoly_list = [round(0.05 + 0.05 * i, 2) for i in range(6)]
    # err_homopoly_list = [0.1]
    genome_file, N, read_length, ids = parse_args()
    Go = load_genome(genome_file)
    err = Error_rate(0)
    Gi, qual_Gi, nb_IDS_Gi = add_IDS_err(Go, err)
    print(Go[:200])
    print(Gi[:200])
    print(qual_Gi[:200])
    print(f"biological distance overall: {nb_IDS_Gi}/{len(Gi)}")
    for err_hom in err_homopoly_list:
        err = Error_rate(err_hom)
        generate_read(Go, Gi, qual_Gi, N, read_length, err)


if __name__ == "__main__":
    main()
