from read import parse_input
from dynamic_prog.local_alignement import LocalED, LocalDTW
import csv
from sys import argv


def evaluate(localM, G, R):
    pos = localM.index_min_last_row()
    origin, alQ, alT = localM.trace_back(pos)
    print(
        f"Alignement (cost {localM.matrix[len(R.sequence)][pos]}) of Q ending at {pos} in T starts at {origin} in T."
    )


def compute_matrices(R, G):
    ldtw = LocalDTW(R.sequence, G)
    ldtw_rev_comp = LocalDTW(reverse(complement(R.sequence)), G)
    led = LocalED(R.sequence, G)
    led_rev_comp = LocalED(reverse(complement(R.sequence)), G)
    return ldtw, ldtw_rev_comp, led, led_rev_comp


def complement(S):
    compl = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    list_comp = [compl[S[i]] for i in range(len(S))]
    return "".join(list_comp)


def reverse(S):
    return S[::-1]


def print_details(R, G, ldtw, ldtw_rev_comp, led, led_rev_comp):
    print("******** Dynamic Time Warp **********")
    print("------------ Forward ----------------")
    evaluate(ldtw, G, R)
    print("------- Reverse complement ----------")
    evaluate(ldtw_rev_comp, G, R)
    print("******** Edit distance **********")
    print("------------ Forward ----------------")
    evaluate(led, G, R)
    print("------- Reverse complement ----------")
    evaluate(led_rev_comp, G, R)


def local_alignement_ED_vs_DTW():
    G, list_read_aligned, list_read_not_aligned = ecoli_badread()

    list_R = list_read_aligned[5:10]
    nb_small = 1  # Number of values to trace back
    for R in list_R:
        if R.end_position > len(G):
            continue
        print_details(R, G)


def stats_local_alignement_ED_vs_DTW(
    G, list_read_aligned, list_read_not_aligned, name_file
):

    output_file = open(name_file + "_score.csv", "w")
    csv_output = csv.writer(output_file)
    csv_output.writerow(["read", "identity percentage", "length", "min DTW", "min ED"])
    list_R = list_read_aligned
    nb_small = 5  # Number of values to trace back
    nb_out_of_range_read = 0
    avg_dist_origin_sam_alignement = 0
    avg_dist_dtw = 0
    avg_dist_ed = 0
    nb_read = 0
    for R in list_R:
        if nb_read % 10 == 0:
            print(f"{nb_read}/{len(list_R)}")
        nb_read += 1

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

        if hasattr(R, "read_identity"):
            csv_output.writerow([R.id, R.read_identity, R.length, score_dtw, score_ed])

        avg_dist_origin_sam_alignement += abs(R.position - R.sam_alignement)
        avg_dist_dtw += abs(R.position - start_dtw)
        avg_dist_ed += abs(R.position - start_ed)
        if abs(R.position - start_dtw) > 10 or abs(R.position - start_ed) > 10:
            print("******** Read information **********")
            print(R.detail_str())
            print(R.sequence)
            print_details(R, G, ldtw, ldtw_rev_comp, led, led_rev_comp)

    avg_dist_origin_sam_alignement /= nb_read
    avg_dist_dtw /= nb_read
    avg_dist_ed /= nb_read
    print(
        f"Average distance between origin and sam alignement found {round(avg_dist_origin_sam_alignement,2)}"
    )
    print(
        f"Average distance between origin and min DTW alignement {round(avg_dist_dtw,2)}"
    )
    print(
        f"Average distance between origin and min ED alignement {round(avg_dist_ed,2)}"
    )
    output_file.close()


def main():
    print("******** Parsing input **********")
    simulator = argv[1]
    genome = argv[2]
    read_pref = argv[3]
    data_prefix = f"../data/{simulator}/"
    G, list_read_aligned, list_read_not_aligned = parse_input(
        data_prefix + genome,
        data_prefix + read_pref,
        simulator,
    )
    if len(argv) == 5:
        print_read(argv[4], list_read_aligned + list_read_not_aligned, G)
    else:
        stats_local_alignement_ED_vs_DTW(
            G, list_read_aligned, list_read_not_aligned, read_pref
        )


def print_read(id, read_list, G):
    for R in read_list:
        if R.id == id:
            print("******** Read information **********")
            print(R.detail_str())
            print(R.sequence)
            ldtw, ldtw_rev_comp, led, led_rev_comp = compute_matrices(R, G)
            print_details(R, G, ldtw, ldtw_rev_comp, led, led_rev_comp)


if __name__ == "__main__":
    main()
