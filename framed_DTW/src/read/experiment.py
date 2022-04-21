from read import Read, parse_input
from dynamic_prog.local_alignement import LocalED, LocalDTW


def evaluate(localM, G, R, nb_small):
    """
    Computes the n_smallest position in the last line and traces them back to
    find the starting position of the local_alignement
    """
    localM.fill()
    print(
        f"Score at the original ending position {localM.matrix[len(R.sequence)][R.end_position]}"
    )
    min_localM = sorted(
        range(len(localM.matrix[len(R.sequence)])),
        key=lambda sub: localM.matrix[len(R.sequence)][sub],
    )[:nb_small]
    for pos in min_localM:
        origin, alQ, alT = localM.trace_back(pos)
        print(
            f"Alignement (cost {localM.matrix[len(R.sequence)][pos]}) of Q ending at {pos} in T starts at {origin} in T."
        )


def compute_origin_min_position(localM, G, R):
    localM.fill()
    min_localM = localM.matrix[len(R.sequence)].index(
        min(localM.matrix[len(R.sequence)])
    )
    origin, alQ, alT = localM.trace_back(min_localM)
    return origin


def local_alignement_ED_vs_DTW():
    print("******** Parsing input **********")
    data_prefix = "../data/badread/"
    G, list_read_aligned, list_read_not_aligned = parse_input(
        data_prefix + "ecoli_10kb.fa",
        data_prefix + "reads_coli.fastq",
        data_prefix + "align_reads_coli.sam",
    )

    list_R = list_read_aligned[0:10]
    nb_small = 10  # Number of values to trace back
    for R in list_R:
        if R.end_position > len(G):
            continue
        print("******** Read information **********")
        print(f"Original position the read was extracted from {R.position}")
        print(f"Sam alignement found (-1 means no alignement found) {R.sam_alignement}")
        print(f"Original ending position the read was extracted from {R.end_position}")
        print("******** Dynamic Time Warp **********")
        ldtw = LocalDTW(R.sequence, G)
        evaluate(ldtw, G, R, nb_small)
        print("******** Edit distance **********")
        led = LocalED(R.sequence, G)
        evaluate(led, G, R, nb_small)


def stats_local_alignement_ED_vs_DTW():
    print("******** Parsing input **********")
    data_prefix = "../data/badread/"
    G, list_read_aligned, list_read_not_aligned = parse_input(
        data_prefix + "ecoli_10kb.fa",
        data_prefix + "reads_coli.fastq",
        data_prefix + "align_reads_coli.sam",
    )

    list_R = list_read_aligned + list_read_not_aligned
    nb_small = 10  # Number of values to trace back
    nb_out_of_range_read = 0
    avg_dist_origin_sam_alignement = 0
    avg_dist_dtw = 0
    avg_dist_ed = 0
    nb_read = 0
    for R in list_R:
        if nb_read % 10 == 0:
            print(f"{nb_read}/{len(list_R)}")
        nb_read += 1
        if R.end_position > len(G):
            nb_out_of_range_read += 1
            continue
        avg_dist_origin_sam_alignement += abs(R.position - R.sam_alignement)
        ldtw = LocalDTW(R.sequence, G)
        avg_dist_dtw += abs(R.position - compute_origin_min_position(ldtw, G, R))
        led = LocalED(R.sequence, G)
        avg_dist_ed += abs(R.position - compute_origin_min_position(led, G, R))
    avg_dist_origin_sam_alignement /= len(list_R) - nb_out_of_range_read
    avg_dist_dtw /= len(list_R) - nb_out_of_range_read
    avg_dist_ed /= len(list_R) - nb_out_of_range_read
    print(f"Genome size: {len(G)}")
    # print(
    #    f"Average distance between origin and sam alignement found {round(avg_dist_origin_sam_alignement,2)}"
    # )
    print(
        f"Average distance between origin and min DTW alignement {round(avg_dist_dtw,2)}"
    )
    print(
        f"Average distance between origin and min ED alignement {round(avg_dist_ed,2)}"
    )


if __name__ == "__main__":
    stats_local_alignement_ED_vs_DTW()
