from dynamic_prog.pattern_matching import PM_DTW
from BlockDTW.DTW_blocks import *
from BlockDTW.DTW_nmk import *
from timer import Timer
from progress_bar import update_progress

import random


def get_random_string(minsize, bound_homopol):
    random_string = ""

    while len(random_string) < minsize:
        # Considering only upper letters
        random_integer = random.randint(97, 100)
        nb_occurrences = random.randint(1, bound_homopol)
        for _ in range(nb_occurrences):
            random_string += chr(random_integer)
    return random_string


def main(nb_tests, size_min, bound_homopol, test_classic=True, max_value=sys.maxsize):
    # nb_tests = 10
    # size_min = 1000
    # bound_homopol = 10
    PT_strings = []
    for _ in range(nb_tests):
        size_p = random.randint(1, size_min)
        PT_strings.append(
            [
                get_random_string(size_p, bound_homopol),
                get_random_string(size_min, bound_homopol),
            ]
        )

    classic_res = []
    # compare times
    if test_classic:
        with Timer() as classical_time:
            for i in range(nb_tests):
                update_progress(i / float(nb_tests))
                P = PT_strings[i][0]
                T = PT_strings[i][1]
                ldtw = PM_DTW(P, T)
                ldtw.fill()
                classic_res.append(ldtw.get_last_value())
            update_progress(1)
        classical_time.print("Durée classique = {} seconds")

    dtw_block_res = []
    nb_blocks = 0
    with Timer() as block_time:
        for i in range(nb_tests):
            update_progress(i / float(nb_tests))
            P = PT_strings[i][0]
            T = PT_strings[i][1]
            dtw = DtwByBlocks(P, T, max_value)
            dtw_block_res.append(dtw.get_br_value())
            nb_blocks += dtw.get_nb_blocks()
        update_progress(1)
    block_time.print("Durée with blocks = {} seconds")
    timeb = block_time.value()
    print(f"Durée per block = {round(timeb,2)}/{nb_blocks} = {timeb/nb_blocks} seconds")

    dtw_border_res = []
    with Timer() as border_time:
        for i in range(nb_tests):
            update_progress(i / float(nb_tests))
            P = PT_strings[i][0]
            T = PT_strings[i][1]
            dtw_borders = DtwByBorders(P, T, max_value)
            dtw_border_res.append(dtw_borders.M_block[-1][-1].Vse)
        update_progress(1)
    border_time.print("Durée with borders compressed = {} seconds")
    timeb = border_time.value()

    # Validations
    for i in range(nb_tests):
        if test_classic:
            assert (
                classic_res[i] == dtw_block_res[i] == dtw_border_res[i]
            ), f"Test failed with {PT_strings[i][0]} {PT_strings[i][1]}: {classic_res[i]} vs {dtw_block_res[i]} vs {dtw_border_res[i]}"
        else:
            assert (
                dtw_block_res[i] == dtw_border_res[i]
            ), f"Test failed with {PT_strings[i][0]} {PT_strings[i][1]}: {dtw_block_res[i]} vs {dtw_border_res[i]}"

    print(f"We performed {nb_tests} tests, all passed !")


if __name__ == "__main__":

    if len(sys.argv) not in [4, 5]:
        sys.stderr.write(
            f"Usage: python {sys.argv[0]} nb_tests min_size_P_T homopol_size_bound (optional_max_value)\n"
        )
    elif len(sys.argv) == 4:
        main(
            nb_tests=int(sys.argv[1]),
            size_min=int(sys.argv[2]),
            bound_homopol=int(sys.argv[3]),
        )
    else:
        main(
            nb_tests=int(sys.argv[1]),
            size_min=int(sys.argv[2]),
            bound_homopol=int(sys.argv[3]),
            test_classic=False,
            max_value=int(sys.argv[4]),
        )
