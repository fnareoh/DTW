from dynamic_prog.pattern_matching import PM_DTW
from BlockDTW.DTW_blocks import *
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


def main(nb_tests, size_min, bound_homopol):
    # nb_tests = 10
    # size_min = 1000
    # bound_homopol = 10
    PT_strings = []
    for _ in range(nb_tests):
        PT_strings.append(
            [
                get_random_string(size_min, bound_homopol),
                get_random_string(size_min, bound_homopol),
            ]
        )

    classic_res = []
    # compare times
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
            dtw = DtwByBlocks(P, T)
            dtw_block_res.append(dtw.get_br_value())
            nb_blocks += dtw.get_nb_blocks()
        update_progress(1)
    block_time.print("Durée with blocks = {} seconds")
    timeb = block_time.value()
    print(f"Durée per block = {round(timeb,2)}/{nb_blocks} = {timeb/nb_blocks} seconds")

    # Validations
    for i in range(nb_tests):
        assert (
            classic_res[i] == dtw_block_res[i]
        ), f"Test failed with {PT_strings[i][0]} {PT_strings[i][1]}: {classic_res[i]} vs {dtw_block_res[i]}"

    print(f"We performed {nb_tests} tests, all passed !")


if __name__ == "__main__":

    if len(sys.argv) != 4:
        sys.stderr.write(
            f"Usage: python {sys.argv[0]} nb_tests min_size_P_T homopol_size_bound\n"
        )
    else:
        main(
            nb_tests=int(sys.argv[1]),
            size_min=int(sys.argv[2]),
            bound_homopol=int(sys.argv[3]),
        )
