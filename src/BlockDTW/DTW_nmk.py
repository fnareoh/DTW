#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Computes a matrix of DTW distances between a pattern P and a text T
    In the matrix a block is defined by a a^height x b^width letters to be compared
    Only distances <= k must be computed
    Computations of a block is in O(height + width).
"""

__author__ = "Garance Gourdel, Pierre Peterlongo"
__email__ = "pierre.peterlongo@inria.fr, garance.gourdel@inria.fr"

import sys, logging

from .borderblock import BorderBlock, unpack


def run_length_compress(T):
    text, length = [], []
    last_char = ""
    nb_occ = 0
    for c in T:
        if c != last_char:
            text.append(last_char)
            length.append(nb_occ)
            last_char, nb_occ = c, 0
        nb_occ += 1
    text.append(last_char)
    length.append(nb_occ)
    return text, length


class DtwByBorders:
    """
    Computes the DTW matrix for a pattern P and a text T in O(#runsP * #runsT * k) time and space.
    """

    def __init__(self, P, T, max_value=sys.maxsize, rle_input=False):
        """Init and compute a DTW matrix

        Computes only borders of blocks in a compressed format taking O(k) space for each blocks.
        Here a match is considered a 0 and a mismatch as 1.

        Args:
            P ([str]): Pattern string -> vertical in the matrix
            T ([str]): target string -> horizontal in the matrix
            max_value (int): maximal value to be computed in the matrix
        """
        self.rle_input = rle_input
        self.P = P
        self.T = T
        self.max_value = max_value
        self.cost = 1  # Here a match is considered a 0 and a mismatch as 1.

        if not rle_input:
            self.rle_P, self.height_P = run_length_compress(self.P)
            self.rle_T, self.width_T = run_length_compress(self.T)

        self.M_block = [[None for _ in self.rle_T] for _ in self.rle_P]

        self.M_block[0][0] = BorderBlock(1, 1, 1, 0, [(0, 0)], [(0, 0)])
        for i, h in enumerate(self.height_P[1:]):
            self.M_block[i + 1][0] = BorderBlock(
                h, 1, 1, self.max_value, [(self.max_value, 0)], [(self.max_value, 0)]
            )
        for j, w in enumerate(self.width_T[1:]):
            self.M_block[0][j + 1] = BorderBlock(
                1, w, 0, 0, [(self.max_value, 0)], [(self.max_value, 0)]
            )

        for i in range(1, len(self.M_block)):
            for j in range(1, len(self.M_block[i])):
                logging.debug(f"i:{i} j:{j}")
                logging.debug(
                    f"P[i]:{self.rle_P[i]} T[j]:{self.rle_T[j]} cost:{int(self.rle_P[i] != self.rle_T[j])}"
                )
                logging.debug(f"Northwest:{self.M_block[i - 1][j - 1].Vse}")

                self.M_block[i][j] = BorderBlock(
                    self.height_P[i],
                    self.width_T[j],
                    int(self.rle_P[i] != self.rle_T[j]),
                    self.M_block[i - 1][j - 1].Vse,
                    self.M_block[i - 1][j].q_bottom,
                    self.M_block[i][j - 1].q_right,
                )
                logging.debug(self.M_block[i][j])

    def get_last_line(self):
        pos = 0
        line = []
        for j, w in enumerate(self.width_T[1:]):
            for (v, l) in self.M_block[-1][j + 1].q_bottom:
                if line == [] or v != line[-1][0]:
                    line.append((v, pos + l))
            pos += w
        return line


def main(P, T, max_value=sys.maxsize):
    """Computes a dtw matrix

    Args:
        P ([str]): pattern
        T ([str]): target
    """
    dtw = DtwByBorders(P, T, max_value)
    print(unpack(dtw.get_last_line(), len(T)))


if __name__ == "__main__":
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    else:
        sys.stderr.write(f"Usage: python {sys.argv[0]} P T <max_value>\n")
