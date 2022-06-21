#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Computes a matrix of DTW distances between a pattern P and a text T
    In the matrix a block is defined by a a^height x b^width letters to be compared
    Only distances <= k must be computed
    Computations of a block is in O(height + width).
"""

__author__ = "Garance Gourdel, Pierre Peterlongo"
__email__ = "pierre.peterlongo@inria.fr, garance.gourdel@inria.fr"

import sys

from .block import Block


class DtwByBlocks:
    """
    Computes the DTW matrix for a pattern P and a text T in O(#runsP |T| + |P| #runsT) time and space.
    """

    def __init__(self, P, T, max_value=sys.maxsize):
        """Init and compute a DTW matrix

        Computes only borders of blocks (a block is computed for each a^height of P versus b^width of T)
        Each block is computed in O(height + width)

        Here a match is considered a 0 and a mismatch as 1.
        Args:
            P ([str]): Pattern string -> vertical in the matrix
            T ([str]): target string -> horizontal in the matrix
            max_value (int): maximal value to be computed in the matrix
        """

        self.P = P
        self.T = T
        self.max_value = max_value
        #      T
        #  --------
        #  |
        # P |
        #  |

        # Determine block positions:
        # list of vertical frontiers:
        self.end_vertical_blocks = []
        for i in range(1, len(self.T)):
            if self.T[i] != self.T[i - 1]:
                self.end_vertical_blocks.append(i - 1)
        self.end_vertical_blocks.append(len(self.T) - 1)

        # list of horizontal frontiers:
        self.end_horizontal_blocks = []
        for i in range(1, len(self.P)):
            if self.P[i] != self.P[i - 1]:
                self.end_horizontal_blocks.append(i - 1)
        self.end_horizontal_blocks.append(len(self.P) - 1)

        # Future array of blocks:
        self.block_matrix = [[] for i in range(len(self.end_horizontal_blocks))]
        for i in range(len(self.end_horizontal_blocks)):
            self.block_matrix[i] = [None for j in range(len(self.end_vertical_blocks))]

        # Compute the blocks:
        self.__compute_blocks__()

    def __repr__(self):
        """Provides a string view of all blocks

        Returns:
            [str]: string view of all blocks
        """
        res = ""
        res += f"P (vertical): {self.P}\n"
        res += f"T (horizontal): {self.T}\n"
        res += f"end_vertical_blocks: {self.end_vertical_blocks}\n"
        res += f"end_horizontal_blocks: {self.end_horizontal_blocks}\n"

        res += f"bottom right value {self.block_matrix[-1][-1].br}\n"

        for h_block_id in range(len(self.end_vertical_blocks)):
            for v_block_id in range(len(self.end_horizontal_blocks)):
                res += f"*** Block {v_block_id} {h_block_id} ***\n{self.block_matrix[v_block_id][h_block_id]}\n"
        return res

    def get_br_value(self):
        return self.block_matrix[-1][-1].br

    def get_nb_blocks(self):
        return len(self.end_horizontal_blocks) * len(self.end_vertical_blocks)

    def __compute_blocks__(self):
        """Compute all blocks from a matrix"""
        # We compute blocks from top to bottom
        ## When we navigate in a line of blocks (h_block_id), blocks are delimited by the end of vertical blocks.
        ## This explain why h_block_id is in end_vertical_blocks
        ## symmetrically explanation for column of blocks
        for h_block_id in range(len(self.end_vertical_blocks)):
            for v_block_id in range(len(self.end_horizontal_blocks)):
                # get 3 top values:
                if v_block_id == 0:  # First column blocks
                    Vnw = 0
                    Vn = 0
                    if h_block_id == 0:  # first block
                        Vw = sys.maxsize
                    else:
                        Vw = self.block_matrix[0][h_block_id - 1].tr
                elif h_block_id == 0:  # First line blocks:
                    # case h_block_id == 0 already performed
                    Vw = sys.maxsize
                    Vnw = sys.maxsize
                    Vn = self.block_matrix[v_block_id - 1][0].bl
                else:  # any other bloc:
                    Vn = self.block_matrix[v_block_id - 1][h_block_id].bl
                    Vw = self.block_matrix[v_block_id][h_block_id - 1].tr
                    Vnw = self.block_matrix[v_block_id - 1][h_block_id - 1].br

                if v_block_id == 0:
                    line_start = 0
                else:
                    line_start = self.end_horizontal_blocks[v_block_id - 1] + 1
                line_end = self.end_horizontal_blocks[v_block_id]
                if h_block_id == 0:
                    col_start = 0
                else:
                    col_start = self.end_vertical_blocks[h_block_id - 1] + 1
                col_end = self.end_vertical_blocks[h_block_id]
                height = line_end - line_start + 1
                width = col_end - col_start + 1

                # get previous cuts
                if v_block_id == 0:
                    h_cuts = []  # we are on the first line, no cut upside
                else:
                    h_cuts = self.block_matrix[v_block_id - 1][h_block_id].bottom_cuts
                if h_block_id == 0:
                    v_cuts = [
                        i for i in range(height)
                    ]  # we are on the first column, only cuts on the left.
                else:
                    v_cuts = self.block_matrix[v_block_id][
                        h_block_id - 1
                    ].rightmost_cuts

                T_letter = self.T[col_end]
                P_letter = self.P[line_end]

                # Block(l, h, False, Vnw, Vw, Vn, h_cuts, v_cuts, line_start, line_end, column_start, column_end)
                current_block = Block(
                    height=height,
                    width=width,
                    equals=(P_letter == T_letter),
                    Vnw=Vnw,
                    Vw=Vw,
                    Vn=Vn,
                    h_cuts=h_cuts,
                    v_cuts=v_cuts,
                    max_value=self.max_value,
                    line_start=line_start,
                    line_end=line_end,
                    column_start=col_start,
                    column_end=col_end,
                )

                self.block_matrix[v_block_id][h_block_id] = current_block


def main(P, T, max_value=sys.maxsize):
    """Computes a dtw matrix

    Args:
        P ([str]): pattern
        T ([str]): target
    """
    dtw = DtwByBlocks(P, T, max_value)
    print(dtw)


if __name__ == "__main__":
    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    else:
        sys.stderr.write(f"Usage: python {sys.argv[0]} P T <max_value>\n")
