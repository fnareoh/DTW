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

from .border_block import BorderBlock


class DtwByBorders:
    """
    Computes the DTW matrix for a pattern P and a text T in O(#runsP * #runsT * k) time and space.
    """

    def __init__(self, P, T, max_value=sys.maxsize):
        """Init and compute a DTW matrix

        Computes only borders of blocks in a compressed format taking O(k) space for each blocks.
        Here a match is considered a 0 and a mismatch as 1.

        Args:
            P ([str]): Pattern string -> vertical in the matrix
            T ([str]): target string -> horizontal in the matrix
            max_value (int): maximal value to be computed in the matrix
        """

        self.P = P
        self.T = T
        self.max_value = max_value
        self.cost = 1  # Here a match is considered a 0 and a mismatch as 1.

        # TODO: RLE of P and T (idealy put an option in the input format)

        # TODO: two for loops for the compressed matrix of blocks
