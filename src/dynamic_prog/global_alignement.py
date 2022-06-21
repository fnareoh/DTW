#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# compatible python3

__author__ = "Garance Gourdel, Pierre Peterlongo"
__email__ = "pierre.peterlongo@inria.fr, garance.gourdel@inria.fr"

import bisect
import sys


class Constant:
    match = 0
    mismatch = 1
    gap = 1


class Alignement:
    """ stores an alignement of a read to a given genome """

    def __init__(self, dist, pos, l, S):
        self.dist = dist
        self.pos = pos
        self.l = l
        self.S = S

    def __lt__(self, other):
        return self.dist < other.dist

    def __str__(self):
        res = f"distance between S and R = {self.dist}"
        res += f"\nposition = {self.pos} len = {self.l}"
        res += f"\nS={self.S}"
        return res


class DynamicMatrix:
    """stores a matrix |S|x|T| (|S|+1 lines and |T|+1columns),
    sequences S and T and the score system (match, mismatch, gap)
    defines some global alignment functions
    """

    def __init__(self, S, T, match, mismatch, gap):
        """ defines and stores initial values"""

        self.S = S
        self.T = T
        self.gap = gap
        self.match = match
        self.mismatch = mismatch

        self.matrix = [0 for i in range(len(S) + 1)]
        for i in range(len(S) + 1):
            self.matrix[i] = [0 for j in range(len(T) + 1)]
        # self.matrix[i][j]
        # i : ith line, in S
        # j : jth column, in T
        # self.matrix[i+1][j+1] in front of S[i] et T[j]

    def __str__(self):
        """ string the matrix"""
        res = ""
        width = 4
        vide = " "
        line = f"{vide:>{2*width}}"
        for j in range(0, len(self.T)):
            line += f"{self.T[j]:>{width}}"
        res += line + "\n"
        line = f"{vide:>{width}}"
        for j in range(0, len(self.T) + 1):
            line += f"{self.matrix[0][j]:>{width}}"
        res += line + "\n"
        for i in range(1, len(self.S) + 1):
            line = f"{self.S[i-1]:>{width}}"
            for j in range(0, len(self.T) + 1):
                line += f"{self.matrix[i][j]:>{width}}"
            res += line + "\n"
        return res

    def score(self, char1, char2):
        """ returns score of a match if char1 == char2, else the score of a mismatch"""
        if char1 == char2:
            return self.match
        else:
            return self.mismatch

    #### Needleman & Wunsch ####

    def initGlobal_NW(self):
        """ initializes first line and first columns for a global alignment"""
        for i in range(1, len(self.matrix)):
            self.matrix[i][0] = self.matrix[i - 1][0] + self.gap
        for j in range(1, len(self.matrix[0])):
            self.matrix[0][j] = self.matrix[0][j - 1] + self.gap

    def fill_NW(self):
        """ fills the matrix for global alignment (Needleman & Wunsch algo)"""
        for i in range(1, len(self.S) + 1):
            # i-th line
            for j in range(1, len(self.T) + 1):
                # j-th column
                self.matrix[i][j] = min(
                    self.matrix[i - 1][j - 1]
                    + self.score(self.S[i - 1], self.T[j - 1]),
                    self.matrix[i][j - 1] + self.gap,
                    self.matrix[i - 1][j] + self.gap,
                )

        return self.matrix[len(self.S)][len(self.T)]

    #### DTW ####

    def initGlobal_DTW(self):
        """ initializes first line and first columns for a DTW"""
        self.initGlobal_NW()

    def fill_DTW(self):
        """ fills the matrix for global alignment (Needleman & Wunsch algo)"""
        for i in range(1, len(self.S) + 1):
            # i-th line
            for j in range(1, len(self.T) + 1):
                # j-th column
                self.matrix[i][j] = (
                    min(
                        self.matrix[i - 1][j - 1],
                        self.matrix[i][j - 1],
                        self.matrix[i - 1][j],
                    )
                    + self.score(self.S[i - 1], self.T[j - 1])
                )

        return self.matrix[len(self.S)][len(self.T)]

    #### GENERIC PRINT ####
    def printGlobalAln(self):
        """prints a global alignment of best score
        and returns the % of id
        """
        i = len(self.S)
        j = len(self.T)
        alS = ""
        alT = ""
        nb_matchs = 0
        while i > 0 and j > 0:
            # test first the diagonal, in order to favor diag in case of ex-aequos
            if (
                self.matrix[i - 1][j - 1] + self.score(self.S[i - 1], self.T[j - 1])
                == self.matrix[i][j]
            ):
                # diag
                alS = self.S[i - 1] + alS
                alT = self.T[j - 1] + alT
                if self.S[i - 1] == self.T[j - 1]:
                    nb_matchs += 1
                i -= 1
                j -= 1
            elif self.matrix[i][j - 1] + self.gap == self.matrix[i][j]:
                # left
                alS = "-" + alS
                alT = self.T[j - 1] + alT
                j -= 1
            else:
                # up
                alT = "-" + alT
                alS = self.S[i - 1] + alS
                i -= 1
        # here either i=0 or j=0

        # dealing with the case of j>0 and i=0
        while j > 0:
            alS = "-" + alS
            alT = self.T[j - 1] + alT
            j -= 1

        # dealing with the case of i>0 and j=0
        while i > 0:
            alT = "-" + alT
            alS = self.S[i - 1] + alS
            i -= 1

        return alT, alS, (float(100 * nb_matchs) / len(alS))


def demo():

    print("******** Needleman Wunsch **********")
    S = "CATGACTAAAG"
    T = "TACTAG"

    print(len(sys.argv))
    if len(sys.argv) == 3:
        S = sys.argv[1]
        T = sys.argv[2]

    dm_NW = DynamicMatrix(S, T, 0, 1, 1)
    dm_NW.initGlobal_NW()
    score = dm_NW.fill_NW()
    print(dm_NW)
    print(f"Best score is {score}")
    alT, alS, pcId = dm_NW.printGlobalAln()
    print(alT, alS, sep="\n")
    print(f"%Id = {pcId:.2f}")

    print("******** Dynamic Time Warp **********")

    dm_DTW = DynamicMatrix(S, T, 0, 1, 1)
    dm_DTW.initGlobal_DTW()
    score = dm_DTW.fill_DTW()
    print(dm_DTW)
    print(f"Best score is {score}")
    alT, alS, pcId = dm_DTW.printGlobalAln()
    print(alT, alS, sep="\n")
    print(f"%Id = {pcId:.2f}")
