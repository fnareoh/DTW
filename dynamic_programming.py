#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# compatible python3

import bisect
import pysam
import sys


class Constant:
    match = 0
    mismatch = 1
    gap = 1


class Read:
    """ stores all information needed about a read """

    def __init__(
        self,
        sequence,
        fastq_text="",
        quality="",
        sam_id=-1,
        sam_alignement=0,
        original_pos=-1,
        original_len=-1,
    ):
        self.sequence = sequence
        self.quality = quality
        if fastq_text != "":
            info = fastq_text[1:].split(" ")
            self.id = info[0]
            assert self.id == sam_id
            self.strand = info[1].split(",")[1]
            self.position = int(info[1].split(",")[2].split("-")[0])
            self.length = int(info[2].split("=")[-1])
            self.error_free_length = int(info[3].split("=")[-1])
            self.read_identity = info[4].split("=")[-1]
        else:
            self.length = len(self.sequence) - 1
            self.strand = "+strand"
            self.position = original_pos
            self.error_free_length = original_len
            self.read_identity = "?%"
        assert self.length == len(self.sequence) - 1
        self.sam_alignement = sam_alignement

    def __str__(self):
        res = f"Read Id = {self.id}"
        return res

    def detail_str(self):
        res = f"Read Id = {self.id}"
        res += f"\nStrand = {self.strand}, Position={self.position}"
        res += f"\nLength = {self.length}, Error-free_length={self.error_free_length}"
        res += f", Read identity = {self.read_identity}"
        res += f"\nSequence = {self.sequence}"
        return res


def parse_input(genome_file, read_fastq_file, read_sam_file):
    file_G = open(genome_file, "r")
    G = file_G.readlines()[1]
    file_G.close()
    fastq_R = open(read_fastq_file, "r")
    sam_R = pysam.AlignmentFile(read_sam_file, "r")
    list_R = []
    nb_read_not_al = 0
    nb_read = 0
    for read in sam_R.fetch():
        nb_read += 1
        text = fastq_R.readline()
        assert text[0] == "@"
        seq = fastq_R.readline()
        optional = fastq_R.readline()
        assert optional[0] == "+"
        qualities = fastq_R.readline()
        R = Read(
            seq,
            fastq_text=text,
            quality=qualities,
            sam_id=read.query_name,
            sam_alignement=read.reference_start,
        )
        if R.sam_alignement == -1:
            nb_read_not_al += 1
        list_R.append(R)
    print(f"Number of read not alligned in the sam file = {nb_read_not_al}")
    print(f"Total number of read = {nb_read}")
    return G, list_R


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
        for i in range(1, len(self.matrix)):
            self.matrix[i][0] = self.matrix[i - 1][0] + self.score(S[0], self.T[j - 1])
        for j in range(1, len(self.matrix[0])):
            self.matrix[0][j] = self.matrix[0][j - 1] + self.score(S[i - 1], self.T[0])

    def fill_DTW(self):
        """ fills the matrix for global alignment (Needleman & Wunsch algo)"""
        for i in range(1, len(self.S) + 1):
            # i-th line
            for j in range(1, len(self.T) + 1):
                # j-th column
                self.matrix[i][j] = min(
                    self.matrix[i - 1][j - 1],
                    self.matrix[i][j - 1],
                    self.matrix[i - 1][j],
                ) + self.score(self.S[i - 1], self.T[j - 1])

        return self.matrix[len(self.S)][len(self.T)]

    #### HOMO EDIT DISTANCE ####
    def homoedit_to_empty(self, x):
        """Homo-edt distance from a string
        to an empty string"""

        n = len(x)
        H = [[0 for i in range(n)] for j in range(n)]
        for i in range(n):
            H[i][i] = 1
        for j in range(2, n + 1):
            for i in range(n - j + 1):
                C = [
                    (H[i][k] + H[k + 1][i + j - 1] - int(x[i] == x[i + j - 1]))
                    for k in range(i, i + j - 1)
                ]
                H[i][i + j - 1] = min(C)
        return H

    def homoedit(self):
        """Homo-edit distance between strings x, y equals
        the minimum number of run insertions and
        deletions required to convert x into y"""
        H_S = self.homoedit_to_empty(self.S)
        H_T = self.homoedit_to_empty(self.T)

        n = len(self.S)
        m = len(self.T)
        # D = [[0 for j in range(m+1)] for i in range(n+1)]

        # fill row 0 and column 0
        for i in range(1, n + 1):
            self.matrix[i][0] = H_S[0][i - 1]
        for j in range(1, m + 1):
            self.matrix[0][j] = H_T[0][j - 1]

        # fill the rest of the table
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                C = []
                if self.S[i - 1] == self.T[j - 1]:
                    C.append(self.matrix[i - 1][j - 1])
                for k in range(i):
                    C.append(self.matrix[k][j] + H_S[k][i - 1])
                for k in range(j):
                    C.append(self.matrix[i][k] + H_T[k][j - 1])
                self.matrix[i][j] = min(C)

        return self.matrix[n][m]

    #### GENERIC PRINT ####
    def printGlobalAln(self):
        """ prints a global alignment of best score
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

    print("******** Homo edit **********")

    dm_HE = DynamicMatrix(S, T, 0, 1, 1)
    score = dm_HE.homoedit()

    print(dm_HE)
    print(f"Best score is {score}")
    alT, alS, pcId = dm_HE.printGlobalAln()
    print(alT, alS, sep="\n")
    print(f"%Id = {pcId:.2f}")


def print_list_alignements(name, original_al, closest_al, R, G):
    print(f"G={G}")
    print(f"R={R} \n")
    for al in closest_al[name]:
        print(al, "\n")
    print("----------")
    print("Alignement to the original substring:", original_al[name], sep="\n\n")


def add_al(name, original_al, closest_al, al, original_len_read, original_pos_read, N):
    if al.l == original_len_read and al.pos == original_pos_read:
        original_al[name] = al
    if len(closest_al[name]) < N:
        bisect.insort(closest_al[name], al)
    elif al.dist < closest_al[name][-1].dist:
        closest_al[name].pop()
        bisect.insort(closest_al[name], al)


def evaluate_all(G, read, k, N):
    """ Computes and prints the N closest alignment of R to all
    substring of G of length len(R) + or -k for all distance """
    R = read.sequence
    original_pos_read = read.position  # original position the read was extracted from
    original_len_read = (
        read.error_free_length
    )  # the original length of the extracted read
    # (before modification)

    cost = Constant()
    if not -k <= original_len_read - len(R) <= k:
        raise ValueError(
            "the original length of the extracted read must be in len(read) + or - k"
        )
    print(
        f"Computing the {N} closest subsequences of G of length"
        + f" {len(R)} +- k={k} "
    )
    closest_al = {"NW": [], "DTW": [], "HE": []}
    original_al = {}
    nb_length_computed = 0
    for l in range(len(R) - k, len(R) + k + 1):
        print(
            f"Number of length of subsequence computed = {nb_length_computed} / {2*k+1}"
        )
        nb_length_computed += 1
        for pos in range(0, len(G) - l):
            S = G[pos : pos + l]
            # Needleman and Wunsch
            dm_NW = DynamicMatrix(S, R, cost.match, cost.mismatch, cost.gap)
            dm_NW.initGlobal_NW()
            score = dm_NW.fill_NW()
            add_al(
                "NW",
                original_al,
                closest_al,
                Alignement(score, pos, l, S),
                original_len_read,
                original_pos_read,
                N,
            )
            # Dynamic time warping distance
            dm_DTW = DynamicMatrix(S, R, cost.match, cost.mismatch, cost.gap)
            dm_DTW.initGlobal_DTW()
            score = dm_DTW.fill_DTW()
            add_al(
                "DTW",
                original_al,
                closest_al,
                Alignement(score, pos, l, S),
                original_len_read,
                original_pos_read,
                N,
            )
            # Homo edit
            dm_HE = DynamicMatrix(S, R, cost.match, cost.mismatch, cost.gap)
            score = dm_HE.homoedit()
            add_al(
                "HE",
                original_al,
                closest_al,
                Alignement(score, pos, l, S),
                original_len_read,
                original_pos_read,
                N,
            )

    print("\n******** Needleman Wunsch **********")
    print_list_alignements("NW", original_al, closest_al, R, G)
    print("\n******** Dynamic Time Warp **********")
    print_list_alignements("DTW", original_al, closest_al, R, G)
    print("\n******** Homo edit **********")
    print_list_alignements("HE", original_al, closest_al, R, G)


def main():
    print("******** Parsing input **********")
    # G, list_read = parse_input(
    #    "data/ecoli_10kb.fa", "data/reads_coli.fastq", "data/align_reads_coli.sam"
    # )
    G = "AAAACCTGGTAATGCTGATTAGCCGCACCGTTTTTACCCGTACGCGGACCTGTATGATGATTTCACCAAGTGCTGACGGGTTGATACCCTGTTGAT"  # genome
    R = Read("CCGCCCCACCGTTTTTAAAACCCGT", original_pos=23, original_len=14)
    k = 15  # the length threshold on the subsequence we consider
    N = 5  # the maximal distance to a subsequence

    list_R = [R]
    for R in list_R:
        evaluate_all(G, R, k, N)


if __name__ == "__main__":
    main()
