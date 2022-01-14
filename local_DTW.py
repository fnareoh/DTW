"""Local alignement according to DTW

This script computes the local alignement matrix of a pattern, given by the
user, against a text, given by the user.

This tool requires two arguments:
- a string T: the text to allign the pattern against.
- a string Q: the pattern to align.
"""

import sys
import os
import random


def run_frontiers(S):
    """ Returns a list of the frontieres of the runs in S """
    frontiers = [0]
    # first column is considered as a special block
    for i in range(1, len(S)):
        if S[i] != S[i - 1]:
            frontiers.append(i)
    frontiers.append(len(S))
    return frontiers


class LocalDTW:
    """stores a matrix |Q|x|T| (|Q|+1 lines and |T|+1columns),
    sequences Q and T and the score system (match, mismatch, gap)
    defines some global alignment functions
    """

    def __init__(self, Q, T, mismatch):
        """
        Defines and computes the following values:
        Q: the pattern
        T: the text
        mismatch: the cost for a mismatch
        match: the cost for a match
        matrix: the distance matrix
        end_vertical_blocs: the list of vertical frontiers
        end_horizontal_blocs: the list of horizontal frontiers
        """
        self.Q = Q  # the pattern
        self.T = T  # the text
        self.mismatch = mismatch  # cost for mismatch
        self.match = 0
        # Shape of the matrix:
        #        T
        #    --------
        #   |
        # Q |
        #   |

        self.matrix = [[0 for j in range(len(T) + 1)] for i in range(len(Q) + 1)]
        # self.matrix[i][j]
        # i : ith line, in Q
        # j : jth column, in T
        # self.matrix[i+1][j+1] in front of Q[i] et T[j]

        # Determine block positions:
        # list of vertical frontiers:
        self.end_vertical_blocs = run_frontiers(T)
        # list of horizontal frontiers:
        self.end_horizontal_blocs = run_frontiers(Q)
        self.init_local_DTW()  # init the row of matrix
        self.fill_DTW()  # fill the matrix according to dynamic programming

    def come_from(self, i, j):
        """
        Returns a string containing three numbers indicating the cells from
        which matrix[i][j] may takes its value, in the following order:
        top, top-left, left
        """
        alpha = self.Q[i - 1]
        beta = self.T[j - 1]
        equals = alpha == beta
        if equals:
            d = self.match
        else:
            d = self.mismatch
        val = self.matrix[i][j]

        top = "0"
        if val == self.matrix[i - 1][j] + d:
            top = "1"
        topleft = "0"
        if val == self.matrix[i - 1][j - 1] + d:
            topleft = "1"
        left = "0"
        if val == self.matrix[i][j - 1] + d:
            left = "1"
        return top + topleft + left

    def latex_str(self):
        res = "\n\n% LATEX VERSION OF THE MATRIX\n"
        res += "\\documentclass{article}\n"
        res += "\\usepackage{graphicx}\n"
        res += "\\usepackage{a4wide}\n"
        res += "\\usepackage[dvipsnames]{xcolor}"
        res += "\\begin{document}\n"
        res += "\\begin{table}[]\n"
        res += "\\begin{center}\n"
        res += "\\resizebox{\\textwidth}{!}{"
        # res += "\\scalebox{0.5}{\n"
        res += "\\begin{tabular}{|cc||"

        for j in range(len(self.T)):
            res += "c"
            if j + 1 in self.end_vertical_blocs:
                res += "|"
        res += "}\n\\hline\n"

        res += " &   "
        for j in range(len(self.T)):
            res += f"& {self.T[j]} "
        res += "\\\\\n"

        res += " & 0 "
        for j in range(len(self.T)):
            res += "& 0 "
        res += "\\\\\n\\hline\n"

        # all other lines
        for i in range(1, len(self.Q) + 1):
            res += f"{self.Q[i-1]} "
            for j in range(0, len(self.T) + 1):
                if self.matrix[i][j] == sys.maxsize:
                    res += " & $\infty$ "
                else:
                    if self.Q[i - 1] == self.T[j - 1]:
                        res += f" &\\textcolor{{red}}{{ {self.matrix[i][j]} }}"  # {self.come_from(i,j)} }}"
                    else:
                        res += f" & {self.matrix[i][j]}"  # {self.come_from(i,j)}"
            res += "\\\\\n"
            if i in self.end_horizontal_blocs:
                res += "\\hline\n"

        res += "\n\\end{tabular} }\n"
        res += "\\end{center}\n"
        res += f"\\caption{{DTW for T={self.T} and Q={self.Q}. 3-value vector indicates the possible coming direction of each value (top, topleft, left) }}\n"
        res += "\\end{table}\n"

        res_q, res_align, res_t = self.get_alignment()
        res += "\\begin{verbatim}\n"
        res += res_t + "\n" + res_align + "\n" + res_q + "\n"
        res += "\\end{verbatim}\n"
        res += "\\end{document}\n"
        res += "% END LATEX VERSION OF THE MATRIX\n\n"

        return res

    def output_latex(self):
        outfilename = "res.tex"
        with open(outfilename, "w") as outfile:
            outfile.write(self.latex_str())
        os.system(f"pdflatex {outfilename}")

    def get_alignment(self):
        i = len(self.Q) - 1
        j = len(self.T) - 1
        min_v = self.matrix[len(self.Q)][len(self.T)]
        min_j = len(self.T) - 1

        res_q = ""
        res_align = ""
        res_t = ""
        while j > 0:
            # print(f"j {j} min_j {min_j} v {self.matrix[len(self.Q)][j+1]} min_v {min_v}")
            if self.matrix[len(self.Q)][j + 1] < min_v:
                min_v = self.matrix[len(self.Q)][j + 1]
                min_j = j
            j -= 1

        for last in range(len(self.T) - 1, min_j, -1):
            res_t = self.T[last] + res_t
            res_q = " " + res_q
            res_align = "*" + res_align

        j = min_j
        while i >= 0 and j >= 0:
            cf = self.come_from(i + 1, j + 1)
            if cf[1] == "1":  # diagonal
                res_q = self.Q[i] + res_q
                res_t = self.T[j] + res_t
                if self.Q[i] == self.T[j]:
                    res_align = "|" + res_align
                else:
                    res_align = "." + res_align
                i -= 1
                j -= 1
                continue
            if (
                cf[0] == "1"
            ):  # top, a letter of the query aligned against nothing (on T)
                res_q = self.Q[i] + res_q
                res_t = "-" + res_t
                res_align = " " + res_align
                i -= 1
                continue
            if (
                cf[2] == "1"
            ):  # left, a letter of the text aligned against nothing (on Q)
                res_q = "-" + res_q
                res_t = self.T[j] + res_t
                res_align = " " + res_align
                j -= 1
                continue
            assert True == False

        while i >= 0:
            res_q = self.Q[i] + res_q
            res_t = "-" + res_t
            res_align = " " + res_align
            i -= 1
        while j >= 0:
            res_q = " " + res_q
            res_t = self.T[j] + res_t
            res_align = "*" + res_align
            j -= 1

        res_q = "Q = " + res_q
        res_align = "    " + res_align
        res_t = "T = " + res_t
        return res_q, res_align, res_t

    def __str__(self):
        """ string the matrix in markdown format"""

        def color(is_green):
            if is_green:
                return "\033[91m"
            else:
                return "\033[92m"

        res = ""
        width = 4
        vide = " "
        inf = "inf"
        line = f"{vide:>{2*width}}"
        is_green = True

        # First line (T)
        for j in range(0, len(self.T)):
            line += f"{self.T[j]:>{width}}"
        res += line + "\n"

        # Qecond line (Os)
        line = f"{vide:>{width}}"
        for j in range(0, len(self.T) + 1):
            line += f"{self.matrix[0][j]:>{width}}"

        # line, normal_color = change_color(normal_color, line)
        res += line + "\n"

        # all other lines
        for i in range(1, len(self.Q) + 1):
            # line = f"\033[00m{self.Q[i-1]:>{width}}"
            line = f"\033[00m{self.Q[i-1]:>{width}}"
            for j in range(0, len(self.T) + 1):
                if self.matrix[i][j] == sys.maxsize:
                    line += f"{inf:>{width}}"
                else:
                    line += f"{color(is_green)}{self.matrix[i][j]:>{width}}"
                if j > 0 and j != len(self.T) and j in self.end_vertical_blocs:
                    is_green = not is_green

            if (
                len(self.end_vertical_blocs) % 2 == 1
            ):  # this line changed the value of the first next block color.
                is_green = not is_green
            if i in self.end_horizontal_blocs:
                is_green = not is_green
            res += line + "\n"

        # res += f"\n horizontal block frontiers: {self.end_horizontal_blocs}"
        # res += f"\n vertical block frontiers: {self.end_vertical_blocs}"
        # for i in range(len(self.end_horizontal_blocs)):
        res += "\033[00m\n"
        res_q, res_align, res_t = self.get_alignment()
        res += res_t + "\n" + res_align + "\n" + res_q + "\n"
        #     res += f" {i}"
        # res += "\n vertical block frontiers:"
        # for i in range(len(self.end_vertical_blocs)):
        #     res += f" {i}"
        return res + "\033[00m\n"

    #### DTW ####
    def dist(self, alpha, beta):
        if alpha == beta:
            return 0
        return self.mismatch

    def init_local_DTW(self):
        """ initializes first line and first columns for a global alignment"""
        for i in range(0, len(self.T) + 1):
            self.matrix[0][i] = 0
        for j in range(1, len(self.Q) + 1):
            self.matrix[j][0] = sys.maxsize

    def fill_DTW(self):
        """ fills the matrix for global alignment (Needleman & Wunsch algo)"""
        for i in range(1, len(self.Q) + 1):
            # i-th line
            for j in range(1, len(self.T) + 1):
                # j-th column
                self.matrix[i][j] = (
                    min(
                        self.matrix[i - 1][j - 1],
                        self.matrix[i][j - 1],
                        self.matrix[i - 1][j],
                    )
                    + self.dist(self.Q[i - 1], self.T[j - 1])
                )
        return self.matrix[len(self.Q)][len(self.T)]

    # notation de pierre est bizarre !
    def check_block_property(self):
        # Check Kuszmaul property
        for bi in range(1, len(self.end_horizontal_blocs)):
            for bj in range(1, len(self.end_vertical_blocs)):
                # positions of the end of the current block
                i = self.end_horizontal_blocs[bi]
                j = self.end_vertical_blocs[bj]
                l_x = self.end_horizontal_blocs[bi] - self.end_horizontal_blocs[bi - 1]
                l_y = self.end_vertical_blocs[bj] - self.end_vertical_blocs[bj - 1]

                if self.Q[i - 1] == self.T[j - 1]:
                    d = 0
                else:
                    d = self.mismatch

                for o_x in range(1, l_x):
                    if o_x < l_y:
                        diag = self.matrix[i - l_x + 1][j - o_x]
                        cost = o_x * d
                    else:
                        diag = self.matrix[i - l_x + 1 + o_x - l_y + 1][j - l_y + 1]
                        cost = (l_y - 1) * d
                    value = self.matrix[i - l_x + o_x + 1][j]
                    top = self.matrix[i - l_x + o_x][j]
                    theoretical_value = min(top + d, diag + cost)
                    if theoretical_value != value:
                        print(
                            f"blocks {bi} {bj}, ending at pos {i} {j}, of length {l_x} {l_y}"
                        )
                        print(
                            f"o_x: position {i-l_x+o_x+1} {j} "
                            + f" value {value}"
                            + f" top {top}"
                            + f" diag {diag} diag+cost: {diag+cost}"
                        )
                        return False
                for o_y in range(1, l_y):
                    if o_y < l_x:
                        diag = self.matrix[i - o_y][j - l_y + 1]
                        cost = o_y * d
                    else:
                        diag = self.matrix[i - l_x + 1][j - l_y + 1 + o_y - l_x + 1]
                        cost = (l_x - 1) * d
                    value = self.matrix[i][j - l_y + o_y + 1]
                    left = self.matrix[i][j - l_y + o_y]
                    theoretical_value = min(left + d, diag + cost)
                    if theoretical_value != value:
                        print(
                            f"blocks {bi} {bj}, ending at pos {i} {j}, of length {l_x} {l_y}"
                        )
                        print(
                            f"o_y: position {i} {j-l_y+o_y+1}"
                            + f" value {value}"
                            + f" left {left}"
                            + f" diag {diag} diag+cost: {diag+cost}"
                        )
                        return False
        return True


def generate_random_string(size, alphabet_size):
    res = ""
    i = 0
    while True:
        cur_chr = chr(random.randint(65, 65 + alphabet_size - 1))
        run_len = random.randint(1, 5)
        if run_len + i > size:
            run_len = size - i
        res += cur_chr * run_len
        i += run_len
        if i >= size:
            break

    return res


def search_counter_example_random_strings():
    max_s = 100
    alphabet_size = 6
    true_property = True
    nb_test = 0
    ldtw = None
    while true_property and nb_test < alphabet_size ** max_s:
        nb_test += 1
        len_T = random.randint(2, max_s)
        len_Q = random.randint(len_T, max_s)
        T = generate_random_string(len_T, alphabet_size)
        Q = generate_random_string(len_Q, alphabet_size)
        ldtw = LocalDTW(Q, T, 1)
        true_property = ldtw.check_block_property()
        if nb_test % 100 == 0:
            print(f"{nb_test} have passed, still no counter example!")
    if nb_test < alphabet_size ** max_s:
        print("A counter example was found!")
        print(ldtw)
        ldtw.output_latex()
        exit()
    print("All test passed succesfuly!")
    exit()


def main():
    T = sys.argv[1]
    Q = sys.argv[2]
    ldtw = LocalDTW(Q, T, 1)
    print(f"{ldtw}")
    search_counter_example_random_strings()
    exit()


if __name__ == "__main__":
    sys.exit(main())
