"""
This is a dynamic programming computation of pattern matching for DTW.
"""

__author__ = "Garance Gourdel, Pierre Peterlongo"
__email__ = "pierre.peterlongo@inria.fr, garance.gourdel@inria.fr"

from json.tool import main
import sys


class Constant:
    match = 0
    mismatch = 1
    gap = 1


class PM_Matrix:
    """
    Stores a matrix |P|x|T| (|P|+1 lines and |T|+1columns),
    strings P and T and the score system (match, mismatch, gap)
    defines some local alignment functions
    """

    def __init__(self, P, T):
        """ defines and stores initial values"""

        consta = Constant()
        self.P = P
        self.T = T
        self.mismatch = consta.mismatch
        self.match = consta.match
        self.gap = consta.gap
        #      T
        #  --------
        #  |
        # P |
        #  |

        self.matrix = [[] for i in range(len(P) + 1)]
        for i in range(len(P) + 1):
            self.matrix[i] = [0 for j in range(len(T) + 1)]

        # initializes first line and first columns for pattern matching
        for j in range(0, len(self.T) + 1):
            self.matrix[0][j] = 0
        for i in range(1, len(self.P) + 1):
            self.matrix[i][0] = sys.maxsize

        # Determine block positions:
        # list of vertical frontiers:
        self.end_vertical_blocs = []
        self.end_vertical_blocs.append(
            0
        )  # first column is considered as a special block
        for i in range(1, len(self.T)):
            if self.T[i] != self.T[i - 1]:
                self.end_vertical_blocs.append(i)
        self.end_vertical_blocs.append(len(self.T))

        # list of horizontal frontiers:
        self.end_horizontal_blocs = []
        self.end_horizontal_blocs.append(
            0
        )  # first line is considered as a special block
        for i in range(1, len(self.P)):
            if self.P[i] != self.P[i - 1]:
                self.end_horizontal_blocs.append(i)
        self.end_horizontal_blocs.append(len(self.P))
        self.fill()

    def fill():
        return self.matrix[len(self.P)][len(self.T)]

    def __str__(self):
        """ string the matrix"""

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

        # Second line (Os)
        line = f"{vide:>{width}}"
        for j in range(0, len(self.T) + 1):
            line += f"{self.matrix[0][j]:>{width}}"

        # line, normal_color = change_color(normal_color, line)
        res += line + "\n"

        # all other lines
        for i in range(1, len(self.P) + 1):
            # line = f"\033[00m{self.P[i-1]:>{width}}"
            line = f"\033[00m{self.P[i-1]:>{width}}"
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

        return res + "\033[00m\n"

    def get_last_value(self):
        return self.matrix[-1][-1]

    def dist(self, alpha, beta):
        if alpha == beta:
            return +self.match
        return self.mismatch

    def from_diag(self, i, j):
        return (
            self.matrix[i - 1][j - 1] + self.dist(self.P[i - 1], self.T[j - 1])
            == self.matrix[i][j]
        )

    def from_left(self, i, j):
        return False

    def from_top(self, i, j):
        return False

    def trace_back(self, pos):
        assert pos < len(self.matrix[0])
        i = len(self.P)
        j = pos
        alP = ""
        alT = ""
        nb_matchs = 0
        while i > 0 and j > 0:
            # test first the diagonal, in order to favor diag in case of ex-aequos
            if self.from_diag(i, j):
                # diag
                alP = self.P[i - 1] + alP
                alT = self.T[j - 1] + alT
                if self.P[i - 1] == self.T[j - 1]:
                    nb_matchs += 1
                i -= 1
                j -= 1
            elif self.from_left(i, j):
                # left
                alP = "-" + alP
                alT = self.T[j - 1] + alT
                j -= 1
            elif self.from_top(i, j):
                # up
                alT = "-" + alT
                alP = self.P[i - 1] + alP
                i -= 1
            else:
                raise ValueError(
                    "This cell does not come from either the top, left or diagonal."
                )
        # here either i=0 or j=0
        assert i == 0
        return j, alP, alT

    def min_last_row_val_index(self):
        m = min(self.matrix[len(self.P)])
        return m, self.matrix[len(self.P)].index(m)

    def index_min_last_row(self):
        return self.matrix[len(self.P)].index(min(self.matrix[len(self.P)]))

    def k_smallest(self, k):
        return sorted(
            range(len(self.matrix[len(self.P)])),
            key=lambda sub: self.matrix[len(self.P)][sub],
        )[:k]

    def compute_origin_min_position(self):
        return self.trace_back(self.index_min_last_row())


class PM_DTW(PM_Matrix):
    #### DTW ####
    def from_left(self, i, j):
        return (
            self.matrix[i][j - 1] + self.dist(self.P[i - 1], self.T[j - 1])
            == self.matrix[i][j]
        )

    def from_top(self, i, j):
        return (
            self.matrix[i - 1][j] + self.dist(self.P[i - 1], self.T[j - 1])
            == self.matrix[i][j]
        )

    def fill(self):
        """ fills the matrix for global alignment (Needleman & Wunsch algo)"""
        for i in range(1, len(self.P) + 1):
            # i-th line
            for j in range(1, len(self.T) + 1):
                # j-th column
                self.matrix[i][j] = (
                    min(
                        self.matrix[i - 1][j - 1],
                        self.matrix[i][j - 1],
                        self.matrix[i - 1][j],
                    )
                    + self.dist(self.P[i - 1], self.T[j - 1])
                )
        return self.matrix[len(self.P)][len(self.T)]


class PM_ED(PM_Matrix):
    #### Edit distance (Needleman & Wunsch) ####
    def from_left(self, i, j):
        return self.matrix[i][j - 1] + self.gap == self.matrix[i][j]

    def from_top(self, i, j):
        return self.matrix[i - 1][j] + self.gap == self.matrix[i][j]

    def fill(self):
        """ fills the matrix for global alignment (Needleman & Wunsch algo)"""
        for i in range(1, len(self.P) + 1):
            # i-th line
            for j in range(1, len(self.T) + 1):
                # j-th column
                self.matrix[i][j] = min(
                    self.matrix[i - 1][j - 1] + self.dist(self.P[i - 1], self.T[j - 1]),
                    self.matrix[i][j - 1] + self.gap,
                    self.matrix[i - 1][j] + self.gap,
                )
        return self.matrix[len(self.P)][len(self.T)]


def main(P, T):
    ldtw = PM_DTW(P, T)
    ldtw.fill()
    print("******** Dynamic Time Warp **********")
    min_dtw, pos_dtw = ldtw.min_last_row_val_index()
    origin, alP, alT = ldtw.trace_back(pos_dtw)
    print(
        f"Alignement (cost {min_dtw}) of P ending at {pos_dtw} in T starts at {origin} in T."
    )
    print(alP)
    print(alT)

    led = PM_ED(P, T)
    led.fill()
    print("******** Edit distance **********")
    min_ed, pos_ed = led.min_last_row_val_index()
    origin, alP, alT = led.trace_back(pos_ed)
    print(
        f"Alignement (cost {min_ed}) of P ending at {pos_ed} in T starts at {origin} in T."
    )
    print(alP)
    print(alT)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
