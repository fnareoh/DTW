"""Computes a matrix of DTW distances between a pattern Q and a text T
    In the matrix a block is defined by a a^height x b^width letters to be compared
    Only distances <= k must be computed
    Computations of a block is in O(height + width).
"""

__author__ = "Garance Gourdel, Pierre Peterlongo"
__email__ = "pierre.peterlongo@inria.fr, garance.gourdel@inria.fr"

import sys  # for maxsize


class Block:
    """ Class for a single block of the DTW matrix """

    def __init__(
        self,
        height: int,
        width: int,
        equals: bool,
        Vnw: int,
        Vw: int,
        Vn: int,
        h_cuts,
        v_cuts,
        max_value=sys.maxsize,
        line_start=None,
        line_end=None,
        column_start=None,
        column_end=None,
    ):
        """initialize and computes a block

        Args:
            l (int): length of the block
            h (int): height of the block
            equals (bool): are letters equal
            Vnw (int): value on top left of the block (pos -1, -1)
            Vn (int): value on top of the block (pos -1, 0)
            Vw (int): value on left of the block (pos 0, -1)
            h_cuts (list<int>): positions where values are increasing on the line on top of the block (outside the block)
            v_cuts (list<int>): positions where values are increasing on the column on the left of the block (outside the block)
            max_value (int): optional : maximal value computed in the matrix. Any value higher than max_value is not computed
            line_start (int): optional : indicates the first line of the block in a whole matrix
            line_end (int): optional : indicates the last line of the block in a whole matrix
            column_start (int): optional : indicates the first column of the block in a whole matrix
            column_end (int): optional : indicates the last column of the block in a whole matrix
        """

        self.width = width
        self.height = height
        self.equals = equals
        self.Vnw = Vnw
        self.Vw = Vw
        self.Vn = Vn
        self.h_cuts = h_cuts  # check that is this is pointer only not a deep copy (else compute directly self.top_cuts here)
        self.v_cuts = v_cuts  # check that is this is pointer only not a deep copy (else compute directly self.leftmost_cuts here)
        self.line_start = line_start
        self.line_end = line_end
        self.column_start = column_start
        self.column_end = column_end
        self.max_value = max_value

        self.bottom_cuts = []
        self.top_cuts = []
        self.leftmost_cuts = []
        self.rightmost_cuts = []
        self.__compute_border__()

    def __repr__(self):
        """Returns a string repr of a block

        Returns:
            [str]: string repr of a block
        """
        res = ""
        if self.line_start != None:
            res += f" lins [{self.line_start},{self.line_end}]\n"
            res += f" cols [{self.column_start},{self.column_end}]\n"
        if self.max_value != sys.maxsize:
            res += f" max value {self.max_value}\n"
        else:
            res += " Not bounded by a maximal value\n"
        res += f" height = {self.height}, width = {self.width}\n"
        res += f" external NW = {self.Vnw}, N = {self.Vn}, W = {self.Vw}\n"
        res += f" external vertical cuts = {self.v_cuts}\n"
        res += f" external horizontal cuts = {self.h_cuts}\n\n"
        res += f" internal left cuts = {self.leftmost_cuts}\n"
        res += f" internal top cuts = {self.top_cuts}\n"
        res += f" internal right cuts = {self.rightmost_cuts}\n"
        res += f" internal bottom cuts = {self.bottom_cuts}\n"

        res += f" top left = {self.tl}\n"
        res += f" top right = {self.tr}\n"
        res += f" bottom left = {self.bl}\n"
        res += f" bottom right = {self.br}\n"

        return res

    def __compute_bottom_cuts__(self):
        """From top and left cuts, computes bottom cuts

        do not exceed max_value
        Returns:
            int: the maximal value obtained in the last line
        """

        if self.height == 1:
            self.bottom_cuts = self.top_cuts
            return self.tr
        current_value = self.bl
        if current_value == self.max_value:
            return current_value
        if self.height < self.width:
            # h-1 positions for the bottom of the first square - reversed order to keep the cut array sorted
            for x in range(self.height - 2, -1, -1):
                if x not in self.leftmost_cuts:
                    self.bottom_cuts.append(self.height - 2 - x)
                    current_value += 1
                    if current_value == self.max_value:
                        return current_value

            # diagonal values
            for y in range(self.width - self.height):
                if y in self.top_cuts:
                    self.bottom_cuts.append(y + self.height - 1)
                    current_value += 1
                    if current_value == self.max_value:
                        return current_value
        else:
            # l-1 positions for the bottom of the first square
            for x in range(self.width - 1):
                if self.height - x - 2 not in self.leftmost_cuts:
                    self.bottom_cuts.append(x)
                    current_value += 1
                    if current_value == self.max_value:
                        return current_value
        return current_value

    def __compute_right_cuts__(self):
        """From top and left cuts, computes right cuts

        do not exceed max_value
        Returns:
            int: the maximal value obtained in the last column
        """
        if self.width == 1:
            self.rightmost_cuts = self.leftmost_cuts
            return self.br

        current_value = self.tr
        if current_value == self.max_value:
            return current_value
        if self.height < self.width:
            # h-1 positions for the right of the last square
            for y in range(self.height - 1):
                if self.width - y - 2 not in self.top_cuts:
                    self.rightmost_cuts.append(y)
                    current_value += 1
                    if current_value == self.max_value:
                        return current_value
        else:
            # l-1 positions for the right of the first square - reversed order to keep the cut array sorted
            for y in range(self.width - 2, -1, -1):
                if y not in self.top_cuts:
                    self.rightmost_cuts.append(self.width - 2 - y)
                    current_value += 1
                    if current_value == self.max_value:
                        return current_value

            # diagonal values
            for x in range(self.height - self.width):
                if x in self.leftmost_cuts:
                    self.rightmost_cuts.append(x + self.width - 1)
                    current_value += 1
                    if current_value == self.max_value:
                        return current_value

        return current_value

    def __compute_left_cuts__(self):
        """Computes internal left cut and bottom left value

            From left cuts and internal Vnw, Vn, and Vw value: create the
            internal rightmost_cuts

        Returns:
            int: the maximal value obtained in the last cell of the first column
        """

        # shortcut in case of a tiny block
        if self.height == 1:
            return self.tl

        current_value = self.tl

        first_impossible_cut = self.height - 1
        last_used_cut = -1  # short cut

        if current_value == self.max_value:
            return current_value

        # special case, the North West value is higher than the North value, by delta
        # This generates delta cuts in the start of the first column.
        # For all following cuts: they come after these first delta cuts
        delta_NW_N = max(0, self.Vnw - self.Vn)
        delta_W_N = max(0, self.Vw - self.Vn)
        delta = min(delta_NW_N, delta_W_N)
        for i in range(min(delta, first_impossible_cut)):
            self.leftmost_cuts.append(i)
            last_used_cut = i
            current_value += 1
            if current_value == self.max_value:
                return current_value

        # special case, corresponds to a cut in -1 position in the external cuts.
        # needs to be performed externally (on should simply add a -1 to external cuts, but this would imply to sort the list)
        if self.Vw == self.Vnw + 1 or self.Vw == sys.maxsize:
            potential_new_cut = max(0, last_used_cut + 1)
            if potential_new_cut < first_impossible_cut:
                self.leftmost_cuts.append(potential_new_cut)
                last_used_cut = potential_new_cut
                current_value += 1
                if current_value == self.max_value:
                    return current_value

        # general case: an external cut position i generates an internal cut position i+1
        # add the delta shift.
        for v_cut in self.v_cuts:
            potential_new_cut = max(v_cut + 1, last_used_cut + 1)
            if potential_new_cut < first_impossible_cut:
                self.leftmost_cuts.append(potential_new_cut)
                last_used_cut = potential_new_cut
                current_value += 1
                if current_value == self.max_value:
                    return current_value

        return current_value

    def __compute_top_cuts__(self):
        """Computes internal top cut and top right value

            From left cuts and internal Vnw, Vn, and Vw value: create the
            internal top_cuts

        Returns:
            int: the maximal value obtained in the last cell of the first row
        """
        # shortcut in case of a tiny block
        if self.width == 1:
            return self.tl

        current_value = self.tl
        first_impossible_cut = self.width - 1
        last_used_cut = -1  # short cut

        if current_value == self.max_value:
            return current_value

        # special case, the North West value is higher than the West value, by delta
        # This generates delta cuts in the start of the first row.
        # For all following cuts: they come after these first delta cuts
        delta_NW_W = max(0, self.Vnw - self.Vw)
        delta_N_W = max(0, self.Vn - self.Vw)
        delta = min(delta_NW_W, delta_N_W)
        for i in range(min(delta, first_impossible_cut)):
            self.top_cuts.append(i)
            last_used_cut = i
            current_value += 1
            if current_value == self.max_value:
                return current_value

        # special case, corresponds to a cut in -1 position in the external cuts
        # simply add a -1 to the self.leftmost_cuts
        if self.Vn == self.Vnw + 1:
            potential_new_cut = max(0, last_used_cut + 1)
            if potential_new_cut < first_impossible_cut:
                self.top_cuts.append(potential_new_cut)
                last_used_cut = potential_new_cut
                current_value += 1
                if current_value == self.max_value:
                    return current_value

        # general case: an external cut position i generates an internal cut position i+1
        # add the delta shift.
        for h_cut in self.h_cuts:
            potential_new_cut = max(h_cut + 1, last_used_cut + 1)
            if potential_new_cut < first_impossible_cut:
                self.top_cuts.append(potential_new_cut)
                last_used_cut = potential_new_cut
                current_value += 1
                if current_value == self.max_value:
                    return current_value

        return current_value

    def __compute_border__(self):
        """Computes the border of the block
        cf https://notability.com/n/M9iqdy50C4dnBZoPSJDaT

        tl (int): top left value (pos 0, 0)
        tr (int): top right value (pos 0, l-1)
        bl (int): bottom left value (pos h-1, 0)
        br (int): bottom right value (pos h-1, l-1)

        bottom_cuts (list<int>): positions were values are increasing on the last line of the block (inside the block)
        rightmost_cuts (list<int>): positions were values are increasing on the last column of the block (inside the block)
        top_cuts (list<int>): positions were values are increasing on the first line of the block (inside the block)
        leftmost_cuts (list<int>): positions were values are increasing on the first column of the block (inside the block)
        """

        if self.equals:
            self.tl = min(self.Vnw, self.Vw, self.Vn)
            self.tr = self.tl
            self.bl = self.tl
            self.br = self.tl
            return

        self.tl = min(self.Vnw, self.Vw, self.Vn) + 1
        if self.tl >= self.max_value:
            self.tl = self.max_value
            self.tr = self.tl
            self.bl = self.tl
            self.br = self.tl
            return

        self.bl = self.__compute_left_cuts__()
        self.tr = self.__compute_top_cuts__()
        self.br = self.__compute_bottom_cuts__()
        br_again = self.__compute_right_cuts__()
        assert (
            br_again == self.br
        ), f"Internal error : {br_again} != {self.br} \n {self}"


def random_tests(nb_tests: int):

    # brut force tests:
    print("Lets make some tests")
    import random

    for _ in range(nb_tests):
        height = random.randint(4, 500)
        width = random.randint(4, 500)

        Vnw = random.randint(0, 10)
        Vw = random.randint(0, Vnw + 1)
        Vn = random.randint(0, Vnw + 1)
        h_cuts = []
        for x in range(height - 1):
            if bool(random.getrandbits(1)):
                h_cuts.append(x)
        v_cuts = []
        for y in range(width - 1):
            if bool(random.getrandbits(1)):
                v_cuts.append(y)
        Block(height, width, False, Vnw, Vw, Vn, h_cuts, v_cuts, max_value=10)
        Block(height, width, False, Vnw, Vw, Vn, h_cuts, v_cuts)

    print(f" Made {nb_tests} tests - every thing is fine")


def main():
    # Manual test
    height = 6
    width = 5
    Vnw = 3
    Vn = 0
    Vw = 0
    v_cuts = {}
    h_cuts = {}
    b = Block(height, width, False, Vnw, Vw, Vn, h_cuts, v_cuts)
    print(b)
    exit()
    random_tests(10000)


if __name__ == "__main__":
    main()
