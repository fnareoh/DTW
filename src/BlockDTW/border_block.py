"""
    Defines a class BorderBlock which compactly represent the values of
    each border cells in a block. Those values are either initialize to a
    specific values or computes from the representation of a neighbouring block.
"""

__author__ = "Garance Gourdel, Pierre Peterlongo"
__email__ = "pierre.peterlongo@inria.fr, garance.gourdel@inria.fr"

import sys  # for maxsize


def unpack(l, width):
    """Returns the uncompressed representation of the values in the list as a list"""
    unpacked = []
    for i, (val, pos) in enumerate(l):
        last_pos = width
        if i + 1 < len(l):
            last_pos = l[i + 1][1]
        unpacked += [val] * (last_pos - pos)
    return unpacked


class BorderBlock:
    """Class for a compact representation of a single block of the DTW matrix"""

    def __init__(
        self,
        height: int,
        width: int,
        cost: int,
        Vnw: int,
        q_North,
        q_West,
        max_value=sys.maxsize,
        coment=None,
    ):
        """initialize and computes a block

        Args:
            l (int): length of the block.
            h (int): height of the block.
            cost (int): cost of matching the letters (0 if equal).
            Vnw (int): value on top left of the block pos + (-1, -1).
            q_North (list<int,int>): sorted list of values in the bottom row of the block above together with the first position of the row this value appears.
            q_West (list<int,int>): sorted list of values in the right column of the block to the left together with the first position of the column this value appears.
            max_value (int): optional : maximal value computed in the matrix. Any value higher than max_value is not computed.
            coment (str): optional : optional coment to give when printing the block.
        """
        self.width = width
        self.height = height
        self.cost = cost
        Vnw = min(Vnw, q_North[0][0], q_West[0][0])
        self.q_top = self.__compute_adjacent_q__(cost, width, Vnw, q_North)
        self.q_left = self.__compute_adjacent_q__(cost, height, Vnw, q_West)
        self.max_value = max_value
        self.coment = coment

        self.__compute_bottom_right__()

    def __compute_adjacent_q__(self, cost, width, Vnw, q_North):
        """From q_North, Vnw, the cost and the width returns q_top.
        From q_west , Vnw, the cost and the height returns q_left."""
        q = []
        for i, (val, pos) in enumerate(q_North):
            last_pos = width
            if i + 1 < len(q_North):
                last_pos = q_North[i + 1][1]

            j = pos
            if j != 0:
                j += 1
            while Vnw + (j + 1) * cost < val + cost and j <= last_pos and j < width:
                q.append((Vnw + (j + 1) * cost, j))
                j += 1
            if j == 0:
                q.append((val + cost, j))
            elif j <= last_pos and j < width:
                q.append((val + cost, j))

        return q

    def __compute_bottom_right__(self):
        """From q_top and q_left height and width, deduce q_right and q_bottom"""
        self.q_bottom = [(self.q_left[-1][0], 0)]
        self.q_right = [(self.q_top[-1][0], 0)]

        # Transfer from left to bottom
        for i, (val, pos) in reversed(list(enumerate(self.q_left[:-1]))):
            nei = self.q_bottom[-1][0]
            d = self.height - pos - 1
            if val + self.cost * d == nei:
                pass
            elif val + self.cost * d > nei:
                self.q_bottom.append((val + self.cost * d, d))
                # BUG: we miss some values in trapezoid shape

        # Transfer from top to right
        # TODO: factorize code
        for i, (val, pos) in reversed(list(enumerate(self.q_top[:-1]))):
            nei = self.q_right[-1][0]
            d = self.width - pos - 1
            print(d, val + self.cost * d)
            if val + self.cost * d > nei:
                print((val - nei) // self.cost + d, d)
                for k in range((val - nei) // self.cost + d + 1, d + 1):
                    self.q_right.append((val + self.cost * k, k))

        # TODO: coppy left to right and top to bottom

        print(self.q_bottom)
        print(self.q_top)
        print(self.q_right)
        # assert self.q_bottom[-1][0] == self.q_left[-1][0]
        self.Vse = self.q_bottom[-1][0]

    def __repr_values__(self):
        """Returns a string containing the value of each of its attribute

        Returns:
            [str]: string repr of a block
        """
        res = ""
        if self.coment != None:
            res += f" comment: [{self.coment}]\n"
        if self.max_value != sys.maxsize:
            res += f" max value {self.max_value}\n"
        else:
            res += " Not bounded by a maximal value\n"
        res += f" height = {self.height}, width = {self.width}, cost = {self.cost}\n"
        res += f" q_top = {self.q_top}\n"
        res += f" q_left = {self.q_left}\n"
        res += f" q_bottom = {self.q_bottom}\n"
        res += f" q_right = {self.q_right}\n"
        res += f" Vse (value in bottom right) = {self.Vse}\n"

        return res

    def __repr__(self):
        """Returns a string containing the value of each of its attribute

        Returns:
            [str]: string repr of a block
        """
        res = ""
        if self.coment != None:
            res += f" comment: [{self.coment}]\n"
        if self.max_value != sys.maxsize:
            res += f" max value {self.max_value}\n"
        else:
            res += " Not bounded by a maximal value\n"

        res += " " + " ".join([str(e) for e in unpack(self.q_top, self.width)]) + "\n"
        left = unpack(self.q_left, self.height)
        right = unpack(self.q_right, self.height)
        for i in range(1, self.height - 1):
            res += f" {left[i]}" + "  " * (self.width - 2) + f" {right[i]}\n"
        res += (
            " " + " ".join([str(e) for e in unpack(self.q_bottom, self.width)]) + "\n"
        )

        return res


def main():
    # Manual test
    height = 6
    width = 5
    Vnw = 3
    q_North = [(2, 0), (3, 2), (4, 3)]
    q_West = [(5, 0), (6, 4), (7, 5)]
    b = BorderBlock(height, width, 1, Vnw, q_North, q_West)
    print(b)
    exit()


if __name__ == "__main__":
    main()
