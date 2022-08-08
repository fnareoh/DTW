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
        """initialize and computes all values of a block that are strictly bellow max_value.

        Args:
            l (int): length of the block.
            h (int): height of the block.
            cost (int): cost of matching the letters (0 if equal).
            Vnw (int): value on top left of the block pos + (-1, -1).
            q_North (list<int,int>): sorted list of values in the bottom row of the block above together with the first position of the row this value appears.
            q_West (list<int,int>): sorted list of values in the right column of the block to the left together with the first position of the column this value appears.
            max_value (int): optional : maximal value computed in the matrix. Any value higher than max_value is not computed and represented by max_value.
            coment (str): optional : optional coment to give when printing the block.
        """
        self.width = width
        self.height = height
        self.cost = cost
        self.max_value = max_value
        self.coment = coment

        Vnw = min(Vnw, q_North[0][0], q_West[0][0])
        self.q_top = self.__compute_adjacent_q__(cost, width, Vnw, q_North)
        self.q_left = self.__compute_adjacent_q__(cost, height, Vnw, q_West)
        assert (
            self.q_top[-1][1] < self.width
        ), f"top row overwlows: {self.q_top[-1][1]} < {self.width} is False"
        assert (
            self.q_left[-1][1] < self.height
        ), f"left column overwlows: {self.q_left[-1][1]} < {self.height} is False"
        assert (
            self.q_top[0][0] == self.q_left[0][0]
        ), f"incoherent value at top left: {self.q_top[0][0]} != {self.q_left[0][0]}"
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
                if (
                    self.max_value is not None
                    and Vnw + (j + 1) * cost >= self.max_value
                ):
                    q.append((self.max_value, j))
                    return q
                q.append((Vnw + (j + 1) * cost, j))
                j += 1
            if j == 0 or (j <= last_pos and j < width):
                if self.max_value is not None and val + cost >= self.max_value:
                    q.append((self.max_value, j))
                    return q
                q.append((val + cost, j))
        return q

    def __transfer_triangle__(self, width, q_top, q_right):
        for i, (val, pos) in reversed(list(enumerate(q_top))):
            last_pos = width
            if i + 1 < len(q_top):
                last_pos = q_top[i + 1][1]
            for k in range(width - last_pos, width - pos):
                if k >= min(self.height, self.width):
                    return
                if len(q_right) == 0 or q_right[-1][0] < val + self.cost * k:
                    if (
                        self.max_value is not None
                        and val + self.cost * k >= self.max_value
                    ):
                        q_right.append((self.max_value, k))
                        return
                    q_right.append((val + self.cost * k, k))
        return

    def __transfer_parralel__(self, width, q_left, q_right):
        for i, (val, pos) in enumerate(q_left[1:]):
            if pos + width - 1 >= max(self.height, self.width):
                return

            nei = q_right[-1][0]
            if val + self.cost * (width - 1) > nei:
                if (
                    self.max_value is not None
                    and val + self.cost * (width - 1) >= self.max_value
                ):
                    q_right.append((self.max_value, pos + width - 1))
                    return
                q_right.append((val + self.cost * (width - 1), pos + width - 1))
        return

    def __compute_bottom_right__(self):
        """From q_top and q_left height and width, deduce q_right and q_bottom"""
        self.q_right, self.q_bottom = [], []
        if self.height == 1:
            self.q_bottom = self.q_top.copy()
            self.q_right = [(self.q_bottom[-1][0], 0)]
        elif self.width == 1:
            self.q_right = self.q_left.copy()
            self.q_bottom = [(self.q_right[-1][0], 0)]
        else:
            # Transfer from left to bottom
            self.__transfer_triangle__(self.height, self.q_left, self.q_bottom)
            # Transfer from top to right
            self.__transfer_triangle__(self.width, self.q_top, self.q_right)

            # Tranfer parralel
            if self.height > self.width:
                self.__transfer_parralel__(self.width, self.q_left, self.q_right)
            elif self.height < self.width:
                self.__transfer_parralel__(self.height, self.q_top, self.q_bottom)

        assert self.q_bottom[-1][1] < self.width
        assert self.q_right[-1][1] < self.height
        assert self.q_bottom[-1][0] == self.q_right[-1][0]
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

        if self.height >= 2:
            res += (
                " " + " ".join([str(e) for e in unpack(self.q_top, self.width)]) + "\n"
            )
        left = unpack(self.q_left, self.height)
        right = unpack(self.q_right, self.height)
        for i in range(1, self.height - 1):
            if self.width >= 2:
                res += f" {left[i]}" + "  " * (self.width - 2) + f" {right[i]}\n"
            else:
                res += f" {right[i]}\n"
        res += " " + " ".join([str(e) for e in unpack(self.q_bottom, self.width)])

        return res


def main():
    # Manual test 1
    height = 6
    width = 5
    Vnw = 3
    q_North = [(2, 0), (3, 2), (4, 3)]
    q_West = [(5, 0), (6, 4), (7, 5)]
    b = BorderBlock(height, width, 1, Vnw, q_North, q_West)
    print(b)

    b = BorderBlock(2, width, 1, 0, [(0, 0)], [(5, 0)])
    print(b)

    b = BorderBlock(2, 1, 1, 0, [(0, 0)], [(5, 0)])
    print(b)

    b = BorderBlock(1, 1, 1, 0, [(0, 0)], [(5, 0)])
    print(b)

    b = BorderBlock(height, width, 1, Vnw, q_North, q_West, max_value=5)
    print(b)

    exit()


if __name__ == "__main__":
    main()
