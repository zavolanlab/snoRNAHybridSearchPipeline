"""
Module with structures
"""

import re


class RNASnoopStructure(object):

    """Representation of a structure from RNASnoop output"""

    def __init__(self, structure):
        """Initialize RNASnoopStructure

        :param structure: @todo

        """
        self.structure = structure
        self.__parse_U_gap()
        self.__parse_gap_right()
        self.__parse_i_t_gap()
        self.__parse_t_i_gap()
        self.__parse_i_b_gap()
        self.__parse_stem_length()

    def __parse_U_gap(self):
        match = re.search(r"<\.\|\.<", self.structure)
        if match:
            self.U_gap = 3
        else:
            self.U_gap = 2

    def __parse_gap_right(self):
        match = re.search(r"\|\.([^\&]+)", self.structure)
        tmp = re.findall(r"<(\.)<", match.group(1))
        if tmp:
            self.gap_right = len([m for m in tmp])
        else:
            self.gap_right = 0

    def __parse_i_t_gap(self):
        match = re.search(r">(\.*)\(", self.structure)
        if match:
            self.i_t_gap = len(match.group(1))
        else:
            self.i_t_gap = 0

    def __parse_t_i_gap(self):
        match = re.search(r"(\)\.*\>)", self.structure)
        if match:
            self.t_i_gap = len(match.group(1)) - 2
        else:
            self.t_i_gap = 0

    def __parse_i_b_gap(self):
        match = re.search(r"(\>[^>]*)$", self.structure)
        tmp = re.search(r"\>([^\)]*)\)", match.group(1))
        if tmp:
            self.i_b_gap = len(tmp.group(1))
        else:
            self.i_b_gap = -1

    def __parse_stem_length(self):
        match = re.search(r">\.*(\(.+\))\.*>", self.structure)
        if match:
            self.stem_length = len(match.group(1))
        else:
            raise Exception("No stem?")

    # def __parse_stem_asymmetry(self):

    def get_components(self, as_string=False):
        """
         :returns tuple: (t_i_gap, U_gap, i_b_gap, i_t_gap, gap_right, stem_length)
        """
        items = (self.t_i_gap,
                 self.U_gap,
                 self.i_b_gap,
                 self.i_t_gap,
                 self.gap_right,
                 self.stem_length)
        if as_string:
            return "\t".join(map(str, items))
        else:
            return items

    def __repr__(self):
        return self.structure

    def __str__(self):
        return self.structure


class RNASnoopResult(object):

    """A container for RNASnoop result"""

    def __init__(self, rnasnoop_list):
        """@todo: to be defined1.

        :param rnasnoop_string: @todo

        """
        items = rnasnoop_list[0].split()
        self.structure = RNASnoopStructure(items[0])
        self.left_duplex_energy = float(items[8]) # LE
        self.right_duplex_energy = float(items[10]) # RE
        self.upper_stem_energy = float(items[12]) # DE
        self.lower_stem_energy = float(items[14]) # TE
        self.total_energy = float(items[6][1:])
        self.sequence = rnasnoop_list[1]
        self.modification_position = int(items[3])

    def get_components(self, as_string=False):
        """
         :returns tuple: (left_duplex_energy, right_duplex_energy, upper_stem_energy, lower_stem_energy)
        """
        items = (self.left_duplex_energy,
                 self.right_duplex_energy,
                 self.upper_stem_energy,
                 self.lower_stem_energy
                )
        if as_string:
            return "\t".join(map(str, items))
        else:
            return items


def parse_RNASnoop_ouput(rna_snoop_string):
    tmp_res = [i for i in rna_snoop_string.split("\n")[2:] if i]
    results = []
    for i in range(0, len(tmp_res), 2):
        results.append(RNASnoopResult(tmp_res[i: i + 2]))
    return results
