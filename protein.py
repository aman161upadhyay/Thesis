from atom import Atom
import numpy


class Protein:
    """
    Model for a protein.
    Reads a PDB file.
    Info https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    """

    def __init__(self, name, data):
        self.name = name
        self.lines = []
        self.lines_ss = []
        self.lines_neighbors = []
        self.resolution_line = []
        for line in data.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("CONECT"):
                self.lines.append(line)
            if line.startswith("HELIX") or line.startswith("SHEET"):
                self.lines_ss.append(line)
            if line.startswith("SEQRES"):
                self.lines_neighbors.append(line)
            if line.startswith("REMARK   2 "):
                self.resolution_line.append(line)

    @staticmethod
    def parse_line_atom_details(line):
        """ Get atom object from line string"""
        atom_name = line[12:16].strip()
        atom_number = int(line[6:11].strip())
        residue_name = (line[17:20].strip())
        residue_number = int(line[22:26].strip())
        chain = line[21:22].strip()
        b_factor = float(line[60:66].strip())
        occupancy = float(line[54:60].strip())
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        return Atom(atom_name, residue_name, residue_number, chain, b_factor, occupancy, atom_number, x, y, z)

    def contains_amino_acid(self, amino_acid_name):
        """ Check if a protein contains an amino acid """
        for line in self.lines:
            if line[17:20].strip() == amino_acid_name:
                return True
        return False

    def get_atoms(self, amino_acid_name, atom_name):
        """ Returns list of atom objects only for SE and MSE """
        results = []
        for line in self.lines:
            if line[12:16].strip() == atom_name and line[17:20].strip() == amino_acid_name:
                atom = Protein.parse_line_atom_details(line)
                results.append(atom)
        return results

    def get_proximal_atoms(self, center_atom, neighbour_atom_type, max_neighbour_distance):
        """
        Get atoms inside a particular radii's sphere
        (exception) Atoms belonging to water residues are excluded
        """
        results = []
        for line in self.lines:
            if line[12:16].strip()[0] == neighbour_atom_type and line[17:20].strip() != "HOH":
                atom = Protein.parse_line_atom_details(line)
                if Atom.distance(center_atom, atom) <= max_neighbour_distance:
                    results.append(atom)
        return results

    def get_connected_hydrogens(self, center_atom):
        """ Return list of connected hydrogen atoms connected to center atom """
        hydrogen_identifier = center_atom.atom_name[1:]

        if hydrogen_identifier == '':
            return []

        results = []
        for line in self.lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):

                if center_atom.residue_number == int(line[22:26].strip()):
                    line_atom_name = line[12:16].strip()
                    if line_atom_name[0] == 'H' and line_atom_name[1:].startswith(hydrogen_identifier):
                        results.append(Protein.parse_line_atom_details(line))
        return results

    def get_connected_atoms(self, center_atom):
        """ Return list of atoms connected to center atom """
        results = []
        for line in self.lines:
            if line.startswith("CONECT"):
                connections = [int(line[start_index: start_index + 5].strip()) for start_index in
                               range(6, len(line) - 1, 5)]
                if connections[0] == center_atom.atom_number:
                    results = [self.get_atom(connection) for connection in connections[1:]]
        return results

    def get_atom(self, atom_number):
        """ Get atom with a specific residue number """
        if atom_number == 0:
            raise ValueError
        for line in self.lines:
            if int(line[6:11].strip()) == atom_number:
                return Protein.parse_line_atom_details(line)

    # def get_ss(self, name, chain_id, residue_number):
    #     """Get the secondary structure of the residue"""
    #     ss = []
    #     for line in self.lines_ss:
    #         if ((line.startswith("HELIX")) and (residue_number >= int(line[21:25].strip())) and
    #                 (residue_number <= int(line[33:37].strip()))):
    #             decl = "Helix"
    #             helix_class = int(line[38:40].strip())
    #             helix_length = int(line[71:76].strip())
    #             ss.append((name, chain_id, residue_number, decl, helix_class, helix_length))
    #             return ss
    #         if ((line.startswith("SHEET")) and (residue_number >= int(line[22:26].strip())) and
    #                 (residue_number <= int(line[33:37].strip()))):
    #             decl = "Sheet"
    #             sheet_sense = line[38:40].strip()
    #             ss.append((name, chain_id, residue_number, decl, sheet_sense))
    #             return ss
    #     ss.append((name, chain_id, residue_number, "Loop"))
    #     return ss

    def get_ss(self, name, chain_id, residue_number):
        """Get the secondary structure of the residue"""

        for line in self.lines_ss:
            if ((line.startswith("HELIX")) and (residue_number >= int(line[21:25].strip())) and
                    (residue_number <= int(line[33:37].strip()))):
                decl = "Helix"
                return decl
            if ((line.startswith("SHEET")) and (residue_number >= int(line[22:26].strip())) and
                    (residue_number <= int(line[33:37].strip()))):
                decl = "Sheet"
                return decl
        decl = "Loop"
        return decl

    def neighbors(self):
        """Extracting the primary sequence from a pdb file"""
        primary_sequence = []
        for line in self.lines_neighbors:
            if line[11:12].strip() == "A":
                primary_sequence.append(line[19:22].strip())
                primary_sequence.append(line[23:26].strip())
                primary_sequence.append(line[27:30].strip())
                primary_sequence.append(line[31:34].strip())
                primary_sequence.append(line[35:38].strip())
                primary_sequence.append(line[39:42].strip())
                primary_sequence.append(line[43:46].strip())
                primary_sequence.append(line[47:50].strip())
                primary_sequence.append(line[51:54].strip())
                primary_sequence.append(line[55:58].strip())
                primary_sequence.append(line[59:62].strip())
                primary_sequence.append(line[63:66].strip())
                primary_sequence.append(line[67:70].strip())

        return primary_sequence

    def get_avg_b_factor(self):
        b_fac = []
        for line in self.lines:
            if line[12:16].strip() == "CA" or line[12:16].strip() == "SE" or line[12:16].strip() == "O"\
                    or line[12:16].strip() == "N":
                b_factor = float(line[60:66].strip())
                b_fac.append(b_factor)
        avg_b_factor = numpy.mean(b_fac)
        return avg_b_factor

    def get_stdev_b_factor(self):
        b_fac = []
        for line in self.lines:
            if line[12:16].strip() == "CA" or line[12:16].strip() == "SE" or line[12:16].strip() == "O"\
                    or line[12:16].strip() == "N":
                b_factor = float(line[60:66].strip())
                b_fac.append(b_factor)
        stdev_b_factor = numpy.std(b_fac)
        return stdev_b_factor

    def get_resolution(self):
        for line in self.resolution_line:
            return line[23:30].strip()
