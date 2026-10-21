import numpy as np
import networkx as nx
import spglib
from ase import Atoms


class AtomsStru(Atoms):

    def __init__(self, symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 cell=None, pbc=None, celldisp=None,
                 constraint=None,
                 calculator=None,
                 info=None,
                 velocities=None):
        super().__init__(symbols, positions, numbers, tags, momenta, masses,
                         magmoms, charges, scaled_positions, cell, pbc,
                         celldisp, constraint, calculator, info, velocities)
        self.l2_norm = self.get_all_distances(mic=True)
        self.adj_matrix = None
        self.graph = None
        self.bond_range = None

    def __set_bond_range__(self, bond_range):
        self.bond_range = bond_range

    def __get_adj_matrix__(self):
        bond_range = self.bond_range
        l2_norm = self.l2_norm
        adj_matrix = np.zeros(l2_norm.shape)
        upper_triangle_indices = np.triu_indices_from(l2_norm, k=1)
        for i, j in zip(*upper_triangle_indices):
            atom_i = self.get_atomic_numbers()[i]
            atom_j = self.get_atomic_numbers()[j]
            bond_ = sorted([atom_i, atom_j])
            bond_name = f'{bond_[0]}-{bond_[1]}'
            if bond_name in bond_range.keys() and bond_range[bond_name][0] <= l2_norm[i, j] <= bond_range[bond_name][1]:
                adj_matrix[i, j] = 1
                adj_matrix[j, i] = 1
        self.adj_matrix = adj_matrix

    def __get_graph__(self):
        G = nx.from_numpy_array(self.adj_matrix)
        for i in range(len(self.numbers)):
            G.nodes[i]["element"] = self.numbers[i]
        self.graph = G

    def get_graph_info(self, bond_range):
        self.bond_range = bond_range
        self.l2_norm = self.get_all_distances(mic=True)
        self.__get_adj_matrix__()
        self.__get_graph__()

    def get_spa_cry(self):
        spagro = spglib.get_spacegroup(self)
        if spagro in range(1, 3):
            crysys = 'triclinic'
        elif spagro in range(3, 16):
            crysys = 'monoclinic'
        elif spagro in range(16, 75):
            crysys = 'orthorhombic'
        elif spagro in range(75, 143):
            crysys = 'tetragonal'
        elif spagro in range(143, 168):
            crysys = 'trigonal'
        elif spagro in range(168, 195):
            crysys = 'hexagonal'
        else:
            crysys = 'cubic'
        return int(spagro.split(' ')[1][1:][:-1]), crysys

    def get_bonds(self):
        bonds = {}
        atoms = self.get_atomic_numbers()
        edges = self.graph.edges
        bond_list = []
        for _e in edges:
            bond_list.append([_e[0], _e[1]])
            bond = sorted([atoms[_e[0]], atoms[_e[1]]])
            bond_str = f'{bond[0]}-{bond[1]}'
            bond_len = round(self.l2_norm[_e[0], _e[1]], 2)
            if bond_str not in bonds:
                bonds[bond_str] = []
            bonds[bond_str].append(bond_len)
        return bond_list, bonds
