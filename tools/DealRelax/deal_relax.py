import math
import networkx as nx
import numpy as np

from ase import Atom, Atoms
from DealRelax.tools.atoms_stru import AtomsStru
from DealRelax.tools.common import greedy_tsp, get_subgraph_unit, \
    find_Zero_Point, get_bond_range, detect_vacuum_direction


class DealRelax:
    """
    处理Atoms对象，得到models等信息
    """

    def __init__(self, atoms: Atoms):
        self.atoms = AtomsStru(atoms)
        self.bond_range = None
        self.__set_bond_range__()

        self.atoms.get_graph_info(self.bond_range)
        self.graph = self.atoms.graph
        self.spagro, self.crysys = self.atoms.get_spa_cry()
        self.bond_list, self.bond_len = self.atoms.get_bonds()
        self.current_tag()
        self.super_vector, self.as_super = self.super_num()
        self.cata_info = self.get_cata_info()  # 催化剂的models信息
        self.mol_info = self.get_mol_info()  # 分子的models信息
        self.as_stru = None
        self.mid_far = False
        self.as_info = self.get_as_info()  # 活性中心的models信息

    def __set_bond_range__(self):
        bond_range = get_bond_range([self.atoms])
        self.bond_range = bond_range

    def get_cata_info(self):
        """
        获取cata的信息
        """
        cata_info = {
            "label": 1,
            "atomic_num": [],
            "bond_list": [],
            "tags": {}
        }
        catas_info = []
        all_tags = self.atoms.get_tags()
        cata_index = np.where(all_tags != 2)[0].tolist()
        # --------------------- 扩包处理------------------------------
        if cata_index:
            cataAtoms = self.atoms[cata_index] * self.super_vector
            cataAtoms.get_graph_info(self.bond_range)
            cata_graph_unit = get_subgraph_unit(cataAtoms.graph, 1)
            for _g in cata_graph_unit:
                g_nodes = list(_g.nodes)
                cata_info['atomic_num'].append(cataAtoms.get_atomic_numbers()[g_nodes].tolist())
                cata_info['bond_list'].append([[g_nodes.index(_b[0]), g_nodes.index(_b[1])] for _b in _g.edges])
                cata_info['tags'] = {
                    'spa_gro': self.spagro,
                    'cry_sys': self.crysys
                }
            catas_info.append(cata_info)
        return catas_info

    def get_mol_info(self):
        """
        获取mol的信息
        """
        mol_info = {
            "label": 2,
            "atomic_num": [],
            "atom_tags": [],
            "bond_list": [],
            "tags": {}
        }
        mols_info = []
        all_tags = self.atoms.get_tags()
        mol_index = np.where(all_tags == 2)[0].tolist()
        if mol_index:
            mol_atoms = self.atoms[mol_index]
            mol_atoms.get_graph_info(self.bond_range)
            mol_comp = list(nx.connected_components(mol_atoms.graph))
            for _mol in mol_comp:
                _mol_info = mol_info.copy()
                _mol_index = [mol_index[_i] for _i in list(_mol)]
                _mol_atoms = self.atoms[_mol_index]
                _mol_atoms.get_graph_info(self.bond_range)
                _mol_info['atomic_num'] = _mol_atoms.get_atomic_numbers().tolist()
                _mol_info['atom_tags'] = _mol_atoms.get_tags().astype(int).tolist()
                _mol_graph = _mol_atoms.graph
                _mol_info['bond_list'] = [[_b[0], _b[1]] for _b in
                                          _mol_graph.edges]
                mols_info.append(_mol_info)
        return mols_info

    def get_as_info(self):
        """
        获取as的信息
        """
        as_info = {
            "label": 3,
            "atomic_num": [],
            "atom_tags": [],
            "bond_list": [],
            "tags": {}
        }
        ass_info = []
        self.__get_AC_conf__(1)
        AC_conf = self.as_stru
        if AC_conf is not None:
            as_info['atomic_num'] = AC_conf.get_atomic_numbers().tolist()
            as_info['atom_tags'] = AC_conf.get_tags().tolist()
            as_info['bond_list'] = [[_b[0], _b[1]] for _b in AC_conf.graph.edges]
            ass_info.append(as_info)
        return ass_info

    def __get_AC_conf__(self, r: int):
        """
        :param r: path radius
        :return: AC构型，AC的mol原子索引
        """
        if 2 in self.atoms.get_tags() and len(set(self.atoms.get_tags())) > 1:
            cata_index = np.where(self.atoms.get_tags() != 2)[0].tolist()
            mid_index = np.where(self.atoms.get_tags() == 2)[0].tolist()
            # 干净表面扩胞
            super_cata_atoms = self.atoms[cata_index] * self.as_super
            # 在扩胞的表面上添加中间体
            mid_atoms = self.atoms[mid_index]
            add_order = greedy_tsp(mid_atoms.get_all_distances(mic=True))
            vector_lst = []
            new_mid_index = set()
            sum_cell_vec = self.atoms.get_cell()[0] + self.atoms.get_cell()[1]
            for _k, _idx in enumerate(add_order):
                if _k == 0:
                    vec = sum_cell_vec
                    sym = mid_atoms[_idx].symbol
                    tag = mid_atoms[_idx].tag
                    vector_lst.append((sym, vec, tag))
                else:
                    vec = mid_atoms.get_distance(add_order[_k - 1], _idx, mic=True, vector=True)
                    sym = mid_atoms[_idx].symbol
                    tag = mid_atoms[_idx].tag
                    vector_lst.append((sym, vec, tag))
            # 添加原子，返回新mid的index
            last_pos = mid_atoms[0].position
            new_stru = super_cata_atoms
            for _item in vector_lst:
                new_stru.append(
                    Atom(symbol=_item[0],
                         position=_item[1] + last_pos,
                         tag=_item[2])
                )
                new_mid_index.add(len(new_stru) - 1)
                last_pos = _item[1] + last_pos
            # 粗取9A
            big_ac_index = set()
            l2_norm = new_stru.get_all_distances(mic=True)
            for _i in list(new_mid_index):
                i_index = np.where(l2_norm[_i] < 9)[0].tolist()
                big_ac_index.update(i_index)
            big_ac = new_stru[list(big_ac_index)]
            # -----耗时步骤----
            big_ac.get_graph_info(self.bond_range)
            actag = big_ac.get_tags()
            big_mid_index = []
            for _i, _v in enumerate(actag):
                if _v == 2:
                    big_mid_index.append(_i)
            # 求解零点位置
            nearestNeighbor = find_Zero_Point(l2_norm=big_ac.l2_norm,
                                              mididx=big_mid_index,
                                              radius=3.5)
            if not set(nearestNeighbor):
                #  中间体距离催化剂太远
                self.mid_far = True
            # # 对于所有的活性位点，取其1邻
            zp_k_ego = set()
            for z in nearestNeighbor:
                # -----耗时步骤----
                z_12_neighbor = list(
                    nx.single_source_shortest_path_length(G=big_ac.graph, source=z, cutoff=r).keys())
                zp_k_ego.update(z_12_neighbor)
            ac_index = list(zp_k_ego.union(big_mid_index))
            AC_configurations = big_ac[ac_index]
            AC_configurations.get_graph_info(self.bond_range)
            self.as_stru = AC_configurations
            self.graph = AC_configurations.graph

    def super_num(self):
        vac_dir = detect_vacuum_direction(self.atoms)
        _a, _b, _c = np.linalg.norm(self.atoms.get_cell(), axis=1)
        n_x = math.ceil(9 / _a)
        n_y = math.ceil(9 / _b)
        n_z = math.ceil(9 / _c)
        super_num = [n_x, n_y, n_z]
        super_num[vac_dir] = 1
        as_super_num = []
        for _i, _s in enumerate(super_num):
            if _i != vac_dir and _s < 3:
                as_super_num.append(3)
            else:
                as_super_num.append(_s)
        return super_num, as_super_num

    def current_tag(self):
        """
        矫正tag，只和tag为2相连的原子，tag为2
        """
        tag2_index = np.where(self.atoms.get_tags() == 2)[0].tolist()
        old_index = []
        self.atoms.get_graph_info(self.bond_range)
        stru = self.atoms
        adj_matrix = stru.adj_matrix
        while list(set(tag2_index) - set(old_index)):
            mid_index = tag2_index[:]
            for i in list(set(tag2_index) - set(old_index)):
                neighbors = np.where(adj_matrix[i] == 1)[0].tolist()
                for j in neighbors:
                    j_neighbors = np.where(adj_matrix[j] == 1)[0].tolist()
                    if not list(set(j_neighbors) - set(tag2_index)):
                        tag2_index.append(j)
            old_index = mid_index[:]
        for _index in tag2_index:
            self.atoms[_index].tag = 2

