import networkx as nx
import numpy as np

from ase import Atoms
from tools.DealRelax.tools.atoms_stru import AtomsStru
from tools.DealRelax.tools.common import (find_Zero_Point,
                                           get_bond_range)

class GetAC:
    """
    处理Atoms对象，得到models等信息
    """

    def __init__(self, atoms: Atoms):
        self.atoms = AtomsStru(atoms)
        self.bond_range = None
        self.__set_bond_range__()
        self.atoms.get_graph_info(self.bond_range)
        self.ac_super = [1, 1, 1]
        self.ac_stru = None  # 活性中心结构
        self.mid_far = False
        self.ac_index = None
        self.__get_AC_conf__(1)

    def __set_bond_range__(self):
        bond_range = get_bond_range([self.atoms])
        self.bond_range = bond_range

    def __get_AC_conf__(self, r: int):
        """
        :param r: path radius
        :return: AC构型，AC的mol原子索引
        """
        if 2 in self.atoms.get_tags() and len(set(self.atoms.get_tags())) > 1:
            cata_index = np.where(self.atoms.get_tags() != 2)[0].tolist()
            mid_index = np.where(self.atoms.get_tags() == 2)[0].tolist()
            mid_atoms = self.atoms[mid_index]
            # 求解零点位置
            nearestNeighbor = find_Zero_Point(l2_norm=self.atoms.l2_norm,
                                              mididx=mid_index,
                                              radius=3.5)
            if not set(nearestNeighbor):
                #  中间体距离催化剂太远
                self.mid_far = True
            # # 对于所有的活性位点，取其1邻
            zp_k_ego = set()
            for z in nearestNeighbor:
                # -----耗时步骤----
                z_12_neighbor = list(
                    nx.single_source_shortest_path_length(G=self.atoms.graph, source=z, cutoff=r).keys())
                zp_k_ego.update(z_12_neighbor)
            ac_index = list(zp_k_ego.union(mid_index))
            AC_configurations = self.atoms[ac_index]
            AC_configurations.get_graph_info(self.bond_range)
            self.ac_index = ac_index
            self.ac_stru = AC_configurations
            self.graph = AC_configurations.graph
