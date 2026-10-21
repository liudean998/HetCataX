import ase
import math
import numpy as np
import networkx as nx

from ..static.element_info import radii
from sklearn.neighbors import KernelDensity
from scipy.signal import find_peaks
from typing import Any, List
from ase import Atoms


def get_bond_range(atom_systems: List[ase.Atoms]) -> dict:
    """
    获取原子体系各种键长的成键范围
    :param atom_systems:
    :return:
    """
    l2_lsts: List[np.ndarray] = [item.get_all_distances(mic=True)
                                 for item in atom_systems]  # 三个体系的二体距离矩阵
    upper_triangle_indices = [np.triu_indices_from(i, k=1) for i in l2_lsts]  # 二体距离矩阵的上三角索引
    atomic_lsts = [item.get_atomic_numbers() for item in atom_systems]  # 三个体系的原子序数列表
    bond_hist = {}
    bond_range = {}
    for l2, upper_tri_idx, atomic_lst in zip(l2_lsts, upper_triangle_indices, atomic_lsts):
        for i, j in zip(*upper_tri_idx):
            atom_i = atomic_lst[i]
            atom_j = atomic_lst[j]
            bond_ = sorted([atom_i, atom_j])
            bond_ = [str(i) for i in bond_]
            bond_name = "-".join(bond_)
            dis_ij_max = (radii[atom_j - 1] + radii[atom_i - 1]) * 1.2
            if 0 < l2[i, j] < dis_ij_max and bond_name not in bond_hist.keys():
                bond_hist[bond_name] = []
        # 统计所有的距离在(0,4)内的所有二体关系
        for i, j in zip(*upper_tri_idx):
            if 4 > l2[i, j] > 0:
                bond_ = sorted([atomic_lst[i], atomic_lst[j]])
                bond_ = [str(i) for i in bond_]
                bond_name = "-".join(bond_)
                if bond_name in bond_hist.keys():
                    bond_hist[bond_name].append(l2[i, j])
    # 高斯核密度对分布进行拟合，寻找峰值和双峰之间的谷值，计算成键距离的范围
    for key, value in bond_hist.items():
        if value:
            kde = KernelDensity(kernel='gaussian', bandwidth=0.08)
            kde.fit(np.array(value).reshape(-1, 1))
            m = np.linspace(0, 4, 4000)
            n = np.exp(kde.score_samples(m.reshape(-1, 1)))
            # from matplotlib import pyplot as plt
            # plt.plot(m, n, label=key)
            # plt.title(f"Gaussian Kernel Density of {key}")
            # plt.show()
            dy_dx_forward = np.gradient(n, m, edge_order=1)
            dy_dx_zero = [i for i in range(len(dy_dx_forward)) if math.isclose(dy_dx_forward[i], 0, abs_tol=0.01)]
            # print(dy_dx_zero)
            peaks = sorted(find_peaks(n)[0])
            # print(peaks)
            if len(peaks) == 0:
                continue
            elif len(peaks) == 1:
                upper_bound_lst = [i for i in dy_dx_zero if i > peaks[0] + 10]
                upper_bound_lst.append(3999)
                upper_bound = min(upper_bound_lst)
                bond_range[key] = (2 * m[peaks[0]] - m[upper_bound], m[upper_bound])
            else:
                # 考虑到一个化学键的长度可能是多峰分布，通过聚类的方式将最有可能是成键的峰位置划为同一簇
                min_peaks_cluse = [peaks[0]]
                for x in peaks[1:]:
                    for y in min_peaks_cluse:
                        if 0.83 * x <= y <= 1.2 * x:
                            min_peaks_cluse.append(x)
                            break
                min_peaks_mean = int(np.mean(min_peaks_cluse))
                if len(min_peaks_cluse) == len(peaks):
                    upper_bound_lst = [i for i in dy_dx_zero if i > peaks[-1] + 10]
                    upper_bound_lst.append(3999)
                    upper_bound = min(upper_bound_lst)
                    bond_range[key] = (2 * m[min_peaks_mean] - m[upper_bound], m[upper_bound])
                else:
                    start, end = peaks[len(min_peaks_cluse) - 1], peaks[len(min_peaks_cluse)]
                    upper_bound_lst = [i for i in range(start + 10, end - 10, 1) if
                                       math.isclose(dy_dx_forward[i], 0, abs_tol=0.01)]
                    upper_bound_lst.extend(
                        [i for i in range(start + 1, end - 2, 1) if dy_dx_forward[i] * dy_dx_forward[i + 1] < 0])
                    upper_bound_lst = list(set(upper_bound_lst))
                    upper_bound = min(upper_bound_lst)
                    bond_range[key] = (2 * m[min_peaks_mean] - m[upper_bound], m[upper_bound])
            # print(f"{key} bond range: {bond_range[key]}")
    return bond_range


def get_subgraph_unit(graph: nx.Graph,
                      k: int = 1,
                      iterable: Any = None) -> List[nx.Graph]:
    """
    以每个点为中心，向外寻k阶邻域，构成邻域子图，获取原图的子图单元列表。如果输入一个可迭代对象则只从这个可迭代对象中的点寻找邻域子图
    :param graph: A NetworkX Graph
    :param k: Include all neighbors of distance <= k from node.
    :param iterable: only find subgraph unit from this iterable (default: None)
    :return:  list
    """
    nodes_lst = iterable if iterable is not None else graph.nodes
    subgraph_unit: list = []
    for node in nodes_lst:
        ego_nodes = nx.single_source_shortest_path_length(graph, node, cutoff=k)
        ego_subgraph = graph.subgraph(ego_nodes)
        iso = 1
        for unit in subgraph_unit:
            if nx.vf2pp_is_isomorphic(ego_subgraph, unit, node_label='element'):
                iso *= 0
                break
        if iso:
            subgraph_unit.append(ego_subgraph)
    return subgraph_unit


def greedy_tsp(distance_matrix) -> list:
    """
    贪心算法求解TSP问题
    :param distance_matrix:
    :return:
    """
    n = len(distance_matrix)
    visited = [False] * n
    tour = [0]
    visited[0] = True
    for _ in range(n - 1):
        current_city = tour[-1]
        min_distance = float('inf')
        next_city = None
        for city in range(n):
            if not visited[city] and distance_matrix[current_city][city] < min_distance:
                min_distance = distance_matrix[current_city][city]
                next_city = city
        visited[next_city] = True
        tour.append(next_city)
    return tour


def find_vector_positions(matrix, vector_list):
    positions = []
    for vector in vector_list:
        position = np.where(np.all(matrix == vector, axis=1))[0]
        if position.size > 0:
            positions.append(position[0])
        else:
            positions.append(None)
    return positions


def find_Zero_Point(l2_norm: np.ndarray, mididx: list, radius: float) -> list:
    nearestNeighbor = set()
    mid = np.abs(mididx)
    for i in mid:
        for j, dis in enumerate(l2_norm[i]):
            if dis <= radius and j not in mid:
                nearestNeighbor.add(j)
    return list(nearestNeighbor)


def detect_vacuum_direction(atoms: Atoms):
    # 判断真空层在哪个轴
    cell_lengths = atoms.cell.cellpar()[:3]  # Deprecated since version 3.21.0: use atoms.cell.cellpar() instead
    atom_positions = atoms.get_positions()
    min_positions = np.min(atom_positions, axis=0)
    max_positions = np.max(atom_positions, axis=0)
    gaps = cell_lengths - (max_positions - min_positions)
    vacuum_direction = np.argmax(gaps)

    return vacuum_direction
