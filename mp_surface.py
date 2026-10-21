from __future__ import annotations

import math
import os
import random
from collections import defaultdict
from pathlib import Path

import numpy as np

from ase import Atoms
from ase.constraints import FixAtoms
from ase.db import connect

from mp_api.client import MPRester

from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.core import Composition, Structure
from pymatgen.core.surface import SlabGenerator, get_symmetrically_distinct_miller_indices
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


# ============================================================
# 从 Materials Project 获取 Cu-X 二元 bulk
# ============================================================

def fetch_cu_binaries(api_key, second_elements):
    with MPRester(api_key) as mpr:
        docs = mpr.materials.summary.search(
            # 必须包含 Cu
            elements=["Cu"],
            # 只能有两种元素
            nelements=(2, 2),
        )
    rows = []
    for doc in docs:
        structure = doc.structure
        # ----------------------------------------------------
        # 跳过无序结构
        # ----------------------------------------------------
        if not structure.is_ordered:
            print(f"[skip disordered] {doc.material_id} {doc.formula_pretty}")
            continue
        # ----------------------------------------------------
        # 找 Cu 之外的第二元素
        # ----------------------------------------------------
        elems = [str(e) for e in structure.composition.elements]
        second = next((e for e in elems if e != "Cu"), None)
        if second is None:
            continue
        # ----------------------------------------------------
        # 只保留随机抽中的元素
        # ----------------------------------------------------
        if second not in second_elements:
            continue
        rows.append(
            {
                "material_id": str(doc.material_id),
                "formula": str(doc.formula_pretty),
                "energy_above_hull": (
                    float(doc.energy_above_hull) if doc.energy_above_hull is not None else None
                ),
                "second_element": second,
                "structure": structure,
            }
        )
    rows.sort(key=lambda r: (r["second_element"], r["material_id"]))
    return rows

# ============================================================
# Pymatgen Structure -> ASE Atoms
# ============================================================
def structure_to_ase(structure: Structure) -> Atoms:
    atoms = AseAtomsAdaptor.get_atoms(structure)
    atoms.set_pbc(True)
    return atoms

# ============================================================
# bulk 标准化
# ============================================================

def standardize_bulk(atoms: Atoms) -> Structure:
    struct = AseAtomsAdaptor.get_structure(atoms)
    return SpacegroupAnalyzer(struct, symprec=SYMPREC).get_conventional_standard_structure()

# ============================================================
# 判断 slab 上下是否等价
# ============================================================
def is_structure_invertible(structure: Structure) -> bool:
    sga = SpacegroupAnalyzer(structure, symprec=SYMPREC)
    for operation in sga.get_symmetry_operations():
        if np.isclose(operation.affine_matrix[2, 2], -1.0):
            return True
    return False
# ============================================================
# 翻转 slab
# ============================================================


def flip_struct(structure: Structure) -> Structure:
    atoms = AseAtomsAdaptor.get_atoms(structure)
    atoms.wrap()
    atoms.rotate(180, "x", rotate_cell=True, center="COM")
    cell = atoms.cell.array.copy()
    if cell[2][2] < 0.0:
        cell[2] = -cell[2]
    if np.cross(cell[0], cell[1])[2] < 0.0:
        cell[1] = -cell[1]
    atoms.set_cell(cell, scale_atoms=False)
    atoms.wrap()
    return AseAtomsAdaptor.get_structure(atoms)


# ============================================================
# 14. 枚举一个 bulk 的全部 surface
# ============================================================


def enumerate_surfaces(bulk_atoms: Atoms, max_miller=MAX_MILLER):
    bulk_struct = standardize_bulk(bulk_atoms)
    all_slabs_info = []
    # --------------------------------------------------------
    # 枚举对称不等价 Miller 面
    # --------------------------------------------------------
    miller_indices = get_symmetrically_distinct_miller_indices(bulk_struct, max_miller)
    for miller in miller_indices:
        slab_gen = SlabGenerator(
            initial_structure=bulk_struct,
            miller_index=miller,
            min_slab_size=MIN_SLAB_SIZE,
            min_vacuum_size=MIN_VACUUM_SIZE,
            lll_reduce=False,
            center_slab=True,
            primitive=True,
            max_normal_search=1,
        )
        # ----------------------------------------------------
        # 获取不同 termination
        # ----------------------------------------------------
        slabs = slab_gen.get_slabs(tol=SLAB_TOL, bonds=None, max_broken_bonds=0, symmetrize=False)
        for slab in slabs:
            # top surface
            all_slabs_info.append((slab, tuple(miller), float(slab.shift), True))
            # 如果上下表面不等价，
            # 再加入翻转后的 bottom surface
            if not is_structure_invertible(slab):
                all_slabs_info.append((flip_struct(slab), tuple(miller), float(slab.shift), False))
    return all_slabs_info


# ============================================================
# 15. xy 扩胞
# ============================================================


def tile_atoms_to_min_xy(atoms: Atoms, min_xy=MIN_XY) -> Atoms:
    x_len = np.linalg.norm(atoms.cell[0])
    y_len = np.linalg.norm(atoms.cell[1])
    nx = int(math.ceil(min_xy / x_len))
    ny = int(math.ceil(min_xy / y_len))
    return atoms.repeat((nx, ny, 1))


# ============================================================
# 16. bulk 配位数
#
# 每个 bulk 只计算一次。
# ============================================================


def bulk_coordination(bulk_atoms: Atoms):
    bulk_struct = AseAtomsAdaptor.get_structure(bulk_atoms)
    sym_struct = SpacegroupAnalyzer(bulk_struct, symprec=SYMPREC).get_symmetrized_structure()
    vnn = VoronoiNN(tol=VORONOI_TOL)
    cn_dict = defaultdict(set)
    for group in sym_struct.equivalent_indices:
        idx = group[0]
        site = sym_struct[idx]
        cn = round(float(vnn.get_cn(sym_struct, idx, use_weights=True)), 5)
        cn_dict[site.species_string].add(cn)
    return cn_dict


# ============================================================
# 17. Voronoi 方法判断表面原子
# ============================================================


def voronoi_surface_tags(bulk_cn, surface_atoms: Atoms):
    surface_struct = AseAtomsAdaptor.get_structure(surface_atoms)
    weights = [float(site.specie.atomic_mass) for site in surface_struct]
    com = np.average(surface_struct.frac_coords, weights=weights, axis=0)
    vnn = VoronoiNN(tol=VORONOI_TOL)
    tags = []
    for idx, site in enumerate(surface_struct):
        # 只判断 slab 上半部分
        if site.frac_coords[2] > com[2]:
            try:
                cn = round(float(vnn.get_cn(surface_struct, idx, use_weights=True)), 5)
                tags.append(1 if cn < min(bulk_cn[site.species_string]) else 0)
            except (RuntimeError, ValueError, KeyError):
                tags.append(1)
        else:
            tags.append(0)
    return tags


# ============================================================
# 18. 高度法判断表面原子
# ============================================================


def height_surface_tags(surface_atoms: Atoms, cutoff=SURFACE_HEIGHT_CUTOFF):
    cell_height = np.linalg.norm(surface_atoms.cell[2])
    scaled = surface_atoms.get_scaled_positions()
    max_z = max(p[2] for p in scaled)
    threshold = max_z - cutoff / cell_height
    return [0 if p[2] < threshold else 1 for p in scaled]


# ============================================================
# 19. 最终处理一个 surface
#
# 1. 转 ASE
# 2. xy 扩胞
# 3. Voronoi + 高度判断 surface atoms
# 4. tag = 1 表面原子
# 5. tag = 0 固定
# ============================================================
def process_surface(bulk_atoms: Atoms, surface_info, bulk_cn):
    slab_struct, miller, shift, top = surface_info
    surface_atoms = AseAtomsAdaptor.get_atoms(slab_struct)
    # --------------------------------------------------------
    # xy 至少扩到 8 × 8 Å
    # --------------------------------------------------------
    surface_atoms = tile_atoms_to_min_xy(surface_atoms)
    # --------------------------------------------------------
    # 检查组成是否发生变化
    # --------------------------------------------------------
    bulk_formula = Composition(bulk_atoms.get_chemical_formula()).reduced_formula
    slab_formula = Composition(surface_atoms.get_chemical_formula()).reduced_formula
    if bulk_formula != slab_formula:
        raise ValueError("Bulk/surface formula mismatch: " f"{bulk_formula} vs {slab_formula}")
    # --------------------------------------------------------
    # Voronoi surface tags
    # --------------------------------------------------------
    vtags = voronoi_surface_tags(bulk_cn, surface_atoms)
    # --------------------------------------------------------
    # Height surface tags
    # --------------------------------------------------------
    htags = height_surface_tags(surface_atoms)
    # --------------------------------------------------------
    # 两种方法只要有一种认为是 surface
    # 就标记为 tag = 1
    # --------------------------------------------------------
    tags = [max(v, h) for v, h in zip(vtags, htags)]
    surface_atoms.set_tags(tags)
    # --------------------------------------------------------
    # tag = 0 固定
    # tag = 1 不固定
    # --------------------------------------------------------
    fixed_mask = [atom.tag == 0 for atom in surface_atoms]
    surface_atoms.set_constraint(FixAtoms(mask=fixed_mask))
    return surface_atoms


# ============================================================
# 20. 主程序
# ============================================================


def main():
    outdir = Path(OUTPUT_DIR)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.SystemRandom()
    # ========================================================
    # 第一步：
    # 指定元素
    # ========================================================
    second_elements = {"Zr", "Al", "Ga", "Sb", "Y", "Pd", "Hf", "Ge", "Zn", "Au", "Sc", "Ru", "Ag"}
    api_key = "LcBedZCVS46Bk9h745vRsP8tsfPewXLf"
    # ========================================================
    # 第二步：
    # 查询元素对应的 Cu-X bulk
    # ========================================================
    print("\n开始从 Materials Project 查询 Cu-X bulk...")
    rows = fetch_cu_binaries(
        api_key=api_key,
        second_elements=second_elements
    )
    # ========================================================
    # 查看元素实际找到了 bulk
    # ========================================================
    found_elements = sorted({row["second_element"] for row in rows})
    missing_elements = sorted(set(second_elements) - set(found_elements))
    print("\n" + "=" * 80)
    print(f"选择元素数：{len(second_elements)}")
    print(f"实际找到 bulk 的元素数：{len(found_elements)}")
    print(f"没有找到 bulk 的元素数：{len(missing_elements)}")
    print(f"bulk 总数：{len(rows)}")
    print("=" * 80)
    print("\n成功找到 bulk 的元素：")
    print(found_elements)
    if missing_elements:
        print("\n没有找到满足条件 bulk 的元素：")
        print(missing_elements)
    if not rows:
        raise RuntimeError("没有找到任何满足条件的 Cu-X bulk。")
    # ========================================================
    # 打印所有 bulk
    # ========================================================
    print("\n" + "=" * 100)
    print(f"一共找到 {len(rows)} 个 Cu 二元 bulk")
    print("=" * 100)
    for i, row in enumerate(rows):
        print(
            f"{i:4d} | {row['material_id']:15s} | {row['formula']:15s} | element = {row['second_element']:3s} | E_hull = {row['energy_above_hull']}"
        )
    print("=" * 100)
    # ========================================================
    # 第三步：
    # 随机 bulk
    #     ↓
    # 枚举这个 bulk 的全部 surface
    #     ↓
    # 从该 bulk 的 surface 中随机一个
    #     ↓
    # 删除这个已经选择过的 surface
    #     ↓
    # 再重新随机 bulk
    # 重复直到得到 MAX_SURFACES 个 surface
    #
    # 注意：
    # 同一个 bulk 可以被多次抽中，
    # 但同一个具体 surface 不会重复。
    # ========================================================
    print("\n" + "=" * 80)
    print("开始按照：随机 bulk -> 枚举 surface -> 随机 surface")
    print("=" * 80)
    surface_cache = {}
    # 每个 bulk 剩余还没有使用过的 surface
    remaining_surfaces = {}
    # --------------------------------------------------------
    # 最终抽中的 (bulk_index, surface_index)
    # --------------------------------------------------------
    selected_candidates = []
    # --------------------------------------------------------
    # 完全没有 surface 的 bulk
    # --------------------------------------------------------
    no_surface_bulks = set()
    total_enumerated = 0
    while len(selected_candidates) < MAX_SURFACES:
        available_bulks = []
        for bulk_index in range(len(rows)):
            # 已知没有 surface
            if bulk_index in no_surface_bulks:
                continue
            if bulk_index not in surface_cache:
                available_bulks.append(bulk_index)
                continue
            if len(remaining_surfaces[bulk_index]) > 0:
                available_bulks.append(bulk_index)
        if not available_bulks:
            print(
                "\nWARNING: 所有 bulk 的 surface 已全部使用。\n"
                f"最终只能得到 {len(selected_candidates)} 个 surface。"
            )
            break

        # ====================================================
        # 3. 先随机选择一个 bulk
        # ====================================================
        bulk_index = rng.choice(available_bulks)
        row = rows[bulk_index]
        mpid = row["material_id"]
        formula = row["formula"]
        element = row["second_element"]
        print(
            f"\n随机 bulk: "
            f"{bulk_index} | {mpid} | {formula} | element={element}"
        )
        # ====================================================
        # 4. 如果这个 bulk 是第一次被抽中
        #    才枚举它所有 surface
        # ====================================================
        if bulk_index not in surface_cache:
            bulk_atoms = structure_to_ase(row["structure"])
            print("    第一次抽到该 bulk，开始枚举全部 surface...")
            info = enumerate_surfaces(bulk_atoms, max_miller=MAX_MILLER)
            print(f"    共枚举到 {len(info)} 个 surface")
            total_enumerated += len(info)
            if len(info) == 0:
                print(
                    f"    WARNING: {mpid} {formula} "
                    "没有生成任何 surface"
                )
                no_surface_bulks.add(bulk_index)
                continue
            # -----------------------------------------------
            # 缓存
            # -----------------------------------------------
            surface_cache[bulk_index] = {
                "row": row,
                "bulk_atoms": bulk_atoms,
                "info": info,
                "bulk_cn": None,
            }
            # -----------------------------------------------
            # 初始化这个 bulk 所有可用 surface
            # -----------------------------------------------
            remaining_surfaces[bulk_index] = list(range(len(info)))
            print()
        # ====================================================
        # 5. 从刚刚选中的 bulk 中随机一个 surface
        # ====================================================
        sidx = rng.choice(remaining_surfaces[bulk_index])
        surface_info = surface_cache[bulk_index]["info"][sidx]
        slab_struct, miller, shift, top = surface_info
        print(
            f"    随机 surface: index={sidx} | "
            f"hkl={miller} | "
            f"{'top' if top else 'bottom'}"
        )
        # ====================================================
        # 6. 保存这个随机结果
        # ====================================================
        selected_candidates.append((bulk_index, sidx))
        # ====================================================
        # 7. 删除已经使用的 surface
        # 保证：
        # 同一个具体 surface 不会重复出现
        # ====================================================
        remaining_surfaces[bulk_index].remove(sidx)
        print(
            f"    当前已选择："
            f"{len(selected_candidates)}/{MAX_SURFACES}"
        )
    print("\n" + "=" * 80)
    print(
        f"最终随机选择了 "
        f"{len(selected_candidates)} 个 surface"
    )
    print(
        f"实际被枚举过的 bulk："
        f"{len(surface_cache)} 个"
    )
    print(
        f"实际枚举出的 surface 总数："
        f"{total_enumerated}"
    )
    print("=" * 80)
    # ========================================================
    # 第五步：
    # 先统计随机后的元素分布
    # ========================================================
    element_surface_count = defaultdict(int)
    bulk_surface_count = defaultdict(int)
    for bulk_index, sidx in selected_candidates:
        row = surface_cache[bulk_index]["row"]
        element_surface_count[row["second_element"]] += 1
        bulk_surface_count[bulk_index] += 1
    print("\n" + "=" * 80)
    print("随机抽取后的元素分布")
    print("=" * 80)
    for element, count in sorted(element_surface_count.items(), key=lambda x: x[1], reverse=True):
        percentage = count / len(selected_candidates) * 100 if selected_candidates else 0.0
        print(f"{element:>3s} : {count:4d} ({percentage:6.2f}%)")
    print("=" * 80)
    print(f"总计：{len(selected_candidates)} 个 surface")
    # ========================================================
    # 第六步：
    # 查看哪些 bulk 没有被抽到
    # ========================================================
    missing_bulks = []
    for i, row in enumerate(rows):
        # 没有成功枚举 surface 的 bulk
        if i not in surface_cache:
            missing_bulks.append(
                (i, row["material_id"], row["formula"], row["second_element"], "no_surface")
            )
            continue
        # 有 surface，但随机时一个也没中
        if bulk_surface_count[i] == 0:
            missing_bulks.append(
                (i, row["material_id"], row["formula"], row["second_element"], "not_selected")
            )
    print("\n" + "=" * 80)
    print("Bulk 抽样情况")
    print("=" * 80)
    print(f"bulk 总数：{len(rows)}")
    print(f"至少抽到一个 surface 的 bulk：{len(rows) - len(missing_bulks)}")
    print(f"没有进入最终数据集的 bulk：{len(missing_bulks)}")
    if missing_bulks:
        print("\n没有进入最终数据集的 bulk：")
        for i, mpid, formula, element, reason in missing_bulks:
            print(f"{i:4d} | {mpid:15s} | {formula:15s} | element = {element:3s} | {reason}")
    # ========================================================
    # 第七步：
    # 最终只建立一个 ASE DB
    # ========================================================
    surface_db_path = outdir / "surfaces.db"
    # --------------------------------------------------------
    # 避免重复运行时向旧 DB 追加
    # --------------------------------------------------------
    if surface_db_path.exists():
        surface_db_path.unlink()
    surface_db = connect(surface_db_path)
    print("\n" + "=" * 80)
    print("开始处理并写入随机抽中的 surface")
    print(f"最终 DB：{surface_db_path}")
    print("=" * 80)
    total_saved = 0
    # ========================================================
    # 第八步：
    # 处理并写入 DB
    # ========================================================
    for bulk_index, sidx in selected_candidates:
        cache = surface_cache[bulk_index]
        row = cache["row"]
        bulk_atoms = cache["bulk_atoms"]
        info = cache["info"]
        mpid = row["material_id"]
        formula = row["formula"]
        surface_info = info[sidx]
        slab_struct, miller, shift, top = surface_info
        # ----------------------------------------------------
        # bulk 配位信息第一次需要时才计算
        # 后面同一 bulk 直接复用
        # ----------------------------------------------------
        if cache["bulk_cn"] is None:
            print(f"    计算 bulk coordination: {mpid}")
            cache["bulk_cn"] = bulk_coordination(bulk_atoms)
        bulk_cn = cache["bulk_cn"]
        # ----------------------------------------------------
        # 真正处理 surface
        # ----------------------------------------------------
        atoms = process_surface(bulk_atoms, surface_info, bulk_cn)
        # ----------------------------------------------------
        # surface 原子数量
        # ----------------------------------------------------
        n_surface_atoms = int(np.sum(atoms.get_tags() == 1))
        # ----------------------------------------------------
        # 写入统一的 surfaces.db
        # 注意：
        # 不使用 "formula" 作为自定义字段名，
        # 避免和 ASE DB 自带字段冲突。
        # ----------------------------------------------------
        surface_db.write(
            atoms,
            # bulk 信息
            bulk_index=int(bulk_index),
            material_id=str(mpid),
            bulk_formula=str(formula),
            second_element=str(row["second_element"]),
            energy_above_hull=float(row["energy_above_hull"]),
            # surface 信息
            surface_index=int(sidx),
            miller_h=int(miller[0]),
            miller_k=int(miller[1]),
            miller_l=int(miller[2]),
            shift=float(shift),
            surface_side=("top" if top else "bottom"),
            # 原子统计
            n_atoms=int(len(atoms)),
            n_surface_atoms=int(n_surface_atoms),
        )
        total_saved += 1
        print(
            f"[{total_saved}/{len(selected_candidates)}] {mpid} {formula} | {row['second_element']} | hkl={miller} | {'top' if top else 'bottom'}"
        )
    # ========================================================
    # 最终结果
    # ========================================================
    print("\n" + "=" * 80)
    print("Done")
    print("=" * 80)
    print(f"找到 bulk 数量：{len(rows)}")
    print(f"枚举 surface 数量：{total_enumerated}")
    print(f"最终保存 surface 数量：{total_saved}")
    print(f"最终 DB：{surface_db_path}")
    print("=" * 80)

if __name__ == "__main__":
    # ============================================================
    # 配置
    # ============================================================
    # 表面生成参数（长度单位：Å）
    MAX_MILLER = 2
    MIN_SLAB_SIZE = 7.0
    MIN_VACUUM_SIZE = 20.0
    MIN_XY = 8.0
    SURFACE_HEIGHT_CUTOFF = 2.0

    # 几何与对称性容差
    SYMPREC = 0.1
    VORONOI_TOL = 0.1
    SLAB_TOL = 0.3

    MAX_SURFACES = 1000
    SAMPLE_PER_ELEMENT = None  # 每种第二元素的 bulk 上限；None 表示不限制3
    MAX_BULKS = None  # 最终 bulk 总数上限；None 表示不限制
    OUTPUT_DIR = r"dir_path"
    main()
