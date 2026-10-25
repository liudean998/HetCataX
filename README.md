# HetCataX — Core Catalysis Data Processing Utilities

**Selected data-processing codes for the HetCataX Catalysis Database**

HetCataX is a catalysis data platform. This repository contains **a subset of the core research and data-processing scripts** used in its workflow, including catalyst-data retrieval, materials structure acquisition, surface generation, adsorption-structure analysis, active-site extraction, bonding analysis, and atomistic descriptor generation.

**Scope:** This is a collection of research utilities, **not** the full HetCataX database, website, backend, or a turnkey data-processing pipeline. Some scripts are standalone; others are helper modules that require project-specific inputs, credentials, or paths.

- **HetCataX platform:** https://hetcatax.com/
- **Source repository:** https://github.com/liudean998/HetCataX

## Overview

The repository supports several stages of a catalysis data workflow:

1. **Retrieve catalyst data** from the HetCataX API.
2. **Obtain bulk materials** from Materials Project and construct candidate catalyst slabs.
3. **Analyze atomic structures**, including atom connectivity, catalyst/adsorbate separation, and local coordination.
4. **Extract active centers** around adsorbates from atomic structures.
5. **Generate structural descriptors**, such as many-body tensor representations (MBTR), for downstream analysis and machine learning.

These components are useful for catalytic materials research, data curation, structure–property analysis, and dataset preparation.

## Repository structure

```text
HetCataX/
├── Readme.md                        # Original repository readme (rename to README.md)
├── requirements.txt                 # Frozen Python environment dependencies
├── get_hetcatax_data.py             # Retrieve catalyst records from HetCataX
├── mp_surface.py                    # Build Cu-based binary-alloy surface structures
├── tools/
│   ├── bond.py                      # Estimate element-pair bonding ranges
│   ├── stru2mbtr.py                 # Generate MBTR structure descriptors
│   ├── DealRelax/
│   │   ├── __init__.py
│   │   ├── deal_relax.py            # Process catalyst/adsorbate structures
│   │   ├── static/
│   │   │   ├── __init__.py
│   │   │   ├── element_info.py      # Element properties and display information
│   │   │   └── vasp_config_default.py  # Default VASP-related settings
│   │   └── tools/
│   │       ├── __init__.py
│   │       ├── atoms_stru.py        # Extended ASE atomic-structure operations
│   │       └── common.py            # Shared bonding and geometry routines
│   └── GetAC/
│       ├── __init__.py
│       ├── get_ac.py                # Extract local catalytic active centers
│       ├── static/
│       │   ├── __init__.py
│       │   ├── element_info.py      # Element properties and display information
│       │   └── vasp_config_default.py  # Default VASP-related settings
│       └── tools/
│           ├── __init__.py
│           ├── atoms_stru.py        # Extended ASE atomic-structure operations
│           └── common.py            # Shared bonding and geometry routines
└── .idea/                           # PyCharm project configuration files
```

The `.idea/` directory contains editor metadata, including `.gitignore`, `HetCataX.iml`, `misc.xml`, `modules.xml`, `vcs.xml`, and inspection profile XML files. These files are not part of the scientific processing workflow.

## File and module descriptions

### Top-level scripts

| File | Purpose |
| --- | --- |
| [`get_hetcatax_data.py`](get_hetcatax_data.py) | Authenticates with the HetCataX API and retrieves catalyst records using filters such as elemental composition, formula, atom-count range, and data source. |
| [`mp_surface.py`](mp_surface.py) | Queries Cu–X binary bulk materials from Materials Project, enumerates surface slabs and terminations, samples candidate surfaces, identifies surface atoms, applies fixed-atom constraints, and writes results to an ASE database. |
| [`requirements.txt`](requirements.txt) | Records packages and versions from a development environment. It should be reviewed before use on a different machine. |
| `Readme.md` / `README.md` | Repository documentation. Keep only one canonical readme filename after updating the repository. |

### General-purpose tools

| File | Purpose |
| --- | --- |
| [`tools/bond.py`](tools/bond.py) | Provides `get_bond_range()` to estimate bonding-distance intervals for element pairs from atomic distances, covalent-radius information, and kernel density estimation. |
| [`tools/stru2mbtr.py`](tools/stru2mbtr.py) | Provides `fp_mbtr()` to calculate two-body (`k2`) and three-body (`k3`) many-body tensor representation features with DScribe. |

### `tools/DealRelax/` — catalyst and adsorption-structure processing

| File | Purpose |
| --- | --- |
| [`deal_relax.py`](tools/DealRelax/deal_relax.py) | Defines `DealRelax` for processing ASE structures into catalyst information, molecular/adsorbate information, bonding relationships, periodic structure information, and local active-site environments. |
| [`tools/atoms_stru.py`](tools/DealRelax/tools/atoms_stru.py) | Defines an ASE-based `AtomsStru` helper with atomic distances, bond/connectivity information, graph representations, and symmetry-related analysis. |
| [`tools/common.py`](tools/DealRelax/tools/common.py) | Contains reusable geometry and graph functions, including bond-range estimation, neighborhood identification, subgraph processing, and vacuum-direction detection. |
| [`static/element_info.py`](tools/DealRelax/static/element_info.py) | Stores chemical-element metadata, such as element symbols, covalent radii, and visualization properties. |
| [`static/vasp_config_default.py`](tools/DealRelax/static/vasp_config_default.py) | Supplies default VASP-related parameter settings and configuration helpers. It does not execute a VASP calculation itself. |
| `__init__.py` files | Mark the `DealRelax`, `static`, and `tools` directories as Python packages. |

### `tools/GetAC/` — active-center extraction

| File | Purpose |
| --- | --- |
| [`get_ac.py`](tools/GetAC/get_ac.py) | Defines `GetAC` to identify local catalytic active-center atoms around an adsorbate and assemble the corresponding local atomic structure and connectivity. |
| [`tools/atoms_stru.py`](tools/GetAC/tools/atoms_stru.py) | Implements ASE-based structural calculations, including interatomic distances, bonding, and graph-related operations. |
| [`tools/common.py`](tools/GetAC/tools/common.py) | Provides geometry, neighborhood, connectivity, and other helper functions used in structural analysis. |
| [`static/element_info.py`](tools/GetAC/static/element_info.py) | Stores element information used by the structure-processing routines. |
| [`static/vasp_config_default.py`](tools/GetAC/static/vasp_config_default.py) | Provides VASP-related default configuration data. |
| `__init__.py` files | Mark the `GetAC`, `static`, and `tools` directories as Python packages. |

**Implementation note:** Some utilities are duplicated between `DealRelax` and `GetAC`, and `GetAC` also imports functions from `DealRelax`. They should not be treated as fully independent packages without checking import paths.

## Key processing workflows

### 1. Retrieving catalyst records from HetCataX

`get_hetcatax_data.py` get_hetcatax_data.py is an example script demonstrating how to retrieve catalyst data from the HetCataX platform using HTTP requests. For detailed API usage and query options, refer to the [RelaxSearchAPI](./RelaxSearchAPI.md) and [ReactSearchAPI](./ReactSearchAPI.md) documentation.

```bash
python get_hetcatax_data.py
```

**Before running:** replace the placeholder credentials in the script with a secure credential-loading mechanism. A valid HetCataX account and API access may be required. The API example is tied to the current platform endpoints and should be adapted if those interfaces change.

### 2. Generating catalyst surfaces from Materials Project

`mp_surface.py` carries out a structure-preparation workflow:

1. Query ordered Cu–X binary bulk structures from Materials Project.
2. Standardize the bulk structures and enumerate symmetry-distinct Miller indices.
3. Generate candidate slab structures and possible surface terminations.
4. Sample candidate surfaces without selecting the same enumerated surface twice.
5. Expand the selected slabs in the lateral directions when required.
6. Identify surface atoms using coordination-based and height-based criteria.
7. Constrain selected non-surface atoms using ASE `FixAtoms`.
8. Save structures and associated metadata in `surfaces.db`.

The script exposes settings such as `MAX_MILLER`, `MIN_SLAB_SIZE`, `MIN_VACUUM_SIZE`, `MIN_XY`, `SURFACE_HEIGHT_CUTOFF`, and `MAX_SURFACES`.

```bash
python mp_surface.py
```

**Before running:** configure the Materials Project API credential and replace the script's placeholder `OUTPUT_DIR` with a valid output directory. The current script contains a hard-coded API key, which **must be revoked/rotated and removed from the repository and its Git history**. A safer code change is to use the already-imported `os` module:

```python
api_key = os.environ["MP_API_KEY"]
```

Then supply the value privately in the runtime environment. **This code change is not yet implemented by the README.**

**Output:** the script creates an ASE SQLite database named `surfaces.db`. It currently recreates that database on execution, so back up existing results before running.

**Surface tagging convention in this script:** `tag == 1` indicates an identified surface atom; `tag == 0` indicates an atom selected for fixation. These tags are specific to this workflow.

### 3. Analyzing catalyst–adsorbate structures

`tools/DealRelax/deal_relax.py` processes atomic configurations to derive catalyst and adsorbate information, including atomic connectivity and localized structural environments. Its helper modules use ASE structures and graph-based representations.

`tools/GetAC/get_ac.py` focuses on active-center extraction: given a labeled atomic structure, it identifies adsorbate-related neighborhoods and builds local active-center configurations. In this module, **`tag == 2` is used to identify adsorbate atoms**.

**Important:** the `mp_surface.py` surface tags and the `GetAC` adsorbate tags serve different purposes. Do not assume they are interchangeable across scripts. These modules are primarily called from Python code and require suitable labeled `ase.Atoms` inputs; they are not documented here as ready-to-run command-line applications.

### 4. Computing bonding ranges and descriptors

- `tools/bond.py` estimates element-pair bonding ranges using interatomic distances and statistical analysis.
- `tools/stru2mbtr.py` produces an MBTR feature vector from an ASE atomic structure, which may be used as an input representation for data analysis or machine-learning experiments.

## Installation

### Create an isolated environment

Python 3.12 is a reasonable starting point for recent `mp-api` and related packages; exact compatibility depends on the pinned versions you install.

```bash
conda create -n hetcatax python=3.12 -y
conda activate hetcatax
```

### Install dependencies

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

**Dependency note:** `requirements.txt` was generated from a particular development environment and may contain platform-specific or local build paths. Check and clean those entries before installation on another workstation or server. For a smaller installation, start from the modules you intend to run.

Major libraries used by the project include:

- `ase` — atomic structures, constraints, and databases
- `pymatgen` and `mp-api` — crystal structures, surface generation, and Materials Project access
- `numpy` and `scipy` — numerical and scientific computing
- `networkx` — atomic connectivity and graph processing
- `spglib` — symmetry-related functionality
- `scikit-learn` — kernel density estimation
- `dscribe` — structure descriptors
- `requests` — HTTP API access

Verify core imports after installation:

```bash
python -c "import ase, pymatgen, networkx, dscribe; from mp_api.client import MPRester; print('Imports OK')"
```

Some utility modules use package-style imports; you may need to run them from an appropriate working directory or adjust `PYTHONPATH`.

## Data handling and reproducibility

- **Credentials:** Never commit passwords, authorization tokens, or API keys. Use environment variables or a secret manager. If a credential has appeared in a public commit, revoke it and address exposure in Git history.
- **Input assumptions:** Several analysis functions expect specific atomic tags, periodic boundary conditions, or valid ASE structure objects.
- **Output management:** Surface-generation runs can overwrite the target database.
- **Random sampling:** The current surface workflow samples candidate materials/surfaces without a configured reproducible random seed.
- **Experimental scope:** These are research scripts and helpers; validation, error handling, and environment portability may need additional work before deployment in automated production pipelines.

## Related resources

- [HetCataX](https://hetcatax.com/)
- [Materials Project](https://materialsproject.org/)
- [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/)
- [Pymatgen](https://pymatgen.org/)
- [DScribe](https://singroup.github.io/dscribe/)

---

*This repository documents selected core data-processing implementations associated with the HetCataX catalysis database. It is not a complete distribution of the HetCataX platform.*
