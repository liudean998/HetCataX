# HetCataX Database API — Structure Search

> **API reference | Intermediate · Molecule · Catalyst**  
> Source: `relax-search.md` (reorganized from the original interface export).  
> This document describes the **documented request/response contract**, not a verified live-service specification.

## 1. Endpoint

| Property | Value |
|---|---|
| Purpose | Search structures in the Intermediate / Molecule / Catalyst database modules |
| Method | `POST` |
| Path | `/self_g/r_relax/` |
| Request format | JSON object |
| Documented responses | `200 OK`, `401 Unauthorized` |
| Base URL | **Not provided in the source document** |
| Authentication | An `Authorization` request header is documented; the precise value format is **not specified** |

A caller must obtain the correct API base URL and authorization value from the HetCataX deployment configuration. Do not assume that the website URL is also the API host.

## 2. Request example

The example below follows the source request structure. `page` has been added to illustrate the documented pagination field; the source example omits it. A value of `1` is **illustrative only**, because the page-numbering convention is not specified.

```json
{
  "base": {
    "model": 0,
    "element": ["C"],
    "formula": {},
    "id": -1000
  },
  "graph": {
    "atom_list": [],
    "atom_tag": [],
    "bond_list": []
  },
  "abstract": {
    "num_of_atoms": [100, 100],
    "source": "",
    "relax_type": 0
  },
  "symmetry": {
    "spa_gro": -1000,
    "cry_sys": "",
    "miller_index": []
  },
  "topo": {
    "mol_mid": -1000,
    "cata_mid": -1000,
    "as_mid": -1000
  },
  "params": {
    "prec": "",
    "encut": -1000,
    "nelect": -1000,
    "ismear": -1000,
    "sigma": -1000,
    "ispin": -1000,
    "nupdown": -1000,
    "lhfcalc": "",
    "hfscreen": -1000,
    "algo": "",
    "gga": "",
    "ldipol": "",
    "idipol": -1000,
    "efield": -1000,
    "ibrion": -1000,
    "isif": -1000,
    "isym": -1000,
    "ldau": "",
    "idvw": -1000,
    "lsol": "",
    "lreal": ""
  },
  "size": 10,
  "page": 1,
  "model_type": "system"
}
```

### 2.1 Top-level fields

| Field | Source type | Description                                                                                                                                                                                                |
|---|---|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `base` | object | Search by chemical elements, formula, or database ID                                                                                                                                                       |
| `graph` | object | Graph-structure query (atoms, tags and bonds)                                                                                                                                                              |
| `abstract` | object | Atom-count, provenance and relaxation-type filters                                                                                                                                                         |
| `symmetry` | object | Space-group, crystal-system and Miller-index filters                                                                                                                                                       |
| `topo` | object | Molecule, catalyst and active-site topological IDs                                                                                                                                                         |
| `params` | object | Electronic-structure / calculation-parameter filters                                                                                                                                                       |
| `size` | integer | Number of table entries per page; source description says default **15**, while its example uses **10**                                                                                                    |
| `page` | integer | Requested page number; source does not define whether counting begins at 0 or 1                                                                                                                            |
| `model_type` | string | Structure module: `system`, `molecule`, or **`catslyst`** as spelled in the original specification (see discrepancies below). Note that `system` corresponds to the `intermediate` module in the database. |

### 2.2 `base` — Basic chemical search

| Field | Type | Description |
|---|---|---|
| `base.model` | integer | `1` when **Only** is selected; `0` for other modes (including the default **At Least** mode) |
| `base.element` | string[] | Element symbols, such as `["C", "H"]` |
| `base.formula` | object | Formula as elemental counts, e.g. `{"C": 2, "H": 4}` |
| `base.id` | integer | Database ID for the ID search mode; sample unset value: `-1000` |

**Search-mode rule:** according to the source, `element`, `formula` and `id` represent alternative search inputs. Populate the input corresponding to the selected mode and leave the other inputs empty/unset using values accepted by the deployed API.

### 2.3 `graph` — Structure-graph filter

| Field | Type | Description |
|---|---|---|
| `graph.atom_list` | string[] | Graph node list |
| `graph.atom_tag` | integer[] | Attributes/tags associated with the nodes |
| `graph.bond_list` | array | Graph edge/bond list |

The precise shape of each bond entry, node-tag semantics and matching algorithm are **not described** in the source document. The empty arrays shown above correspond to an unused graph filter.

### 2.4 `abstract`, `symmetry`, and `topo`

| Field | Type | Description |
|---|---|---|
| `abstract.num_of_atoms` | integer[2] | Inclusive/exclusive boundary semantics are unspecified; source describes it as `[minimum, maximum]` |
| `abstract.source` | string | Data source or provenance filter |
| `abstract.relax_type` | integer in example | Relaxation type; source mentions `"slab"` and `"bulk"` but does **not** define the mapping to numeric values |
| `symmetry.spa_gro` | integer | Space-group filter |
| `symmetry.cry_sys` | string | Crystal-system filter |
| `symmetry.miller_index` | array in example | Miller-index filter; the parameter table inconsistently calls it a string |
| `topo.mol_mid` | integer | Molecule topological ID |
| `topo.cata_mid` | integer | Catalyst topological ID |
| `topo.as_mid` | integer | Active-site topological ID |

### 2.5 `params` — Calculation-parameter filters

These names are preserved from the original request example. The example generally uses `-1000` for unset **numeric** parameters and `""` for unset **text** parameters; the exact sentinel rules should be confirmed with the deployed API.

| Fields | Types indicated by source | Notes |
|---|---|---|
| `prec`, `algo`, `gga` | string | Calculation settings |
| `encut`, `nelect`, `ismear`, `ispin`, `ibrion`, `isif`, `isym` | numeric / integer | Source example uses `-1000` for unset values |
| `sigma`, `nupdown`, `hfscreen`, `efield` | number | Some are labeled integer in the parameter table; validate if non-integer values are needed |
| `lhfcalc`, `ldipol`, `ldau`, `lsol`, `lreal` | string | Source example uses empty strings |
| `idipol` | integer | Calculation parameter |
| `idvw` | integer | **Spelling in request example**; parameter/response documentation also uses `ivdw` |

## 3. Example request

Save the JSON above as `relax-request.json`. Then use the correct base URL and credentials for your deployment:

```bash
curl -X POST "${HETCATAX_API_BASE_URL}/self_g/r_relax/" \
  -H "Content-Type: application/json" \
  -H "Authorization: ${HETCATAX_AUTH_VALUE}" \
  --data-binary @relax-request.json
```

> The source only documents a header named `Authorization`. It does **not** specify `Bearer`, token prefix, cookie use or a login endpoint. Set `HETCATAX_AUTH_VALUE` to the value required by your deployment.

## 4. Response

### `200 OK` — Search results

**Abbreviated example** adapted from an actual record in the source export; other documented fields and additional returned records are omitted for readability.

```json
{
  "code": 200,
  "data": {
    "infos": [
      {
        "id": 30375,
        "model_id": 34916,
        "num_of_atoms": 100,
        "source": "OC20",
        "formula": "CH<sub>2</sub>OV<sub>48</sub>Cu<sub>16</sub>Ta<sub>32</sub>",
        "encut": 350,
        "gga": "RP",
        "as_mid": 34916,
        "cata_mid": [34915],
        "mol_mid": [34914],
        "spa_gro": 1,
        "miller": [2, 2, 1],
        "model_type": "System"
      }
    ],
    "ids": [30375]
  }
}
```

The original exported sample includes many more entries in both `infos` and `ids`. Its schema additionally documents `data.count`, but that field is **absent from the shown response example**.

### Response envelope

| Field | Type | Meaning |
|---|---|---|
| `code` | integer | Application code (`200` on documented success) |
| `data` | object | Search response body |
| `data.infos` | object[] | Information for the displayed result records |
| `data.ids` | integer[] | List of matching data IDs shown in the response example |
| `data.count` | integer | Total matching records, according to the response-schema table; **not present in example** |

### Result fields in `data.infos`

| Fields | Meaning / documented presentation |
|---|---|
| `id`, `model_id` | Record ID and structure/model ID |
| `num_of_atoms` | Atom count |
| `formula` | Chemical formula; sample contains HTML `<sub>` notation |
| `model_type`, `source` | Record type and data source |
| `cry_sys`, `spa_gro`, `miller` | Crystal system, space group and Miller index |
| `as_mid`, `cata_mid`, `mol_mid` | Active-site, catalyst and molecule topological IDs |
| `fr_energy`, `zero_energy`, `wo_entrp` | Thermodynamic / energy fields listed in the original schema; specific units and definitions **not supplied** |
| `prec`, `encut`, `nelect`, `ismear`, `sigma`, `ispin`, `nupdown` | Calculation parameters |
| `lhfcalc`, `hfscreen`, `algo`, `gga`, `ldipol`, `idipol`, `efield` | Calculation parameters |
| `ibrion`, `isif`, `isym`, `ldau`, `ivdw`, `lsol`, `lreal` | Calculation parameters |

Many calculation fields can be `null` in records. The original output example represents `cata_mid` and `mol_mid` as **arrays**, and `miller` as a **three-element array**, even where the schema table gives other types.

### `401 Unauthorized`

The source provides the following response body:

```json
{}
```

The documentation does not define other HTTP response codes or a structured error-message schema.

## 5. Source-contract differences to confirm

The following differences are present in the uploaded `relax-search.md` and should be clarified **before generating typed SDKs or enforcing strict validation**:

1. **`model_type` spelling:** the request-table text lists `catslyst` (not `catalyst`); output examples use `System` (capitalized). Actual accepted values and case sensitivity need verification.
2. **Pagination:** the request sample uses `size: 10`, while the table says default `15`. The sample omits `page` although the table marks it required.
3. **`abstract.relax_type`:** sample value `0` is numeric, but description names `"slab"` and `"bulk"` without defining a numeric mapping.
4. **`symmetry.miller_index`:** example is an array; table states string. Output `miller` is an array in the example.
5. **`idvw` vs `ivdw`:** request uses `idvw`; the response and some parameter descriptions use `ivdw`.
6. **Topological-ID types:** example output uses arrays for `cata_mid`/`mol_mid` and integer for `as_mid`, although the schema table gives string/null.
7. **`data.count`:** documented in the response schema, omitted from the exported example.
8. **`params` completeness:** the sample includes fields that are not uniformly listed in the parameter table; use the actual backend schema as the final authority.
9. **Authentication and host:** authorization value format, base URL, and access-control policy are not specified in the source.

---

**Documentation scope:** Structure-search API only (`POST /self_g/r_relax/`). No structure-detail, export, authentication or mutation endpoints are defined by this source.
