# HetCataX Database API — Reaction Search

> **API reference | Reaction**  
> Source: `react-search.md` (reorganized from the original interface export).  
> This document describes the **documented request/response contract**, not a verified live-service specification.

## 1. Endpoint

| Property | Value |
|---|---|
| Purpose | Search reaction records using IS / TS / FS structure filters, catalyst composition, reaction properties, symmetry, and calculation parameters |
| Method | `POST` |
| Path | `/self_g/r_react/` |
| Request format | JSON object |
| Documented responses | `200 OK`, `401 Unauthorized` |
| Base URL | **Not provided in the source document** |
| Authentication | An `Authorization` request header is documented; its exact value format is **not specified** |

Terminology retained from the source: **IS** (initial state), **TS** (transition state), **FS** (final state), and **Cata** (catalyst).

## 2. Request example

The example follows the original payload and adds `size`/`page`, which are described in its parameter table but absent from its request example. `page: 1` is illustrative; the indexing convention is not documented.

```json
{
  "is": [
    {
      "element": ["C"],
      "model": "at_least",
      "formula": {},
      "stru": {},
      "topo": -1000
    }
  ],
  "ts": [],
  "fs": [],
  "cata": {
    "element": ["Cu"],
    "model": "at_least",
    "topo": {}
  },
  "abstract": {
    "reactant_num": [0, 5],
    "product_num": [],
    "barrier": [0, 8],
    "enthalpy": [-10, 10]
  },
  "symmetry": {
    "cry_sys": "",
    "spa_gro": -1000,
    "miller_index": []
  },
  "params": {
    "prec": "",
    "encut": 350,
    "nelect": -1000,
    "ismear": -1000,
    "sigma": -1000,
    "ispin": -1000,
    "nupdown": -1000,
    "lhfcalc": "",
    "hfscreen": -1000,
    "algo": "",
    "gga": "",
    "ldipol": -1000,
    "efield": -1000,
    "ibrion": -1000,
    "isif": -1000,
    "isym": -1000,
    "ldau": "",
    "idvw": -1000,
    "lsol": "",
    "lreal": ""
  },
  "size": 15,
  "page": 1
}
```

**Warning:** `cata.topo: {}` above is kept **exactly as shown in the original example**, while the parameter table calls `cata.topo` an **integer**. Check the backend-supported unset value instead of assuming `{}` is valid.

### 2.1 Top-level fields

| Field | Source type | Description |
|---|---|---|
| `is` | object[] | Initial-state (IS) search conditions |
| `ts` | object[] | Transition-state (TS) search conditions |
| `fs` | object[] | Final-state (FS) search conditions |
| `cata` | object | Catalyst search conditions |
| `abstract` | object | Reactant/product counts, barrier and reaction enthalpy |
| `symmetry` | object | Crystal system, space group and Miller index |
| `params` | object | Calculation-parameter filters |
| `size` | integer | Number of table rows per page; documented default is `15` |
| `page` | integer | Requested page; page index origin is not documented |

### 2.2 `is`, `ts`, `fs` — Reaction-state search

Each array may contain state-search objects with the same set of fields.

| Field (inside each state object) | Source type | Meaning |
|---|---|---|
| `element` | string[] | Chemical element symbols; example `["C", "H"]` |
| `model` | string | Element match mode: `at_least` or `only` |
| `formula` | object | Formula-count mapping, e.g. `{"C": 2, "H": 4}` |
| `stru` | object | Structure-drawing / graph filter |
| `stru.atom_list` | string[] | Graph nodes / atoms |
| `stru.bond_list` | array | Graph edges / bonds |
| `topo` | integer | Topological ID; the source shows `-1000` for an unset value |

**Mutual exclusivity:** the source states that `element`, `formula`, `stru`, and `topo` are mutually exclusive search alternatives within a state filter. Populate the chosen criterion and leave the others in their documented unset/empty state.

The structure of individual `bond_list` entries and the query-matching behavior for multiple state objects are not specified.

### 2.3 `cata` — Catalyst search

| Field | Source type | Meaning |
|---|---|---|
| `cata.element` | string[] | Catalyst elemental composition |
| `cata.model` | string | Matching mode, as for state element filters |
| `cata.topo` | integer in field table | Catalyst topological ID; original example contains `{}` instead |

### 2.4 `abstract` and `symmetry` — Scientific filters

| Field | Source type | Meaning |
|---|---|---|
| `abstract.reactant_num` | integer[] | Number-of-reactants filter; original example `[0, 5]` |
| `abstract.product_num` | integer[] | Number-of-products filter; original example empty |
| `abstract.barrier` | number[2] | Barrier-energy interval; original description gives default `[-100, 100]` |
| `abstract.enthalpy` | number[2] | Reaction-enthalpy interval; original description gives default `[-100, 100]` |
| `symmetry.cry_sys` | string | Crystal system |
| `symmetry.spa_gro` | integer | Space group |
| `symmetry.miller_index` | array | Miller-index filter |

The original export does **not** define energy units, exact range-boundary inclusion, or whether an empty count filter disables the condition.

### 2.5 `params` — Calculation-parameter filters

The source example generally uses `-1000` as an unset numeric marker, and `""` for unused text fields. This convention should be verified for each field against the deployed API.

| Fields | Source type / example | Notes |
|---|---|---|
| `prec`, `algo`, `gga`, `lhfcalc`, `ldau`, `lsol`, `lreal` | string | Calculation settings; unused values shown as `""` |
| `encut`, `nelect`, `ismear`, `ispin`, `nupdown`, `ibrion`, `isif`, `isym` | numeric / integer | Source example uses `-1000` when unused |
| `sigma`, `hfscreen`, `efield` | number | Source table and example are not fully consistent about numeric types |
| `ldipol` | integer in example | Field table also lists `ldipol`, without a clear consistent type across references |
| `idipol` | string in request-field table | Mentioned by the table but **absent from the source request example** |
| `idvw` | integer | **Spelling in source request example**; the response uses `ivdw` |

## 3. Example request

Save the JSON example above as `react-request.json` and set your deployment's actual host and authorization value:

```bash
curl -X POST "${HETCATAX_API_BASE_URL}/self_g/r_react/" \
  -H "Content-Type: application/json" \
  -H "Authorization: ${HETCATAX_AUTH_VALUE}" \
  --data-binary @react-request.json
```

> The uploaded source does not specify whether `Authorization` should contain a bare token, `Bearer <token>`, or another scheme. Do not invent an authorization prefix.

## 4. Response

### `200 OK` — Reaction search results

**Abbreviated example** based on the first returned reaction in the source export:

```json
{
  "code": 200,
  "data": {
    "ids": [2192],
    "infos": [
      {
        "id": 2192,
        "reaction": "C<sub>2</sub> -> C<sub>2</sub>",
        "catalyst": "H<sub>112</sub>K<sub>48</sub>Cu<sub>32</sub>",
        "params": null,
        "type": "ar",
        "source": "OC20",
        "cry_sys": null,
        "spa_gro": 1,
        "miller": null,
        "encut": 350,
        "gga": "RP"
      }
    ]
  }
}
```

This is **not** a complete server response: additional records and fields from the original example are omitted. The response-schema table lists `data.count`, but the exported example does **not** include it.

### Response envelope

| Field | Type | Meaning |
|---|---|---|
| `code` | integer | Application status code, shown as `200` on success |
| `data` | object | Search results |
| `data.ids` | integer[] | IDs included in the search response |
| `data.infos` | object[] | Reaction metadata for displayed rows |
| `data.count` | integer | Count of matching results in the response schema; not present in the sample payload |

### Fields in `data.infos`

| Field / group | Meaning |
|---|---|
| `id` | Reaction record ID |
| `type` | Reaction type |
| `reaction` | Reaction text; sample contains HTML `<sub>` markup |
| `catalyst` | Catalyst formula/text; sample contains HTML `<sub>` markup |
| `source` | Dataset/source of record |
| `h`, `ea` | Fields listed in the response-schema table; **the source does not define their physical meaning or units here** |
| `cry_sys`, `spa_gro`, `miller` | Crystallographic / symmetry metadata |
| `params` | Present in the exported sample (e.g. `null`), though not separately detailed by the schema |
| `prec`, `encut`, `ismear`, `sigma`, `ispin` | Calculation parameters |
| `lhfcalc`, `hfscreen`, `algo`, `gga`, `ldipol`, `idipol`, `efield`, `ldau`, `ivdw` | Additional calculation parameters |

Parameter fields can be `null`. The source lists additional parameter names in the request that do not all appear in its response-field table.

### `401 Unauthorized`

The uploaded document shows the response body as:

```json
{}
```

No other error response bodies are documented.

## 5. Source-contract differences to confirm

The uploaded `react-search.md` contains these inconsistencies or omissions:

1. **Pagination:** `size` and `page` appear in the parameter table but **not** in the original JSON request example.
2. **Catalyst topological ID:** `cata.topo` is an object (`{}`) in the sample and an integer in the table.
3. **Mutually exclusive state filters:** `element`, `formula`, `stru` and `topo` should not be populated as independent active filters in the same state-search object.
4. **`params` field-name mismatch:** request example uses `idvw`, response uses `ivdw`; `idipol` appears in the table but not in the example.
5. **Data count:** `data.count` is documented in the response model, but the example omits it.
6. **`h` and `ea`:** appear in the response schema but not the exported example; no unit definitions are provided.
7. **Type ambiguities:** numeric parameter types and whether an empty array means “no filter” are not fully defined.
8. **Host and authentication:** neither the base URL nor the required `Authorization` format is specified.
9. **Energy ranges:** default `[-100, 100]` is stated for barrier and enthalpy, but units and interval inclusivity are unspecified.

---

**Documentation scope:** Reaction-search API only (`POST /self_g/r_react/`). The source does not define reaction-detail, trajectory-download, authentication or mutation endpoints.
