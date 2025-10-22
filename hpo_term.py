#!/usr/bin/env python3
"""
Map HPO terms to system-level categories (direct children of HP:0000118).

Usage examples:
  # Plain text list (one per line: HPO ID or exact term name)
  python hpo_system_mapper.py --in terms.txt --out mapped.csv

  # CSV with a specific column holding terms/IDs
  python hpo_system_mapper.py --in phenos.csv --col phenotype --out mapped.csv

  # Provide a local hp.obo (otherwise it will download)
  python hpo_system_mapper.py --in terms.txt --obo hp.obo --out mapped.csv
"""

import argparse
import csv
import os
import sys
import re
from typing import Dict, List, Set, Tuple, Optional

try:
    import obonet  # pip install obonet networkx
    import networkx as nx
except ImportError as e:
    sys.stderr.write("This script requires 'obonet' and 'networkx'. Install with:\n"
                     "  pip install obonet networkx\n")
    raise

import urllib.request

HP_OBO_URL = "https://purl.obolibrary.org/obo/hp.obo"
PHENOTYPIC_ABN = "HP:0000118"

def fetch_hp_obo(obo_path: str) -> str:
    """Download hp.obo if not present."""
    if os.path.exists(obo_path):
        return obo_path
    sys.stderr.write(f"Downloading HPO ontology to '{obo_path}'...\n")
    urllib.request.urlretrieve(HP_OBO_URL, obo_path)
    return obo_path

def load_hpo_graph(obo_path: str) -> nx.MultiDiGraph:
    """Load HPO as a NetworkX graph."""
    return obonet.read_obo(obo_path)

def build_name_index(G: nx.MultiDiGraph) -> Dict[str, str]:
    """
    Build a case-insensitive index from exact primary term name -> HPO ID.
    (You can extend this to alt_ids or synonyms if you want fuzzier matching.)
    """
    name_to_id = {}
    for node, data in G.nodes(data=True):
        name = data.get("name")
        if name:
            name_to_id[name.lower()] = node
    return name_to_id

def get_system_roots(G: nx.MultiDiGraph, restrict_to: Optional[Set[str]] = None) -> Dict[str, str]:
    """
    Get the direct children of HP:0000118 (Phenotypic abnormality) as system-level roots.
    Returns dict {id: name}. Optionally supply 'restrict_to' as a set of IDs to keep.
    """
    if PHENOTYPIC_ABN not in G:
        raise ValueError("HP:0000118 not found in the ontology graph.")

    children = set()
    # Edges are from child -> parent in OBO graphs produced by obonet
    for child, parent, k in G.out_edges(PHENOTYPIC_ABN, keys=True):
        # Careful: obonet uses parent <- child edges; we want children of PHENOTYPIC_ABN,
        # which means nodes with an edge (child -> PHENOTYPIC_ABN). The code below
        # should traverse incoming edges to PHENOTYPIC_ABN.
        pass

    # The above note indicates we should actually look at in_edges of PHENOTYPIC_ABN
    children = set(u for (u, v, k) in G.in_edges(PHENOTYPIC_ABN, keys=True))

    systems = {}
    for cid in children:
        if restrict_to and cid not in restrict_to:
            continue
        nm = G.nodes[cid].get("name", "")
        systems[cid] = nm
    return systems

def resolve_input_to_id(term: str, name_to_id: Dict[str, str]) -> Optional[str]:
    """Resolve an input string to an HPO ID, using either ID recognition or exact name match."""
    t = term.strip()
    if not t:
        return None
    # Looks like an HPO ID?
    if re.match(r"^HP:\d{7}$", t, re.IGNORECASE):
        # Normalize uppercase prefix
        return "HP:" + t.split(":")[1].zfill(7)
    # Else try exact name match (case-insensitive)
    return name_to_id.get(t.lower())

def ancestors(G: nx.MultiDiGraph, node: str) -> Set[str]:
    """Return all ancestors (all superclasses) of a node following is_a/part_of relationships."""
    if node not in G:
        return set()
    # Move upward: follow edges from child -> parent.
    # In obonet graphs, edges are directed from child to parent.
    anc = nx.ancestors(G, node)
    return anc

def map_term_to_systems(G: nx.MultiDiGraph,
                        systems: Dict[str, str],
                        node_id: str) -> List[Tuple[str, str]]:
    """Return list of (system_id, system_name) for a given HPO node."""
    if node_id not in G:
        return []
    ancs = ancestors(G, node_id)
    hits = [(sid, systems[sid]) for sid in systems.keys() if sid in ancs or sid == node_id]
    # Sort by name for stable output
    hits.sort(key=lambda x: x[1].lower())
    return hits

def read_inputs(path: str, column: Optional[str]) -> List[str]:
    """
    Read input terms/IDs. If column is None, treat as plain text, one per line.
    If CSV/TSV with column name, auto-detect delimiter by sniffing.
    """
    if column is None:
        with open(path, 'r', encoding='utf-8') as f:
            return [line.strip() for line in f if line.strip()]
    # CSV/TSV with named column
    with open(path, 'r', encoding='utf-8') as f:
        sample = f.read(2048)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=',\t;')
        reader = csv.DictReader(f, dialect=dialect)
        if column not in reader.fieldnames:
            raise ValueError(f"Column '{column}' not found. Available: {reader.fieldnames}")
        return [row[column].strip() for row in reader if row.get(column, "").strip()]

def write_output(out_path: str,
                 rows: List[Dict[str, str]]) -> None:
    fieldnames = ["input", "resolved_hpo_id", "resolved_hpo_name",
                  "system_level_ids", "system_level_names"]
    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)

def main():
    ap = argparse.ArgumentParser(description="Map HPO terms/IDs to system-level categories.")
    ap.add_argument("--in", dest="in_path", required=True, help="Input file (txt/CSV/TSV).")
    ap.add_argument("--col", dest="column", default=None,
                    help="Column name if input is CSV/TSV. If omitted, expects one item per line.")
    ap.add_argument("--obo", dest="obo_path", default="hp.obo",
                    help="Path to hp.obo (will download if missing).")
    ap.add_argument("--out", dest="out_path", required=True, help="Output CSV path.")
    ap.add_argument("--restrict_ids", dest="restrict_ids", default=None,
                    help="Optional comma-separated list of system root IDs to keep (e.g., "
                         "\"HP:0000707,HP:0001626,...\").")
    args = ap.parse_args()

    obo_path = fetch_hp_obo(args.obo_path)
    G = load_hpo_graph(obo_path)
    name_to_id = build_name_index(G)

    restrict: Optional[Set[str]] = None
    if args.restrict_ids:
        restrict = set([x.strip() for x in args.restrict_ids.split(",") if x.strip()])

    system_roots = get_system_roots(G, restrict_to=restrict)

    inputs = read_inputs(args.in_path, args.column)

    results = []
    for raw in inputs:
        resolved_id = resolve_input_to_id(raw, name_to_id)
        resolved_name = G.nodes[resolved_id]["name"] if resolved_id and resolved_id in G else None

        systems = map_term_to_systems(G, system_roots, resolved_id) if resolved_id else []
        sys_ids = "|".join([sid for sid, _ in systems]) if systems else ""
        sys_names = "|".join([sn for _, sn in systems]) if systems else ""

        results.append({
            "input": raw,
            "resolved_hpo_id": resolved_id or "",
            "resolved_hpo_name": resolved_name or "",
            "system_level_ids": sys_ids,
            "system_level_names": sys_names
        })

    write_output(args.out_path, results)
    sys.stderr.write(f"Wrote {len(results)} rows to {args.out_path}\n")

if __name__ == "__main__":
    main()
