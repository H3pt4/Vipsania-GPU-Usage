#!/usr/bin/env python3
"""
Genome lengths of a clade, from NCBI Datasets.

Given a taxon ID, a taxon level and a minimum assembly level, this script:

  1. pulls all assembly reports below the taxon (NCBI `datasets` CLI)
  2. filters them  (atypical removed, min. assembly level, min. length,
                    max. number of scaffolds)
  3. keeps ONE assembly per species (reference > representative > highest
     level > RefSeq > annotated > newest)
  4. counts the unique taxa at the chosen level (species/genus/family/...)
  5. returns / saves the list of genome lengths (one per species)

Filtering happens BEFORE choosing the assembly per species, so a species is only
lost if none of its assemblies passes the filters.

Command line:
    python clade_genome_lengths.py 9443 -l family -a chromosome
    python clade_genome_lengths.py 7215 -l genus -a "complete genome" -o dros

Outputs (prefix defaults to clade_<taxid>):
    <prefix>_lengths.txt   one genome length (bp) per line
    <prefix>_table.tsv     table of the selected assemblies + metadata
    + a summary printed to stdout

From Python:
    from clade_genome_lengths import analyse_clade
    res = analyse_clade("9443", taxon_level="family", min_assembly_level="chromosome")
    res["n_taxa"]      # number of unique families
    res["lengths"]     # list of genome lengths (bp)
    res["rows"]        # list of dicts (the table)

Requirements: NCBI `datasets` CLI in PATH, Python >= 3.8, no other packages.
"""

import argparse
import csv
import json
import shutil
import statistics
import subprocess
import sys
from typing import Any, Dict, Iterable, List, Optional

LEVEL_RANK = {"contig": 1, "scaffold": 2, "chromosome": 3, "complete genome": 4}
CATEGORY_RANK = {"reference genome": 2, "representative genome": 1}
TAXON_LEVELS = ["species", "genus", "family", "order", "class", "phylum"]
LINEAGE_COLUMNS = ["genus", "family", "order", "class", "phylum"]

# COLUMNS = (
#     ["species_taxid", "species_name"] + LINEAGE_COLUMNS + [
#         "organism_name", "n_assemblies_passing_filters", "accession",
#         "source_database", "refseq_category", "assembly_name", "assembly_level",
#         "release_date", "bioproject", "biosample",
#         "total_sequence_length", "gc_percent", "contig_n50", "scaffold_n50",
#         "number_of_contigs", "number_of_scaffolds", "annotated",
#     ]
# )

COLUMNS = (
    ["species_taxid", "species_name"] + LINEAGE_COLUMNS + [
        "organism_name", "n_assemblies_passing_filters", "accession", "refseq_category", "assembly_name", "assembly_level",
        "release_date", "total_sequence_length", "annotated",
    ]
)


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def get(d: Any, path: str, default: Any = "") -> Any:
    for key in path.split("."):
        if not isinstance(d, dict) or key not in d:
            return default
        d = d[key]
    return d if d is not None else default


def num(value: Any) -> Optional[float]:
    try:
        return float(str(value).replace(",", "").strip())
    except ValueError:
        return None


def run_datasets(args: List[str], api_key: Optional[str]) -> Iterable[dict]:
    cmd = ["datasets"] + args + ["--as-json-lines"]
    if api_key:
        cmd += ["--api-key", api_key]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        line = line.strip()
        if line:
            yield json.loads(line)
    err = proc.stderr.read() if proc.stderr else ""
    if proc.wait() != 0:
        raise RuntimeError(f"`{' '.join(cmd)}` failed:\n{err}")


def lineage_lookup(taxids: List[str], api_key: Optional[str]) -> Dict[str, Dict[str, Dict[str, str]]]:
    """taxid -> {rank: {"id":..., "name":...}} for every rank in the classification."""
    out: Dict[str, Dict[str, Dict[str, str]]] = {}
    for i in range(0, len(taxids), 200):
        batch = taxids[i:i + 200]
        for rep in run_datasets(["summary", "taxonomy", "taxon"] + batch, api_key):
            tax = rep.get("taxonomy", {})
            tid = str(tax.get("tax_id", ""))
            if not tid:
                continue
            cl = tax.get("classification", {})
            out[tid] = {rank: {"id": str(v.get("id", "")), "name": v.get("name", "")}
                        for rank, v in cl.items() if isinstance(v, dict)}
    return out


def passes_filters(rec: dict, min_rank: int, min_length: float, max_scaffolds: int,
                   keep_atypical: bool) -> bool:
    ai = rec.get("assembly_info", {})
    if not keep_atypical and get(ai, "atypical.is_atypical") is True:
        return False
    if LEVEL_RANK.get(str(ai.get("assembly_level", "")).lower(), 0) < min_rank:
        return False
    length = num(get(rec, "assembly_stats.total_sequence_length"))
    if length is None or length < min_length:
        return False
    n_scaf = num(get(rec, "assembly_stats.number_of_scaffolds"))
    if n_scaf is None:
        n_scaf = num(get(rec, "assembly_stats.number_of_contigs"))
    if n_scaf is None or n_scaf > max_scaffolds:
        return False
    return True


def sort_key(rec: dict):
    ai = rec.get("assembly_info", {})
    return (
        CATEGORY_RANK.get(str(ai.get("refseq_category", "")).lower(), 0),
        LEVEL_RANK.get(str(ai.get("assembly_level", "")).lower(), 0),
        1 if rec.get("accession", "").startswith("GCF_") else 0,
        1 if rec.get("annotation_info") else 0,
        str(ai.get("release_date", "")),
    )


def to_row(rec: dict, sp_id: str, sp_name: str, lineage: Dict[str, Dict[str, str]],
           n: int) -> Dict[str, Any]:
    ai, st = "assembly_info", "assembly_stats"
    row = {
        "species_taxid": sp_id,
        "species_name": sp_name,
        "organism_name": get(rec, "organism.organism_name"),
        "n_assemblies_passing_filters": n,
        "accession": rec.get("accession", ""),
        "refseq_category": get(rec, f"{ai}.refseq_category"),
        "assembly_name": get(rec, f"{ai}.assembly_name"),
        "assembly_level": get(rec, f"{ai}.assembly_level"),
        #"release_date": get(rec, f"{ai}.release_date"),
        "total_sequence_length": int(num(get(rec, f"{st}.total_sequence_length")) or 0),
        #"gc_percent": get(rec, f"{st}.gc_percent"),
        #"contig_n50": get(rec, f"{st}.contig_n50"),
        #"scaffold_n50": get(rec, f"{st}.scaffold_n50"),
        #"number_of_contigs": get(rec, f"{st}.number_of_contigs"),
        #"number_of_scaffolds": get(rec, f"{st}.number_of_scaffolds"),
        "annotated": "yes" if rec.get("annotation_info") else "no",
    }
    for rank in LINEAGE_COLUMNS:
        row[rank] = lineage.get(rank, {}).get("name", "")
    return row


# --------------------------------------------------------------------------- #
# main logic (importable)
# --------------------------------------------------------------------------- #
def analyse_clade(
    taxon: str,
    taxon_level: str = "species",
    min_assembly_level: str = "chromosome",
    min_length: float = 2_000_000,
    max_scaffolds: int = 20_000,
    keep_atypical: bool = False,
    refseq_only: bool = False,
    api_key: Optional[str] = None,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Return {"n_taxa", "taxa", "lengths", "rows"} for the clade."""
    taxon_level = taxon_level.lower()
    min_assembly_level = min_assembly_level.lower()
    if taxon_level not in TAXON_LEVELS:
        raise ValueError(f"taxon_level must be one of {TAXON_LEVELS}")
    if min_assembly_level not in LEVEL_RANK:
        raise ValueError(f"min_assembly_level must be one of {list(LEVEL_RANK)}")
    if shutil.which("datasets") is None:
        raise RuntimeError("NCBI `datasets` CLI not found in PATH.")

    log = (lambda m: print(m, file=sys.stderr)) if verbose else (lambda m: None)

    # 1) all assembly reports
    cmd = ["summary", "genome", "taxon", str(taxon)]
    if refseq_only:
        cmd += ["--assembly-source", "RefSeq"]
    log(f"Fetching assembly reports for taxon {taxon} ...")
    records = list(run_datasets(cmd, api_key))
    log(f"  {len(records)} assemblies retrieved")

    # 2) filters
    min_rank = LEVEL_RANK[min_assembly_level]
    records = [r for r in records
               if passes_filters(r, min_rank, min_length, max_scaffolds, keep_atypical)]
    log(f"  {len(records)} assemblies pass the filters")
    if not records:
        return {"n_taxa": 0, "taxa": [], "lengths": [], "rows": []}

    # 3) resolve lineage, group by species, choose best assembly
    org_ids = sorted({str(get(r, "organism.tax_id")) for r in records
                      if get(r, "organism.tax_id")})
    log(f"Resolving lineage for {len(org_ids)} taxids ...")
    lineages = lineage_lookup(org_ids, api_key)

    groups: Dict[str, List[dict]] = {}
    info: Dict[str, Dict[str, Any]] = {}
    for r in records:
        tid = str(get(r, "organism.tax_id"))
        lin = lineages.get(tid, {})
        sp = lin.get("species", {"id": tid, "name": get(r, "organism.organism_name")})
        groups.setdefault(sp["id"], []).append(r)
        info[sp["id"]] = {"name": sp["name"], "lineage": lin}

    rows: List[Dict[str, Any]] = []
    for sp_id, recs in groups.items():
        best = max(recs, key=sort_key)
        rows.append(to_row(best, sp_id, info[sp_id]["name"], info[sp_id]["lineage"], len(recs)))
    rows.sort(key=lambda r: r["species_name"])

    # 4) taxa at the chosen level
    def level_name(row: Dict[str, Any]) -> str:
        return row["species_name"] if taxon_level == "species" else row.get(taxon_level, "")

    taxa = sorted({level_name(r) for r in rows if level_name(r)})
    n_missing = sum(1 for r in rows if not level_name(r))
    if n_missing:
        log(f"  note: {n_missing} species have no '{taxon_level}' in their lineage")

    # 5) genome lengths (one per species)
    lengths = [r["total_sequence_length"] for r in rows]
    return {"n_taxa": len(taxa), "taxa": taxa, "lengths": lengths, "rows": rows}


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("taxon", help="NCBI taxonomy ID (or scientific name) of the clade")
    ap.add_argument("-l", "--taxon-level", default="species", choices=TAXON_LEVELS,
                    help="level at which unique taxa are counted (default: species)")
    ap.add_argument("-a", "--min-assembly-level", default="chromosome",
                    choices=list(LEVEL_RANK), help="minimum assembly level (default: chromosome)")
    ap.add_argument("--min-length", type=float, default=2_000_000,
                    help="minimum genome length in bp (default: 2000000)")
    ap.add_argument("--max-scaffolds", type=int, default=20_000,
                    help="maximum number of scaffolds (default: 20000)")
    ap.add_argument("--keep-atypical", action="store_true", help="do not discard atypical assemblies")
    ap.add_argument("--refseq-only", action="store_true", help="only consider RefSeq assemblies")
    ap.add_argument("--api-key", help="NCBI API key")
    ap.add_argument("-o", "--prefix", help="output file prefix (default: clade_<taxon>)")
    args = ap.parse_args()

    try:
        res = analyse_clade(args.taxon, args.taxon_level, args.min_assembly_level,
                            args.min_length, args.max_scaffolds, args.keep_atypical,
                            args.refseq_only, args.api_key)
    except (RuntimeError, ValueError) as e:
        sys.exit(f"ERROR: {e}")

    if not res["rows"]:
        sys.exit("No assemblies left after filtering.")

    prefix = args.prefix or f"clade_{str(args.taxon).replace(' ', '_')}"
    with open(f"{prefix}_lengths.txt", "w") as fh:
        fh.write("\n".join(str(x) for x in res["lengths"]) + "\n")
    with open(f"{prefix}_table.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(res["rows"])

    L = res["lengths"]
    print(f"Clade                : {args.taxon}")
    print(f"Unique {args.taxon_level:<14}: {res['n_taxa']}")
    print(f"Genomes (1/species)  : {len(L)}")
    print(f"Length min/median/max: {min(L):,} / {int(statistics.median(L)):,} / {max(L):,} bp")
    print(f"Saved: {prefix}_lengths.txt, {prefix}_table.tsv")


if __name__ == "__main__":
    main()