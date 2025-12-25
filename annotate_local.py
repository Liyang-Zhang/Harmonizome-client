#!/usr/bin/env python3
"""
Annotate genes using local Harmonizome edge files defined in a YAML/JSON config.

Key features
- Reads a config listing datasets (edge file paths, gene column, fields to keep).
- Supports input genes via CLI or a TSV file containing a column named `gene` or `基因`
  (other columns are preserved, annotations are appended to the right).
- Each dataset outputs one JSON column + one count column by default to avoid column explosion.

Usage examples
--------------
python annotate_local.py --genes EGFR TP53 KRAS ERBB2 \
  --config config/local_datasets.yml --out annotations.tsv

python annotate_local.py \
  --genes-file genes.tsv \  # TSV with a header containing gene/基因列
  --config config/local_datasets.yml \
  --out annotations.tsv

Config format (YAML or JSON, simplified)
----------------------------------------
{
  "datasets": [
    {
      "name": "GTEx_Tissue",
      "path": "downloaded_edges/GTEx_Tissue_Gene_Expression_Profiles/gene_attribute_edges.txt.gz",
      "format": "tsv",
      "gene_field": "source",
      "fields": ["target", "weight"],
      "output": {
        "mode": "json",          # json | join
        "count_field": "count"   # optional count column name
      }
    }
  ]
}
"""
from __future__ import annotations
import argparse
import csv
import gzip
import json
from pathlib import Path
from typing import Dict, List, Any, Tuple

try:
    import yaml  # type: ignore
except Exception:  # pragma: no cover - fallback if PyYAML missing
    yaml = None


def load_config(path: Path) -> Dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() in {".yaml", ".yml"}:
        if yaml is None:
            raise RuntimeError("PyYAML not installed; cannot read YAML config")
        return yaml.safe_load(text)
    return json.loads(text)


def read_edges(path: Path, gene_field: str) -> Dict[str, List[Dict[str, Any]]]:
    gene_to_rows: Dict[str, List[Dict[str, Any]]] = {}
    opener = gzip.open if path.suffix.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            g = row.get(gene_field)
            if not g:
                continue
            gene_to_rows.setdefault(g.strip(), []).append(row)
    return gene_to_rows


def aggregate_rows(rows: List[Dict[str, Any]], fields: List[str]) -> List[Dict[str, Any]]:
    agg: List[Dict[str, Any]] = []
    for r in rows:
        agg.append({k: r.get(k, "") for k in fields})
    return agg


def annotate(genes: List[str], config: Dict[str, Any]) -> List[Dict[str, Any]]:
    datasets = config.get("datasets", [])
    # Preload all datasets into memory (ok for current test sizes)
    loaded = []
    for ds in datasets:
        name = ds["name"]
        path = Path(ds["path"])
        gene_field = ds["gene_field"]
        fields = ds.get("fields") or []
        output_cfg = ds.get("output", {})
        rows_map = read_edges(path, gene_field)
        loaded.append((name, rows_map, fields, output_cfg))

    results: List[Dict[str, Any]] = []
    for g in genes:
        row: Dict[str, Any] = {"gene": g}
        for name, rows_map, fields, output_cfg in loaded:
            hits = rows_map.get(g, [])
            mode = (output_cfg.get("mode") or "json").lower()
            count_field = output_cfg.get("count_field") or f"{name}_count"
            if mode == "json":
                values = aggregate_rows(hits, fields) if hits else []
                row[f"{name}_json"] = json.dumps(values, ensure_ascii=False)
                row[count_field] = len(hits)
            elif mode == "join":
                join_field = output_cfg.get("join_field") or (fields[0] if fields else "")
                sep = output_cfg.get("sep", "|")
                row[f"{name}_joined"] = sep.join([h.get(join_field, "") for h in hits]) if hits else ""
                row[count_field] = len(hits)
            else:
                row[f"{name}_json"] = json.dumps([], ensure_ascii=False)
                row[count_field] = 0
        results.append(row)
    return results


def write_tsv(rows: List[Dict[str, Any]], path: Path) -> None:
    if not rows:
        return
    fields = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def read_input_table(path: Path) -> Tuple[List[Dict[str, Any]], str, List[str]]:
    """Read TSV with a gene column named 'gene' or '基因' (case-insensitive)."""
    rows: List[Dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError("Input file has no header")
        # detect gene column
        gene_col = None
        for col in reader.fieldnames:
            if col.lower() == "gene" or col == "基因":
                gene_col = col
                break
        if gene_col is None:
            raise RuntimeError("Input file must have a 'gene' or '基因' column")
        for row in reader:
            if not row.get(gene_col):
                continue
            rows.append(row)
    return rows, gene_col, reader.fieldnames or []


def main() -> int:
    ap = argparse.ArgumentParser(description="Annotate genes using local Harmonizome edge files")
    ap.add_argument("--genes", nargs="*", help="Gene symbols (space-separated)")
    ap.add_argument(
        "--genes-file",
        help="TSV file with a header containing 'gene' or '基因'; other columns are preserved",
    )
    ap.add_argument("--config", default="config/local_datasets.yml", help="Config YAML/JSON path")
    ap.add_argument("--out", default="annotations.tsv", help="Output TSV path")
    args = ap.parse_args()

    base_rows: List[Dict[str, Any]] = []
    gene_col = "gene"
    base_fields: List[str] = [gene_col]
    if args.genes_file:
        base_rows, gene_col, base_fields = read_input_table(Path(args.genes_file))
    genes: List[str] = []
    if args.genes:
        genes.extend(args.genes)
    if not genes and base_rows:
        genes = [r[gene_col] for r in base_rows]
    if args.genes and base_rows:
        # append ad-hoc genes as extra rows
        for g in args.genes:
            base_rows.append({gene_col: g})
    if not genes:
        raise SystemExit("No genes provided")

    cfg = load_config(Path(args.config))
    ann_rows = annotate(genes, cfg)
    ann_by_gene = {r["gene"]: r for r in ann_rows}
    # determine annotation fields to append
    ann_fields: List[str] = []
    if ann_rows:
        for k in ann_rows[0].keys():
            if k != "gene":
                ann_fields.append(k)

    merged: List[Dict[str, Any]] = []
    if not base_rows:
        base_rows = [{gene_col: g} for g in genes]
    for r in base_rows:
        gene_val = r.get(gene_col) or r.get("gene") or ""
        out_row = dict(r)
        ann = ann_by_gene.get(gene_val, {})
        for f in ann_fields:
            out_row[f] = ann.get(f, "")
        merged.append(out_row)

    # final header: original columns + new annotation columns (avoid duplicates)
    header = []
    for col in base_fields:
        if col not in header:
            header.append(col)
    for col in merged[0].keys():
        if col not in header:
            header.append(col)
    # rewrite rows with ordered keys
    ordered_rows = [{k: row.get(k, "") for k in header} for row in merged]
    with Path(args.out).open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        for row in ordered_rows:
            w.writerow(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
