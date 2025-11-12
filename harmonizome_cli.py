#!/usr/bin/env python3
"""
Harmonizome CLI — a tiny, dependency‑light downloader for Harmonizome datasets.

Features
- Select datasets by short key (e.g. "gtextissue") or by Name=path mapping
- Read dataset list from a TSV/CSV/plain file
- Choose which downloadable files to fetch (or try "all")
- Optional on-the-fly decompression of *.gz to *.txt (streaming)
- Concurrent downloads with robust retries
- Skips files that already exist (unless --force)

Examples
---------
# 1) Download two datasets' matrices
harmonizome_cli.py download \
  -d gtextissue -d biogrid \
  -t gene_attribute_matrix.txt.gz -t gene_attribute_edges.txt.gz \
  -o ./harmonizome_data

# 2) Try all downloadables for GTEx and HPA, auto-decompress
harmonizome_cli.py download -d gtextissue -d hpatissuesmrna -t all --decompress

# 3) Provide datasets via a text file (one per line; supports "Name\tpath" or just "path")
harmonizome_cli.py download --datasets-file my_datasets.txt -t gene_attribute_matrix.txt.gz

# 4) List common dataset keys & downloadables the tool knows about
harmonizome_cli.py list

"""
from __future__ import annotations
import argparse
import concurrent.futures as cf
import csv
import gzip
import io
import os
import sys
import textwrap
from pathlib import Path
from typing import Iterable, List, Tuple

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE_URL = "https://maayanlab.cloud/static/hdfs/harmonizome/data"

# A compact starter map of common datasets (add more as needed)
DATASET_MAP = {
    # short_key: human_readable_name (path is the key itself)
    "gtextissue": "GTEx Tissue Gene Expression Profiles",
    "gtexsample": "GTEx Tissue Sample Gene Expression Profiles",
    "gtexeqtl": "GTEx eQTL",
    "biogrid": "BioGRID Protein-Protein Interactions",
    "encodetf": "ENCODE Transcription Factor Binding Site Profiles",
    "encodehm": "ENCODE Histone Modification Site Profiles",
    "reactome": "Reactome Pathways",
    "reactomeppi": "Reactome Biomolecular Interactions",
    "kegg": "KEGG Pathways",
    "keggppi": "KEGG Biomolecular Interactions",
    "hpatissuesmrna": "HPA Tissue Gene Expression Profiles",
    "hpasamples": "HPA Tissue Sample Gene Expression Profiles",
    "hprd": "HPRD Protein-Protein Interactions",
    "omim": "OMIM Gene-Disease Associations",
    "hpo": "HPO Gene-Disease Associations",
    "gdsc": "GDSC Cell Line Gene Expression Profiles",
    "tcga": "TCGA Signatures of Differentially Expressed Genes for Tumors",
    "jasparpwm": "JASPAR Predicted Transcription Factor Targets",
    "targetscan": "TargetScan Predicted Conserved microRNA Targets",
}

# Common downloadables. You can pass -t all to attempt them all.
ALL_DOWNLOADABLES = [
    "gene_attribute_matrix.txt.gz",
    "gene_attribute_edges.txt.gz",
    "gene_set_library_crisp.txt.gz",
    "gene_set_library_up_crisp.txt.gz",
    "gene_set_library_dn_crisp.txt.gz",
    "attribute_set_library_crisp.txt.gz",
    "attribute_set_library_up_crisp.txt.gz",
    "attribute_set_library_dn_crisp.txt.gz",
    "gene_similarity_matrix_cosine.txt.gz",
    "attribute_similarity_matrix_cosine.txt.gz",
    "gene_list_terms.txt.gz",
    "attribute_list_entries.txt.gz",
    "processing_script.m",
]

# ----------------------- HTTP session with retries ----------------------- #

def make_session(total_retries: int = 5, backoff: float = 0.5) -> requests.Session:
    s = requests.Session()
    retry = Retry(
        total=total_retries,
        read=total_retries,
        connect=total_retries,
        respect_retry_after_header=True,
        backoff_factor=backoff,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=16, pool_maxsize=16)
    s.mount("http://", adapter)
    s.mount("https://", adapter)
    return s

# ----------------------------- Utilities -------------------------------- #

def parse_dataset_item(s: str) -> Tuple[str, str]:
    """Accept forms like
    - "gtextissue" -> ("GTEx Tissue Gene Expression Profiles", "gtextissue") if in map
    - "Name=path" -> ("Name", "path")
    - "path" -> ("path", "path")
    - "Name\tpath" -> ("Name", "path")
    """
    if "\t" in s:
        name, path = s.split("\t", 1)
        return name.strip(), path.strip()
    if "=" in s:
        name, path = s.split("=", 1)
        return name.strip(), path.strip()
    # single token: try known map; else use token as both name and path
    token = s.strip()
    if token in DATASET_MAP:
        return DATASET_MAP[token], token
    return token, token


def read_datasets_file(fp: Path) -> List[Tuple[str, str]]:
    items: List[Tuple[str, str]] = []
    with fp.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            items.append(parse_dataset_item(line))
    return items


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


# ----------------------------- Downloaders ------------------------------ #

def build_url(dataset_path: str, downloadable: str) -> str:
    return f"{BASE_URL}/{dataset_path}/{downloadable}"


def download_one(
    session: requests.Session,
    dataset_name: str,
    dataset_path: str,
    downloadable: str,
    outdir: Path,
    decompress: bool,
    force: bool,
    timeout: int,
) -> Tuple[str, str, str, bool, str]:
    """Return (dataset_name, dataset_path, downloadable, success, message)."""
    dataset_dir = outdir / dataset_name
    ensure_dir(dataset_dir)

    target = dataset_dir / downloadable
    if decompress and downloadable.endswith(".gz"):
        target = dataset_dir / downloadable[:-3]  # strip .gz -> .txt

    if target.exists() and not force:
        return (dataset_name, dataset_path, downloadable, True, "exists, skipped")

    url = build_url(dataset_path, downloadable)
    try:
        with session.get(url, stream=True, timeout=timeout) as resp:
            if resp.status_code != 200:
                return (dataset_name, dataset_path, downloadable, False, f"HTTP {resp.status_code}")

            if decompress and downloadable.endswith(".gz"):
                # stream-decompress gzip to text
                with resp.raw as raw:
                    gzf = gzip.GzipFile(fileobj=raw, mode="rb")
                    # text wrapper with utf-8 decode
                    with io.TextIOWrapper(gzf, encoding="utf-8", newline="") as reader, target.open(
                        "w", encoding="utf-8", newline=""
                    ) as out:
                        for chunk in iter(lambda: reader.read(1024 * 64), ""):
                            out.write(chunk)
            else:
                # write binary as-is
                with target.open("wb") as out:
                    for chunk in resp.iter_content(chunk_size=1024 * 64):
                        if chunk:
                            out.write(chunk)

        return (dataset_name, dataset_path, downloadable, True, "ok")
    except Exception as e:
        return (dataset_name, dataset_path, downloadable, False, str(e))


# ------------------------------- CLI ----------------------------------- #

def cmd_list(args: argparse.Namespace) -> int:
    print("Known dataset keys → names (you can still pass Name=path or custom path):\n")
    for k, v in sorted(DATASET_MAP.items()):
        print(f"  {k:<18} {v}")
    print("\nCommon downloadables (use -t all to attempt all of these):\n")
    for x in ALL_DOWNLOADABLES:
        print(f"  {x}")
    return 0


def flatten_downloadables(tokens: List[str]) -> List[str]:
    result: List[str] = []
    for t in tokens:
        if t.lower() == "all":
            result.extend(ALL_DOWNLOADABLES)
        else:
            result.append(t)
    # de-duplicate while preserving order
    seen = set()
    uniq: List[str] = []
    for x in result:
        if x not in seen:
            seen.add(x)
            uniq.append(x)
    return uniq


def cmd_download(args: argparse.Namespace) -> int:
    # Collect datasets
    datasets: List[Tuple[str, str]] = []
    for item in args.dataset or []:
        datasets.append(parse_dataset_item(item))
    if args.datasets_file:
        datasets.extend(read_datasets_file(Path(args.datasets_file)))
    if not datasets:
        print("[ERROR] No datasets provided. Use -d/--dataset or --datasets-file.", file=sys.stderr)
        return 2

    downloadables = flatten_downloadables(args.downloadable or [])
    if not downloadables:
        print("[ERROR] No downloadables provided. Use -t/--downloadable or -t all.", file=sys.stderr)
        return 2

    outdir = Path(args.outdir or ".").resolve()
    ensure_dir(outdir)

    # Prepare work list
    jobs = []
    for name, path in datasets:
        for dl in downloadables:
            jobs.append((name, path, dl))

    print(f"→ Datasets: {len(datasets)}; files to attempt: {len(jobs)}; workers: {args.workers}")

    session = make_session(total_retries=args.retries, backoff=args.backoff)

    results = []
    with cf.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [
            ex.submit(
                download_one,
                session,
                name,
                path,
                dl,
                outdir,
                args.decompress,
                args.force,
                args.timeout,
            )
            for (name, path, dl) in jobs
        ]
        for fut in cf.as_completed(futs):
            results.append(fut.result())

    # Summarize
    ok = [r for r in results if r[3]]
    fail = [r for r in results if not r[3]]

    for (name, path, dl, success, msg) in results:
        status = "OK" if success else "FAIL"
        print(f"[{status}] {name} ({path}) :: {dl} -> {msg}")

    print(f"\nDone. Success: {len(ok)}  Fail: {len(fail)}  (outdir={outdir})")
    return 0 if not fail else 1


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="harmonizome_cli",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__ or "Harmonizome CLI"),
    )

    sub = p.add_subparsers(dest="cmd", required=True)

    sp_list = sub.add_parser("list", help="Show known dataset keys and common downloadables")
    sp_list.set_defaults(func=cmd_list)

    sp_dl = sub.add_parser("download", help="Download selected datasets/files")
    sp_dl.add_argument(
        "-d",
        "--dataset",
        action="append",
        help=(
            "Dataset selector; supports forms: 'gtextissue' (known key), 'Name=path', 'path', or 'Name\\tpath'. "
            "Use multiple -d to add more."
        ),
    )
    sp_dl.add_argument(
        "--datasets-file",
        help=(
            "Text/CSV/TSV file with one dataset per line. Each line can be 'path', 'Name=path', or 'Name\tpath'."
        ),
    )
    sp_dl.add_argument(
        "-t",
        "--downloadable",
        action="append",
        help=(
            "Which files to fetch, e.g. 'gene_attribute_matrix.txt.gz'. Use -t all to attempt common files. "
            "Repeat -t for multiple types."
        ),
    )
    sp_dl.add_argument("-o", "--outdir", default=".", help="Output directory (default: current directory)")
    sp_dl.add_argument("--decompress", action="store_true", help="Stream-decompress *.gz to *.txt while downloading")
    sp_dl.add_argument("--force", action="store_true", help="Re-download even if target file already exists")
    sp_dl.add_argument("-w", "--workers", type=int, default=4, help="Parallel download workers (default: 4)")
    sp_dl.add_argument("--retries", type=int, default=5, help="HTTP retry attempts (default: 5)")
    sp_dl.add_argument("--backoff", type=float, default=0.5, help="Retry backoff factor (default: 0.5)")
    sp_dl.add_argument("--timeout", type=int, default=60, help="HTTP timeout seconds per request (default: 60)")
    sp_dl.set_defaults(func=cmd_download)

    return p


def main(argv: List[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
