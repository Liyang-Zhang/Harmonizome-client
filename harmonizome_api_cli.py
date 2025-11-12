#!/usr/bin/env python3
"""
Harmonizome API CLI — Query Harmonizome REST & Knowledge Graph APIs from the terminal.

Features
- Subcommands for REST entities: gene, dataset, attribute, protein, resource, gene-set, hgnc-family, naming-authority
- Get single entity or list (with cursor pagination & --all to exhaust)
- Optional --show-associations for gene / gene-set
- Knowledge Graph (KG) queries: neighbors, path, constrained relations, remove nodes
- Robust HTTP: retries, timeouts; pretty or JSON output
- Zero non‑stdlib deps except 'requests'

Examples
--------
# 1) REST: list first 100 genes (JSON by default)
harmonizome_api_cli.py gene list

# 2) REST: list genes with cursor and auto-pagination to fetch all
harmonizome_api_cli.py gene list --all

# 3) REST: get gene NANOG base info
harmonizome_api_cli.py gene get NANOG

# 4) REST: get gene NANOG with associations
harmonizome_api_cli.py gene get NANOG --show-associations

# 5) REST: list datasets and render a TSV summary to stdout
harmonizome_api_cli.py dataset list --format tsv

# 6) KG: immediate neighbors of STAT3 (limit 10)
harmonizome_api_cli.py kg neighbors --start Gene --term STAT3 --limit 10

# 7) KG: STAT3 participates_in_(GO Bio Process 2023) only
harmonizome_api_cli.py kg neighbors --start Gene --term STAT3 \
  --relation "participates_in_(GO Bio Process 2023)" --limit 5

# 8) KG: shortest path STAT3 -> MAPK1
harmonizome_api_cli.py kg path --start Gene --term STAT3 --end Gene --end-term MAPK1

"""
from __future__ import annotations
import argparse
import json
import os
import sys
import time
from typing import Any, Dict, Iterable, List, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BASE = "https://maayanlab.cloud/Harmonizome/api/1.0"
KG_BASE = "https://harmonizome-kg.maayanlab.cloud/api/knowledge_graph"

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
    adapter = HTTPAdapter(max_retries=retry, pool_connections=8, pool_maxsize=8)
    s.mount("http://", adapter)
    s.mount("https://", adapter)
    s.headers.update({"User-Agent": "harmonizome-api-cli/1.0"})
    return s

# ----------------------------- Utilities -------------------------------- #

def eprint(*a: Any, **k: Any) -> None:
    print(*a, file=sys.stderr, **k)


def dumps(obj: Any, pretty: bool) -> str:
    return json.dumps(obj, ensure_ascii=False, indent=2 if pretty else None)


def print_tsv(rows: List[Dict[str, Any]], fields: Optional[List[str]] = None) -> None:
    if not rows:
        return
    if fields is None:
        # union of keys in order of first row, then append unseen keys
        fields = list(rows[0].keys())
        seen = set(fields)
        for r in rows[1:]:
            for k in r.keys():
                if k not in seen:
                    seen.add(k)
                    fields.append(k)
    print("\t".join(fields))
    for r in rows:
        line = []
        for f in fields:
            v = r.get(f, "")
            if isinstance(v, (dict, list)):
                v = json.dumps(v, ensure_ascii=False)
            line.append(str(v))
        print("\t".join(line))


def safe_get(session: requests.Session, url: str, params: Optional[Dict[str, Any]] = None, timeout: int = 30) -> Dict[str, Any]:
    r = session.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.json()


# ------------------------------ REST layer ------------------------------- #
ENTITY_PATHS = {
    "attribute": "attribute",
    "dataset": "dataset",
    "gene": "gene",
    "gene-set": "gene_set",
    "hgnc-family": "hgnc_family",
    "naming-authority": "naming_authority",
    "protein": "protein",
    "resource": "resource",
}


def rest_list(session: requests.Session, entity: str, cursor: Optional[int], all_pages: bool, delay: float, timeout: int) -> Dict[str, Any]:
    path = ENTITY_PATHS[entity]
    url = f"{BASE}/{path}"
    params: Dict[str, Any] = {}
    if cursor is not None:
        params["cursor"] = cursor
    page = safe_get(session, url, params=params, timeout=timeout)
    if not all_pages:
        return page

    # Accumulate pages
    entities: List[Any] = []
    total = page.get("count")
    next_url = page.get("next")
    entities.extend(page.get("entities", []))
    while next_url:
        full = f"https://maayanlab.cloud{next_url}" if next_url.startswith("/") else next_url
        time.sleep(delay)
        p = safe_get(session, full, timeout=timeout)
        entities.extend(p.get("entities", []))
        next_url = p.get("next")
    return {"count": total if total is not None else len(entities), "entities": entities}


def rest_get(session: requests.Session, entity: str, name: str, show_associations: bool, timeout: int) -> Dict[str, Any]:
    path = ENTITY_PATHS[entity]
    url = f"{BASE}/{path}/{name}"
    params = {"showAssociations": str(show_associations).lower()} if show_associations else None
    return safe_get(session, url, params=params, timeout=timeout)


# --------------------------- Knowledge Graph ---------------------------- #

def kg_query(session: requests.Session, payload: Dict[str, Any], timeout: int) -> Dict[str, Any]:
    params = {"filter": json.dumps(payload, ensure_ascii=False)}
    return safe_get(session, KG_BASE, params=params, timeout=timeout)


# ------------------------------- CLI ----------------------------------- #

def add_common_rest_args(sp: argparse.ArgumentParser) -> None:
    sp.add_argument("--cursor", type=int, help="Start index for pagination (default API page=100)")
    sp.add_argument("--all", action="store_true", help="Fetch all pages by following 'next'")
    sp.add_argument("--delay", type=float, default=0.1, help="Delay between pages when --all (seconds)")
    sp.add_argument("--timeout", type=int, default=30, help="HTTP timeout seconds (default: 30)")
    sp.add_argument("--format", choices=["json", "tsv"], default="json", help="Output format for list")
    sp.add_argument("--pretty", action="store_true", help="Pretty JSON output")


def cmd_rest_list(args: argparse.Namespace) -> int:
    session = make_session(args.retries, args.backoff)
    data = rest_list(session, args.entity, args.cursor, args.all, args.delay, args.timeout)
    if args.format == "json":
        print(dumps(data, args.pretty))
    else:
        rows = data.get("entities", [])
        # best-effort flatten: keep top-level scalar fields
        flat_rows: List[Dict[str, Any]] = []
        for r in rows:
            flat = {k: v for k, v in r.items() if not isinstance(v, (dict, list))}
            # keep href/symbol/name if nested
            for k in ("symbol", "name", "href"):
                if k not in flat and k in r:
                    flat[k] = r[k]
            flat_rows.append(flat)
        print_tsv(flat_rows)
    return 0


def cmd_rest_get(args: argparse.Namespace) -> int:
    session = make_session(args.retries, args.backoff)
    data = rest_get(session, args.entity, args.name, args.show_associations, args.timeout)
    print(dumps(data, args.pretty))
    return 0


def build_rest_subparser(sub: argparse._SubParsersAction, entity_key: str, title: str) -> None:
    sp = sub.add_parser(entity_key, help=f"REST entity: {title}")
    sp.add_argument("--retries", type=int, default=5)
    sp.add_argument("--backoff", type=float, default=0.5)
    sp_sub = sp.add_subparsers(dest="action", required=True)

    sp_list = sp_sub.add_parser("list", help=f"List {title}s (paginated)")
    add_common_rest_args(sp_list)
    sp_list.set_defaults(func=cmd_rest_list, entity=entity_key)

    sp_get = sp_sub.add_parser("get", help=f"Get one {title}")
    sp_get.add_argument("name", help="Entity name/symbol (e.g., NANOG or dataset name)")
    sp_get.add_argument("--show-associations", action="store_true", help="Include associations when supported")
    sp_get.add_argument("--timeout", type=int, default=30)
    sp_get.add_argument("--retries", type=int, default=5)
    sp_get.add_argument("--backoff", type=float, default=0.5)
    sp_get.add_argument("--pretty", action="store_true")
    sp_get.set_defaults(func=cmd_rest_get, entity=entity_key)


# ------------------------------- KG CLI -------------------------------- #

def cmd_kg_neighbors(args: argparse.Namespace) -> int:
    session = make_session(args.retries, args.backoff)
    payload: Dict[str, Any] = {
        "start": args.start,
        "start_field": args.start_field,
        "start_term": args.term,
    }
    if args.limit is not None:
        payload["limit"] = args.limit
    if args.relation:
        payload["relation"] = [{"name": r, "limit": args.relation_limit} for r in args.relation]
    if args.remove:
        payload["remove"] = args.remove
    data = kg_query(session, payload, args.timeout)
    print(dumps(data, args.pretty))
    return 0


def cmd_kg_path(args: argparse.Namespace) -> int:
    session = make_session(args.retries, args.backoff)
    payload: Dict[str, Any] = {
        "start": args.start,
        "start_field": args.start_field,
        "start_term": args.term,
        "end": args.end,
        "end_field": args.end_field,
        "end_term": args.end_term,
    }
    if args.path_length is not None:
        payload["path_length"] = args.path_length
    if args.limit is not None:
        payload["limit"] = args.limit
    if args.relation:
        payload["relation"] = [{"name": r, "limit": args.relation_limit} for r in args.relation]
    if args.remove:
        payload["remove"] = args.remove
    data = kg_query(session, payload, args.timeout)
    print(dumps(data, args.pretty))
    return 0


def build_kg_subparser(sub: argparse._SubParsersAction) -> None:
    sp = sub.add_parser("kg", help="Knowledge Graph queries")
    sp.add_argument("--retries", type=int, default=5)
    sp.add_argument("--backoff", type=float, default=0.5)
    sp_sub = sp.add_subparsers(dest="action", required=True)

    # neighbors
    n = sp_sub.add_parser("neighbors", help="Immediate neighbors of a start node (optionally restricted by relation)")
    n.add_argument("--start", required=True, help="Start node type (e.g., Gene, HPO, Disease, Pathway)")
    n.add_argument("--start-field", default="label", help="Field to query in metadata (default: label)")
    n.add_argument("--term", required=True, help="Term value (e.g., STAT3)")
    n.add_argument("--relation", action="append", help="Restrict to relation name(s); repeat to add multiple")
    n.add_argument("--relation-limit", type=int, default=5, help="Per-relation limit (default: 5)")
    n.add_argument("--remove", action="append", help="IDs to remove from the subnetwork (repeatable)")
    n.add_argument("--limit", type=int, default=10, help="Neighbor limit (default: 10)")
    n.add_argument("--timeout", type=int, default=60)
    n.add_argument("--pretty", action="store_true")
    n.set_defaults(func=cmd_kg_neighbors)

    # path
    p = sp_sub.add_parser("path", help="Find shortest path between two nodes")
    p.add_argument("--start", required=True, help="Start node type (e.g., Gene)")
    p.add_argument("--start-field", default="label")
    p.add_argument("--term", required=True, help="Start term (e.g., STAT3)")
    p.add_argument("--end", required=True, help="End node type (e.g., Gene, HPO)")
    p.add_argument("--end-field", default="label")
    p.add_argument("--end-term", required=True, help="End term (e.g., MAPK1)")
    p.add_argument("--path-length", type=int, help="Max path length (default: server side)")
    p.add_argument("--relation", action="append", help="Restrict to relation name(s); repeat to add multiple")
    p.add_argument("--relation-limit", type=int, default=5)
    p.add_argument("--remove", action="append")
    p.add_argument("--limit", type=int, help="Overall limit for neighbors considered")
    p.add_argument("--timeout", type=int, default=60)
    p.add_argument("--pretty", action="store_true")
    p.set_defaults(func=cmd_kg_path)


# ------------------------------- Parser -------------------------------- #

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="harmonizome_api_cli", description="Harmonizome REST & KG API CLI")
    sub = p.add_subparsers(dest="cmd", required=True)

    # REST entity subcommands
    for key, title in [
        ("gene", "Gene"),
        ("dataset", "Dataset"),
        ("attribute", "Attribute"),
        ("protein", "Protein"),
        ("resource", "Resource"),
        ("gene-set", "Gene Set"),
        ("hgnc-family", "HGNC Family"),
        ("naming-authority", "Naming Authority"),
    ]:
        build_rest_subparser(sub, key, title)

    # KG
    build_kg_subparser(sub)

    return p


# -------------------------------- Main --------------------------------- #

def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except requests.HTTPError as e:
        eprint(f"HTTP error: {e}")
        if e.response is not None:
            try:
                eprint(e.response.text)
            except Exception:
                pass
        return 1
    except Exception as e:
        eprint(f"Error: {e}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
