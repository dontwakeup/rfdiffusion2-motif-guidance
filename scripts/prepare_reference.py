#!/usr/bin/env python3
"""
Deterministic extraction of motif + boundary reference atoms

Inputs:
  - inputs/motif.json
  - inputs/reference.pdb

Output:
  - inputs/reference.npz
  - inputs/reference_summary.json

The .npz contains:
  - motif_resnums: (L,) int
  - boundary_resnums: (B,) int
  - motif_bb: (L, 3, 3) float32  (atoms ordered N,CA,C)
  - boundary_bb: (B, 3, 3) float32
  - meta_json: (1,) str  (JSON string)
"""

import argparse
import json
import os
import platform
import sys
from dataclasses import dataclass
from datetime import datetime, timezone

import numpy as np
from Bio import __version__ as biopython_version
from Bio.PDB import PDBParser


BB_ATOMS = ("N", "CA", "C")


@dataclass(frozen=True)
class MotifSpec:
    pdb_id: str
    source: str
    chain: str
    motif_residues: str  # "45-62"
    boundary_flank: int
    notes: str = ""


def _parse_range(s: str) -> tuple[int, int]:
    """
    Parse a contiguous range like "45-62" -> (45, 62)
    """
    if not s or "-" not in s:
        raise ValueError(f"motif_residues must be like '45-62', got: {s!r}")
    a, b = s.split("-", 1)
    start = int(a.strip())
    end = int(b.strip())
    if end < start:
        raise ValueError(f"Invalid range: {s!r} (end < start)")
    return start, end


def _load_motif_json(path: str) -> MotifSpec:
    with open(path, "r", encoding="utf-8") as f:
        d = json.load(f)
    required = ["pdb_id", "source", "chain", "motif_residues", "boundary_flank"]
    missing = [k for k in required if k not in d]
    if missing:
        raise ValueError(f"Missing keys in motif.json: {missing}")
    return MotifSpec(
        pdb_id=str(d["pdb_id"]),
        source=str(d["source"]),
        chain=str(d["chain"]),
        motif_residues=str(d["motif_residues"]),
        boundary_flank=int(d["boundary_flank"]),
        notes=str(d.get("notes", "")),
    )


def _resnum_from_residue(residue) -> int | None:
    """
    Biopython residue id: (hetero_flag, resseq, icode)
    We ignore hetero residues and require standard resseq int
    """
    hetflag, resseq, icode = residue.id
    if hetflag.strip():  # HETATM or water
        return None
    return int(resseq)


def _collect_chain_residues(chain) -> dict[int, object]:
    """
    Return mapping resseq -> residue (standard residues only)
    If duplicate resseq occurs (rare), last one wins
    """
    m: dict[int, object] = {}
    for res in chain.get_residues():
        resnum = _resnum_from_residue(res)
        if resnum is None:
            continue
        m[resnum] = res
    return m


def _extract_backbone_xyz(residue) -> np.ndarray | None:
    """
    Return (3,3) array for N,CA,C or None if any atom missing
    """
    coords = []
    for name in BB_ATOMS:
        if name not in residue:
            return None
        atom = residue[name]
        coords.append(atom.coord.astype(np.float32))
    return np.stack(coords, axis=0)  # (3,3)


def _extract_block(
    res_map: dict[int, object], resnums: list[int]
) -> tuple[np.ndarray, list[int], list[int]]:
    """
    Extract backbond coords for a list of residue numbers

    Returns:
        bb: (K,3,3) float32 for residues that were successfully extracted
        ok_resnums: residue numbers included
        missing_resnums: residue numbers missing or incomplete
    """
    bb_list: list[np.ndarray] = []
    ok: list[int] = []
    missing: list[int] = []

    for r in resnums:
        res = res_map.get(r)
        if res is None:
            missing.append(r)
            continue
        xyz = _extract_backbone_xyz(res)
        if xyz is None:
            missing.append(r)
            continue
        bb_list.append(xyz)
        ok.append(r)

    if bb_list:
        bb = np.stack(bb_list, axis=0).astype(np.float32)
    else:
        bb = np.zeros((0, 3, 3), dtype=np.float32)
    return bb, ok, missing


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--motif-json", required=True)
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    args = ap.parse_args()

    spec = _load_motif_json(args.motif_json)
    start, end = _parse_range(spec.motif_residues)
    motif_resnums = list(range(start, end + 1))
    boundary_resnums = list(
        range(start - spec.boundary_flank, end + spec.boundary_flank)
    )

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(spec.pdb_id, args.pdb)

    # Choose first model deterministically
    model = next(structure.get_models())

    if spec.chain not in model:
        raise KeyError(f"Chain {spec.chain!r} not found in PDB: {args.pdb}")

    chain = model[spec.chain]
    res_map = _collect_chain_residues(chain)

    motif_bb, motif_ok, motif_missing = _extract_block(res_map, motif_resnums)
    boundary_bb, boundary_ok, boundary_missing = _extract_block(
        res_map, boundary_resnums
    )

    # boundary should include motif; if boundry trimmed by missing residues, that's fine, but record it
    meta = {
        "pdb_id": spec.pdb_id,
        "pdb_source_path": os.path.normpath(args.pdb),
        "chain": spec.chain,
        "motif_residues_input": spec.motif_residues,
        "boundary_flank": spec.boundary_flank,
        "notes": spec.notes,
        "motif_resnums_requested": motif_resnums,
        "motif_resnums_extracted": motif_ok,
        "motif_missing": motif_missing,
        "boundary_resnums_requested": boundary_resnums,
        "boundary_resnums_extracted": boundary_ok,
        "boundary_missing": boundary_missing,
        "bb_atoms_order": list(BB_ATOMS),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "platform": {
            "python": sys.version.split()[0],
            "biopython": biopython_version,
            "numpy": np.__version__,
            "os": platform.platform(),
        },
    }

    # Save .npz deterministically
    np.savez_compressed(
        args.out,
        motif_resnums=np.asarray(motif_ok, dtype=np.int32),
        boundary_resnums=np.asarray(boundary_ok, dtype=np.int32),
        motif_bb=motif_bb,
        boundary_bb=boundary_bb,
        meta_json=np.asarray([json.dumps(meta, sort_keys=True)], dtype=object),
    )

    with open(args.summary, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, sort_keys=True)

    print("=== Reference extraction summary ===")
    print(f"PDB: {spec.pdb_id}  chain: {spec.chain}")
    print(f"Motif requested: {spec.motif_residues}  flank: ±{spec.boundary_flank}")
    print(f"Motif extracted residues: {len(motif_ok)} / {len(motif_resnums)}")
    if motif_missing:
        print(f"Motif missing/incomplete: {motif_missing}")
    print(f"Boundary extracted residues: {len(boundary_ok)} / {len(boundary_resnums)}")
    if boundary_missing:
        print(f"Boundary missing/incomplete: {boundary_missing}")
    print(f"Wrote: {args.out}")
    print(f"Wrote: {args.summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
