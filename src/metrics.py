"""Structural evaluation metrics (RMSD, TM-score, GDT_TS, lDDT) and CSV export."""

import json
import re
import subprocess
import tempfile
from pathlib import Path
import csv

from structure_io import read_structure_as_pdb
from pdb_validation import validate_pdb

TOOLS_DIR = Path(__file__).resolve().parent.parent / "tools"
TMSCORE = TOOLS_DIR / "TMscore"
OST_PYTHON = TOOLS_DIR / "ost-env" / "bin" / "python"
OST_WORKER = Path(__file__).resolve().parent / "ost_lddt_worker.py"

_TMSCORE_PATTERNS = {
    "rmsd": r"RMSD of\s+the common residues=\s*([\d.]+)",
    "tm_score": r"TM-score\s*=\s*([\d.]+)",
    "gdt_ts": r"GDT-TS-score=\s*([\d.]+)",
    "n_common": r"Number of residues in common=\s*(\d+)",
}


def _validate_pair(pred_pdb, native_pdb):
    """Validate both structures; return the number of residues in the native one."""
    n_native = len(validate_pdb(native_pdb))
    validate_pdb(pred_pdb, expected_length=n_native)
    return n_native


def tm_score(pred_pdb, native_pdb, timeout=60):
    """
    Compares a predicted structure with the ground-truth structure using TMscore.

    Args:
        pred_pdb: PDB-format text (CA-only) of the prediction.
        native_pdb: PDB-format text (CA-only) of the ground truth.

    Returns:
        dict with rmsd, tm_score (0-1), gdt_ts (0-1) and n_common (residues compared).
    """
    if not TMSCORE.exists():
        raise FileNotFoundError(f"{TMSCORE} not found, run scripts/install_tmscore.sh")
    n_native = _validate_pair(pred_pdb, native_pdb)

    # TMscore reads files only, so the text goes through a temporary directory
    with tempfile.TemporaryDirectory() as tmp:
        pred_path, native_path = Path(tmp) / "pred.pdb", Path(tmp) / "native.pdb"
        pred_path.write_text(pred_pdb)
        native_path.write_text(native_pdb)
        result = subprocess.run(
            [str(TMSCORE), str(pred_path), str(native_path)],
            capture_output=True,
            text=True,
            timeout=timeout,
        )

    metrics = {}
    for key, pattern in _TMSCORE_PATTERNS.items():
        m = re.search(pattern, result.stdout)
        if m is None:
            raise RuntimeError(
                f"'{key}' not found in TMscore output (exit {result.returncode}):\n{result.stdout}{result.stderr}"
            )
        metrics[key] = int(m.group(1)) if key == "n_common" else float(m.group(1))

    if metrics["n_common"] != n_native:
        raise RuntimeError(
            f"TMscore matched {metrics['n_common']} of {n_native} residues"
        )
    return metrics


def lddt_score(pred_pdb, native_pdb, timeout=300):
    """
    Computes CA-only lDDT of a predicted structure against the ground truth
    using OpenStructure's lDDTScorer (bb_only=True).
    """
    if not OST_PYTHON.exists():
        raise FileNotFoundError(f"{OST_PYTHON} not found, run scripts/install_ost.sh")
    _validate_pair(pred_pdb, native_pdb)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        pred_path = tmp / "pred.pdb"
        native_path = tmp / "native.pdb"
        pred_path.write_text(pred_pdb)
        native_path.write_text(native_pdb)

        result = subprocess.run(
            [str(OST_PYTHON), str(OST_WORKER), str(pred_path), str(native_path)],
            capture_output=True,
            text=True,
            timeout=timeout,
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"OpenStructure lDDTScorer failed (exit {result.returncode}):\n{result.stderr}"
        )

    lines = [line for line in result.stdout.splitlines() if line.strip()]
    if not lines:
        raise RuntimeError(f"OpenStructure produced no output:\n{result.stderr}")

    data = json.loads(lines[-1])
    raw_lddt = data.get("lddt")
    if isinstance(raw_lddt, list):
        raw_lddt = raw_lddt[0]
    if raw_lddt is None:
        raise RuntimeError("OpenStructure lDDTScorer returned None (too few CA contacts)")
    return float(raw_lddt)


def evaluate_pair(pred_pdb: str, native_pdb: str) -> dict[str, float]:
    """Compute RMSD, TM_score, GDT_TS, and lDDT for a single structure pair."""
    tm = tm_score(pred_pdb, native_pdb)
    lddt = lddt_score(pred_pdb, native_pdb)
    return {
        "RMSD": tm["rmsd"],
        "TM_score": tm["tm_score"],
        "GDT_TS": tm["gdt_ts"],
        "lDDT": lddt,
    }


def evaluate_predictions(
    pred_dir: str | Path,
    ref_dir: str | Path,
    output_csv: str | Path,
) -> list[dict[str, object]]:
    """Compare all predicted PDBs in `pred_dir` against `ref_dir` and write results to CSV."""
    pred_dir = Path(pred_dir)
    ref_dir = Path(ref_dir)
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    pred_files = sorted(pred_dir.glob("*.pdb"))
    if not pred_files:
        raise FileNotFoundError(f"No predicted .pdb files found in {pred_dir}")

    rows: list[dict[str, object]] = []
    for pred_path in pred_files:
        protein_id = pred_path.stem
        ref_path = ref_dir / f"{protein_id}.pdb"
        if not ref_path.exists():
            ref_path = ref_dir / f"{protein_id}.cif"
        if not ref_path.exists():
            raise FileNotFoundError(f"Missing reference for {protein_id} in {ref_dir}")

        try:
            scores = evaluate_pair(
                read_structure_as_pdb(pred_path),
                read_structure_as_pdb(ref_path),
            )
            rows.append({"proteinID": protein_id, **scores})
        except Exception as e:
            print(f"[FAILED] {protein_id}: {e}")

    fieldnames = ["proteinID", "RMSD", "TM_score", "GDT_TS", "lDDT"]
    with output_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved {len(rows)} evaluation rows to {output_csv}")
    return rows
