import json
import re
import subprocess
import tempfile
from pathlib import Path

from pdb_validation import validate_pdb

TOOLS_DIR = Path(__file__).resolve().parent.parent / "tools"
TMSCORE = TOOLS_DIR / "TMscore"
OST = TOOLS_DIR / "ost-env" / "bin" / "ost"

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
    Computes lDDT of a predicted structure against the ground truth with OpenStructure.
    Only CA atoms are used (--bb-lddt).

    Args:
        pred_pdb: PDB-format text (CA-only) of the prediction.
        native_pdb: PDB-format text (CA-only) of the ground truth.

    Returns:
        lDDT in the range [0, 1]. Higher is better.
    """
    if not OST.exists():
        raise FileNotFoundError(f"{OST} not found, run scripts/install_ost.sh")
    _validate_pair(pred_pdb, native_pdb)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "pred.pdb").write_text(pred_pdb)
        (tmp / "native.pdb").write_text(native_pdb)
        result = subprocess.run(
            [
                str(OST),
                "compare-structures",
                "-m",
                str(tmp / "pred.pdb"),
                "-r",
                str(tmp / "native.pdb"),
                "-o",
                str(tmp / "out.json"),
                "--residue-number-alignment",  # residues are numbered 1..L in both files
                "--bb-lddt",  # lDDT on CA atoms only
                "--lddt-no-stereochecks",  # CA-only files have no side chains to check
            ],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        out_path = tmp / "out.json"
        if not out_path.exists():
            raise RuntimeError(
                f"OpenStructure produced no output (exit {result.returncode}):\n{result.stderr}"
            )
        data = json.loads(out_path.read_text())

    if data.get("status") != "SUCCESS":
        raise RuntimeError(f"OpenStructure failed:\n{data.get('traceback')}")
    return data["bb_lddt"]
