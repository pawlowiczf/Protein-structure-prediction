"""Validation of CA-only PDB text, independent of how it is produced or used."""

from io import StringIO

import numpy as np
from Bio.PDB import PDBParser


class PDBValidationError(ValueError):
    """Raised when a PDB string is not a usable CA-only structure."""


def validate_pdb(pdb_text, expected_length=None):
    """Check that `pdb_text` is a well-formed CA-only structure.

    Checks: parses without errors, one model, one chain, only CA atoms,
    no duplicate residues, all coordinates finite, and optionally the
    number of residues. Returns the CA coordinates as an (L, 3) array.
    Raises PDBValidationError otherwise.
    """
    if not pdb_text or not pdb_text.strip():
        raise PDBValidationError("empty PDB string")
    if not any(line.startswith(("ATOM", "HETATM")) for line in pdb_text.splitlines()):
        raise PDBValidationError("no ATOM records")

    try:
        structure = PDBParser(PERMISSIVE=False, QUIET=True).get_structure(
            "x", StringIO(pdb_text)
        )
    except Exception as e:  # Biopython raises several exception types
        raise PDBValidationError(f"cannot parse PDB: {e}") from e

    models = list(structure)
    if len(models) != 1:
        raise PDBValidationError(f"expected 1 model, found {len(models)}")
    chains = list(models[0])
    if len(chains) != 1:
        raise PDBValidationError(f"expected 1 chain, found {len(chains)}")

    coords = []
    for residue in chains[0]:
        atoms = [a.get_id() for a in residue]
        if atoms != ["CA"]:
            raise PDBValidationError(
                f"residue {residue.get_id()[1]} {residue.get_resname()} has atoms {atoms}, expected only CA"
            )
        coords.append(residue["CA"].get_coord())

    coords = np.asarray(coords, dtype=float)
    if len(coords) == 0:
        raise PDBValidationError("no residues")
    if not np.isfinite(coords).all():
        raise PDBValidationError("non-finite coordinates (nan/inf)")
    if expected_length is not None and len(coords) != expected_length:
        raise PDBValidationError(
            f"expected {expected_length} residues, found {len(coords)}"
        )
    return coords
