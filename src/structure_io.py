"""Reading and converting protein structures (PDB / mmCIF)."""

from io import StringIO
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, Select


class _CAOnly(Select):
    def accept_atom(self, atom):
        return atom.get_id() == "CA"


def cif_to_pdb_string(cif_path, ca_only=True):
    """Convert an mmCIF file to a PDB-format string (CA atoms only by default)."""
    structure = MMCIFParser(QUIET=True).get_structure("x", str(cif_path))
    io = PDBIO()
    io.set_structure(structure)
    buf = StringIO()
    io.save(buf, select=_CAOnly() if ca_only else Select())
    return buf.getvalue()


def read_structure_as_pdb(path, ca_only=True):
    """Return PDB text for a .pdb or .cif file (CIF is converted on the fly)."""
    path = Path(path)
    if path.suffix.lower() in {".cif", ".mmcif"}:
        return cif_to_pdb_string(path, ca_only=ca_only)
    return path.read_text()
