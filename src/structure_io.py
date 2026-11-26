"""Reading, writing, and converting protein structures and sequences (PDB / mmCIF / FASTA)."""

from io import StringIO
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.SeqUtils import seq1
from skbio.io import read

from pdb_validation import PDBValidationError, validate_pdb

THREE_LETTER = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
    "X": "UNK",
}


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


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Read sequences from a FASTA file as a list of (seq_id, sequence) tuples."""
    records = []
    for seq in read(str(path), format="fasta"):
        seq_id = seq.metadata.get("id") or f"sequence_{len(records) + 1}"
        records.append((str(seq_id), str(seq)))
    return records


def write_pdb(sequence: str, coords: np.ndarray, path: str | Path) -> Path:
    """Write one CA atom per residue in wwPDB ATOM format."""
    coords = np.asarray(coords)
    lines = []
    for i, (aa, xyz) in enumerate(zip(sequence, coords), start=1):
        res_name = THREE_LETTER.get(aa, "UNK")
        line = (
            f"ATOM  {i:5d}  CA  {res_name} A{i:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00           C  "
        )
        lines.append(line)
    lines.append("END")

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")
    return out_path


def prepare_reference_from_cif(
    cif_path: str | Path,
    pdb_dir: str | Path,
    fasta_dir: str | Path | None = None,
) -> str:
    """Convert a single .cif file into validated CA-only .pdb and optional .fasta files."""
    cif_path = Path(cif_path)
    pdb_text = cif_to_pdb_string(cif_path, ca_only=True)
    validate_pdb(pdb_text)

    resnames = [
        line[17:20].strip() for line in pdb_text.splitlines() if line.startswith("ATOM")
    ]
    sequence = seq1("".join(resnames))

    pdb_path = Path(pdb_dir) / f"{cif_path.stem}.pdb"
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(pdb_text)

    if fasta_dir is not None:
        fasta_path = Path(fasta_dir) / f"{cif_path.stem}.fasta"
        fasta_path.parent.mkdir(parents=True, exist_ok=True)
        fasta_path.write_text(f">{cif_path.stem}\n{sequence}\n")

    return sequence


def prepare_split_subset(
    cif_dir: str | Path,
    pdb_dir: str | Path,
    fasta_dir: str | Path,
    limit: int | None = None,
) -> None:
    """Batch-convert .cif files from a split directory, skipping already generated pairs."""
    cif_dir = Path(cif_dir)
    pdb_dir = Path(pdb_dir)
    fasta_dir = Path(fasta_dir)
    pdb_dir.mkdir(parents=True, exist_ok=True)
    fasta_dir.mkdir(parents=True, exist_ok=True)

    cif_files = sorted(cif_dir.rglob("*.cif"))
    if not cif_files:
        raise FileNotFoundError(f"No .cif files found in {cif_dir}")

    ready = converted = skipped = failed = 0
    for cif_path in cif_files:
        if limit is not None and ready >= limit:
            break

        pdb_out = pdb_dir / f"{cif_path.stem}.pdb"
        fasta_out = fasta_dir / f"{cif_path.stem}.fasta"

        if (
            pdb_out.exists()
            and pdb_out.stat().st_size > 0
            and fasta_out.exists()
            and fasta_out.stat().st_size > 0
        ):
            skipped += 1
            ready += 1
            continue

        try:
            prepare_reference_from_cif(cif_path, pdb_dir, fasta_dir)
            converted += 1
            ready += 1
        except PDBValidationError as e:
            failed += 1
            print(f"[REJECTED] {cif_path.name}: {e}")

    print(
        f"Summary ({cif_dir.name}): "
        f"converted={converted}, skipped={skipped}, failed={failed}"
    )
