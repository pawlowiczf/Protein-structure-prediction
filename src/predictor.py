"""Sequence encoding and structure prediction pipeline."""

from collections.abc import Callable
from pathlib import Path

import torch

from structure_io import read_fasta, write_pdb

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY") + ["X", "PAD"]
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
MAX_LENGTH = 1024
BATCH_SIZE = 16


def encode_batch(
    records: list[tuple[str, str]],
    max_length: int | None = MAX_LENGTH,
) -> tuple[torch.Tensor, torch.Tensor, list[str]]:
    """One-hot encode a batch of (seq_id, sequence) pairs."""
    standard = set("ACDEFGHIKLMNPQRSTVWY")
    n_aa = len(AMINO_ACIDS)

    if max_length is None:
        target_len = max(len(seq) for _, seq in records)
    else:
        target_len = max_length

    tokens = torch.zeros(len(records), target_len, n_aa)
    mask = torch.zeros(len(records), target_len, dtype=torch.bool)
    sequences = []

    for i, (_, sequence) in enumerate(records):
        letters = [aa if aa in standard else "X" for aa in sequence.upper()]
        letters = letters[:target_len]
        sequences.append("".join(letters))

        for j, aa in enumerate(letters):
            tokens[i, j, AA_TO_IDX[aa]] = 1
            mask[i, j] = True
        for j in range(len(letters), target_len):
            tokens[i, j, AA_TO_IDX["PAD"]] = 1

    return tokens, mask, sequences


def dummy_model(tokens: torch.Tensor, mask: torch.Tensor) -> list[torch.Tensor]:
    """Baseline random predictor returning (L, 3) CA coordinates per protein."""
    batch_size, length, _ = tokens.shape
    coords = torch.randn(batch_size, length, 3) * 10

    return [coords[i][mask[i]] for i in range(batch_size)]


def predict_directory(
    fasta_dir: str | Path,
    out_dir: str | Path,
    model_fn: Callable[[torch.Tensor, torch.Tensor], list[torch.Tensor]] = dummy_model,
    batch_size: int = BATCH_SIZE,
    max_length: int | None = None,
    limit: int | None = None,
    seed: int | None = 0,
) -> list[Path]:
    """Run batch structure prediction for all FASTA files in `fasta_dir`."""
    fasta_dir = Path(fasta_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(fasta_dir.glob("*.fasta"))[:limit]
    if not fasta_files:
        raise FileNotFoundError(f"No .fasta files in {fasta_dir}")

    records = []
    for fp in fasta_files:
        records.extend(read_fasta(fp))

    if seed is not None:
        torch.manual_seed(seed)

    written = []
    for start in range(0, len(records), batch_size):
        batch = records[start : start + batch_size]
        tokens, mask, sequences = encode_batch(batch, max_length=max_length)

        with torch.no_grad():
            predictions = model_fn(tokens, mask)

        for (seq_id, _), sequence, coords in zip(batch, sequences, predictions):
            pred_path = out_dir / f"{seq_id}.pdb"
            ca_coords = coords[: len(sequence)].detach().cpu().numpy()
            written.append(write_pdb(sequence, ca_coords, pred_path))

    print(f"Generated {len(written)} predictions in {out_dir}")
    return written
