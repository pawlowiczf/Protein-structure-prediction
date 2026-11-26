"""Worker script executed inside tools/ost-env to compute CA-only lDDT."""

import json
import sys

import ost
from ost import io
from ost.mol.alg import lddt


def compute_ca_lddt(pred_path: str, native_path: str) -> float | None:
    ost.PushVerbosityLevel(0)
    model = io.LoadPDB(pred_path)
    target = io.LoadPDB(native_path)
    mapping = {target.chains[0].name: model.chains[0].name}
    scorer = lddt.lDDTScorer(target, bb_only=True)
    return scorer.lDDT(model, chain_mapping=mapping)


if __name__ == "__main__":
    score = compute_ca_lddt(sys.argv[1], sys.argv[2])
    print(json.dumps({"lddt": score}))
