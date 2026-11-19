"""
Run PlaToLoCo on a large FASTA file in batches, one job at a time, and merge
the results. Uses the functions from the group's run_platoloco.py, which must
be in the same folder.

Usage:
    python3 run_platoloco_batches.py input.fasta OUTPUT_DIR [batch_size]

Why batches: one failed job no longer loses the whole run, and an interrupted
run can be restarted with the same command - finished batches are skipped.
Jobs are never run in parallel (the PlaToLoCo tools share temporary files).

Outputs in OUTPUT_DIR:
  batches/batch_NNN/      results of each batch (CSV + FASTA per method)
  <method>.csv / .fasta   merged results for SEG_strict, SEG_intermediate,
                          CAST and fLPS_strict
  platoloco_request.json  parameters used (identical for every batch)
  batches.tsv             proteins and regions per batch
"""
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_platoloco import (make_payload, submit_job, wait_for_job,  # noqa: E402
                           get_rows, write_outputs)

METHODS = ["SEG_strict", "SEG_intermediate", "CAST", "fLPS_strict"]


def read_records(path):
    """Split a FASTA file into a list of record strings."""
    text = Path(path).read_text().strip()
    return [">" + rec.strip() + "\n" for rec in text.split(">") if rec.strip()]


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit("Usage: python3 run_platoloco_batches.py input.fasta OUTPUT_DIR [batch_size]")
    records = read_records(sys.argv[1])
    out = Path(sys.argv[2])
    size = int(sys.argv[3]) if len(sys.argv) == 4 else 1000
    batch_root = out / "batches"
    batch_root.mkdir(parents=True, exist_ok=True)

    n_batches = (len(records) + size - 1) // size
    print(f"{len(records)} proteins in {n_batches} batches of up to {size}")
    (out / "platoloco_request.json").write_text(
        json.dumps({**make_payload(">x\nA\n"), "sequences": "<per batch>"}, indent=2))

    summary = []
    for i in range(n_batches):
        bdir = batch_root / f"batch_{i + 1:03d}"
        done = bdir / "DONE"
        chunk = records[i * size:(i + 1) * size]
        if not done.exists():
            print(f"\n=== Batch {i + 1}/{n_batches} ({len(chunk)} proteins) ===")
            bdir.mkdir(exist_ok=True)
            (bdir / "input.fasta").write_text("".join(chunk))
            token = submit_job(make_payload("".join(chunk)))
            print(f"Submitted job: {token}")
            wait_for_job(token)
            rows = get_rows(token)
            write_outputs(rows, bdir)
            done.write_text(f"{token}\n")
        else:
            print(f"Batch {i + 1}/{n_batches} already done, skipped")
        n_reg = sum(len(pd.read_csv(bdir / f"{m}.csv")) for m in METHODS)
        summary.append((i + 1, len(chunk), n_reg))

    # Merge batches
    for m in METHODS:
        parts = [pd.read_csv(d / f"{m}.csv") for d in sorted(batch_root.glob("batch_*"))]
        pd.concat(parts, ignore_index=True).to_csv(out / f"{m}.csv", index=False)
        with open(out / f"{m}.fasta", "w") as fh:
            for d in sorted(batch_root.glob("batch_*")):
                fh.write((d / f"{m}.fasta").read_text())

    pd.DataFrame(summary, columns=["batch", "proteins", "regions_all_methods"]) \
        .to_csv(out / "batches.tsv", sep="\t", index=False)
    print(f"\nMerged results written to {out}")
    for m in METHODS:
        t = pd.read_csv(out / f"{m}.csv")
        print(f"  {m:17s} {len(t):7d} regions in {t['protein_id'].nunique():6d} proteins")


if __name__ == "__main__":
    main()
