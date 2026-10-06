"""Execute the additive direct-PARI amendment without replacing v1 artifacts."""
from pathlib import Path
from time import perf_counter
import hashlib
import json
import os
import statistics
import subprocess
from direct_reference import MODULI

HERE = Path(__file__).resolve().parent
PYTHON = "/home/cx/miniforge3/envs/ecs-qcb-baselines/bin/python"


def main():
    old = [json.loads(line) for line in (HERE / "raw.jsonl").read_text().splitlines()]
    assert all(row["status"] == "PASS" for row in old)
    expected = {(row["N"], row["L"], row["a"]): row["value"] for row in old}
    validation = json.loads((HERE / "validation.json").read_text())
    direct = {(row["N"], row["L"], row["a"]): row["direct_sum"]
              for row in validation["direct_checks"]}
    identity_paths = [HERE / name for name in ("AMENDMENT_DIRECT_PARI.md", "direct_worker.py",
                                              "direct_run.py", "raw.jsonl", "environment.json")]
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in identity_paths}
    (HERE / "direct_identity.json").write_text(json.dumps(hashes, indent=2) + "\n")
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    rows = []
    started = perf_counter()
    with (HERE / "direct_raw.jsonl").open("w") as output:
        for N in MODULI:
            for L in (2, 3):
                for a in (1, 2, 3):
                    for repetition in range(3):
                        start = perf_counter()
                        command = [PYTHON, str(HERE / "direct_worker.py"), str(N), str(L), str(a)]
                        try:
                            run = subprocess.run(command, capture_output=True, text=True,
                                                 env=env, timeout=120, check=False)
                            row = (json.loads(run.stdout) if run.returncode == 0 else
                                   {"status": "ERROR", "stdout": run.stdout, "stderr": run.stderr,
                                    "returncode": run.returncode})
                        except subprocess.TimeoutExpired as error:
                            row = {"status": "TIMEOUT", "stdout": str(error.stdout), "stderr": str(error.stderr)}
                        row.update(N=N, L=L, a=a, algorithm="pari_direct", repetition=repetition,
                                   process_seconds=perf_counter() - start)
                        if row["status"] == "PASS":
                            key = (N, L, a)
                            row["matches_original_algorithms"] = row["value"] == expected[key]
                            row["matches_direct"] = None if key not in direct else row["value"] == direct[key]
                            if not row["matches_original_algorithms"] or row["matches_direct"] is False:
                                row["status"] = "VALUE_MISMATCH"
                        output.write(json.dumps(row) + "\n")
                        output.flush()
                        rows.append(row)
                        if row["status"] != "PASS":
                            raise RuntimeError(f"Direct baseline failure: {row}")
                    print(json.dumps({"finished_case": [N, L, a], "runs": len(rows)}), flush=True)
    groups = []
    for N in MODULI:
        for L in (2, 3):
            for a in (1, 2, 3):
                selected = [row for row in rows if (row["N"], row["L"], row["a"]) == (N, L, a)]
                group = {"N": N, "L": L, "a": a, "algorithm": "pari_direct",
                         "value": selected[0]["value"], "passed_repetitions": len(selected)}
                for key in ("compute_seconds", "process_seconds", "max_rss_kib"):
                    metric = [row[key] for row in selected]
                    group[key] = {"median": statistics.median(metric), "min": min(metric), "max": max(metric)}
                groups.append(group)
    summary = {"status": "PASS", "workers": len(rows), "case_count": len(expected),
               "seconds": perf_counter() - started, "groups": groups}
    (HERE / "direct_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "groups"}))


if __name__ == "__main__":
    main()
