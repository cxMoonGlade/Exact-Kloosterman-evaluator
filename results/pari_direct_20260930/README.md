# Saved low-rank Direct PARI measurements

These files preserve the completed 30 September 2026 point-counting experiment. They are separate from the four-query CSV datasets in the parent directory and are not counted in those datasets' worker/output totals.

| File | Contents |
| --- | --- |
| `direct_raw.jsonl` | 90 successful Direct PARI workers: 30 `(N,L,a)` cases, three repetitions each, one exact integer per worker |
| `direct_summary.json` | Original per-case compute/process/RSS medians, minima and maxima |
| `raw.jsonl` | 180 preceding historical workers: 90 Python AH research-reference and 90 Sage/PARI runs on the same cases |
| `validation.json` | Saved field checks, original AH-adapter checks and 18 independent direct-sum checks at `N=2,4,6` |
| `direct_identity.json` | Original hashes of the successful worker, runner, comparison records and environment; also an amendment-document hash |
| `environment.json`, `software_versions.json` | Original environment and dependency/source identities |
| `direct_failed_attempt.jsonl`, `direct_failed_identity.json` | One earlier metadata-assertion failure and its source identities, excluded from successful timings |

All listed files are copied byte-for-byte. Machine-local paths in the historical environment and failure receipt identify the original execution, not prerequisites for the new query launcher. The original amendment document and historical AH implementation archive are not distributed in this focused addition; their recorded hashes remain provenance references. No third-party paper text or library source is copied here. The failed identity's `direct_worker.py` hash corresponds to the archived `core/pari_direct_20260930/direct_worker_failed.py`.

The Direct PARI clock includes field/parameter/curve construction and exact recovery, but excludes imports and PARI runtime initialization. See the [implementation and timing documentation](../../core/pari_direct_20260930/README.md). Original AH records are comparison evidence for that older Python implementation; no matched current-native-AH/Bessel curve is implied.

Run `python3 verify_results.py` from the repository root to recheck all 90 Direct PARI records against the 180 historical outputs and the 18 saved direct-sum values, reconstruct the summary and check source identities. This verifies archived evidence; it does not rerun validation kernels, execute the old AH implementation, or regenerate timing measurements.
