# Direct PARI point counting at ranks two and three

This is the actual direct-cypari2 worker used in the completed 30 September 2026 low-rank experiment. It is separate from the [full-table PARI cyclic-polynomial method](../pari_variable_rank_20260930/code/methods.py).

## Task and measured scope

The recorded input consists of a supplied binary-field polynomial, a nonzero parameter label `a`, and rank `L`. The worker constructs the field and parameter, calls `ellinit([1, 0, 0, a, 0], generator)` and `ellcard()` for `E_a: y^2 + xy = x^3 + ax`, and returns `#E_a - 2^N - 1` at rank two. At rank three it applies `Kl_3(a) = Kl_2(a)^2 - 2^N`. It computes one integer for the requested parameter without constructing a field-wide Kloosterman table.

The measured grid is `N in {2,4,6,8,12}`, `L in {2,3}`, `a in {1,2,3}`, with three fresh processes per case: **30 cases and 90 successful Direct PARI workers**. The field polynomials are in `direct_reference.MODULI`. The nonzero labels are polynomial-basis coordinates; `a=2` denotes the supplied field generator, not the integer two embedded in characteristic two.

The original comparison used a Python AH research reference and a Sage/PARI wrapper. These are historical measurements, not a matched comparison with this repository's current native AH/Bessel implementations. The Direct PARI refinement followed the original 180-worker comparison and reused its inputs. The archive supports exact-output and finite-implementation checks at these low ranks; it supplies no variable-rank speedup claim.

## Files and provenance

- [direct_worker.py](direct_worker.py) and [direct_reference.py](direct_reference.py) are copied byte-for-byte from the measured implementation and its field/reference dependency.
- [direct_run.py](direct_run.py) is the unchanged historical orchestration source. It binds the original Conda Python executable, expects original results alongside the script, and rewrites result files. It is retained for provenance; use the command below to run a fresh individual query.
- [direct_worker_failed.py](direct_worker_failed.py) preserves an earlier metadata-assertion error after computation. Its single failed attempt is retained separately and is not part of the 90 successful timings.
- [query.py](query.py) is a new input-checking launcher around the unchanged worker. It selects the current Python interpreter and rejects inputs outside the recorded grid. It does not rewrite archived results and was not used for the saved timings. The frozen worker itself does not reject unsupported ranks, so use this launcher for new queries.
- Original and distributed source hashes are recorded in [SOURCE_MANIFEST.json](../SOURCE_MANIFEST.json). Measurements, summary statistics, saved validation, software versions and original identity records are in [results/pari_direct_20260930](../../results/pari_direct_20260930/).

## Run a fresh query

Use a Python environment with cypari2 and PARI installed. The recorded environment used Python 3.13.15, cypari2 2.2.4 and PARI 2.17.3. Sage is not required by the Direct PARI worker. The worker uses Unix `resource` accounting; the recorded platform was Linux.

From the repository root:

```sh
python core/pari_direct_20260930/query.py 8 2 3
python core/pari_direct_20260930/query.py 8 3 3
```

Each command starts one fresh worker and writes one JSON result to standard output. Its value can be checked against the archived records. Newly produced timings are separate observations and do not replace the saved campaign.

## Timing boundary

`compute_seconds` includes polynomial/field/parameter construction, curve construction, point counting, conversion to the output integer and the rank-three identity. It excludes interpreter imports, PARI runtime initialization, metadata assertions and JSON output. Irreducibility of the supplied polynomials was checked in the earlier validation; finding a polynomial and repeating that validation are outside this clock.

`entry_import_runtime_seconds` records the in-worker import/runtime interval. `process_seconds` in the archive is the parent-observed process duration, including startup, completion and JSON parsing. `max_rss_kib` is whole-process peak resident memory, including interpreter and libraries. These definitions differ from the separate four-query full-table campaign. Short timings are descriptive measurements; scheduling and host activity were not controlled enough to infer a universal ordering.

The standard-library command `python3 verify_results.py` checks the saved outputs, repetitions, stage totals, reconstructed statistics and distributed hashes. It does not rerun the mathematical kernels or a benchmark.
