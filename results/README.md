# Retained experimental results

This directory contains compact result tables, exact integer outputs and figures for the completed experiments. It omits the full local receipt archive, experiment protocols, historical development audits and runtime installations.

| Dataset label | Worker appearances | Exact outputs | Scope |
| --- | ---: | ---: | --- |
| `six_case` | 36 | 144 | N=8,16; L=16,24,32; AH and Bessel |
| `aligned_reference` | 150 | 600 | N,L=4,6,8,12,16; retained AH and Bessel |
| `aligned_baseline` | 375 | 1500 | Same 25 inputs; five full-table baselines including PARI |
| `odd_even` | 240 | 960 | N=8,16; ten odd/even ranks; four methods |
| `even_serial` | 252 | 1008 | N=8,16,24,32,48,64; seven even ranks; AH and Bessel |
| `even_n36` | 42 | 168 | N=36; seven even ranks; AH and Bessel |
| `figure1_extension` | 120 | 480 | N=24,32,48,64; L=4,6,8,12,16; AH and Bessel |

These are 1,215 dataset appearances of 1,107 original worker streams. Exactly 108 original streams occur in more than one dataset. Repeated `sample_id` values identify reused measurements; do not count them as independent repetitions. PARI contributes 75 successful workers and 300 exact outputs in `aligned_baseline`. Method labels retain the original code names: `ordinary_int` and `ordinary_new` denote Bessel matrix recurrence; `ah_block_second` and `ah_new` denote the corresponding AH implementation.

`exact_outputs.csv` contains decimal signed integers and their field polynomial `f`, parameter `a`, query position, rank, degree, method, repeat, dataset, original attempt ID and distributed source hash. The gzip copy contains the same bytes. Parse integer fields with arbitrary precision: spreadsheet or floating-point conversions can lose digits. Comparisons during extraction were restricted to the same `(dataset, N, L, f, a)`. Agreement on large fields is not an independent full-table validation.

`timing_samples.csv` records every parent process wall clock in seconds, including startup, acquisition, four queries and process completion. `receipt_sha256` identifies the distributed parent receipt (the six-case campaign stores its receipts in one aggregate stream). `timing_summary.csv` recomputes each three-repeat median, minimum and maximum from these clocks. Ranges are observed ranges, not confidence intervals. The separately named original summary tables retain their existing units; `six_case_times_ms.csv` uses milliseconds.

`output_sources.csv` maps compact sample IDs to the full original and distributed raw-stream hashes and logical source paths. `EXTRACTION_CHECK.json` gives counts, source-manifest identities and copied result-file hashes. Other `CHECK.json` files are retained outputs of earlier complete local statistics verification; their referenced full source archives are not included here. Timing batches have different scheduling conditions; old serial workers did not retain a CPU-model field.

The three SVG figures are frozen outputs of the original complete supplement. The aligned figure shows 175 method/input cells; the later 40-cell extension is supplied separately in CSV form. The wider odd-rank parity-scaling campaign was incomplete at export and contributes no results here.

To reconstruct this compact extraction from separately retained complete local snapshots, using only the Python standard library:

```sh
python3 extract_results.py --artifact /path/to/artifact \
  --extension /path/to/artifact_extension \
  --pari-results /path/to/pari_results --output /path/to/new_results
```

The source directories need the original distributed manifests and saved statistics. The extractor verifies every used receipt/raw-stream hash and exact receipt/raw agreement, verifies fixed completed counts and within-dataset integer consistency, regenerates timing summaries, and copies the already verified result tables and figures. It imports no scientific implementation, performs no mathematical experiment, and leaves the source snapshots unchanged.
