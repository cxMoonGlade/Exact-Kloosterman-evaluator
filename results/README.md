# Arithmetic benchmark data

This directory contains timings, exact hyper-Kloosterman integers and figures for the arithmetic evaluators, together with low-rank Direct PARI measurements. The rank field `L` corresponds to $r=\ell-1$ in the [half-period quantum-walk task](../docs/quantum_task.md).

| Dataset label | Worker appearances | Exact outputs | Scope |
| --- | ---: | ---: | --- |
| `six_case` | 36 | 144 | N=8,16; L=16,24,32; AH and Bessel |
| `aligned_reference` | 150 | 600 | N,L=4,6,8,12,16; retained AH and Bessel |
| `aligned_baseline` | 375 | 1500 | Same 25 inputs; five full-table baselines including PARI |
| `odd_even` | 240 | 960 | N=8,16; ten odd/even ranks; four methods |
| `even_serial` | 252 | 1008 | N=8,16,24,32,48,64; seven even ranks; AH and Bessel |
| `even_n36` | 42 | 168 | N=36; seven even ranks; AH and Bessel |
| `figure1_extension` | 120 | 480 | N=24,32,48,64; L=4,6,8,12,16; AH and Bessel |

These four-query datasets contain 1,215 recorded appearances of 1,107 original worker runs. The 108 runs shared between datasets carry the same `sample_id`; use this field to deduplicate measurements when combining tables. PARI contributes 75 successful workers and 300 exact outputs in `aligned_baseline`.

Method labels match the source code. `ordinary_int` and `ordinary_new` denote Bessel matrix recurrence; `ah_block_second` and `ah_new` denote the corresponding AH implementation.

## Low-rank Direct PARI measurements

[`pari_direct_20260930/`](pari_direct_20260930/) contains 90 Direct PARI single-query runs, their summaries and source identities, plus 180 earlier Python-AH/Sage-PARI comparison records. These measurements and their validation records are stored separately from the four-query datasets above. The AH comparator in this experiment is the earlier Python implementation.

Direct PARI computation time includes field, parameter and curve construction, point counting and integer recovery. Import time, runtime initialization and elapsed time observed by the launching process are recorded separately.

## Four-query dataset formats

`exact_outputs.csv` stores each result as a signed decimal integer, identified by its field polynomial `f`, parameter `a`, query position, rank and field degree. It also records the method, repeat, dataset, original attempt ID and distributed source hash. The gzip copy contains the same bytes. Read integer fields with arbitrary precision to preserve every output digit.

Output comparisons match the full `(dataset, N, L, f, a)` identity. At large field degrees, they check consistency between the recorded implementations.

`timing_samples.csv` records elapsed seconds from process launch through completion, including startup, coefficient acquisition and four queries. The `receipt_sha256` field identifies the distributed parent receipt; the six-case dataset stores its receipts in one aggregate stream. `timing_summary.csv` gives the median, minimum and maximum of each set of three repetitions. Its ranges are the observed extrema. Original summary tables keep their stated units, including milliseconds in `six_case_times_ms.csv`.

`output_sources.csv` maps sample IDs to original and distributed raw-stream hashes and logical source paths. `EXTRACTION_CHECK.json` records counts, source-manifest identities and result-file hashes. Other `CHECK.json` files summarize earlier verification against separately retained local source archives. Measurements retain their original scheduling labels; older serial runs have no recorded CPU model.

The three SVG figures reproduce the original supplement. The aligned figure contains 175 method/input cells, with the later 40-cell extension supplied separately as CSV.

## Rebuild the tables from source records

With a local copy of the original source archives, regenerate the compact data using the Python standard library:

```sh
python3 extract_results.py --artifact /path/to/artifact \
  --extension /path/to/artifact_extension \
  --pari-results /path/to/pari_results --output /path/to/new_results
```

The input directories must contain the original distributed manifests and saved statistics. The extractor checks source hashes, receipt/output agreement, dataset counts and integer consistency, then regenerates timing summaries and copies the result tables and figures. It reads the saved measurements and writes the extraction to `--output`.
