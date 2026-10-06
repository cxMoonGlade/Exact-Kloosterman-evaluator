# Exact Kloosterman evaluation over binary fields with variable rank

Core code and experimental results for **Exact Kloosterman evaluation over binary fields with variable rank** (6 October 2026).

本仓库收录 AH、Bessel 矩阵递推、PARI 及其他全表基线的核心代码，以及已完成实验的逐次计时、精确整数输出、统计表与图表。

## Core code

The input specifies a binary field representation, a nonzero field parameter and a rank; the output is an exact signed Kloosterman integer. `N` is the binary-field degree and `L` is the rank.

| Component | Entry point |
| --- | --- |
| Shared Artin–Hasse coefficient construction | [block_acquisition.py](core/optimized_pair_20261004/reconstruction/block_second/block_acquisition.py), [block.c](core/optimized_pair_20261004/reconstruction/block_second/block.c) |
| Bessel matrix recurrence | [ordinary_int.py](core/optimized_pair_20261004/reconstruction/carrier_coverage/ordinary_int.py) |
| Common parameter construction | [parameter.py](core/optimized_pair_20261004/reconstruction/parameter_minpoly/parameter.py) |
| Specialization and ordered-matrix readout | [rank_query.py](core/optimized_pair_20261004/reconstruction/rank_extension/rank_query.py), [rank_resident.c](core/optimized_pair_20261004/reconstruction/rank_extension/rank_resident.c) |
| Expanded-grid readout | [grid_query.py](core/optimized_pair_20261004/reconstruction/grid_12x12_20261004/grid_query.py), [grid_resident.c](core/optimized_pair_20261004/reconstruction/grid_12x12_20261004/grid_resident.c) |
| PARI exact cyclic-polynomial baseline | [methods.py](core/pari_variable_rank_20260930/code/methods.py), function `run_pari` |
| Integer convolution | [convolution.py](core/rank_cost_optimized_20260930/reference/code/convolution.py) |
| Packed convolution and integer FWHT | [fast_convolution.py](core/rank_cost_optimized_20260930/code/fast_convolution.py), [fwht.py](core/rank_cost_optimized_20260930/code/fwht.py) |

[`core/`](core/) preserves the measured mathematical source files and their local dependencies in the original relative layout. Source hashes are recorded in [`core/SOURCE_MANIFEST.json`](core/SOURCE_MANIFEST.json). These are frozen research implementations: fresh native runs need installed dependencies, relocated paths and new local build identities. The original machine-specific review/build gates are not a portable installer. Native binaries, historical review documents and unrelated research files are not distributed.

The recorded matrix-method environment used python-flint 0.9.0 and FLINT 3.6.0, with one FLINT thread. PARI used cypari2 2.2.4, cysignals 1.12.6 and PARI 2.17.3. Third-party libraries must be installed separately.

## Experimental results

[`results/`](results/) contains individual measurements and exact output integers, together with the recomputed CSV summaries and result figures. Every successful timed worker performs initialization and four exact queries; the main wall clock includes process launch through completion.

| Dataset | Inputs | Successful workers |
| --- | --- | ---: |
| Six-case AH/Bessel comparison | `N=8,16`; `L=16,24,32`; 2 methods × 3 repeats | 36 |
| Seven-method comparison, including PARI | `N,L ∈ {4,6,8,12,16}`; 7 methods × 3 repeats | 525 |
| Adjacent odd/even ranks | `N=8,16`; `L=3,4,7,8,15,16,23,24,31,32`; 4 methods × 3 repeats | 240 |
| Larger-field even ranks | `N=8,16,24,32,48,64`; `L=4,8,16,24,32,48,64`; 2 methods × 3 repeats | 252 |
| Separate `N=36` batch | `L=4,8,16,24,32,48,64`; 2 methods × 3 repeats | 42 |
| Larger-field Figure 1 extension | `N=24,32,48,64`; `L=4,6,8,12,16`; 2 methods × 3 repeats | 120 |

Some datasets reuse the same original observations; these counts must not be added as independent repetitions. Timings from the serial and parallel scheduling batches remain labelled separately. The broader odd-rank scaling campaign and full 12×12 grid were incomplete at export and are not represented as completed experiments.

AH has the lower median in the six initial cases. In the separate 25-point comparison, AH has the lower median at 9 points and Bessel at 16. The data do not establish a uniform AH advantage. Large-field output consistency does not constitute an independent full-table oracle.

## PARI

PARI is explicitly included: [core implementation](core/pari_variable_rank_20260930/code/methods.py) and [25-point result table](results/pari_25_cells.csv), with **75 successful workers and 300 exact integer outputs**.

`run_pari` uses cypari2 and an inline GP character-polynomial function. It constructs the supplied binary field, finds a primitive element, computes a full integer cyclic-polynomial power and extracts the requested coefficients using discrete logarithms. No external `.gp` file is required. The registered adapter accepts ranks `L ∈ {4,6,8,12,16}`.

The measured aligned campaign initialized an 8 MiB PARI stack with a 72 GiB maximum, under an 80 GiB per-worker address-space limit. Its wrapper bypassed the older `setup_pari()` helper, which has different resource settings. The optional `run_shared()` helper in that source belongs to a different historical method and is not part of the PARI baseline.

## Verify the exported results

From the repository root, using Python 3.9 or later:

```sh
python3 verify_results.py
```

This standard-library command checks the exported file hashes, exact-output consistency and timing summaries using only the distributed measurements. It does not execute mathematical kernels or start a new benchmark. Source receipt hashes provide provenance back to the original records; the full historical research archive is not included in this focused code-and-results release.

## License

Project code is distributed under the existing [MIT license](LICENSE). FLINT, python-flint, PARI and their other dependencies are external software with their own licenses.
