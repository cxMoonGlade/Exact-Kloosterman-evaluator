# Exact Kloosterman evaluation over binary fields with variable rank

This repository contains the code and experimental results for the manuscript of the same title (6 October 2026). Given a binary-field representation, a nonzero field parameter and a rank, the implementations return an exact signed Kloosterman integer. `N` denotes the binary-field degree and `L` the rank.

The main implementations use Artin-Hasse (AH) and Bessel matrix recurrences. The repository also includes a full-table PARI cyclic-polynomial baseline, an independent PARI elliptic-curve point-counting baseline for low ranks, and the timings, exact integer outputs, statistical summaries and figures from completed experiments.

## Verify the results

From the repository root, with Python 3.9 or later:

```sh
python3 verify_results.py
```

This command uses only the Python standard library. It checks file hashes and exact-output consistency, and recomputes timing summaries from the saved measurements. Verification reads the distributed data; it does not run mathematical kernels or launch a benchmark. Receipt hashes link the exported data to the original records; the complete research archive is outside the scope of this release.

## Core code

| Component | Entry point |
| --- | --- |
| Shared Artin-Hasse coefficient construction | [block_acquisition.py](core/optimized_pair_20261004/reconstruction/block_second/block_acquisition.py), [block.c](core/optimized_pair_20261004/reconstruction/block_second/block.c) |
| Bessel matrix recurrence | [ordinary_int.py](core/optimized_pair_20261004/reconstruction/carrier_coverage/ordinary_int.py) |
| Common parameter construction | [parameter.py](core/optimized_pair_20261004/reconstruction/parameter_minpoly/parameter.py) |
| Specialization and ordered-matrix readout | [rank_query.py](core/optimized_pair_20261004/reconstruction/rank_extension/rank_query.py), [rank_resident.c](core/optimized_pair_20261004/reconstruction/rank_extension/rank_resident.c) |
| Expanded-grid readout | [grid_query.py](core/optimized_pair_20261004/reconstruction/grid_12x12_20261004/grid_query.py), [grid_resident.c](core/optimized_pair_20261004/reconstruction/grid_12x12_20261004/grid_resident.c) |
| Direct PARI point counting, ranks 2 and 3 | [direct_worker.py](core/pari_direct_20260930/direct_worker.py), [checked query launcher](core/pari_direct_20260930/query.py), [documentation](core/pari_direct_20260930/README.md) |
| PARI exact cyclic-polynomial baseline | [methods.py](core/pari_variable_rank_20260930/code/methods.py), function `run_pari` |
| Integer convolution | [convolution.py](core/rank_cost_optimized_20260930/reference/code/convolution.py) |
| Packed convolution and integer FWHT | [fast_convolution.py](core/rank_cost_optimized_20260930/code/fast_convolution.py), [fwht.py](core/rank_cost_optimized_20260930/code/fwht.py) |

[`core/`](core/) preserves the source files used in the experiments, their local dependencies and their original relative layout. [`core/SOURCE_MANIFEST.json`](core/SOURCE_MANIFEST.json) records the source hashes.

### Running the research implementations

The code is a frozen research snapshot. Fresh native runs require separately installed dependencies, updated paths and new local build identities. The original review and build checks are specific to the measurement machine, so they need adaptation for a new environment. The release includes source code; native binaries, historical review documents and unrelated research files are omitted.

The matrix methods were measured with python-flint 0.9.0 and FLINT 3.6.0, using one FLINT thread. The PARI environment used cypari2 2.2.4, cysignals 1.12.6 and PARI 2.17.3.

## Experimental results

[`results/`](results/) contains individual measurements, exact output integers, recomputed CSV summaries and figures. AH has the lower median runtime in all six initial cases. In the separate 25-point comparison, AH has the lower median at 9 points and Bessel at 16. Neither method has a uniform advantage across that comparison.

The main datasets time initialization and four exact queries per successful worker, from process launch through completion. The low-rank Direct PARI experiment uses one query per worker and records computation time and total process duration separately; see [PARI](#pari) for the timing definitions.

| Dataset | Inputs | Successful workers |
| --- | --- | ---: |
| Six-case AH/Bessel comparison | `N=8,16`; `L=16,24,32`; 2 methods × 3 repeats | 36 |
| Seven-method comparison, including PARI | `N,L ∈ {4,6,8,12,16}`; 7 methods × 3 repeats | 525 |
| Adjacent odd/even ranks | `N=8,16`; `L=3,4,7,8,15,16,23,24,31,32`; 4 methods × 3 repeats | 240 |
| Larger-field even ranks | `N=8,16,24,32,48,64`; `L=4,8,16,24,32,48,64`; 2 methods × 3 repeats | 252 |
| Separate `N=36` batch | `L=4,8,16,24,32,48,64`; 2 methods × 3 repeats | 42 |
| Larger-field Figure 1 extension | `N=24,32,48,64`; `L=4,6,8,12,16`; 2 methods × 3 repeats | 120 |
| Historical Direct PARI point counting | `N=2,4,6,8,12`; `L=2,3`; `a=1,2,3`; 3 repeats | 90 |

Some datasets share original observations, so the counts above cannot be summed as independent repetitions. Serial and parallel scheduling batches are labelled separately. The broader odd-rank scaling experiment and full 12×12 grid were incomplete at export and are excluded from the completed results. Agreement between large-field outputs is a consistency check, not an independent full-table reference.

## PARI

The two PARI baselines perform different computations and use different timing conventions.

### Full-table cyclic powering

The [full-table implementation](core/pari_variable_rank_20260930/code/methods.py) has a [25-point result table](results/pari_25_cells.csv) covering 75 successful workers and 300 exact integer outputs.

`run_pari` constructs the supplied binary field, finds a primitive element, computes a full integer cyclic-polynomial power and uses discrete logarithms to extract the requested coefficients. It uses cypari2 with an inline GP character-polynomial function, so no external `.gp` file is needed. The adapter accepts ranks `L ∈ {4,6,8,12,16}`.

The aligned comparison used an 8 MiB initial PARI stack, a 72 GiB stack maximum and an 80 GiB address-space limit per worker. Its wrapper bypassed `setup_pari()`, whose resource settings differ from those used in the measurements. The optional `run_shared()` helper implements a separate historical method and is outside this baseline.

### Low-rank point counting

The [Direct PARI worker](core/pari_direct_20260930/direct_worker.py) calls `ellinit` and `ellcard` for a single supplied parameter. At rank two, it returns the elliptic-curve point count minus `2^N + 1`; at rank three, it uses `Kl_3(a) = Kl_2(a)^2 - 2^N`. This computation requires no full cyclic-polynomial table.

The [saved results](results/pari_direct_20260930/) contain 90 successful runs on 30 cases, the original source identities and 180 records from the preceding AH/Sage-PARI comparison. The inputs were `N=2,4,6,8,12`, `L=2,3` and three nonzero parameter labels. The AH comparison used an older Python research reference; it does not compare Direct PARI with the current native implementations. An initial run that failed a metadata assertion is archived separately.

The computation time includes field, parameter and curve construction, point counting and integer recovery. Imports, PARI runtime initialization and the parent-observed process duration are recorded separately. This single-query computation time is not directly comparable with the four-query process wall time used for cyclic powering. See the [implementation, input limits and fresh-query commands](core/pari_direct_20260930/README.md) for details.

## Quantum simulation motivation

The computed integers also determine observables of the manuscript's specific ideal $N$-qubit evolution, which is built from a binary-field Gauss kernel. After $L$ steps, the amplitude from a nonzero computational-basis label $y$ to $x$ follows from the rank-$L$ Kloosterman integer at $a=x/y$, with a known normalization and additive correction.

The value at $a=1$ determines the full trace, including the fixed zero-basis state. Its squared modulus gives the unnormalized spectral form factor. These algebraic relations give individual amplitudes and spectral quantities exactly without constructing the full state vector.

## License

Project code is distributed under the [MIT license](LICENSE). FLINT, python-flint, PARI and their dependencies are separate software with their own licenses.
