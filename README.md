# Exact Kloosterman evaluation over binary fields with variable rank

Core code and experimental results for **Exact Kloosterman evaluation over binary fields with variable rank** (6 October 2026).

This repository contains Artin–Hasse (AH) and Bessel matrix recurrence implementations, a full-table PARI cyclic-polynomial baseline, an independent low-rank PARI elliptic-curve point-counting baseline, and per-run timings, exact integer outputs, statistical summaries, and figures from completed experiments.

## Quantum simulation motivation

The integers computed here also determine observables of the manuscript's specific ideal $N$-qubit evolution built from a binary-field Gauss kernel. After $L$ steps, the amplitude from a nonzero computational-basis label $y$ to $x$ is obtained from the rank-$L$ Kloosterman integer at $a=x/y$ by a known normalization and additive correction. The value at $a=1$ determines the full trace, including the fixed zero-basis state, and hence the unnormalized spectral form factor (the squared modulus of that trace). This algebraic post-processing gives individual amplitudes and spectral quantities exactly without constructing the full state vector.

## Core code

The input specifies a binary field representation, a nonzero field parameter and a rank; the output is an exact signed Kloosterman integer. `N` is the binary-field degree and `L` is the rank.

| Component | Entry point |
| --- | --- |
| Shared Artin–Hasse coefficient construction | [block_acquisition.py](core/optimized_pair_20261004/reconstruction/block_second/block_acquisition.py), [block.c](core/optimized_pair_20261004/reconstruction/block_second/block.c) |
| Bessel matrix recurrence | [ordinary_int.py](core/optimized_pair_20261004/reconstruction/carrier_coverage/ordinary_int.py) |
| Common parameter construction | [parameter.py](core/optimized_pair_20261004/reconstruction/parameter_minpoly/parameter.py) |
| Specialization and ordered-matrix readout | [rank_query.py](core/optimized_pair_20261004/reconstruction/rank_extension/rank_query.py), [rank_resident.c](core/optimized_pair_20261004/reconstruction/rank_extension/rank_resident.c) |
| Expanded-grid readout | [grid_query.py](core/optimized_pair_20261004/reconstruction/grid_12x12_20261004/grid_query.py), [grid_resident.c](core/optimized_pair_20261004/reconstruction/grid_12x12_20261004/grid_resident.c) |
| Direct PARI point counting, ranks 2 and 3 | [direct_worker.py](core/pari_direct_20260930/direct_worker.py), [checked query launcher](core/pari_direct_20260930/query.py), [documentation](core/pari_direct_20260930/README.md) |
| PARI exact cyclic-polynomial baseline | [methods.py](core/pari_variable_rank_20260930/code/methods.py), function `run_pari` |
| Integer convolution | [convolution.py](core/rank_cost_optimized_20260930/reference/code/convolution.py) |
| Packed convolution and integer FWHT | [fast_convolution.py](core/rank_cost_optimized_20260930/code/fast_convolution.py), [fwht.py](core/rank_cost_optimized_20260930/code/fwht.py) |

[`core/`](core/) preserves the measured mathematical source files and their local dependencies in the original relative layout. Source hashes are recorded in [`core/SOURCE_MANIFEST.json`](core/SOURCE_MANIFEST.json). These are frozen research implementations: fresh native runs need installed dependencies, relocated paths and new local build identities. The original machine-specific review/build gates are not a portable installer. Native binaries, historical review documents and unrelated research files are not distributed.

The recorded matrix-method environment used python-flint 0.9.0 and FLINT 3.6.0, with one FLINT thread. PARI used cypari2 2.2.4, cysignals 1.12.6 and PARI 2.17.3. Third-party libraries must be installed separately.

## Experimental results

[`results/`](results/) contains individual measurements and exact output integers, together with the recomputed CSV summaries and result figures. The main datasets below use initialization and four exact queries per successful timed worker, with a wall clock from process launch through completion. The separately archived low-rank Direct PARI campaign uses one query per worker and records both mathematical computation and process duration; its timing definitions are described under [PARI](#pari).

| Dataset | Inputs | Successful workers |
| --- | --- | ---: |
| Six-case AH/Bessel comparison | `N=8,16`; `L=16,24,32`; 2 methods × 3 repeats | 36 |
| Seven-method comparison, including PARI | `N,L ∈ {4,6,8,12,16}`; 7 methods × 3 repeats | 525 |
| Adjacent odd/even ranks | `N=8,16`; `L=3,4,7,8,15,16,23,24,31,32`; 4 methods × 3 repeats | 240 |
| Larger-field even ranks | `N=8,16,24,32,48,64`; `L=4,8,16,24,32,48,64`; 2 methods × 3 repeats | 252 |
| Separate `N=36` batch | `L=4,8,16,24,32,48,64`; 2 methods × 3 repeats | 42 |
| Larger-field Figure 1 extension | `N=24,32,48,64`; `L=4,6,8,12,16`; 2 methods × 3 repeats | 120 |
| Separate historical Direct PARI point counting | `N=2,4,6,8,12`; `L=2,3`; `a=1,2,3`; 3 repeats | 90 |

Some datasets reuse the same original observations; these counts must not be added as independent repetitions. Timings from the serial and parallel scheduling batches remain labelled separately. The broader odd-rank scaling campaign and full 12×12 grid were incomplete at export and are not represented as completed experiments.

AH has the lower median in the six initial cases. In the separate 25-point comparison, AH has the lower median at 9 points and Bessel at 16. The data do not establish a uniform AH advantage. Large-field output consistency does not constitute an independent full-table oracle.

## PARI

Two distinct PARI baselines are included.

### Full-table cyclic powering

The full-table baseline has [core implementation](core/pari_variable_rank_20260930/code/methods.py) and [25-point result table](results/pari_25_cells.csv), with **75 successful workers and 300 exact integer outputs**.

`run_pari` uses cypari2 and an inline GP character-polynomial function. It constructs the supplied binary field, finds a primitive element, computes a full integer cyclic-polynomial power and extracts the requested coefficients using discrete logarithms. No external `.gp` file is required. The registered adapter accepts ranks `L ∈ {4,6,8,12,16}`.

The measured aligned campaign initialized an 8 MiB PARI stack with a 72 GiB maximum, under an 80 GiB per-worker address-space limit. Its wrapper bypassed the older `setup_pari()` helper, which has different resource settings. The optional `run_shared()` helper in that source belongs to a different historical method and is not part of the PARI baseline.

### Low-rank point counting

The actual [Direct PARI worker](core/pari_direct_20260930/direct_worker.py) calls `ellinit` and `ellcard` for one supplied parameter. At rank two it returns the elliptic-curve point count minus `2^N + 1`; rank three uses `Kl_3(a) = Kl_2(a)^2 - 2^N`. It does not construct the full cyclic-polynomial table.

The [saved results](results/pari_direct_20260930/) contain **90 successful runs on 30 cases**, the original source identities, and the preceding 180 historical AH/Sage-PARI comparison records. This campaign used `N=2,4,6,8,12`, `L=2,3`, and three nonzero parameter labels. Its AH comparison is with an older Python research reference, not the current native implementations. The initial failed metadata assertion is also retained separately.

The low-rank computation clock includes field, parameter and curve construction, point counting and integer recovery; it excludes imports and PARI runtime initialization. Those costs and the parent process duration are recorded separately. Do not mix this single-query clock with the four-query cyclic-powering wall clock. [Implementation, input limits and fresh-query commands](core/pari_direct_20260930/README.md) document the distinction.

## Verify the exported results

From the repository root, using Python 3.9 or later:

```sh
python3 verify_results.py
```

This standard-library command checks the exported file hashes, exact-output consistency and timing summaries using only the distributed measurements. It does not execute mathematical kernels or start a new benchmark. Source receipt hashes provide provenance back to the original records; the full historical research archive is not included in this focused code-and-results release.

## License

Project code is distributed under the existing [MIT license](LICENSE). FLINT, python-flint, PARI and their other dependencies are external software with their own licenses.
