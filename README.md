# Exact hyper-Kloosterman evaluation for half-period quantum walks on norm graphs

AH (Artin–Hasse) computes hyper-Kloosterman sums exactly over binary fields, allowing both the field degree and rank to grow. These integers determine position probabilities for quantum walks on norm graphs at half their minimum operator period.

This repository contains the arithmetic evaluator, comparison methods and saved benchmark results. Given a binary-field representation, a nonzero field parameter and a rank, the evaluator returns a signed integer. The [quantum task description](docs/quantum_task.md) explains how these evaluations give exact probabilities and samples of complete vertex labels.

## Quantum task

Let $q=2^N$ and $E=\mathbb F_{q^\ell}$, with $\ell\ge2$. The norm graph joins $u,v\in E$ when $N_{E/F_q}(u+v)=1$. We consider its continuous-time Laplacian walk from any initial vertex at the half-period $t_*=\pi/2^\ell$.

Each nonzero-displacement amplitude is a normalized hyper-Kloosterman sum of rank $r=\ell-1$, evaluated at the displacement's norm. Exact evaluation therefore gives individual output probabilities. Combining these values with Deligne's bound gives a rejection sampler for the complete position distribution.

Including preprocessing, probability queries take deterministic polynomial bit time and sampling takes expected polynomial bit time, jointly in $N$ and $\ell$. The cases $N=1$ and $\ell=2$ use direct samplers. The [task description](docs/quantum_task.md) gives the formulas, expected query count and total costs.

In the code and data, `N` denotes the binary-field degree and `L` the arithmetic rank: $L=r=\ell-1$. Each sampled position has $n=N\ell$ bits.

## Verify the results

From the repository root, with Python 3.9 or later:

```sh
python3 verify_results.py
```

The command uses the Python standard library to check file hashes, compare saved exact outputs and recompute timing summaries. It reads the distributed measurements without running the arithmetic kernels.

## Core code

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

[`core/`](core/) contains the benchmarked source files and their local dependencies in the original directory layout. [`core/SOURCE_MANIFEST.json`](core/SOURCE_MANIFEST.json) records their hashes.

### Running the research implementations

To run the native implementations, install the dependencies and compile the included sources. The build scripts retain paths, build identities and review-file checks from the original environment; adapt these settings to your machine before compiling.

The matrix methods were measured with python-flint 0.9.0 and FLINT 3.6.0, using one FLINT thread. The PARI environment used cypari2 2.2.4, cysignals 1.12.6 and PARI 2.17.3.

## Arithmetic benchmark results

The [arithmetic benchmarks](results/README.md) compare AH with Bessel matrix recurrences, convolution and Walsh–Hadamard tables, PARI cyclic powering, and low-rank point counting. AH has the lower median runtime in all six initial cases. In the separate 25-point AH/Bessel comparison, AH has the lower median at 9 points and Bessel at 16.

The main datasets measure process time from launch through initialization, four exact queries and completion. Direct PARI uses one query per process, with computation time and total process time recorded separately.

| Dataset | Inputs | Completed runs |
| --- | --- | ---: |
| Six-case AH/Bessel comparison | `N=8,16`; `L=16,24,32`; 2 methods × 3 repeats | 36 |
| Seven-method comparison, including PARI | `N,L ∈ {4,6,8,12,16}`; 7 methods × 3 repeats | 525 |
| Adjacent odd/even ranks | `N=8,16`; `L=3,4,7,8,15,16,23,24,31,32`; 4 methods × 3 repeats | 240 |
| Larger-field even ranks | `N=8,16,24,32,48,64`; `L=4,8,16,24,32,48,64`; 2 methods × 3 repeats | 252 |
| Even ranks at `N=36` | `L=4,8,16,24,32,48,64`; 2 methods × 3 repeats | 42 |
| Larger-field comparison | `N=24,32,48,64`; `L=4,6,8,12,16`; 2 methods × 3 repeats | 120 |
| Direct PARI point counting | `N=2,4,6,8,12`; `L=2,3`; `a=1,2,3`; 3 repeats | 90 |

Some datasets reuse measurements. The [data guide](results/README.md) identifies this overlap and defines the output and timing fields. The saved results include exact integer comparisons across methods; large-field agreement checks consistency where a full-table reference is unavailable.

## PARI

PARI supplies two baselines: full-table cyclic powering and single-parameter elliptic-curve point counting.

### Full-table cyclic powering

The [full-table implementation](core/pari_variable_rank_20260930/code/methods.py) has a [25-point result table](results/pari_25_cells.csv) covering 75 runs and 300 exact integer outputs.

`run_pari` constructs the supplied binary field and finds a primitive element. It then computes an integer cyclic-polynomial power and extracts the requested coefficients using discrete logarithms. The GP function is embedded in the cypari2 adapter, which accepts ranks `L ∈ {4,6,8,12,16}`.

The comparison used an 8 MiB initial PARI stack, a 72 GiB stack maximum and an 80 GiB address-space limit per process. Its wrapper set these values directly, bypassing `setup_pari()`. The separate `run_shared()` helper was not used in this comparison.

### Low-rank point counting

The [Direct PARI implementation](core/pari_direct_20260930/direct_worker.py) calls `ellinit` and `ellcard` for a single supplied parameter. At rank two, it returns the elliptic-curve point count minus `2^N + 1`. At rank three, it uses `Kl_3(a) = Kl_2(a)^2 - 2^N`.

The [saved results](results/pari_direct_20260930/) contain 90 Direct PARI runs on 30 cases and 180 runs comparing a Python AH reference with Sage/PARI. The inputs were `N=2,4,6,8,12`, `L=2,3` and three nonzero parameter labels.

Computation time covers field, parameter and curve construction, point counting and integer recovery. Imports, PARI initialization and total process time have separate fields. See the [implementation guide](core/pari_direct_20260930/README.md) for input limits and commands.

## License

Project code is distributed under the [MIT license](LICENSE). FLINT, python-flint, PARI and their dependencies are separate software with their own licenses.
