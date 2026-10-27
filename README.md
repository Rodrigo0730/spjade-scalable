# spJADE variants

Alternative algorithms and kernel choices for spJADE, the spatial JADE estimator for independent component analysis of heteroskedastic random fields. It builds on spFOBI and spJADE from Morales-Martínez, Nordhausen and Ruiz (2026), *Statistics & Probability Letters* 238, 110858 (doi:10.1016/j.spl.2026.110858). The code for that paper is at [Rodrigo0730/spBSS](https://github.com/Rodrigo0730/spBSS).

## Methods (`R/analysis.R`)
| Function | Idea |
|---|---|
| `spFOBI`, `spJADE` | Reference implementations; `spJADE` takes a kernel normalization scheme (`none`, `frobenius`, `max_abs`) |
| `nn_spjade` | Joint diagonalization with Adam, without forming the cumulant matrices |
| `riemann_nn_spjade` | Riemannian optimization on the orthogonal group |
| `implicit_nn_spjade`, `rff_nn_spjade`, `hybrid_nn_spjade` | Implicit cumulants, random Fourier features, and a hybrid of the two approaches |
| `iterative_spjade` | Learns the spatial kernel and the unmixing matrix jointly, blending geographic and source-space kernels |
| `spJADE_subsample` | Stochastic subsampling of locations |

`R/datagen.R` simulates the source fields (settings 1–4).

## Results so far (setting 3)
**Kernel normalization and f0** (`scripts/efficiency.R`, `efficiency/setting_3/`): compares the normalization schemes with and without the f0 (diagonal) kernel, for d = 5, 10, 20, 40.

**Computation time** (`figures/time_vs_n.pdf`, n = d² up to 1600):
- Classical `spJADE` is the fastest method at every size tested (0.02–0.03 s per fit).
- The optimization-based variants are 10–80 times slower. They are written in plain R.
- `riemann_nn_spjade` comes closest. Its time is nearly flat in n (0.31 to 0.38 s), and its accuracy is close to `spJADE`'s (median MD 0.131 vs 0.111 at n = 1600).
- `nn_spjade` matches `spJADE`'s accuracy.
- `implicit` and `hybrid` lose accuracy (median MD around 0.7).

**Subsampling** (`figures/subsample_tradeoff.pdf`): the time saved grows with n. At n = 6400, using 30% of the locations cuts the time from 0.153 to 0.035 s, while the median MD rises from 0.074 to 0.139. At n ≤ 1600 the fit is already fast, so subsampling saves little.

**Open questions:**
- At what n, if any, does the O(n) scaling of the Riemannian variant overtake `spJADE`? That needs n well beyond 1600, or compiled versions of the optimization loops.
- How well does iterative kernel refinement work (`scripts/benchmark_iterative.R`)?

## Scripts (`scripts/`)
- `dataMPI.R`, `example_mpi.sh`: generate the simulated datasets in `data/setting_X/data_d.rds`
- `benchmark_debug.R`: checks that `nn_spjade` matches `spJADE` on a small case
- `benchmark_nn_spjade.R`, `benchmark_all.R`: accuracy and timing of the variants against `spJADE`
- `benchmark_d80.R`: timing at d = 80 (n = 6400)
- `benchmark_subsample.R`: accuracy and time against the subsampling fraction (d = 20, 40, 80)
- `benchmark_iterative.R`: iterative kernel refinement against geographic kernels
- `efficiency.R`: kernel normalization and f0 comparison
- `op_count_analysis.R`: operation counts of the methods
- `plot_nn_benchmark.R`, `plot_scaling.R`: figures

Benchmark outputs are written to `results/`, which is not tracked. Run the scripts from the repository root. R packages are listed in `requirements.txt`.

## License
GNU GPLv3
