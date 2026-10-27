# Scalable spatial JADE

Faster and more flexible variants of spJADE for spatial independent component analysis of heteroskedastic random fields. It builds on spFOBI and spJADE from Morales-Martínez, Nordhausen and Ruiz (2026), *Statistics & Probability Letters* 238, 110858 (doi:10.1016/j.spl.2026.110858). The code for that paper is at [Rodrigo0730/spBSS](https://github.com/Rodrigo0730/spBSS).

## Methods (`R/analysis.R`)
| Function | Idea |
|---|---|
| `spFOBI`, `spJADE` | Reference implementations |
| `nn_spjade` | Joint diagonalization with Adam, without forming the cumulant matrices; gives the same result as `spJADE` |
| `implicit_nn_spjade`, `riemann_nn_spjade`, `rff_nn_spjade`, `hybrid_nn_spjade` | Speed-oriented variants: implicit cumulants, Riemannian optimization on the orthogonal group, random Fourier features, and a hybrid |
| `iterative_spjade` | Learns the spatial kernel and the unmixing matrix jointly, blending geographic and source-space kernels |
| `spJADE_subsample` | Stochastic subsampling of locations |

`R/datagen.R` simulates the source fields (settings 1–4).

## Scripts (`scripts/`)
- `dataMPI.R`, `example_mpi.sh`: generate the simulated datasets in `data/setting_X/data_d.rds`
- `benchmark_debug.R`: checks that `nn_spjade` matches `spJADE` exactly on a small case
- `benchmark_nn_spjade.R`, `benchmark_all.R`: accuracy and timing of all variants against `spJADE` (setting 3, d = 5 to 60)
- `benchmark_d80.R`: timing at d = 80 (n = 6400)
- `benchmark_subsample.R`: accuracy against subsampling fraction (d = 20, 40, 80)
- `benchmark_iterative.R`: iterative kernel refinement against geographic kernels
- `efficiency.R`, `op_count_analysis.R`, `plot_nn_benchmark.R`: efficiency study, operation counts, figures

Run the scripts from the repository root. R packages are listed in `requirements.txt`.

## License
GNU GPLv3
