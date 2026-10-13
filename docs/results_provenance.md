# Gallery result provenance

Fresh executions of selected archived final-revision Python blocks. Plots restyled; solvers and numerical parameters retained. These are individual example runs, not aggregate benchmark scores.

Rebuild with `python tools/build_showcase.py`; use `--case cylinder` (or another case ID) for one result. Dependencies: NumPy, SciPy, Matplotlib, pandas, scikit-learn, Pillow.

All source files are preserved. The builder extracts their fenced Python blocks and executes them in `.showcase_work/<case>/`. Python plotting uses a noninteractive backend. The cylinder adds read-only snapshot collection in the original time loop; the GIF shows those computed states at fixed color limits. No physical or numerical parameter was changed.

The gallery reports solver diagnostics and in-sample regression fits honestly. It does not certify discretization accuracy; execution success and benchmark accuracy are distinct.

## Cylinder flow

Source: `src/ipynb_ver/sample_results/_GPT5mini_/cylinder_for_Reproduction/sample1/solu_rev_2.dat`
SHA-256: `d6a03f4cdf9e5db922756bf1846a7fd5bbfeb8b70c32d458a6d303196daa0744`

A fresh execution of the archived final revision. The animation contains actual velocity snapshots from the time integration; the original grid and final time are retained.

- Solver execution: 9.62 s
- Grid: 220 × 120
- Final time: 20
- Reynolds number: 100
- Snapshots: 66
- Maximum final speed: 1.5221
- Max divergence, all cells: 4.309e+01
- Max divergence, outside cylinder: 1.977e-01

The unchanged archived solver has nonzero discrete divergence (see diagnostics). This is an illustrative reproduction, not a validated incompressible reference solution. Snapshot arrays are stored as float32 for display; the solver runs in its original precision.

Log: `assets/results/logs/cylinder.txt`

## Sod shock tube

Source: `src/ipynb_ver/sample_results/_3Programmers_for_Reproduction/PDEbench/SodShock/GPT_ans/sample1/solu_rev_2.dat`
SHA-256: `8a6ad3b10180b898ee265df73e7296fcd7c6858b6bfb8ba3c6109b09531c705a`

Density, velocity and pressure from one archived GPT final-revision solver, rerun on its original 400-cell grid.

- Solver execution: 2.04 s
- Grid: 400 cells
- Final time: 0.200
- Density range: 0.125–1.000



Log: `assets/results/logs/sod.txt`

## Lid-driven cavity

Source: `src/ipynb_ver/sample_results/_3Programmers_for_Reproduction/PDEbench/Lid_Driven/GPT_ans/sample1/solu_rev_2.dat`
SHA-256: `296da20009f58c047632eb77b72f025420fda6bd1f01eef0c7dcdafc50ff6a89`

Velocity magnitude and streamlines from the archived artificial-compressibility implementation. The run's own stopping criterion and diagnostics are preserved.

- Solver execution: 8.66 s
- Grid: 81 × 81
- Iterations: 22974
- Final max update: 9.983e-07
- Tolerance: 1e-06



Log: `assets/results/logs/cavity.txt`

## Unsteady Navier–Stokes

Source: `src/ipynb_ver/sample_results/_3Programmers_for_Reproduction/PDEbench/Unsteady_NS/GPT_ans/sample1/solu_rev_2.dat`
SHA-256: `6b168740b62c55858eac8fe4cc2de7d837c4833b94ff80b404e81eec95141ce0`

The two velocity components and pressure of the archived generated solver at its original final time. Displayed fields are illustrative outputs, not a claim of validation against the reference.

- Solver execution: 7.24 s
- Grid: 81 × 41
- Final time: 0.5
- Time step: 0.002
- Final divergence max: 1.100e-01

A successful execution is not an accuracy certification. Divergence diagnostics are reported from this run; no solver correction or reference-error claim is added.

Log: `assets/results/logs/unsteady_ns.txt`

## Viscous Burgers equation

Source: `src/ipynb_ver/sample_results/_3Programmers_for_Reproduction/PDEbench/Burgers_equ/GPT_ans/sample1/solu_rev_2.dat`
SHA-256: `78928b0c181e263cd736492c2b66ca68d8fbb2e2b6831a8797e7bbb8f93c51fe`

Space–time evolution and the final velocity profile from an archived MUSCL-type limited solver with implicit diffusion.

- Solver execution: 0.88 s
- Grid: 201 points
- Final time: 1.000
- Saved states: 230



Log: `assets/results/logs/burgers.txt`

## Poisson / Laplace equation

Source: `src/ipynb_ver/sample_results/_3Programmers_for_Reproduction/PDEbench/Poisson_equ/GPT_ans/solu_rev_2.dat`
SHA-256: `55c94d63f9e9c82618d4a9e1eebdff44b51f99a7d9184fc4e68062a543eef3f6`

The archived Poisson task is a zero-source Laplace solve on a square with four circular holes. Hole interiors are masked in the presentation.

- Solver execution: 1.06 s
- Grid: 201 × 201
- Field range: 0.0000–1.0000
- Linear residual L2: 5.255e-14



Log: `assets/results/logs/poisson.txt`

## Helmholtz-type equation

Source: `src/ipynb_ver/sample_results/_3Programmers_for_Reproduction/PDEbench/Hemholtz/GPT_ans/solu_rev_2.dat`
SHA-256: `d2be84a4a7943244eb661db4c4bc3be00ee847f4579c1d67d44539eff58ad174`

A screened elliptic solve from the archive (named Poisson–Boltzmann in the generated code), on a square with four differently sized holes.

- Solver execution: 1.12 s
- Grid: 201 × 201
- Field range: -0.7107–1.3410
- Linear residual L2: 3.411e-03



Log: `assets/results/logs/helmholtz.txt`

## Hilbert linear systems

Source: `src/ipynb_ver/sample_results/_3Programmers_for_Reproduction/Hilbert/GPT_ans/sample1/solu_rev_2.dat`
SHA-256: `5441d76b04166a1fc344622a4293b332632b5494158eba15da5d5c6edd0f4151`

Forward errors and condition numbers from one generated implementation. Failed decompositions remain visible as missing values rather than being silently repaired.

- Solver execution: 1.19 s
- Matrix sizes: 5–50
- Regularization: 1e-12
- Failed original Cholesky: 8



Log: `assets/results/logs/hilbert.txt`

## Keyhole scaling

Source: `src/ipynb_ver/sample_results/_3Programmers_for_Reproduction/Keyhole/GPT_ans/solu_rev_2.dat`
SHA-256: `cdbd80a03f745ed6227c167ef4f7ad0b7a951dffa56adef57e3f5f2aebaee074`

The archived rational-exponent search and log–log regression rerun on dataset_keyhole.csv. R² describes the fitted dataset; it is not a held-out score.

- Solver execution: 3.12 s
- Samples: 90
- In-sample R²: 0.9510



Log: `assets/results/logs/keyhole.txt`

## Rayleigh–Bénard scaling

Source: `src/ipynb_ver/sample_results/_3Programmers_for_Reproduction/Rayleigh/GPT_ans/solu_rev_2.dat`
SHA-256: `2836a30bb7238b6eea58b6f3390e1be5667223fdb896bcd6b4cb0a97d8ec8df2`

The archived dimensionless-group candidate search rerun on dataset_rayben.csv. R² describes the fitted dataset; it is not a held-out score.

- Solver execution: 2.94 s
- Samples: 182
- Candidates: 128
- In-sample R²: 0.9923



Log: `assets/results/logs/rayleigh.txt`
