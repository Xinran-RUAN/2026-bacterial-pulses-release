# Bacterial traveling pulses

**An asymptotic preserving scheme for kinetic equations describing bacterial traveling pulses**  
Canan Akkoyunlu, Jose A. Carrillo, Gissell Estrada-Rodriguez, and Xinran Ruan.  
[Preprint, arXiv:2611.11764v1](https://arxiv.org/html/2611.11764v1)

## Code and data

- MATLAB code: [GitHub repository](https://github.com/Xinran-RUAN/2026-bacterial-pulses-release)
- Simulation data: [Zenodo, DOI: 10.5281/zenodo.24129973](https://doi.org/10.5281/zenodo.24129973)

## Setup

Use MATLAB R2023a (tested with 9.14.0.2206163). No additional toolbox is
required. Download the code from GitHub and `data.zip` from Zenodo.
Extract the archive so that `data/` is alongside `code/` and this README.
Place the MAT files in the `data/figureNN/` folders listed below, with the
shared `baseline.mat` in `data/figure01/`.

Before running any commands below, set the MATLAB Current Folder to
`2026-bacterial-pulses-release`, which contains `code/`, `data/`, and
this README.

## Redraw from saved data

From the release directory in MATLAB, run:

```matlab
addpath(fullfile(pwd, 'code'));
out = fullfile(tempdir, 'bacterial-pulse-figures');
redraw_paper('all', out);
```

Use the command in the last column to redraw one result. All file paths
in the table are relative to `data/`. Each case name shown without an
extension is a `.mat` file. Shared cases are stored in one place, so
Figure 5 and Table 1 use files from other figure folders.

| Result | Data files and plotted settings | Source or simulation entry | Saved-data plotting command |
| --- | --- | --- | --- |
| Fig. 1(a,b) | `figure01/sign_macro`, density at `t=30,50` | `run_paper_case` and a peak-aligned, mass-matched analytical profile from Eqs. (5.1)–(5.2) | `redraw_paper(1,out)` |
| Fig. 2(a) | Eq. (5.2) evaluated on `chiS,chiN=0:0.02:2` | `analytic_speed_surface` | `redraw_paper(2,out)` |
| Fig. 2(b) | `figure02/historical_centroids.mat`, `macro_x` on `chiS,chiN=0:0.1:2` | Supplied centroid data | `redraw_paper(2,out)` |
| Fig. 2(c) | Same file, `kinetic_x` on the same grid with `epsilon=0.01` | Supplied centroid data | `redraw_paper(2,out)` |
| Fig. 3 | `figure01/sign_macro` and `figure03/` cases `sign_kinetic_e0p1`, `sign_kinetic_e0p05`, `sign_kinetic_e0p01`. Density at `t=1,10,20,30,40,50` | `run_paper_case` | `redraw_paper(3,out)` |
| Fig. 4 | All four Fig. 3 cases plus `figure04/` cases `sign_kinetic_e0p02`, `sign_kinetic_e0p005`, `sign_kinetic_e0p002`, `sign_kinetic_e0p001`. Full-domain discrete density `L2` error at `T=50` against `sign_macro` | `run_paper_case` | `redraw_paper(4,out)` |
| Fig. 5 | Same four files as Fig. 3. Saved center of mass and analytical speed from Eq. (5.2) | `run_paper_case` | `redraw_paper(5,out)` |
| Fig. 6 | `figure06/exp_a0p001_kinetic_e0p01`. Density at `t=1,10,20,30,40,50` | `run_paper_case` | `redraw_paper(6,out)` |
| Fig. 7(a) | Fig. 6 file and `figure07/` cases `exp_a0p001_macro`, `exp_a0p001_kinetic_e0p1`, `exp_a0p001_kinetic_e0p05`. Final density for `a=0.001` | `run_paper_case` | `redraw_paper(7,out)` |
| Fig. 7(b) | Fig. 6 file, `figure07/` cases `exp_a0p0005_kinetic_e0p01`, `exp_a0p0001_kinetic_e0p01`, and `figure03/sign_kinetic_e0p01`. Final density for `epsilon=0.01` | `run_paper_case` | `redraw_paper(7,out)` |
| Table 1 | Same four files as Fig. 7(b). Least-squares centroid slopes on `[30,50]` give `0.690, 0.707, 0.720, 0.724` in table order | `run_paper_case` | `redraw_paper('table',out)` |

The redraw writes `figure01a`–`figure01b`, `figure02a`–`figure02c`, and
`figure03`–`figure07` as PNG and PDF under `out`. Table 1 is
`table_1_wave_speeds.csv` in the same directory. Figure 2(a) is computed
from Eq. (5.2), and Figures 2(b,c) are plotted from the supplied centroid
data. Figure 5 uses PCHIP interpolation on `0:0.1:50` and a local linear fit within
`+/-0.5` of each plotting time.

## Rerun from initial conditions

For a fine-grid result, pass each case ID from the table to
`run_paper_case`. Use a separate output directory so the supplied data
stay intact. For example, to recompute and draw Figure 1:

```matlab
addpath(fullfile(pwd, 'code'));
out = fullfile(tempdir, 'bacterial-pulse-rerun-figures');
runs = tempname;
assert(run_paper_case('sign_macro', runs));
redraw_paper(1, out, runs);
```

For all 14 fine-grid cases, use the catalog returned by `paper_cases`:

```matlab
addpath(fullfile(pwd, 'code'));
out = fullfile(tempdir, 'bacterial-pulse-rerun-figures');
runs = tempname;
cases = paper_cases;
for k = 1:numel(cases)
    assert(run_paper_case(cases(k).id, runs));
end
for n = [1, 3:7]
    redraw_paper(n, out, runs);
end
redraw_paper('table', out, runs);
```

New simulation data are saved in `runs`, and the regenerated figures
and table are saved in `out`. The third argument to `redraw_paper`
selects the newly computed data for Figures 1 and 3–7 and Table 1.

## Data and parameters

- A case name gives the response (`sign` or `exp`), the solver
  (`macro` or `kinetic`), the exponential parameter `a` if applicable,
  and the kinetic `epsilon` after `_e`. The letter `p` marks a decimal
  point. For example, `exp_a0p001_kinetic_e0p01` has `a=0.001` and
  `epsilon=0.01`. The macroscopic solver uses `epsilon=0` in the response.
- The shared baseline MAT stores parameter struct `p` and grid struct
  `g`. The fine-grid cases use `(0,150)`, `dx=0.025`, `dt=0.001`,
  `T=50`, and 101 velocity nodes on `[-1,1]` with `dv=0.02`.
  The initial fields are `rho=exp(-3*x)`, `S=0`, `N=1000`, and
  kinetic `g=0`. The sensitivities are `chiS=0.5` and `chiN=1.1`.
  The chemical parameters are `DS=2`, `DN=0`, `alpha=0.05`,
  and `beta=gamma=1`. Variables use the manuscript's model units.
- The fine-grid v7.3 MAT files store `out.grid.x` (6,000 cell centers),
  `out.snap` (501 snapshots on `0:0.1:50`), and `out.cm_t` and
  `out.cm_x` (centroid times and positions). Each snapshot has density
  `rho`, chemoattractant `S`, and nutrient `N`, each a 1 × 6,000 array.
  Kinetic microstates `g` have 101 velocity rows and 6,001 edge columns.
  The `a=0.0005` and `0.0001` cases use complete time-50 checkpoint
  files. `read_completed_case` reads both layouts.
- The Figure 2 v7 MAT stores `macro_x` (21 × 21 × 100) and `kinetic_x`
  (21 × 21 × 60). Dimensions are `chiS`, `chiN`, and time, with unit-time
  samples in `macro_t` and `kinetic_t`. Panels (b) and (c) average the
  final 10 and 6 unit-time centroid differences, respectively.
