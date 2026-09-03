# Beyond-BLAST

This code stems from the BLAST code presented [here](https://github.com/sofiachiarenza/blast_paper) and whose results are described in the paper (https://arxiv.org/abs/2410.03632).

This repository tries to go beyond that (hence the name) and to compute the 3D angular power spectrum.

# beyond_BLAST notebook guide

This README explains what `beyond_BLAST.ipynb` does, how it relates to the Julia code in `blast_code/src/`, and what the main numerical objects mean.

## Purpose

The notebook extends the BLAST workflow to build a quantity of the form `S_l(k1,k2)` rather than only standard projected angular spectra. In practice, it loads cosmological background data, constructs galaxy and shear kernels, expands them in Chebyshev polynomials, evaluates Bessel/Chebyshev mixed integrals, combines them with a matter power spectrum, and produces diagnostic plots.

## Main files

| File | Role |
|---|---|
| `beyond_BLAST.ipynb` | Main driver notebook. Sets parameters, runs the pipeline, saves arrays and plots. |
| `blast_code/src/Blast.jl` | Core module. Loads BLAST internals and precomputed `T_tilde` artifacts. |
| `blast_code/src/background.jl` | Background cosmology utilities such as `H(z)` and `chi(z)`. |
| `blast_code/src/galaxy_galaxy.jl` | Galaxy kernel ingredients and Chebyshev coefficients. |
| `blast_code/src/shear_shear.jl` | Shear kernel ingredients and Chebyshev coefficients. |
| `blast_code/src/my_funcs.jl` | Clenshaw-Curtis grids/weights and the custom `W_tilde_computed` routine. |
| `blast_code/data/` | Input arrays such as background tables, `dN/dz`, and matter power spectra. |

## Pipeline overview

- `bb.jl` is the entry point included by the notebook; its `include` statements load the 11 source modules shown above.
- `config.jl` reads the three packaged input files used to build the cosmology and redshift grids.
- The `W_tilde.npy` reuse path is hard-coded as an absolute path in the notebook and is only read when `reuse = true`.
- The run directory is timestamped, so its exact name is created at execution time. Plot files are produced by the plotting functions and may vary with the executed cells.
- `notebooks/bb_code/src/blast_tutorials.jl`, `Blast.jl`, and `shear_shear.jl` exist in the source directory but are not included by `bb.jl` for this notebook.


```mermaid
flowchart LR
    notebook["notebooks/main_nb.ipynb"] --> env["notebooks/bb_code/Project.toml"]
    notebook --> entry["notebooks/bb_code/src/bb.jl"]

    entry --> source["src modules:<br/>cosmo · background · projected_matter · chebcoefs · integrals · funcs · plots · gg · paths · config · plot_config"]
    config --> data["data/background/z.npy<br/>data/background/chi.npy<br/>data/dNdzs_fullwidth.npz"]

    notebook -. reuse .-> reuse["out/.../W_tilde.npy"]
    notebook --> generated["out/runs/run_<timestamp>/<br/>output_info.txt · quantities/ · plots/"]

    env --> packages["Julia packages:<br/>NPZ · HDF5 · Plots · QuadGK · FFTW · Cosmology · ..."]
    notebook -. cites .-> paper["arXiv:2410.03632"]

