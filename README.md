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

```mermaid
flowchart TD
<<<<<<< HEAD
    A["Start notebook and activate Julia environment"] --> B["using Pkg Pkg.activate(blast_code) Pkg.resolve() Pkg.instantiate()"]
    A --> C["Include BLAST and custom source files"]
=======
    A["Start notebook"] --> B["Activate Julia environment"]
    B --> C["Include BLAST and custom source files"]
>>>>>>> refs/remotes/origin/main
    C --> D["Load background arrays z and chi"]
    D --> E["Build interpolators z(chi) and chi(z)"]
    E --> F["Set grids in chi, z, ell, k"]
    F --> G["Build shear prefactor"]
    F --> H["Build galaxy prefactor"]
    G --> I["Compute Chebyshev coefficients for shear"]
    H --> J["Compute Chebyshev coefficients for galaxy"]
    I --> K["Load or compute W_tilde"]
    J --> K
    K --> L["Contract W_tilde with kernel coefficients"]
    L --> M["Build W_final"]
    M --> N["Load matter power spectrum P(k,z)"]
    N --> O["Build k-space weights"]
    O --> P["Assemble S_l(k1,k2)"]
    P --> Q["Save arrays and diagnostic plots"]
<<<<<<< HEAD

=======
>>>>>>> refs/remotes/origin/main

