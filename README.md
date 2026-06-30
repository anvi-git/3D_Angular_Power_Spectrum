# Beyond-BLAST

This code stems from the BLAST code presented [here](https://github.com/sofiachiarenza/blast_paper) and whose results are described in the paper (https://arxiv.org/abs/2410.03632).

This repository tries to go beyond that (hence the name) and to compute the 3D angular power spectrum.

```mermaid
flowchart TD
    A["Open beyond_BLAST notebook"] --> B["Precompile / import Julia packages"]
    B --> C["Load cosmology inputs and data files"]
    C --> D["Set cosmology and pipeline parameters"]
    D --> E["Preprocess inputs: interpolation & grids"]
    E --> F["Construct Chebyshev grids and operators"]
    F --> G["Compute Chebyshev coefficients (FastChebInterp)"]
    G --> H["Run FFT / FastTransforms steps"]
    H --> I["Evaluate large-scale integrals and kernels"]
    I --> J["Assemble power spectrum / transfer functions"]
    J --> K["Compare with BLAST baseline and diagnostics"]
    K --> L["Plot figures and save intermediate outputs"]
    L --> M["Run validation checks and quick unit tests"]
    M --> N["Refine parameters and re-run segments"]
    N --> G
    L --> O["Save notebook, figures, and tables"]
    O --> P["End"]

