# Beyond-BLAST

This code stems from the BLAST code presented [here](https://github.com/sofiachiarenza/blast_paper) and whose results are described in the paper (https://arxiv.org/abs/2410.03632).

This repository tries to go beyond that (hence the name) and to compute the 3D angular power spectrum.

```mermaid
flowchart TD
    A["Start notebook"] --> B["Activate Julia environment"]
    B --> C["Include BLAST and custom source files"]
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

