# LaTeX Wiki: 3D Power Spectrum & Weak Lensing Theory

## Structure & Maintenance Guide

This wiki is organized as a **modular LaTeX project** with cross-linked concepts, equations, examples, and explanations. Each section is self-contained but interconnected via labels and hyperlinks.

### Folder Organization

```
sheet/
├── wiki_main.tex           # Master document (compile this)
├── index_concepts.tex      # Concept index with cross-links
├── README_WIKI.md          # This file
│
├── 00_frontmatter/         # Introduction, notation, guide
│   ├── intro.tex
│   ├── notation.tex
│   ├── reading_guide.tex
│   └── quick_ref.tex       # Quick reference formulas
│
├── 01_background/          # Cosmological background theory
│   ├── cosmological_model.tex
│   ├── background_evolution.tex
│   └── distance_measures.tex
│
├── 02_power_spectrum/      # 3D power spectrum foundations
│   ├── cartesian_fourier.tex
│   ├── spherical_decomposition.tex
│   ├── statistical_framework.tex
│   └── matter_perturbations.tex
│
├── 03_spherical_fb/        # Spherical Fourier-Bessel framework
│   ├── bessel_functions.tex
│   ├── sfb_basis.tex
│   ├── sfb_expansion.tex
│   └── orthonormality.tex
│
├── 04_weak_lensing/        # Weak lensing theory
│   ├── lensing_geometry.tex
│   ├── convergence_shear.tex
│   ├── spin_harmonics.tex
│   └── lensing_power_spectrum.tex
│
├── 05_3d_observables/      # 3D observables & applications
│   ├── angular_power_spectrum.tex
│   ├── redshift_binning.tex
│   ├── sfb_kernels.tex
│   └── cross_spectra.tex
│
├── 06_numerical/           # Numerical methods
│   ├── integration_schemes.tex
│   ├── quadrature_strategies.tex
│   ├── performance.tex
│   └── benchmarking.tex
│
├── _common/                # Shared resources
│   ├── preamble.tex        # Package definitions & settings
│   ├── commands.tex        # Custom LaTeX commands
│   ├── notation_defs.tex   # Notation definitions
│   └── bib.bib             # Bibliography (BibTeX)
│
└── old/                    # Archive of previous work
    └── [original files]
```

### How to Use This Wiki

1. **Compile**: `pdflatex wiki_main.tex` or use your LaTeX editor
2. **Navigation**: Use table of contents and clickable cross-references
3. **Search**: In PDF viewer, search for concept names (all are indexed)
4. **Adding content**: Create new `.tex` files in appropriate folders, include in `wiki_main.tex`

### Cross-Referencing Convention

Every concept has a **label** for internal linking:
```latex
\section{My Concept}\label{sec:my_concept}
\subsection{Details}\label{subsec:my_concept_details}
```

Reference via: `\nameref{sec:my_concept}` or `\autoref{sec:my_concept}`

### Maintenance Workflow

**Adding a new concept**:
1. Create file: `02_power_spectrum/my_new_concept.tex`
2. Write with structure: `\section{...}\label{sec:...}`, `\subsection{}`, etc.
3. Include in `wiki_main.tex`: `\input{02_power_spectrum/my_new_concept.tex}`
4. Add entry to `index_concepts.tex` (one line per concept)
5. Link from related concepts using `\autoref{}`

**Updating existing content**:
- Edit the file directly
- Recompile to check cross-links
- No need to modify master file unless reorganizing

**Version control**:
- Commit `.tex` files only (not `.pdf`, `.aux`, `.log`)
- Add to `.gitignore`: `*.pdf *.aux *.log *.out *.toc *.bbl *.blg`

### Key Features

✅ **Modular**: Each concept in separate file  
✅ **Cross-linked**: Labels & auto-references throughout  
✅ **Searchable**: PDF index and concept listing  
✅ **Equation-heavy**: Full mathematical exposition  
✅ **Maintained**: Clear structure for daily updates  
✅ **Exhaustive**: Examples, explanations, derivations  

---

**Start reading at**: [01_background/cosmological_model.tex](01_background/cosmological_model.tex) or jump to [index_concepts.tex](index_concepts.tex)
