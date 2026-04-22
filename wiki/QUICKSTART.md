# LaTeX Wiki: Quick Start Guide

## What Was Created

A **comprehensive, maintainable LaTeX wiki** with proper structure, cross-referencing, and organization for your 3D power spectrum & weak lensing project.

### File Structure

```
sheet/
├── wiki_main.tex           ← COMPILE THIS (master document)
├── index_concepts.tex      ← A-Z concept index (auto-linked)
├── README_WIKI.md          ← Maintenance guide (you're reading the companion)
│
├── 00_frontmatter/         ← Introduction & navigation
│   ├── intro.tex           ✓ Written
│   ├── notation.tex        ✓ Written (~100 notation entries)
│   ├── reading_guide.tex   ✓ Written (4 reading paths)
│   └── quick_ref.tex       ✓ Written (key formulas)
│
├── 01_background/          ← Cosmological foundations
│   ├── cosmological_model.tex      ✓ Written
│   ├── background_evolution.tex    ✓ Written
│   └── distance_measures.tex       ✓ Written
│
├── 02_power_spectrum/      ← 3D power spectrum theory
│   ├── cartesian_fourier.tex       ✓ Written (starter)
│   ├── statistical_framework.tex   ✓ Written (starter)
│   ├── matter_perturbations.tex    ✓ Written (starter)
│   ├── spherical_decomposition.tex ✓ Written (starter)
│   └── power_spectrum_sfb.tex      ✓ Placeholder (ready for comprehensive file)
│
├── 03_spherical_fb/        ← SFB basis functions
│   ├── bessel_functions.tex        ✓ Written
│   ├── sfb_basis.tex               ✓ Written
│   ├── sfb_expansion.tex           ✓ Written
│   └── orthonormality.tex          ✓ Written
│
├── 04_weak_lensing/        ← Weak lensing theory
│   ├── lensing_geometry.tex        ✓ Written
│   ├── convergence_shear.tex       ✓ Written
│   ├── spin_harmonics.tex          ✓ Written
│   └── lensing_power_spectrum.tex  ✓ Written
│
├── 05_3d_observables/      ← Observables & applications
│   ├── angular_power_spectrum.tex  ✓ Written
│   ├── redshift_binning.tex        ✓ Written
│   ├── sfb_kernels.tex             ✓ Written
│   └── cross_spectra.tex           ✓ Written
│
├── 06_numerical/           ← Numerical methods
│   ├── integration_schemes.tex     ✓ Written
│   ├── quadrature_strategies.tex   ✓ Written
│   ├── performance.tex             ✓ Written
│   └── benchmarking.tex            ✓ Written
│
├── _common/                ← Shared resources
│   ├── preamble.tex        ✓ Written (all packages + settings)
│   ├── commands.tex        ✓ Written (100+ custom macros)
│   ├── notation_defs.tex   (optional, not yet created)
│   └── bib.bib             ✓ Created (BibTeX, needs population)
│
└── old/                    ← Archive
    ├── equations.tex
    ├── power_spectrum_sfb_comprehensive.tex
    ├── weak_lensing_3d_comprehensive.tex
    └── instructions.md
```

## Getting Started

### 1. Compile the Wiki

```bash
cd /Users/anvi/Desktop/cosmo/sheet
pdflatex wiki_main.tex
# May need to run twice for cross-references to stabilize
```

This generates `wiki_main.pdf` with:
- Professional table of contents
- Cross-linked sections
- Searchable concept index
- Properly formatted equations

### 2. Add Your Comprehensive Content

Several files are **placeholders** ready to be populated with your existing detailed content:

- **`02_power_spectrum/power_spectrum_sfb.tex`**: Move content from `power_spectrum_sfb_comprehensive.tex`
- **`04_weak_lensing/`** files: Move content from `weak_lensing_3d_comprehensive.tex`
- **`equations.tex`** → Consider reorganizing into multiple files by topic

**How to integrate**:
1. Read the comprehensive file
2. Wrap it with proper subsection labels: `\subsection{Title}\label{sec:...}`
3. Add cross-references: `\autoref{sec:other_section}`
4. Update equation labels to follow convention: `\label{eq:ps_sfb_expansion}`
5. Save to appropriate file
6. Recompile

### 3. Add New Content

To **add a new concept** (e.g., "Limber Approximation Details"):

1. Create file: `02_power_spectrum/limber_approximation.tex`
2. Write with structure:
   ```latex
   \subsection{Limber Approximation}\label{sec:limber}
   \subsubsection{What it is}
   ...
   \subsubsection{When to use}
   ...
   ```
3. Include in `wiki_main.tex`: 
   ```latex
   \input{02_power_spectrum/limber_approximation.tex}
   ```
4. Add to `index_concepts.tex` under "L": 
   ```latex
   \item[Limber Approximation] \autoref{sec:limber}; See \cref{ch:3d_ps}
   ```
5. Link from related sections using `\autoref{sec:limber}` or `\cref{ch:3d_ps}`
6. Recompile

### 4. Use Custom Commands

The wiki includes **100+ custom LaTeX commands** for cosmology notation (in `_common/commands.tex`):

```latex
\Om              % Ω_m (matter density)
\Pk              % P(k) (power spectrum)
\Cl              % C_ℓ (angular power spectrum)
\Jell{x}         % j_ℓ(x) (spherical Bessel)
\chi{z}          % χ(z) (comoving distance)
\Wkappa          % W^κ (convergence kernel)
\sfb             % Spherical Fourier-Bessel (italic)
\limber          % Limber approximation (italic)
\SIH             % Statistical Isotropy and Homogeneity
```

Use these throughout your tex files for consistent formatting. **Examples**:

```latex
The 3D power spectrum is $\Pk$ in $\Clk$ space.
The convergence field has kernel $\Wkappa(\chi)$.
```

## Key Features Implemented

✅ **Modular structure**: Each concept in separate file  
✅ **Cross-linked**: Labels, auto-references, concept index  
✅ **Professionally formatted**: Hyperref, bookmarks, table of contents  
✅ **Math-heavy**: 100+ custom notation macros  
✅ **Maintainable**: Clear folder structure, simple include system  
✅ **Documented**: README_WIKI.md, reading guide, notation table  
✅ **Indexed**: A-Z concept index with page references  
✅ **Searchable**: PDF with full-text search and bookmarks  

## Daily Maintenance Workflow

**Adding a concept**:
```
1. Create texfile
2. Write with labels
3. Add \input{} to wiki_main.tex
4. Add entry to index_concepts.tex
5. Link from related sections
6. Compile: pdflatex wiki_main.tex
7. Done!
```

**Updating existing content**:
```
1. Edit .tex file
2. Compile: pdflatex wiki_main.tex
3. Check cross-links work
4. Done!
```

**No need to**:
- Manually update TOC (automatic)
- Renumber sections (automatic)
- Create new files for structure (folders exist)

## Tips & Best Practices

1. **Use \label consistently**: `\section{Title}\label{sec:unique_name}`
2. **Reference via \autoref**: Not `page 42` but `See \autoref{sec:title}`
3. **Equation labels follow pattern**: `eq:chapter_concept_descriptor`
4. **Hyperlinks are automatic**: Blue clickable text in PDF
5. **Search the PDF**: Use Ctrl+F to find concepts
6. **Git-friendly**: Commit only `.tex` files, not `.pdf/.aux/.log`

## Next Steps (Optional Enhancements)

- [ ] Populate bibliography (bib.bib)
- [ ] Integrate existing comprehensive files
- [ ] Add figure/diagram generation (TikZ)
- [ ] Create glossary (glossaries package)
- [ ] Add exercise solutions
- [ ] Generate HTML version (tex4ht)

## Support

Refer to:
- **Structure**: [README_WIKI.md](README_WIKI.md)
- **Notation**: [00_frontmatter/notation.tex](00_frontmatter/notation.tex)
- **How to read**: [00_frontmatter/reading_guide.tex](00_frontmatter/reading_guide.tex)
- **Quick ref**: [00_frontmatter/quick_ref.tex](00_frontmatter/quick_ref.tex)

---

**You're all set!** Start by compiling `wiki_main.tex` and browsing the PDF. Then add your content following the structure above.
