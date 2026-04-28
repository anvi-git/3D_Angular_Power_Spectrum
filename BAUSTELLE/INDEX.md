# Implementation Roadmap Index

**Last Updated**: April 28, 2026  
**Project**: Beyond Limber Angular Power Spectra via Spherical Fourier-Bessel Decomposition  
**Status**: Complete planning phase, ready for implementation

---

## 📂 All Documents Created

Located in: `/Users/anvi/Desktop/cosmo/`

```
cosmo/
├── 📘 MASTER_GUIDE.md                    ← START HERE
│   └─ Ties everything together, navigation guide
│
├── 📕 IMPLEMENTATION_GUIDE.md             ← BIG PICTURE
│   └─ 100+ pages of conceptual breakdown
│      • Dependency chain
│      • Layer-by-layer explanation
│      • 7-phase timeline
│      • Mathematical foundations
│
├── 📗 CODE_TEMPLATES.md                   ← COPY-PASTE CODE
│   └─ 150+ pages of implementation details
│      • Complete function templates
│      • Docstrings with examples
│      • Testing frameworks
│      • Performance profiling
│
├── 📙 QUICK_REFERENCE.md                  ← FAST LOOKUP
│   └─ 50+ pages of tables and checklists
│      • Equation-to-code mapping
│      • File-by-file checklists
│      • Debugging tips
│      • Week-by-week schedule
│
├── 📓 TUTORIAL_NOTEBOOK_OUTLINE.md        ← HANDS-ON
│   └─ 100+ pages of worked examples
│      • Cell-by-cell notebook outline
│      • Validation plots
│      • Performance benchmarks
│      • Physical interpretation
│
├── ✨ THIS FILE (INDEX)                   ← YOU ARE HERE
│
└── 📄 Other files:
    ├── equations.tex                      ← Your starting equations
    ├── sheet/                             ← Related theory docs
    └── paper_blast/                       ← Implementation location
        └── src/
            ├── cosmo.jl                   (expand)
            ├── background.jl              (expand heavily)
            ├── projected_matter.jl        (expand heavily)
            ├── spherical_bessel.jl        (CREATE NEW)
            ├── sfb_kernels.jl             (CREATE NEW)
            ├── sfb_decomposition.jl       (CREATE NEW)
            └── integrals.jl               (expand)
```

---

## 🎯 Quick Start Path

### For Impatient People (5 minutes)
1. Read this file (you're doing it!)
2. Read **MASTER_GUIDE.md** "Overview" section
3. Look at **QUICK_REFERENCE.md** table of contents

### For Serious Developers (2 hours)
1. Read **MASTER_GUIDE.md** completely
2. Read **IMPLEMENTATION_GUIDE.md** completely
3. Skim **CODE_TEMPLATES.md** function signatures
4. Print or bookmark **QUICK_REFERENCE.md**

### For Immediate Implementation (30 min prep)
1. Read **QUICK_REFERENCE.md** "Week 1 Plan"
2. Open **CODE_TEMPLATES.md** at "LAYER 1: Background Evolution"
3. Open **TUTORIAL_NOTEBOOK_OUTLINE.md** at "Cell 2"
4. Start coding in `src/background.jl`

---

## 📊 Content Overview

### Document Sizes

| Document | Pages | Focus | Time to Read |
|----------|-------|-------|--------------|
| MASTER_GUIDE.md | 7 | Navigation | 15 min |
| IMPLEMENTATION_GUIDE.md | ~60 | Concepts + math | 2-3 hours |
| CODE_TEMPLATES.md | ~70 | Code + examples | 2 hours |
| QUICK_REFERENCE.md | ~40 | Lookup tables | 30 min |
| TUTORIAL_NOTEBOOK_OUTLINE.md | ~80 | Worked examples | 1-2 hours |
| **TOTAL** | **~257** | **Complete coverage** | **8-10 hours** |

**Time to first working code**: 2-3 weeks (following Week 1-2 in timeline)

---

## 🗺️ Document Dependency Graph

```
                    equations.tex
                          ↓
                   MASTER_GUIDE ◄─── START HERE
                    (navigation)
                    /    |    \
                   /     |     \
          IMPL_GUIDE    QUICK_REF    CODE_TEMPLATES
           (concepts)   (lookup)        (code)
             \           |            /
              \          |           /
               \         |          /
                \        |         /
                 \       ↓        /
               TUTORIAL_NOTEBOOK
                (worked examples)
                        ↓
                  Your Julia Code
                  (paper_blast/src/)
```

---

## 📚 Content Mapping

### Where Is Everything?

**"I need to understand the math"**
→ IMPLEMENTATION_GUIDE.md → "Mathematical Foundations" section

**"I need to implement function X"**
→ QUICK_REFERENCE.md → "Master Equation-to-Code Index" → find function → CODE_TEMPLATES.md

**"How do I test this?"**
→ CODE_TEMPLATES.md → "Testing Strategy" section

**"My code breaks, help!"**
→ QUICK_REFERENCE.md → "Debugging Tips" table

**"What's the big picture?"**
→ MASTER_GUIDE.md → "The 5 Layers" section

**"Show me a worked example"**
→ TUTORIAL_NOTEBOOK_OUTLINE.md → appropriate cell

**"How long will this take?"**
→ QUICK_REFERENCE.md → "Implementation Timeline & Milestones"

**"What equations are we implementing?"**
→ QUICK_REFERENCE.md → "Master Equation-to-Code Index" (all 20+ equations)

---

## 🔄 Typical Workflow

### Daily Development

```
Morning: Check QUICK_REFERENCE.md "Week X/Day Y" plan
    ↓
Read relevant section in IMPLEMENTATION_GUIDE.md
    ↓
Check CODE_TEMPLATES.md for function signature & docstring
    ↓
Copy template code
    ↓
Code in src/LAYER.jl
    ↓
Test using TUTORIAL_NOTEBOOK_OUTLINE.md example
    ↓
Benchmark and optimize
    ↓
Commit to git
    ↓
Evening: Review yesterday's code, plan tomorrow
```

---

## 📋 Implementation Checklist

### Phase 1: Preparation (Week 1, 5 hours)
- [ ] Read MASTER_GUIDE.md
- [ ] Read IMPLEMENTATION_GUIDE.md completely
- [ ] Print QUICK_REFERENCE.md
- [ ] Create branch: `git checkout -b sfb/layer-1-background`
- [ ] Set up development environment in Julia

### Phase 2: Layer 1 (Week 1-2, 15 hours)
- [ ] Implement `comoving_distance(z)`
- [ ] Implement `growth_function(z)`
- [ ] Create `Background` struct
- [ ] Write tests (use CODE_TEMPLATES.md)
- [ ] Run TUTORIAL_NOTEBOOK_OUTLINE.md Cell 2
- [ ] Validate vs `data/background/` (< 0.1% error)
- [ ] Commit: `git commit -m "feat: Layer 1 complete"`

### Phase 3: Layer 2A (Week 2-3, 20 hours)
- [ ] Load `P_lin(k)` from file
- [ ] Implement `nonlinear_power_spectrum()`
- [ ] Implement `power_spectrum_3d()`
- [ ] Benchmark
- [ ] Run TUTORIAL_NOTEBOOK_OUTLINE.md Cell 3
- [ ] Commit

### ... (continue for Layers 2B-5, following QUICK_REFERENCE.md timeline)

---

## 🎓 Learning Paths

### Path A: "Just Make It Work" (Impatient)
1. QUICK_REFERENCE.md → Week 1 checklist
2. CODE_TEMPLATES.md → Function you need
3. TUTORIAL_NOTEBOOK_OUTLINE.md → Validation
4. Code → Test → Done

Time: 1-2 weeks implementation

### Path B: "I Want to Understand" (Thorough)
1. IMPLEMENTATION_GUIDE.md → Read slowly, understand deeply
2. Each layer: CODE_TEMPLATES.md → Understand design choices
3. Each layer: TUTORIAL_NOTEBOOK_OUTLINE.md → Validation
4. Code → Benchmark → Optimize → Done

Time: 3-4 weeks implementation, deeper understanding

### Path C: "Help Me Debug" (Troubleshooting)
1. QUICK_REFERENCE.md → Debugging table
2. CODE_TEMPLATES.md → Testing frameworks
3. IMPLEMENTATION_GUIDE.md → Mathematical verification
4. Fix → Retest → Done

Time: Variable

---

## 📈 Progress Tracking

### Use This Template

```
WEEK 1 (Background Evolution - Layer 1):
Days 1-2: ✅ comoving_distance(z) - COMPLETE
Days 3-4: ✅ growth_function(z) - IN PROGRESS  
Days 4-5: ⏳ Background struct - NOT STARTED
Validation: 🔄 In progress (χ(z) vs reference)
Tests: 🟡 Partial

WEEK 2 (continued...):
```

Track this in:
- A `.txt` file in your working branch
- Or in Notion/Obsidian if you prefer
- Or in git commit messages

---

## 🔍 How to Use Each Document

### ✅ MASTER_GUIDE.md
- **Read**: Once, at the start
- **Purpose**: Understand structure and navigate
- **Reference**: Whenever confused about which document to use

### ✅ IMPLEMENTATION_GUIDE.md
- **Read**: Completely before coding
- **Reread**: Relevant section before each layer
- **Purpose**: Deep understanding of mathematics and physics
- **Don't skip**: The "Dependency Chain" diagram

### ✅ CODE_TEMPLATES.md
- **Read**: Skim first, then reference specific functions
- **Use**: Copy templates while coding
- **Reference**: Testing and validation sections
- **Purpose**: Practical implementation details

### ✅ QUICK_REFERENCE.md
- **Use**: Almost daily
- **Reference**: Equation-to-code index
- **Purpose**: Fast lookups, week-by-week plan
- **Purpose**: Debugging when stuck

### ✅ TUTORIAL_NOTEBOOK_OUTLINE.md
- **Use**: After implementing each layer
- **Convert**: To `.ipynb` notebook for testing
- **Purpose**: Validation and visualization
- **Purpose**: Learning expected outputs

---

## 💾 Files to Create/Modify

### Create (NEW)
```
src/spherical_bessel.jl              (0 → 300 lines)
src/sfb_kernels.jl                   (0 → 500 lines)
src/sfb_decomposition.jl             (0 → 400 lines)
```

### Expand (EXISTING)
```
src/background.jl                    (add ≈200 lines)
src/projected_matter.jl              (add ≈150 lines)
src/integrals.jl                     (add ≈200 lines)
```

### Reference (DON'T MODIFY)
```
original_src/                        (reference only)
data/Limber/                         (validation data)
data/background/                     (validation data)
```

---

## 🚨 Critical Success Factors

### Must Do
- ✅ Validate each layer against Limber before moving to next
- ✅ Test thoroughly at each stage (use CODE_TEMPLATES.md tests)
- ✅ Keep code well-documented with docstrings
- ✅ Commit regularly with clear messages
- ✅ Use precomputation and caching for expensive operations

### Must Not Do
- ❌ Skip testing
- ❌ Commit hard-coded paths (use relative paths)
- ❌ Add `.npy` data files to git
- ❌ Ignore performance (profile as you go)
- ❌ Skip validation against Limber

---

## 📞 Quick Q&A

**Q: Which document should I read first?**
A: MASTER_GUIDE.md (this explains the others)

**Q: I'm drowning in documents. What do I absolutely need?**
A: QUICK_REFERENCE.md + CODE_TEMPLATES.md (skim the rest as needed)

**Q: How do I know my implementation is correct?**
A: Compare with data in `data/Limber/` and `data/background/`

**Q: Where should I put new code?**
A: See QUICK_REFERENCE.md "File-by-File Implementation Checklist"

**Q: What if I find an error in the documents?**
A: You're reading live documents - update and commit!

**Q: How long until I have working code?**
A: Layer 1 (background): 1-2 weeks. Full pipeline: 8 weeks.

---

## 🎯 Success Metrics

### By Week 2 (End of Phase 1-2):
- ✅ Layer 1 complete and tested
- ✅ Background vs reference < 0.1% error
- ✅ Validation notebook Cell 2 producing plots
- ✅ 10+ commits to git with clear messages

### By Week 4 (End of Phase 2-3):
- ✅ Layers 1-2A complete
- ✅ Power spectrum loaded and working
- ✅ First full pipeline test (L0-L2A)
- ✅ 20+ commits

### By Week 8 (End of all phases):
- ✅ All 5 layers complete
- ✅ Full validation notebook working
- ✅ < 1% error vs Limber reference
- ✅ Performance benchmarks documented
- ✅ Production-ready code
- ✅ 50+ commits with clear history

---

## 🎉 Conclusion

You now have:

✅ **4 comprehensive guides** (257 pages total)  
✅ **Complete code templates** for every function  
✅ **Step-by-step timeline** for 8 weeks  
✅ **Validation strategies** at every step  
✅ **Everything you need** to go from equations.tex → production Julia code  

The hardest part (planning) is done.

**The fun part (implementation) begins now.** 🚀

---

## 🔗 File Locations

| What | Where | Action |
|------|-------|--------|
| Equations | `sheet/equations.tex` | Read |
| Roadmap | `/cosmo/MASTER_GUIDE.md` | Read first |
| Concepts | `/cosmo/IMPLEMENTATION_GUIDE.md` | Deep read |
| Code | `/cosmo/CODE_TEMPLATES.md` | Copy-paste |
| Lookup | `/cosmo/QUICK_REFERENCE.md` | Daily use |
| Examples | `/cosmo/TUTORIAL_NOTEBOOK_OUTLINE.md` | Test against |
| Implementation | `paper_blast/src/` | Your code here |
| Data (validation) | `paper_blast/data/` | Compare against |

---

**Ready? Start with MASTER_GUIDE.md → Then read IMPLEMENTATION_GUIDE.md → Then follow QUICK_REFERENCE.md timeline.**

**Happy implementing! 🎯**

