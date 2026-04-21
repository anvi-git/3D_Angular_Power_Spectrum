# Normalization and Computation Check

## Issues Found

### 1. **Spherical Harmonic Integral** ❌ CRITICAL ERROR
**Line**: Angular integral step

**Current statement**:
$$\int d \Omega_{A} d \Omega_B Y^*_{l,m}(\theta_A, \phi_A) Y_{l',m'}(\theta_B, \phi_B) = \delta_{ll'} \delta_{mm'}$$

**Problem**: This is **incorrect**. You're integrating over two *independent* angular coordinates. These should separate:

$$\int d \Omega_{A} d \Omega_B Y^*_{l,m}(\theta_A, \phi_A) Y_{l',m'}(\theta_B, \phi_B) = \left[\int d \Omega_A Y^*_{l,m}(\theta_A, \phi_A)\right] \left[\int d \Omega_B Y_{l',m'}(\theta_B, \phi_B)\right] = 0$$

**Correct approach**: You need to keep both $(\theta_A, \phi_A)$ and $(\theta_B, \phi_B)$ integrated with the plane wave exponential, not separate them prematurely.

---

### 2. **Plane Wave Expansion Step** ⚠️ AMBIGUOUS/INCOMPLETE

**Line**: Angular integral of exponential

**Current statement**:
$$\int d\Omega e^{i\vec{\textbf{k}''\cdot(\vec{\textbf{r}}_A-\vec{\textbf{r}}_B)}} = 4\pi i^{-l} j_l(k''r)Y_{lm}^{*}(\theta'',\phi'')$$

**Problems**:
- Notation is unclear: what is $(l,m)$ and what are $(\theta'', \phi'')$?
- What exactly is being integrated? Both $\Omega_A$ and $\Omega_B$, or just one?
- If integrating both angular coordinates, the result should involve a spherical harmonic addition theorem or be of the form $j_l(k''r_A)j_l(k''r_B)$, not a single term.

**Standard reference**: The plane wave expansion is:
$$e^{i\vec{k}·\vec{r}} = 4\pi \sum_{l,m} i^l j_l(kr) Y^*_{l,m}(\hat{k}) Y_{l,m}(\hat{r})$$

---

### 3. **Final Coefficient Simplification** ⚠️ CHECK DIMENSIONS

**Line**: Final equation for $C_l^{AB,obs}(k,k')$

**Coefficient calculation**:
$$\left(\frac{2}{\pi}\right)^2 \cdot \frac{4\pi}{(2\pi)^3} = \frac{4}{\pi^2} \cdot \frac{4\pi}{8\pi^3} = \frac{4}{\pi^2} \cdot \frac{1}{2\pi^2} = \frac{2}{\pi^4}$$

**Issue**: This seems overly complex. Typically after all normalizations cancel, you'd expect something like:
- $\frac{1}{2\pi^2}$ or 
- $\frac{2}{\pi^2}$ or similar

**Simplified form**: The factor should probably be:
$$\left(\frac{2}{\pi}\right)^2 \cdot \frac{4\pi}{(2\pi)^3} = \frac{4}{\pi^2} \cdot \frac{1}{2\pi^2} = \frac{1}{\pi^4}$$

Wait, let me recalculate more carefully...

---

## What Needs Fixing

1. **Restart the angular integration step** with both $(\theta_A, \phi_A)$ and $(\theta_B, \phi_B)$ participating in the plane wave integral, not separated beforehand.

2. **Clarify the plane wave expansion**: Show explicitly which angles are being integrated and how the Bessel functions combine.

3. **Verify the final coefficient** by working through the full derivation step-by-step (can compare against standard SFB literature if available).

4. **Check against limits**: For a delta function selection function $\phi_A(r) = \delta(r - r_0)$, the formula should reduce to something recognizable.

---

## Suggested Reference

Compare with standard SFB derivations in cosmology (e.g., FLAME papers, Szapudi et al., or similar). The $W_l^A(k,k'')$ and $W_l^B(k',k'')$ window functions are correct conceptually—the issue is in the path getting there.

