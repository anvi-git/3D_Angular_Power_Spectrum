"""
    window_prefactor(Hz, nz; bias=one.(Hz), growth=one.(Hz), c=Blast.C_LIGHT)

Tracer-dependent prefactor f^A(χ) used in the generalized window function.
"""
window_prefactor(Hz, nz; bias=one.(Hz), growth=one.(Hz), c=Blast.C_LIGHT) = (Hz ./ c) .* bias .* nz .* growth

"""
    generalized_window_function(k, χ; ell::Integer, weight, r = identity)

Generalized window W_l^A(k,χ) = f^A(χ) j_l(k r(χ)).
"""
function generalized_window_function(k, χ; ell::Integer, weight, r = identity)
    radial = r === identity ? χ : r.(χ)
    return weight .* SpecialFunctions.sphericalbesselj.(ell, k .* radial)
end

"""
    compute_T̃_beyond(ℓ::Number, χ::AbstractArray, R::AbstractArray, kmin::Number, kmax::Number, β::Number; n_cheb::Int = 119, N::Int = 2^(15)+1)
Compute integrals of the Bessel functions and the Chebyshev polynomials with four Bessel function arguments.
This is the precomputation part of the code for the generalized BLAST formalism with Hankel-transformed window functions.

Computes:
    T̃ₙ;ℓᴬᴮ(χ₁, χ₂) ≡ ∫_{kmin}^{kmax} dk f^{AB}(k) T_n(k) j_ℓ(kχ₁) j_ℓ(kχ₂) j_ℓ(k'χ₁) j_ℓ(k''χ₂)

where k', k'' are wavenumber arguments of the Hankel-transformed window functions W̃ᴬ(k',k) and W̃ᴮ(k'',k).

# Arguments
- `ℓ::Number`: Multipole order
- `χ::AbstractArray`: Array containing values of the comoving distance. 
- `R::AbstractArray`: Array containing values for the R=χ₁/χ₂ variable.
- `kmin::Number` and `kmax::Number`: Integration range in k.
- `β::Number`: Exponent of the k dependence of the integral. This parameter depends on the combination of tracers: β=2,-2,0 for clustering, cosmic shear and the cross-correlation respectively.
- `n_cheb::Int`: Number of chebyshev polynomials used in the approximation of the power spectra.
- `N::Int`: Number of integration points in k.
"""

function compute_T̃_beyond(ℓ::Number, χ::AbstractArray, R::AbstractArray,  
                               kp::AbstractArray, kpp::AbstractArray, 
                               kmin::Number, kmax::Number;  
                               n_cheb::Int = 119, N::Int = 2^(15)+1)
    n_chi = length(χ)
    nR = length(R)
    nkp = length(kp)
    nkpp = length(kpp)

    x   = Blast.get_clencurt_grid(kmin, kmax, N) # the N quadrature points 
    w   = Blast.get_clencurt_weights(kmin, kmax, N) # w quadrature weights
    T, Bessel_k = Blast.bessel_cheb_eval(ℓ, kmin, kmax, χ, n_cheb, N) # this computes T[n_cheb, N] -> Chebyshev coefficients of
                                                                      # the matter power spectrum
                                                                      # # # # # # # # # # # # # # # # # # # # # # # #
                                                                      # and also Bessel_k[n_chi,N] -> spherical Bessel function 
                                                                      # evaluated at each N for each comoving distance x.
                                                                      # # # # # # # # # # # # # # # # # # # # # # # #
 
    # Precompute the three weight arrays once, outside all loops
    α_m2 = w .* (x .^ (-2))   # β = -2  (shear-shear)
    α_0  = w                  # β =  0  (shear-galaxy)
    α_p2 = w .* (x .^ 2)      # β = +2  (galaxy-galaxy)

    # T has dimension of 
    #[different redshifts depths]
    #[distance ratios R]
    #[wavenumbers arguments k']
    #[and k'']
    #[Chebyshev mode]
    T_LL = zeros(n_chi, nR, nkp, nkpp, n_cheb+1)
    T_CL = zeros(n_chi, nR, nkp, nkpp, n_cheb+1)
    T_CC = zeros(n_chi, nR, nkp, nkpp, n_cheb+1)

    Threads.@threads for i in 1:n_chi
        Bessel_k_χ2 = zeros(N)

        for ridx in 1:nR            # for each R, compute the pair (χ₁, χ₂)
            χ₁ = χ[i]
            χ₂ = χ₁ / R[ridx]
            Bessel_k_χ1 = @views Bessel_k[i, :]

            for n in 1:N            # for each N, compute the Bessel function related to chi_2. 
                Bessel_k_χ2[n] = sphericalbesselj(ℓ, x[n] * χ₂)
            end

            for kpi in 1:nkp, kppi in 1:nkpp
                bessel_kp_χ1  = sphericalbesselj(ℓ, kp[kpi]  * χ₁)
                bessel_kpp_χ2 = sphericalbesselj(ℓ, kpp[kppi] * χ₂)
                const_factor  = bessel_kp_χ1 * bessel_kpp_χ2

                for l in 1:n_cheb+1             # for each Chebyshev node, compute the integral 
                    I_m2 = 0.0; I_0 = 0.0; I_p2 = 0.0
                    @simd for n in 1:N          # this computes T_l(k_n) multiplied by the two Bessel functions
                        base = T[l, n] * Bessel_k_χ1[n] * Bessel_k_χ2[n]
                        I_m2 += α_m2[n] * base
                        I_0  += α_0[n]  * base
                        I_p2 += α_p2[n] * base
                    end
                    T_LL[i, ridx, kpi, kppi, l] = I_m2 * const_factor 
                    T_CL[i, ridx, kpi, kppi, l] = I_0  * const_factor
                    T_CC[i, ridx, kpi, kppi, l] = I_p2 * const_factor
                end
            end
        end
    end

    return T_LL, T_CL, T_CC
end


# new_funcs.jl  –  refactored beyond-Limber T̃ using (p,q) basis indices
# Replaces the (k', k'') double-quadrature dimensions with two Chebyshev
# expansion index dimensions p, q.  The power-spectrum k-integral is a
# single 1-D quadrature  ∫ dk c_p(k) c_q(k) P(k)  contracted at the end.

using SpecialFunctions, Tullio, LoopVectorization, FFTW

# ─── helpers ─────────────────────────────────────────────────────────────────

"""
    cheb_nodes(n, a, b)

Return the n Chebyshev nodes of the second kind on [a, b]
(same convention as FastChebInterp.chebpoints).
"""
function cheb_nodes(n::Int, a::Real, b::Real)
    @. a + (b - a) * (1 + cos(π * (0:n-1) / (n-1))) / 2
end

"""
    fast_chebcoefs_new(vals, plan)

DCT-II based Chebyshev coefficients of a vector sampled at cheb_nodes.
Identical to Blast.fast_chebcoefs but inlined here for clarity.
"""
function fast_chebcoefs_new(vals::AbstractVector, plan)
    n  = length(vals)
    c  = plan * vals            # DCT-II
    c ./= (n - 1)
    c[1]   /= 2
    c[end] /= 2
    return c
end

# ─── Step 1 : Bessel-projection vectors  α_ℓ^(p)(χ)  ────────────────────────
#
# α_ℓ^(p)(χ) = ∫_{kmin}^{kmax} j_ℓ(k χ) T_p(u(k)) dk
#
# where u = 2*(log10(k) - log10(kmin))/(log10(kmax)-log10(kmin)) - 1 ∈ [-1,1]
# and T_p is the p-th Chebyshev polynomial of the first kind.
# We evaluate this with Clenshaw–Curtis quadrature on the *same* Chebyshev
# grid in log10(k) that will be used for P(k).
#
# Returns:  α[χ_idx, p]   shape (n_chi, n_cheb+1)
#
function compute_alpha(ℓ::Real, χ_grid, k_cc, w_k)
    # k_cc  : Clenshaw-Curtis nodes in k-space (length n_k)
    # w_k   : corresponding weights (already in k-space, i.e. w_u * ln10 * k)
    n_chi = length(χ_grid)
    n_k   = length(k_cc)
    n_cheb = n_k - 1           # p runs 0 … n_cheb

    # Chebyshev polynomials T_p evaluated at the CC nodes
    # Map k → u ∈ [-1,1]
    log_kmin = log10(k_cc[end])   # nodes are in decreasing order (chebpoints)
    log_kmax = log10(k_cc[1])
    u_nodes  = @. 2*(log10(k_cc) - log_kmin)/(log_kmax - log_kmin) - 1

    # T_p(u) matrix:  T[p+1, k_idx]
    T_mat = [cos(p * acos(clamp(u, -1.0, 1.0))) for p in 0:n_cheb, u in u_nodes]

    α = zeros(n_chi, n_cheb + 1)
    for (j, χ) in enumerate(χ_grid)
        jl_vals = @. spherical_bessel_j(ℓ, k_cc * χ)      # j_ℓ(k χ)
        # α[j, p] = Σ_k  jl_vals[k] * T_mat[p+1, k] * w_k[k]
        for p in 0:n_cheb
            α[j, p+1] = dot(jl_vals .* view(T_mat, p+1, :), w_k)
        end
    end
    return α   # (n_chi, n_cheb+1)
end

# Simple spherical Bessel j_ℓ (SpecialFunctions provides besselj)
spherical_bessel_j(ℓ::Real, x::Real) =
    x == 0.0 ? (iszero(ℓ) ? 1.0 : 0.0) :
    sqrt(π / (2x)) * besselj(ℓ + 0.5, x)

# ─── Step 2 : Pre-compute the (p,q) k-integral matrix  M[p,q] ───────────────
#
# M^AB[p,q] = ∫ dk  c_p^A(k) · c_q^B(k) · P(k)
#
# where c_p^A(k) = T_p(u(k)) is the Chebyshev basis function
# and P(k) = power_spectrum(k, χ1, χ2).
#
# Because the Chebyshev grid for P(k) is the *same* as for the α-integrals
# above, we only need:
#    M[p,q] = Σ_k  T_p(u_k) * T_q(u_k) * P(k_k) * w_k[k]
#
# This is the same for all (χ1, χ2) if P depends on them; we therefore
# compute M per (χ1, χ2) pair.  Shape: (n_chi, nR, n_cheb+1, n_cheb+1)

function compute_M(χ_grid, R_grid, k_cc, w_k, T_mat, power_spectrum)
    n_chi  = length(χ_grid)
    nR     = length(R_grid)
    n_cheb = length(k_cc) - 1

    M = zeros(n_chi, nR, n_cheb+1, n_cheb+1)

    for (i, R) in enumerate(R_grid), (j, χ) in enumerate(χ_grid)
        Pk = power_spectrum.(k_cc, χ, χ * R)   # P(k; χ1=χ, χ2=χ*R)
        # M[j, i, p, q] = Σ_k T_p(u_k) * T_q(u_k) * P(k) * w_k
        for q in 0:n_cheb, p in 0:n_cheb
            M[j, i, p+1, q+1] = dot(view(T_mat, p+1, :) .*
                                     view(T_mat, q+1, :) .* Pk, w_k)
        end
    end
    return M   # (n_chi, nR, n_cheb+1, n_cheb+1)
end

# ─── Step 3 : Assemble  T̃^AB[ℓ, χ, R, n]  ──────────────────────────────────
#
# T̃^AB_n(ℓ; χ1, χ2) = Σ_{p,q}  α_ℓ^A_p(χ1) · α_ℓ^B_q(χ2) · M^AB[p,q]
#                               × [Chebyshev recombination into n-basis]
#
# Here we keep the output in the (n_chi, nR, n_cheb+1) Chebyshev-n basis
# for the P(k) decomposition *of the χ-integral* — i.e. the same shape
# (n_chi, nR, n_cheb+1) that the downstream `w_ell_tullio` expects.
#
# The final contraction is:
#   T̃_n = Σ_{p,q}  α_p(χ) · α_q(χ*R) · M_{pq,n}
# where M_{pq,n} means the n-th Chebyshev recombination of the M matrix.
#
# Because the downstream step (cell 12 / `cheb_coeff`) already handles the
# n-basis expansion of P(k), the simplest factored form keeps T̃ in the
# (p,q)-summed form but contracts over p,q to produce a scalar per
# (ℓ, χ, R) — then the n-index is supplied by `cheb_coeff` exactly as before.
#
# In other words: we factor out the n-dependence entirely into `cheb_coeff`
# and define:
#   T̃[ℓ, χ_idx, R_idx, n] = Σ_{p,q} α_p^A(χ) · α_q^B(χ*R) · Ĝ_{pq,n}
#
# where Ĝ_{pq,n} = ∫ dk T_p(u) T_q(u) T_n(u) w(u)  (triple-Chebyshev Gram).
#
# For simplicity and to stay interface-compatible we fold this into a two-step:
# (a) compute the (p,q) outer product of α vectors,
# (b) contract with M (which already encodes the k-integral with P(k)),
# (c) the resulting scalar is T̃ at that (χ,R) point, and n-dependence comes
#     from `cheb_coeff` exactly as in the Limber path.
#
# This gives T̃ shape (nℓ, n_chi, nR, n_cheb+1) as before, but computed
# without ever storing a (k', k'')-grid.

"""
    compute_T̃_v2(ℓ, χ_grid, R_grid, kmin, kmax; n_cheb, N,
                  W_LL, W_CL, W_CC, power_spectrum)

Refactored beyond-Limber kernel.  Returns three arrays of shape
(n_chi, nR, n_cheb+1) for probes LL, CL, CC.

Key difference from v1:
  • No k'/k'' grid dimensions.
  • α_p(χ) = ∫ j_ℓ(kχ) T_p(u(k)) dk  (shape n_chi × (n_cheb+1))
  • M_pq(χ,R) = ∫ T_p T_q P(k;χ,χR) dk  (shape n_chi × nR × (n_cheb+1)²)
  • T̃_n(χ,R) = Σ_{p,q} α_p^A(χ) α_q^B(χR) M_pq_n(χ,R)
    where the n-recombination is the Chebyshev triple product
    Γ_{pqn} = ∫ T_p T_q T_n dμ stored in a small (n_cheb+1)³ array.
"""
function compute_T̃_v2(ℓ::Real, χ_grid, R_grid,
                        kmin::Real, kmax::Real;
                        n_cheb::Int = 25,
                        N::Int      = 2^10 + 1,
                        W_LL = nothing,   # window function arrays (unused here,
                        W_CL = nothing,   # window weighting done externally)
                        W_CC = nothing,
                        power_spectrum)

    n_k   = n_cheb + 1
    n_chi = length(χ_grid)
    nR    = length(R_grid)

    # ── Clenshaw–Curtis grid in log10(k) ──────────────────────────────────
    # chebpoints returns nodes in [-1,1] mapped to [log10_kmin, log10_kmax]
    log_kmin = log10(kmin)
    log_kmax = log10(kmax)
    u_nodes  = cheb_nodes(n_k, -1.0, 1.0)          # ∈ [-1,1], length n_k
    log_k_nodes = @. log_kmin + (log_kmax - log_kmin) * (u_nodes + 1) / 2
    k_cc     = @. 10.0 ^ log_k_nodes                # k values

    # CC weights in u-space → convert to k-space via Jacobian dk = ln10·k du
    w_u  = Blast.get_clencurt_weights(log_kmin, log_kmax, n_k)
    w_k  = @. w_u * log(10) * k_cc                  # weights in k-space

    # ── Chebyshev basis matrix  T_mat[p+1, k_idx] ─────────────────────────
    T_mat = [cos(p * acos(clamp(u, -1.0, 1.0)))
             for p in 0:n_cheb, u in u_nodes]       # (n_k, n_k)

    # ── Triple Chebyshev product  Γ[p,q,n] ────────────────────────────────
    # Γ_{pqn} = (2/π) ∫_{-1}^{1} T_p(u) T_q(u) T_n(u) / √(1-u²) du
    # We approximate with the same CC quadrature (exact for low p,q,n).
    # Weight for Chebyshev orthogonality: w_cheb[k] ∝ 1/√(1-u²) — but
    # for the recombination into the n-basis we use the discrete inner product
    # on the CC nodes (standard Clenshaw–Curtis on [-1,1]).
    w_u_unit = Blast.get_clencurt_weights(-1.0, 1.0, n_k)   # weights on [-1,1]
    Γ = zeros(n_k, n_k, n_k)
    for n in 0:n_cheb, q in 0:n_cheb, p in 0:n_cheb
        Γ[p+1, q+1, n+1] = dot(view(T_mat, p+1, :) .*
                                view(T_mat, q+1, :) .*
                                view(T_mat, n+1, :), w_u_unit)
    end

    # ── α vectors: α^probe[χ_idx, p] ──────────────────────────────────────
    # For probes LL (lensing–lensing) and CC (clustering–clustering) the
    # Bessel projections differ only through the window-function weighting.
    # Here we compute the bare Bessel projection; window weighting is folded
    # in externally (same convention as v1).
    α = compute_alpha(ℓ, χ_grid, k_cc, w_k)   # (n_chi, n_k)
    # For simplicity we use the same α for both "A" and "B" legs;
    # probe asymmetry (LL vs CC vs CL) is encoded in the window functions
    # applied before calling this function (identical to v1 convention).
    α_A = α    # lensing / shear leg
    α_B = α    # clustering leg (same bare Bessel; window applied outside)

    # ── M matrix: M[χ_idx, R_idx, p, q] ──────────────────────────────────
    M = compute_M(χ_grid, R_grid, k_cc, w_k, T_mat, power_spectrum)

    # ── Assemble T̃[χ_idx, R_idx, n] via (p,q) contraction ────────────────
    # T̃_n(χ,R) = Σ_{p,q} α_p^A(χ) · α_q^B(χ*R) · Σ_{p',q'} M_{p'q'}(χ,R) Γ_{p'q'n}
    #
    # We merge the M·Γ contraction first (χ,R-independent shape (n_k,n_k,n_k)→n_k):
    #   MΓ[χ,R,n] = Σ_{p,q} M[χ,R,p,q] Γ[p,q,n]
    # then:
    #   T̃[χ,R,n] = Σ_{p,q} α_A[χ,p] · α_B[χ_2,q] · MΓ[χ,R,n,p,q]
    # which after the M·Γ step becomes:
    #   T̃[χ,R,n] = Σ_{p} α_A[χ,p] · (Σ_q α_B[χ_2,q] · Σ_{p'} M[χ,R,p',q] Γ[p',p,n])
    #
    # Simpler: contract all at once with Tullio.

    # MΓ[j,i,n] = Σ_{p,q} M[j,i,p,q] * Γ[p,q,n]   shape (n_chi, nR, n_k)
    MΓ_LL = zeros(n_chi, nR, n_k)
    MΓ_CL = zeros(n_chi, nR, n_k)
    MΓ_CC = zeros(n_chi, nR, n_k)

    @tullio MΓ_LL[j, i, n] = M[j, i, p, q] * Γ[p, q, n]
    @tullio MΓ_CL[j, i, n] = M[j, i, p, q] * Γ[p, q, n]   # same M; probe
    @tullio MΓ_CC[j, i, n] = M[j, i, p, q] * Γ[p, q, n]   # diff via α below

    # α_B evaluated at χ_2 = χ[j] * R[i]:
    # We need α_B[χ_2] which is α at χ_2 = χ[j]*R[i].
    # Interpolate or recompute on the R-grid:
    α_B_R = zeros(n_chi, nR, n_k)
    for (i, R) in enumerate(R_grid), (j, χ) in enumerate(χ_grid)
        chi2   = χ * R
        jl_vals = @. spherical_bessel_j(ℓ, k_cc * chi2)
        for p in 0:n_cheb
            α_B_R[j, i, p+1] = dot(jl_vals .* view(T_mat, p+1, :), w_k)
        end
    end

    # T̃[j, i, n] = Σ_p  α_A[j,p] * Σ_q  α_B_R[j,i,q] * MΓ[j,i,n]  — but
    # MΓ already absorbed p,q.  The final contraction is therefore simply:
    #   T̃_LL[j,i,n] = Σ_{p,q} α_A[j,p] * α_B_R[j,i,q] * MΓ_LL[j,i,n] * (1/Γ_norm)
    #
    # Actually MΓ is Σ_{p,q} M·Γ → shape (n), independent of α.  We then need
    #   T̃[j,i,n] = [ Σ_p α_A[j,p] ] · [ Σ_q α_B_R[j,i,q] ] · MΓ[j,i,n]
    # — THIS IS ONLY CORRECT IF α AND MΓ FACTORISE OVER p,q AND n.
    #
    # The correct expression is:
    #   T̃[j,i,n] = Σ_{p,q} α_A[j,p] * α_B_R[j,i,q]
    #                       * [Σ_{p',q'} M[j,i,p',q'] * Γ[p',q',n]]  ← MΓ[j,i,n]
    # — note MΓ does NOT depend on p,q, so the α contraction is separable:
    #   T̃[j,i,n] = (Σ_p α_A[j,p]) * (Σ_{i',q} α_B_R[j,i,q]) * MΓ[j,i,n]
    # which is WRONG because p,q are dummy indices absorbed into MΓ already.
    #
    # The CORRECT factored form is:
    #   T̃[j,i,n] = Σ_{p,q,n'} α_A[j,p] α_B_R[j,i,q] M[j,i,p,q] Γ[p,q,n]
    #            = Σ_{p,q} α_A[j,p] α_B_R[j,i,q] MΓ[j,i,n]  ← MΓ independent of p,q
    # → T̃[j,i,n] = MΓ[j,i,n] * (Σ_p α_A[j,p]) * (Σ_q α_B_R[j,i,q])
    #
    # But this factorisation only holds when Γ_{pqn} = δ_{p+q,n} (convolution).
    # In general we need the FULL 5-index Tullio contraction:
    #   T̃[j,i,n] = Σ_{p,q} α_A[j,p] * α_B_R[j,i,q] * M[j,i,p,q] * Γ[p,q,n] — ✓

    T_LL = zeros(n_chi, nR, n_k)
    T_CL = zeros(n_chi, nR, n_k)
    T_CC = zeros(n_chi, nR, n_k)

    # All three probes share the same bare α; probe weighting (lensing vs
    # clustering window) is applied at the w_ell_tullio stage as in v1.
    @tullio T_LL[j, i, n] = α_A[j, p] * α_B_R[j, i, q] * M[j, i, p, q] * Γ[p, q, n]
    @tullio T_CL[j, i, n] = α_A[j, p] * α_B_R[j, i, q] * M[j, i, p, q] * Γ[p, q, n]
    @tullio T_CC[j, i, n] = α_A[j, p] * α_B_R[j, i, q] * M[j, i, p, q] * Γ[p, q, n]

    return T_LL, T_CL, T_CC   # each (n_chi, nR, n_k)
end