# Table

|Name|Type|Size|Description|
|---|---|---|---|
|z_b|Array{Float64,1}|(N_z,)|redshift array from background|
|χ_b|Array{Float64,1}|(N_z,)|comoving distance array from background|
|z_of_χ|AkimaInterpolation|callable, not an array|interpolates redshift as a function of comoving distance|
|χ_of_z|AkimaInterpolation|callable, not an array|interpolates comoving distance as a function of redshift|
|n5k_bins|Dict{String,Any}|not an array|N5K redshift-distribution archive loaded from .npz|
|z_n5k|Array{Float64,1}|(N_n5k,)|clustering redshift grid|
|n_z_matrix|Array{Float64,2}|(N_z,N_bin)|clustering number-density matrix|
|norms|Vector{Float64}|(N_bin,)|normalization integrals for each tomographic bin|
|n_z_norm|Array{Float64,2}|(N_z,N_bin)|normalized clustering number-density matrix|
|z_array|LinRange|length(z_n5k)|redshift grid used to evaluate n(z)|
|chi_array|Array{Float64,1}|(N_n5k,)|comoving distances sampled from z_array|
|cosmo|Blast.FlatΛCDM|not an array|flat-ΛCDM cosmology object|
|n_chi|Int64|scalar|number of comoving-distance sample points|
|x_array|LinRange|length(n_chi)|comoving-distance grid|
|z_range|Array{Float64,1}|(n_chi,)|redshift values corresponding to x_array|
|b_0|Float64|scalar|bias normalization|
|b_z_array|Array{Float64,1}|(n_chi,)|bias evaluated on the redshift grid|
|heath_integral|Function|callable|unnormalized Heath growth integral|
|compute_growth_factor|Function|callable|normalized growth-factor helper|
|growth_array|Array{Float64,1}|(n_chi,)|growth factor values on z_range|
|n_bins|Int64|scalar|number of tomographic bins|
|nz|Array{Float64,2}|(n_chi,n_bins)|interpolated tomographic n(z) samples|
|nz_func|AkimaInterpolation|callable, not an array|per-bin n(z) interpolant used in the loop|
|nz_norm|Float64|scalar|integral of one bin's n(z) interpolant|
|prefac|Array{Float64,1}|(n_chi,)|χ^2 b(z) D(z) prefactor|
|kernel_composed|Array{Float64,2}|(n_chi,n_bins)|kernel prefactor multiplied by n(z)|
|ℓ|Any|not an array|multipole collection imported from Blast|
|x_min|Int64|scalar|minimum comoving distance|
|x_max|Int64|scalar|maximum comoving distance|
|χ|LinRange|length(n_chi)|comoving-distance grid for later transforms|
|R|Array{Float64,1}|(nR,)|positive Chebyshev ratio nodes|
|nR|Int64|scalar|number of retained ratio nodes|
|kmax|Float64|scalar|maximum wavenumber|
|kmin|Float64|scalar|minimum wavenumber|
|n_cheb|Int64|scalar|number of Chebyshev nodes|
|β|Int64|scalar|probe-dependent exponent|
|k_cheb|Array{Float64,1}|(n_cheb,)|Chebyshev nodes in log10(k)|
|N|Int64|scalar|number of integration points|
|k_array|Array{Float64,1}|(n_cheb,)|wavenumber grid on Chebyshev nodes|
|zmin|Float64|scalar|minimum redshift used for Chebyshev mapping|
|zmax|Float64|scalar|maximum redshift used for Chebyshev mapping|
|t_of_z|Function|callable|maps z to Chebyshev coordinate t|
|z_of_t|Function|callable|inverse map from t to z|
|j|UnitRange{Int64}|length(n_cheb)|index range for Chebyshev node construction|
|t_nodes|Array{Float64,1}|(n_cheb,)|Chebyshev-Lobatto nodes on [-1,1]|
|z_nodes|Array{Float64,1}|(n_cheb,)|Chebyshev nodes mapped to redshift|
|chi_interp|AkimaInterpolation|callable, not an array|interpolates χ(z) from sampled values|
|bias_interp|AkimaInterpolation|callable, not an array|interpolates bias as a function of redshift|
|growth_interp|AkimaInterpolation|callable, not an array|interpolates growth factor as a function of redshift|
|nz_interp|Vector{Any}|(n_bins,)|one n(z) interpolant per tomographic bin|
|K_vals|Array{Float64,2}|(n_cheb,n_bins)|kernel values at Chebyshev nodes|
|plan_fft_1d|Function|callable|builds a one-dimensional FFTW plan|
|fast_chebcoefs_1d|Function|callable|converts nodal values to Chebyshev coefficients|
|plan|FFTW.r2rFFTWPlan|not an array|cached FFTW plan reused for all bins|
|cheb_coeff_K|Array{Float64,2}|(n_bins,n_cheb)|Chebyshev coefficients for K_vals|
|inv_Hz_array|Array{Float64,1}|(n_chi,)|inverse Hubble factor sampled on z_range|
|inv_Hz_interp|AkimaInterpolation|callable, not an array|interpolates 1/H(z) as a function of redshift|
|K_tilde_vals|Array{Float64,2}|(n_cheb,n_bins)|kernel values including the H(z) factor|
|cheb_coeff_K_tilde|Array{Float64,2}|(n_bins,n_cheb)|Chebyshev coefficients for K_tilde_vals|
|nk|Int64|scalar|length of the final k grid|
|nk1|Int64|scalar|length of k1_array|
|nk2|Int64|scalar|length of k2_array|
|W_base_k1|Array{Float64,4}|(length(ℓ),nk1,nk,n_cheb)|base W modes for k1|
|W_base_k2|Array{Float64,4}|(length(ℓ),nk2,nk,n_cheb)|base W modes for k2|
|computation_time|Float64|scalar|elapsed time for W_base construction|
|W_tilde_k1|Array{Float64,4}|(length(ℓ),n_bins,nk1,nk)|Chebyshev-composed W modes for k1|
|W_tilde_k2|Array{Float64,4}|(length(ℓ),n_bins,nk2,nk)|Chebyshev-composed W modes for k2|
|pk_dict|Dict{String,Any}|not an array|power-spectrum archive loaded from .npz|
|Pklin|Array{Float64,2}|(length(z),length(k))|linear matter power spectrum|
|Pknonlin|Array{Float64,2}|(length(z),length(k))|nonlinear matter power spectrum|
|k|Array{Float64,1}|(length(k),)|wavenumber grid from the power-spectrum archive|
|z|Array{Float64,1}|(length(z),)|redshift grid from the power-spectrum archive|
|y|LinRange|length(k)|log10(k) grid for interpolation|
|x|LinRange|length(z)|redshift grid for interpolation|
|InterpPmm|Any|not an array|interpolant for the linear power spectrum|
|InterpPmm_nl|Any|not an array|interpolant for the nonlinear power spectrum|
|power_spectrum|Function|callable|linear power-spectrum wrapper|
|power_spectrum_nl|Function|callable|nonlinear power-spectrum wrapper|
|dtype|DataType|not an array|floating-point type used for benchmarks|
|benchmark_gg|Dict{String,Any}|not an array|galaxy-galaxy benchmark archive|
|benchmark_ll|Dict{String,Any}|not an array|lensing-lensing benchmark archive|
|benchmark_gl|Dict{String,Any}|not an array|galaxy-lensing benchmark archive|
|gg|Array{Float64,2}|(N_ell,N_gg)|galaxy-galaxy C_ℓ samples|
|ll|Array{Float64,2}|(N_ell,N_ll)|lensing-lensing C_ℓ samples|
|gl|Array{Float64,2}|(N_ell,N_gl)|galaxy-lensing C_ℓ samples|
|ell|Array{Float64,1}|(N_ell,)|multipole grid for benchmarks|
|gg_reshaped|Array{Float64,3}|(length(ell),10,10)|galaxy-galaxy benchmark reshaped to matrix form|
|ll_reshaped|Array{Float64,3}|(length(ell),5,5)|lensing-lensing benchmark reshaped to matrix form|
|gl_reshaped|Array{Float64,3}|(length(ell),10,5)|galaxy-lensing benchmark reshaped to matrix form|
|S_l|Array{Float64,4}|(length(ℓ),n_bins,nk,nk)|spectral object returned by compute_Sℓ|
|w0|Array{Float64,1}|(n_cheb,)|single-mode contribution for diagnostic plotting|
|w_full|Array{Float64,1}|(nk1,)|full W_tilde slice for diagnostic plotting|
|ik1|Int64|scalar|selected k1 index for a diagnostic slice|
|ik2|Int64|scalar|selected k2 index for a diagnostic slice|
|nl|Int64|scalar|number of multipoles|
|ncols|Int64|scalar|number of subplot columns|
|nrows|Int64|scalar|number of subplot rows|
|contract_Sℓ|Function|callable|contracts S_l with projection vectors|
|proj_k1|Array{Float64,1}|(size(S_l,3),)|projection vector for the first k axis|
|proj_k2|Array{Float64,1}|(size(S_l,4),)|projection vector for the second k axis|
|C_l_like|Array{Float64,2}|(length(ℓ),n_bins)|contracted C_ℓ-like result|
