# Line-by-line explanation of `new_beyond_BLAST.ipynb`

## Reading notes

When a line defines a scalar, the explanation says so explicitly. When a line defines a vector, matrix, interpolation object, FFT plan, or plotting object, that is also stated explicitly. Complexity statements are rough asymptotic descriptions such as linear time, `O(N log N)`, or adaptive-quadrature cost when no fixed polynomial bound is the right description. [file:1]

```julia
#Background quantities
z_b = npzread("blast_code/data/background/z.npy") # array 
χ_b = npzread("blast_code/data/background/chi.npy") # array 
# using Akima interpolation
z_of_χ = DataInterpolations.AkimaInterpolation(z_b, χ_b); # z(χ)
χ_of_z = DataInterpolations.AkimaInterpolation(χ_b, z_b); # χ(z)
```
This cell reads external data files into memory. It also creates interpolation objects, which are callable approximations built from tabulated arrays. [file:1]

- `z_b = npzread("blast_code/data/background/z.npy") # array `: Reads a NumPy background redshift file. `z_b` is a one-dimensional array, likely a `Vector{Float64}`, with shape `(N_z,)`. Reading cost is linear in the number of stored entries.
- `χ_b = npzread("blast_code/data/background/chi.npy") # array `: Reads a NumPy background comoving-distance file. `χ_b` is another one-dimensional array, expected to align elementwise with `z_b`, so its shape is also `(N_z,)`.
- `z_of_χ = DataInterpolations.AkimaInterpolation(z_b, χ_b); # z(χ)`: Constructs an Akima interpolant intended to represent `z(χ)`. The result is a callable interpolation object. Precomputation cost is roughly linear in the tabulated sample count.
- `χ_of_z = DataInterpolations.AkimaInterpolation(χ_b, z_b); # χ(z)`: Constructs an Akima interpolant intended to represent `χ(z)`. The result is another callable interpolation object.


```julia
# Load N5K n(z), needed to compute lensing and clustering kernels
n5k_bins = npzread("blast_code/data/dNdzs_fullwidth.npz")
z_n5k = n5k_bins["z_cl"]
z_n5k = sort(z_n5k) # sort the redshift array in ascending order
n_z_matrix = n5k_bins["dNdz_cl"];
```

- `n5k_bins = npzread("blast_code/data/dNdzs_fullwidth.npz")`: Reads an `.npz` archive of N5K redshift-distribution data. The result is a dictionary-like container keyed by strings.
- `z_n5k = n5k_bins["z_cl"]`: Extracts the clustering redshift grid. `z_n5k` is a one-dimensional array with shape `(N_n5k,)`.
- `z_n5k = sort(z_n5k) # sort the redshift array in ascending order`: Sorts the redshift grid into ascending order. Sorting a vector of length `N` costs `O(N log N)`.
- `n_z_matrix = n5k_bins["dNdz_cl"];`: Extracts the clustering number-density matrix. This is a two-dimensional array with shape `(N_z, N_bin)`, where rows are redshift samples and columns are tomographic bins. 

```julia
p = plot()
for i in 1:size(n_z_matrix, 2)
    plot!(p, z_n5k, n_z_matrix[:, i], label="bin $i")
end
xlabel!(p, "z")
ylabel!(p, "dN/dz")
display(p)
```

- `p = plot()`: Creates a plot object configured for the kernel display. `p` is a plotting object rather than a scientific array. [file:1]
- `for i in 1:size(n_z_matrix, 2)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    plot!(p, z_n5k, n_z_matrix[:, i], label="bin $i")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `xlabel!(p, "z")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `ylabel!(p, "dN/dz")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `display(p)`: Renders the plot in the notebook output area. [file:1]

```julia
norms = [quadgk(x -> DataInterpolations.AkimaInterpolation(n_z_matrix[:, i], z_n5k, extrapolation=ExtrapolationType.Linear)(x),
                minimum(z_n5k), maximum(z_n5k))[1] for i in 1:size(n_z_matrix, 2)]
n_z_norm = n_z_matrix ./ reshape(norms, 1, :);
#### PLOTTING the results
p = plot()
for i in 1:size(n_z_norm, 2)
    plot!(p, z_n5k, n_z_norm[:, i], label="bin $i")
end
xlabel!(p, "z")
ylabel!(p, "dN/dz")
title!(p, "N5K Redshift Distribution (Clustering) - Normalized")
display(p)
```
- `norms = [quadgk(x -> DataInterpolations.AkimaInterpolation(n_z_matrix[:, i], z_n5k, extrapolation=ExtrapolationType.Linear)(x),`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content
- `minimum(z_n5k), maximum(z_n5k))[1] for i in 1:size(n_z_matrix, 2)]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `n_z_norm = n_z_matrix ./ reshape(norms, 1, :);`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `for i in 1:size(n_z_norm, 2)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot!(p, z_n5k, n_z_norm[:, i], label="bin $i")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

```julia
#### REDSHIFT ARRAY
z_array = LinRange(minimum(z_n5k), maximum(z_n5k), size(z_n5k)[1]) # array of redshifts to evaluate n(z) on
#### COMOVING DISTANCES ARRAY
chi_array = χ_of_z.(z_array); # array of comoving distances corresponding to z_range
```

### Cell role

Several lines define grids or allocate arrays whose dimensions set the scale of later computations. [file:1]

### Line-by-line

- `#### REDSHIFT ARRAY`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `z_array = LinRange(minimum(z_n5k), maximum(z_n5k), size(z_n5k)[1]) # array of redshifts to evaluate n(z) on`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `#### COMOVING DISTANCES ARRAY`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `chi_array = χ_of_z.(z_array); # array of comoving distances corresponding to z_range`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 13

```julia
open("dimensions.txt", "a") do f
    write(f, "------------------------ \n")
    write(f, "REDSHIFT ARRAY\n")
    write(f, "Redshift z_array is a vector or dimension: $(size(z_array))\n")
    write(f, "MIN = $(minimum(z_array)) || MAX = $(maximum(z_array))\n")
    write(f, "------------------------ \n")
    write(f, "COMOVING DISTANCE ARRAY\n")
    write(f, "Comoving distance chi_array is a vector or dimension: $(size(chi_array))\n")
    write(f, "MIN = $(minimum(chi_array)) || MAX = $(maximum(chi_array))\n")
    write(f, "------------------------ \n")
    write(f, "Values of normalized N5K dN/dz (clustering) for the first bin are: \n$(n_z_norm[:, 1])\n")
    write(f, "------------------------ \n")
end;
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `open("dimensions.txt", "a") do f`: Opens `dimensions.txt` in append mode so new diagnostics are added to the existing file. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "REDSHIFT ARRAY\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Redshift z_array is a vector or dimension: $(size(z_array))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN = $(minimum(z_array)) || MAX = $(maximum(z_array))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "COMOVING DISTANCE ARRAY\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Comoving distance chi_array is a vector or dimension: $(size(chi_array))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN = $(minimum(chi_array)) || MAX = $(maximum(chi_array))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Values of normalized N5K dN/dz (clustering) for the first bin are: \n$(n_z_norm[:, 1])\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `end;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 14

```julia
cosmo = Blast.FlatΛCDM()
n_chi = 10;                  # Number of comoving distance points
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `cosmo = Blast.FlatΛCDM()`: Constructs a flat-ΛCDM cosmology object using the local `Blast` definitions. This is a structured object with constant-size parameter fields. [file:1]
- `n_chi = 10;                  # Number of comoving distance points`: Defines the number of comoving-distance points as an integer scalar. It sets the length of several later vectors. [file:1]

## Code cell 15

```julia
open("dimensions.txt", "a") do f
    write(f, "------------------------ \n")
    write(f, "NUMBER OF COMOVING DISTANCE POINTS: $n_chi\n")
end;
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `open("dimensions.txt", "a") do f`: Opens `dimensions.txt` in append mode so new diagnostics are added to the existing file. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "NUMBER OF COMOVING DISTANCE POINTS: $n_chi\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `end;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 16

```julia
# Array of comoving distances
x_array = LinRange(26, 7000, n_chi)
# Array of redshifts corresponding to the comoving distances
z_range = z_of_χ.(x_array);
```

### Cell role

Several lines define grids or allocate arrays whose dimensions set the scale of later computations. [file:1]

### Line-by-line

- `# Array of comoving distances`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `x_array = LinRange(26, 7000, n_chi)`: Creates a linearly spaced comoving-distance grid from 26 to 7000 with `n_chi` points. It behaves like a one-dimensional vector of length `n_chi`. [file:1]
- `# Array of redshifts corresponding to the comoving distances`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `z_range = z_of_χ.(x_array);`: Evaluates the redshift interpolant elementwise on `x_array`, producing a one-dimensional redshift vector of shape `(n_chi,)`. Complexity is `O(n_chi)`. [file:1]

## Code cell 17

```julia
open("dimensions.txt", "a") do f
    write(f, "------------------------ \n")
    write(f, "COMOVING DISTANCES AND CORRESPONDING REDSHIFTS\n")
    write(f, "Array of comoving distances x_array has dimensions: $(size(x_array))\n")
    write(f, "MIN = $(minimum(x_array)) || MAX = $(maximum(x_array)) \n")
    write(f, "Array of redshifts z_range has dimensions: $(size(z_range))\n")
    write(f, "MIN = $(minimum(z_range)) || MAX = $(maximum(z_range)) \n")
    write(f, "------------------------ \n")
end;
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `open("dimensions.txt", "a") do f`: Opens `dimensions.txt` in append mode so new diagnostics are added to the existing file. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "COMOVING DISTANCES AND CORRESPONDING REDSHIFTS\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Array of comoving distances x_array has dimensions: $(size(x_array))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN = $(minimum(x_array)) || MAX = $(maximum(x_array)) \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Array of redshifts z_range has dimensions: $(size(z_range))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN = $(minimum(z_range)) || MAX = $(maximum(z_range)) \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `end;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 18

```julia
###### BIAS ######
# compute the bias. the equation is b(z) = b_0 * sqrt(1+z). we set b_0 = 1.0.!
b_0 = 1.0
b_z_array = b_0 .* sqrt.(1 .+ z_range); # bias as a function of redshift
```

### Cell role

Some lines define reusable functions whose cost is mostly paid when they are evaluated later. [file:1]

### Line-by-line

- `###### BIAS ######`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `# compute the bias. the equation is b(z) = b_0 * sqrt(1+z). we set b_0 = 1.0.!`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `b_0 = 1.0`: Defines the bias normalization as a scalar floating-point value. [file:1]
- `b_z_array = b_0 .* sqrt.(1 .+ z_range); # bias as a function of redshift`: Builds the bias vector from the model `b(z)=b_0 sqrt(1+z)`. The result has shape `(n_chi,)` and is computed in linear time. [file:1]

## Code cell 19

```julia
open("dimensions.txt", "a") do f
    write(f, "------------------------ \n")
    write(f, "BIAS ARRAY\n")
    write(f, "Array of bias values b_z_array has dimensions: $(size(b_z_array))\n")
    write(f, "MIN = $(minimum(b_z_array)) || MAX = $(maximum(b_z_array)) \n")
  end;
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `open("dimensions.txt", "a") do f`: Opens `dimensions.txt` in append mode so new diagnostics are added to the existing file. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "BIAS ARRAY\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Array of bias values b_z_array has dimensions: $(size(b_z_array))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN = $(minimum(b_z_array)) || MAX = $(maximum(b_z_array)) \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `  end;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 20

```julia
###### GROWTH FACTOR ######
#i want to obtain the growth factor D(z) from the background quantities. i can use the cosmology package to compute this.
function heath_integral(cosmo, z)
    integrand(zp) = (1.0 + zp) / (Blast.compute_adimensional_hubble_factor(zp, cosmo)^3)
    integral_val, _ = quadgk(integrand, z, Inf, rtol=1e-8)
    return Blast.compute_adimensional_hubble_factor(z, cosmo) * integral_val
end
function compute_growth_factor(cosmo, z_range)
    D_unnorm_z0 = heath_integral(cosmo, 0.0)
    D_array = [heath_integral(cosmo, zz) / D_unnorm_z0 for zz in z_range]
    return D_array
end
# Array of growth factor values
growth_array = compute_growth_factor(cosmo, z_range);
```

### Cell role

Some lines define reusable functions whose cost is mostly paid when they are evaluated later. [file:1]

### Line-by-line

- `###### GROWTH FACTOR ######`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `#i want to obtain the growth factor D(z) from the background quantities. i can use the cosmology package to compute this.`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `function heath_integral(cosmo, z)`: Begins a function definition for the unnormalized Heath growth integral. The function maps a cosmology object and a scalar redshift to a scalar result. [file:1]
- `    integrand(zp) = (1.0 + zp) / (Blast.compute_adimensional_hubble_factor(zp, cosmo)^3)`: Defines the scalar quadrature integrand used inside the growth calculation. [file:1]
- `    integral_val, _ = quadgk(integrand, z, Inf, rtol=1e-8)`: Computes an adaptive Gauss–Kronrod integral from `z` to infinity. Runtime is data-dependent because quadrature refines adaptively to the requested tolerance. [file:1]
- `    return Blast.compute_adimensional_hubble_factor(z, cosmo) * integral_val`: Returns the unnormalized growth quantity at redshift `z`. The result is a scalar. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `function compute_growth_factor(cosmo, z_range)`: Begins a helper function that computes a whole growth-factor vector over the redshift grid. [file:1]
- `    D_unnorm_z0 = heath_integral(cosmo, 0.0)`: Computes the normalization at redshift zero as a scalar. [file:1]
- `    D_array = [heath_integral(cosmo, zz) / D_unnorm_z0 for zz in z_range]`: Builds the normalized growth-factor vector with a comprehension. The output shape is `(length(z_range),)`, so here `(n_chi,)`. The cost is approximately one adaptive quadrature per grid point. [file:1]
- `    return D_array`: Returns the one-dimensional growth-factor array. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `# Array of growth factor values`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `growth_array = compute_growth_factor(cosmo, z_range);`: Evaluates the growth-factor function on the working grid. `growth_array` is a vector of length `n_chi`. [file:1]

## Code cell 21

```julia
open("dimensions.txt", "a") do f
    write(f, "------------------------ \n")
    write(f, "GROWTH FACTOR ARRAY\n")
    write(f, "Array of growth factor values D_z_array has dimensions: $(size(growth_array))\n")
    write(f, "MIN = $(minimum(growth_array)) || MAX = $(maximum(growth_array)) \n")
  end;
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `open("dimensions.txt", "a") do f`: Opens `dimensions.txt` in append mode so new diagnostics are added to the existing file. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "GROWTH FACTOR ARRAY\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Array of growth factor values D_z_array has dimensions: $(size(growth_array))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN = $(minimum(growth_array)) || MAX = $(maximum(growth_array)) \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `  end;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 22

```julia
p1 = plot(z_range, 
          b_z_array, 
          label=L"$b(z)$", 
          xlabel=L"$z$", 
          ylabel=L"$b(z)$", 
          title=L"$b(z)$")
p2 = plot(z_range, 
          x_array, 
          label=L"z and $\chi$",
          title = L"$\chi$",
          xlabel=L"$z$", 
          ylabel=L"$\chi$")
p3 = plot(z_range, 
          growth_array, 
          label="Growth Factor D(z)", 
          xlabel=L"$z$", 
          ylabel=L"$D(z)$", 
          title=L"$D(z)$")
plot(p1, p2, p3, layout=(1,3), titlefontsize=15, labelfontsize=15, legendfontsize=15, size=(1200,400), legend=false)
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `p1 = plot(z_range, `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          b_z_array, `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          label=L"$b(z)$", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          xlabel=L"$z$", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          ylabel=L"$b(z)$", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          title=L"$b(z)$")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `p2 = plot(z_range, `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          x_array, `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          label=L"z and $\chi$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          title = L"$\chi$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          xlabel=L"$z$", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          ylabel=L"$\chi$")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `p3 = plot(z_range, `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          growth_array, `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          label="Growth Factor D(z)", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          xlabel=L"$z$", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          ylabel=L"$D(z)$", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `          title=L"$D(z)$")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot(p1, p2, p3, layout=(1,3), titlefontsize=15, labelfontsize=15, legendfontsize=15, size=(1200,400), legend=false)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 23

```julia
open("dimensions.txt", "a") do f
    write(f, "------------------------ \n")
    write(f, "Array of redshifts: z_range has dimensions: $(size(z_range)) || z_array has dimensions: $(size(z_array))\n")
  end;
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `open("dimensions.txt", "a") do f`: Opens `dimensions.txt` in append mode so new diagnostics are added to the existing file. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Array of redshifts: z_range has dimensions: $(size(z_range)) || z_array has dimensions: $(size(z_array))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `  end;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 24

```julia
n_bins = size(n_z_norm, 2) # number of tomographic bins
nz = zeros(n_chi, n_bins) # size nx per number of bins
for b in 1:n_bins
     nz_func = DataInterpolations.AkimaInterpolation(n_z_norm[:, b], z_array, extrapolation=ExtrapolationType.Linear)
     nz_norm, _ = quadgk(nz_func, minimum(z_range), maximum(z_range), rtol=1e-8)
     #println("Bin $b: Normalized n(z) integral: ", nz_norm)
     nz[:, b] = nz_func.(z_range)
     #println("Bin $b: n(z) has size: ", size(nz))
end
```

### Cell role

It also creates interpolation objects, which are callable approximations built from tabulated arrays. Several lines define grids or allocate arrays whose dimensions set the scale of later computations. [file:1]

### Line-by-line

- `n_bins = size(n_z_norm, 2) # number of tomographic bins`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `nz = zeros(n_chi, n_bins) # size nx per number of bins`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `for b in 1:n_bins`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     nz_func = DataInterpolations.AkimaInterpolation(n_z_norm[:, b], z_array, extrapolation=ExtrapolationType.Linear)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     nz_norm, _ = quadgk(nz_func, minimum(z_range), maximum(z_range), rtol=1e-8)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     #println("Bin $b: Normalized n(z) integral: ", nz_norm)`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `     nz[:, b] = nz_func.(z_range)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     #println("Bin $b: n(z) has size: ", size(nz))`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 25

```julia
open("dimensions.txt", "a") do f
    write(f, "------------------------ \n")
    write(f, "Number of tomographic bins: $(n_bins)\n")
    write(f, "Size of n(z) array nz: $(size(nz))\n")
  end;
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `open("dimensions.txt", "a") do f`: Opens `dimensions.txt` in append mode so new diagnostics are added to the existing file. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Number of tomographic bins: $(n_bins)\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Size of n(z) array nz: $(size(nz))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `  end;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 26

```julia
prefac = x_array.^2 .* b_z_array .* growth_array
kernel_composed = reshape(prefac, :, 1) .* nz;
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `prefac = x_array.^2 .* b_z_array .* growth_array`: Forms the prefactor `χ^2 b(z) D(z)` elementwise. The result is a one-dimensional array of length `n_chi`, computed in `O(n_chi)`. [file:1]
- `kernel_composed = reshape(prefac, :, 1) .* nz;`: Reshapes `prefac` into a column and multiplies it against the tomography matrix `nz`. If `nz` has shape `(n_chi,n_bins)`, then `kernel_composed` has shape `(n_chi,n_bins)`. Complexity is `O(n_chi n_bins)`. [file:1]

## Code cell 27

```julia
open("dimensions.txt", "a") do f
    write(f, "------------------------ \n")
    write(f, "Number of tomographic bins n_bins: $(n_bins)\n")
    write(f, "Size of n(z) array nz: $(size(nz))\n")
    write(f, "Size of prefactor array (chi^2 * b(z) * D(z)): $(size(prefac))\n")
    write(f, "Size of kernel_composed array: $(size(kernel_composed))\n")
    write(f, "The prefactor array is the product of chi^2, b(z), and D(z) evaluated at the redshifts corresponding to the comoving distances in x_array.\n")
    write(f, "Size of array z_range: $(size(z_range))\n")
    write(f, "Size of array x_array: $(size(x_array))\n")
    write(f, "Size of array b_z_array: $(size(b_z_array))\n")
    write(f, "Size of array growth_array: $(size(growth_array))\n")
    write(f, "Size of array nz: $(size(nz))\n")
  end;
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `open("dimensions.txt", "a") do f`: Opens `dimensions.txt` in append mode so new diagnostics are added to the existing file. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Number of tomographic bins n_bins: $(n_bins)\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Size of n(z) array nz: $(size(nz))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Size of prefactor array (chi^2 * b(z) * D(z)): $(size(prefac))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Size of kernel_composed array: $(size(kernel_composed))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "The prefactor array is the product of chi^2, b(z), and D(z) evaluated at the redshifts corresponding to the comoving distances in x_array.\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Size of array z_range: $(size(z_range))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Size of array x_array: $(size(x_array))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Size of array b_z_array: $(size(b_z_array))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Size of array growth_array: $(size(growth_array))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Size of array nz: $(size(nz))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `  end;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 28

```julia
p = plot(
    xlabel=L"$z$",
    ylabel=L"\mathcal{K}(z)",
    title=L"\mathcal{K}(z)\ \mathrm{for\ all\ bins}",
    xlabelfontsize=15,
    ylabelfontsize=15,
    titlefontsize=20,
    legendfontsize=8,
    size=(900, 600)
)
for b in 1:size(kernel_composed, 2)
    plot!(p, z_range, kernel_composed[:, b], label="bin $b")
end
display(p)
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `p = plot(`: Creates a plot object configured for the kernel display. `p` is a plotting object rather than a scientific array. [file:1]
- `    xlabel=L"$z$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    ylabel=L"\mathcal{K}(z)",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    title=L"\mathcal{K}(z)\ \mathrm{for\ all\ bins}",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xlabelfontsize=15,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    ylabelfontsize=15,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    titlefontsize=20,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    legendfontsize=8,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    size=(900, 600)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `for b in 1:size(kernel_composed, 2)`: Loops over tomographic bins, so the iteration count is `n_bins`. [file:1]
- `    plot!(p, z_range, kernel_composed[:, b], label="bin $b")`: Adds one curve per bin to the plot. The sliced object `kernel_composed[:, b]` is a vector of length `n_chi`. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `display(p)`: Renders the plot in the notebook output area. [file:1]

## Code cell 29

```julia
ℓ = Blast.ℓ
x_min = 26                                            # minimum comoving distance
x_max = 7000                                          # maximum comoving distance
χ = LinRange(x_min, x_max, n_chi)                     # array of comoving distances
R = chebpoints(n_chi, -1, 1)                          # array of ratios between comoving distances
R = reverse(R[R.>0])
nR = length(R)                                        # number of ratios R
kmax = 200/13                                         # maximum wavenumber (small scales)
kmin = 2.5/x_max                                      # minimum wavenumber (large scales)
n_cheb = 120                                          # number of Chebyshev nodes 
β = 2                                                 # exponent depending on the probe: 0 for Galaxy - Shear, 2 for Galaxy - Galaxy, -2 for Shear - Shear
k_cheb = chebpoints(n_cheb, log10(kmin), log10(kmax)) # number of Chebyshev "points"
N = 2^(15)+1;                                         # number of integration points
k_array = k1_array = k2_array = 10 .^ chebpoints(n_cheb, log10(kmin), log10(kmax)); # wavenumbers k' and k'' on Chebyshev grid

# Define k' and k'' on Chebyshev grid
#k1 = k2 = 10 .^ chebpoints(n_cheb, log10(kmin), log10(kmax)) # wavenumbers k' and k''
```

### Cell role

Several lines define grids or allocate arrays whose dimensions set the scale of later computations. [file:1]

### Line-by-line

- `ℓ = Blast.ℓ`: Copies the multipole collection or multipole-related constant from `Blast` into the local variable `ℓ`. The exact type is defined in the local source module. [file:1]
- `x_min = 26                                            # minimum comoving distance`: Defines the minimum comoving distance as a scalar. [file:1]
- `x_max = 7000                                          # maximum comoving distance`: Defines the maximum comoving distance as a scalar. [file:1]
- `χ = LinRange(x_min, x_max, n_chi)                     # array of comoving distances`: Creates another comoving-distance grid of length `n_chi`. [file:1]
- `R = chebpoints(n_chi, -1, 1)                          # array of ratios between comoving distances`: Creates Chebyshev points on `[-1,1]`. The result is a one-dimensional array-like object of length `n_chi`. [file:1]
- `R = reverse(R[R.>0])`: Keeps only positive Chebyshev-ratio values and reverses their order. The result is a shorter one-dimensional array of length `nR`. [file:1]
- `nR = length(R)                                        # number of ratios R`: Stores the number of retained positive ratio nodes as an integer scalar. [file:1]
- `kmax = 200/13                                         # maximum wavenumber (small scales)`: Defines the maximum wavenumber as a scalar. [file:1]
- `kmin = 2.5/x_max                                      # minimum wavenumber (large scales)`: Defines the minimum wavenumber as a scalar. [file:1]
- `n_cheb = 120                                          # number of Chebyshev nodes `: Defines the number of Chebyshev nodes as an integer scalar. It controls several later matrix dimensions. [file:1]
- `β = 2                                                 # exponent depending on the probe: 0 for Galaxy - Shear, 2 for Galaxy - Galaxy, -2 for Shear - Shear`: Defines the probe-dependent exponent `β` as an integer scalar. [file:1]
- `k_cheb = chebpoints(n_cheb, log10(kmin), log10(kmax)) # number of Chebyshev "points"`: Creates Chebyshev nodes in logarithmic wavenumber space. The result is a vector of length `n_cheb`. [file:1]
- `N = 2^(15)+1;                                         # number of integration points`: Defines the number of integration points as the integer `32769`. [file:1]
- `k_array = k1_array = k2_array = 10 .^ chebpoints(n_cheb, log10(kmin), log10(kmax)); # wavenumbers k' and k'' on Chebyshev grid`: Creates one logarithmic wavenumber grid and assigns it to three variable names. Each behaves as a one-dimensional array of length `n_cheb`. [file:1]
- `# Define k' and k'' on Chebyshev grid`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `#k1 = k2 = 10 .^ chebpoints(n_cheb, log10(kmin), log10(kmax)) # wavenumbers k' and k''`: Comment line. It documents intent but does not execute numerical work. [file:1]

## Code cell 30

```julia
open("dimensions.txt", "a") do f
    write(f, "------------------------ \n")
    write(f, "INPUT PARAMETERS \n")
    write(f, "MIN COMOVING DISTANCE = $x_min, MAX COMOVING DISTANCE = $x_max, \n")
    write(f, "NUMBER OF COMOVING DISTANCE POINTS = $n_chi, \n")
    write(f, "NUMBER OF RATIOS R = $nR, \n")
    write(f, "MIN RATIO Rmin = $(minimum(R)), MAX RATIO Rmax = $(maximum(R)), \n")
    write(f, "MIN WAVENUMBER kmin = $kmin, MAX WAVENUMBER kmax = $kmax, \n")
    write(f, "NUMBER OF CHEBYSHEV NODES = $n_cheb, \n")
    write(f, "NUMBER OF INTEGRATION POINTS = $N, \n")
    write(f, "EXPONENT β = -2 for shear-shear, 0 for galaxy-shear and 2 for galaxy-galaxy \n")
  end;
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `open("dimensions.txt", "a") do f`: Opens `dimensions.txt` in append mode so new diagnostics are added to the existing file. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "INPUT PARAMETERS \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN COMOVING DISTANCE = $x_min, MAX COMOVING DISTANCE = $x_max, \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "NUMBER OF COMOVING DISTANCE POINTS = $n_chi, \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "NUMBER OF RATIOS R = $nR, \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN RATIO Rmin = $(minimum(R)), MAX RATIO Rmax = $(maximum(R)), \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN WAVENUMBER kmin = $kmin, MAX WAVENUMBER kmax = $kmax, \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "NUMBER OF CHEBYSHEV NODES = $n_cheb, \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "NUMBER OF INTEGRATION POINTS = $N, \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "EXPONENT β = -2 for shear-shear, 0 for galaxy-shear and 2 for galaxy-galaxy \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `  end;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 31

```julia
# I choose a grid of Chebyshev nodes in z, which I will use to evaluate 
# the kernels and the power spectrum.
zmin = minimum(z_range)
zmax = maximum(z_range)

t_of_z(z) = (2*z - (zmax + zmin)) / (zmax - zmin)
z_of_t(t) = (zmin + zmax)/2 + (zmax - zmin)/2 * t

j = 1:n_cheb
#t_nodes = cos.((2 .* j .- 1) .* pi ./ (2*n_cheb))     # size N, in [-1,1]
t_nodes = [cos(π * k / (n_cheb - 1)) for k in 0:n_cheb-1]
z_nodes = @. zmin + (zmax - zmin) / 2 * (1 + t_nodes);
#z_nodes = z_of_t.(t_nodes);                            # Chebyshev nodes in z
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `# I choose a grid of Chebyshev nodes in z, which I will use to evaluate `: Comment line. It documents intent but does not execute numerical work. [file:1]
- `# the kernels and the power spectrum.`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `zmin = minimum(z_range)`: Computes and stores the minimum of the redshift grid as a scalar. [file:1]
- `zmax = maximum(z_range)`: Computes and stores the maximum of the redshift grid as a scalar. [file:1]
- `t_of_z(z) = (2*z - (zmax + zmin)) / (zmax - zmin)`: Defines the affine map from redshift `z` to the Chebyshev coordinate `t` in `[-1,1]`. [file:1]
- `z_of_t(t) = (zmin + zmax)/2 + (zmax - zmin)/2 * t`: Defines the inverse affine map from `t` back to redshift `z`. [file:1]
- `j = 1:n_cheb`: Creates a unit-range index object from 1 to `n_cheb`. [file:1]
- `#t_nodes = cos.((2 .* j .- 1) .* pi ./ (2*n_cheb))     # size N, in [-1,1]`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `t_nodes = [cos(π * k / (n_cheb - 1)) for k in 0:n_cheb-1]`: Explicitly constructs the Chebyshev–Lobatto nodes. `t_nodes` is a one-dimensional array of length `n_cheb`, built in `O(n_cheb)`. [file:1]
- `z_nodes = @. zmin + (zmax - zmin) / 2 * (1 + t_nodes);`: Maps the Chebyshev nodes from `t`-space to redshift space. `z_nodes` is a one-dimensional array of length `n_cheb`. [file:1]
- `#z_nodes = z_of_t.(t_nodes);                            # Chebyshev nodes in z`: Comment line. It documents intent but does not execute numerical work. [file:1]

## Code cell 32

```julia
open("dimensions.txt", "a") do f
    write(f, "------------------------ \n")
    write(f, "REDSHIFT GRID")
    write(f, "MIN REDSHIFT = $zmin, MAX REDSHIFT = $zmax, \n")
    write(f, "Dimension of z_nodes = $n_cheb \n")
    write(f, "Dimension of t_nodes = $n_cheb \n")
    write(f, "MIN Chebyshev node in z: $(z_nodes[1]), MAX Chebyshev node in z: $(z_nodes[end]) \n")
    write(f, "MIN Chebyshev node in t: $(t_nodes[1]), MAX Chebyshev node in t: $(t_nodes[end]) \n")
  end;
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `open("dimensions.txt", "a") do f`: Opens `dimensions.txt` in append mode so new diagnostics are added to the existing file. [file:1]
- `    write(f, "------------------------ \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "REDSHIFT GRID")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN REDSHIFT = $zmin, MAX REDSHIFT = $zmax, \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Dimension of z_nodes = $n_cheb \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Dimension of t_nodes = $n_cheb \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN Chebyshev node in z: $(z_nodes[1]), MAX Chebyshev node in z: $(z_nodes[end]) \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN Chebyshev node in t: $(t_nodes[1]), MAX Chebyshev node in t: $(t_nodes[end]) \n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `  end;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 33

```julia
h1 = histogram(z_nodes, bins=20, 
               label="Chebyshev nodes in z", 
               xlabel="Redshift z", 
               ylabel="Frequency", 
               title="Distribution of Chebyshev Nodes in z",
               legendfontsize=10,
               legend=:bottomright)
h2 = histogram(t_nodes, bins=20, 
               label="Chebyshev nodes in t", 
               xlabel="t", 
               ylabel="Frequency", 
               title="Distribution of Chebyshev Nodes in t",
               legendfontsize=10,
               legend=:bottomright)
plot(h1, h2, layout=(1,2), size=(1200,400))
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `h1 = histogram(z_nodes, bins=20, `: Creates a histogram plot object for `z_nodes`. [file:1]
- `               label="Chebyshev nodes in z", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `               xlabel="Redshift z", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `               ylabel="Frequency", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `               title="Distribution of Chebyshev Nodes in z",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `               legendfontsize=10,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `               legend=:bottomright)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `h2 = histogram(t_nodes, bins=20, `: Creates a histogram plot object for `t_nodes`. [file:1]
- `               label="Chebyshev nodes in t", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `               xlabel="t", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `               ylabel="Frequency", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `               title="Distribution of Chebyshev Nodes in t",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `               legendfontsize=10,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `               legend=:bottomright)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot(h1, h2, layout=(1,2), size=(1200,400))`: Combines the two histograms into one two-panel plot. [file:1]

## Code cell 34

```julia
chi_interp    = AkimaInterpolation(x_array, z_range, extrapolation=ExtrapolationType.Linear)
bias_interp   = AkimaInterpolation(b_z_array, z_range, extrapolation=ExtrapolationType.Linear)
growth_interp = AkimaInterpolation(growth_array, z_range, extrapolation=ExtrapolationType.Linear)
# nz has shape (96, 10); build one interpolator per tomographic bin
nz_interp = [AkimaInterpolation(nz[:, i], z_range,
                                 extrapolation=ExtrapolationType.Linear)
             for i in 1:n_bins];
K_vals = zeros(n_cheb, n_bins)
for i in 1:n_bins
    for k in 1:n_cheb
        z  = z_nodes[k]
        χ  = chi_interp(z)
        b  = bias_interp(z)
        D  = growth_interp(z)
        ni = max(0.0, nz_interp[i](z))   # clamp fisico
        K_vals[k, i] = χ^2 * b * D * ni
    end
end
# ── 6. Plan the DCT-I (REDFT00) once, reuse for all bins ────────────────────
#       plan_fft mirrors the BLAST helper exactly
function plan_fft_1d(vals::AbstractVector)
    kind = size(vals, 1) > 1 ? FFTW.REDFT00 : FFTW.DHT
    p = FFTW.plan_r2r(deepcopy(vals), kind, [1];
                      flags=FFTW.PATIENT, timelimit=Inf)
    return p
end
function fast_chebcoefs_1d(vals::AbstractVector, plan::FFTW.r2rFFTWPlan)
    coefs = plan * vals
    N = length(coefs)
    coefs ./= 2 * (N - 1)          # normalise
    coefs[2:end-1] .*= 2           # double interior coefficients
    return coefs
end

plan = plan_fft_1d(K_vals[:, 1]);   # plan built on a single column (n_cheb elements)

# ── 7. Compute Chebyshev coefficients  c[bin, node_index] ───────────────────
#       Shape: (n_bins, n_cheb)  →  c_n^(i) for bin i, order n
cheb_coeff_K = zeros(n_bins, n_cheb)

for i in 1:n_bins
    cheb_coeff_K[i, :] = fast_chebcoefs_1d(K_vals[:, i], plan)
end
```

### Cell role

It also creates interpolation objects, which are callable approximations built from tabulated arrays. Some lines define reusable functions whose cost is mostly paid when they are evaluated later. Several lines define grids or allocate arrays whose dimensions set the scale of later computations. The main algorithmic idea here is a spectral decomposition using FFTW-backed real transforms. [file:1]

### Line-by-line

- `chi_interp    = AkimaInterpolation(x_array, z_range, extrapolation=ExtrapolationType.Linear)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `bias_interp   = AkimaInterpolation(b_z_array, z_range, extrapolation=ExtrapolationType.Linear)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `growth_interp = AkimaInterpolation(growth_array, z_range, extrapolation=ExtrapolationType.Linear)`: Creates an interpolation object intended to return `D(z)` from redshift input. [file:1]
- `# nz has shape (96, 10); build one interpolator per tomographic bin`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `nz_interp = [AkimaInterpolation(nz[:, i], z_range,`: Builds one interpolation object per tomographic bin and stores them in a one-dimensional container of length `n_bins`. [file:1]
- `                                 extrapolation=ExtrapolationType.Linear)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `             for i in 1:n_bins];`: Starts the outer loop over tomographic bins. [file:1]
- `K_vals = zeros(n_cheb, n_bins)`: Allocates the kernel-value matrix with shape `(n_cheb,n_bins)`. Memory and fill cost scale as `O(n_cheb n_bins)`. [file:1]
- `for i in 1:n_bins`: Starts the outer loop over tomographic bins. [file:1]
- `    for k in 1:n_cheb`: Starts the inner loop over Chebyshev nodes. The nested structure creates `n_bins * n_cheb` scalar evaluations. [file:1]
- `        z  = z_nodes[k]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        χ  = chi_interp(z)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        b  = bias_interp(z)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        D  = growth_interp(z)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        ni = max(0.0, nz_interp[i](z))   # clamp fisico`: Evaluates the `i`-th number-density interpolation and clamps it below by zero to avoid negative interpolation artifacts. `ni` is a scalar. [file:1]
- `        K_vals[k, i] = χ^2 * b * D * ni`: Stores one kernel value in the matrix `K_vals`. Filling the entire matrix is `O(n_cheb n_bins)`. [file:1]
- `    end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `# ── 6. Plan the DCT-I (REDFT00) once, reuse for all bins ────────────────────`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `#       plan_fft mirrors the BLAST helper exactly`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `function plan_fft_1d(vals::AbstractVector)`: Defines a helper to create an FFTW real-to-real transform plan for one-dimensional vectors. [file:1]
- `    kind = size(vals, 1) > 1 ? FFTW.REDFT00 : FFTW.DHT`: Chooses DCT-I for vectors longer than one element and a Hartley transform otherwise. [file:1]
- `    p = FFTW.plan_r2r(deepcopy(vals), kind, [1];`: Constructs an FFTW plan object on a copied vector. Planning may be expensive once, but it is amortized over repeated transforms. [file:1]
- `                      flags=FFTW.PATIENT, timelimit=Inf)`: Requests a thorough FFTW planning strategy with no time limit. This increases setup time to reduce later execution time. [file:1]
- `    return p`: Returns the FFTW plan object. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `function fast_chebcoefs_1d(vals::AbstractVector, plan::FFTW.r2rFFTWPlan)`: Defines a helper that converts nodal values into Chebyshev coefficients using the supplied FFTW plan. [file:1]
- `    coefs = plan * vals`: Applies the planned transform, returning a coefficient vector of the same length as `vals`. The transform cost is approximately `O(n_cheb log n_cheb)`. [file:1]
- `    N = length(coefs)`: Stores the transform length as an integer scalar. [file:1]
- `    coefs ./= 2 * (N - 1)          # normalise`: Normalizes all coefficients in place. This is a linear-time vector operation. [file:1]
- `    coefs[2:end-1] .*= 2           # double interior coefficients`: Doubles the interior coefficients in place, leaving endpoint coefficients unchanged. [file:1]
- `    return coefs`: Returns the coefficient vector. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plan = plan_fft_1d(K_vals[:, 1]);   # plan built on a single column (n_cheb elements)`: Builds the FFTW plan from the first column of `K_vals`, which is a vector of length `n_cheb`. The same plan can be reused for every bin because all columns share the same size. [file:1]
- `# ── 7. Compute Chebyshev coefficients  c[bin, node_index] ───────────────────`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `#       Shape: (n_bins, n_cheb)  →  c_n^(i) for bin i, order n`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `cheb_coeff_K = zeros(n_bins, n_cheb)`: Allocates the Chebyshev-coefficient matrix with shape `(n_bins,n_cheb)`. [file:1]
- `for i in 1:n_bins`: Starts the outer loop over tomographic bins. [file:1]
- `    cheb_coeff_K[i, :] = fast_chebcoefs_1d(K_vals[:, i], plan)`: Computes the Chebyshev coefficients for the `i`-th tomographic bin and stores them in row `i`. Across all bins, this costs about `O(n_bins n_cheb log n_cheb)`. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 35

```julia
open("dimensions.txt", "a") do f
    write(f, "------------------------\n")
    write(f, "CHEBYSHEV DECOMPOSITION OF WINDOW KERNEL K(z)\n")
    write(f, "K(z) = chi^2(z) * b(z) * D(z) * n(z)\n")
    write(f, "Number of Chebyshev nodes n_cheb = $n_cheb\n")
    write(f, "Chebyshev nodes in t: t_nodes has dimension $(size(t_nodes))\n")
    write(f, "MIN t = $(minimum(t_nodes)) || MAX t = $(maximum(t_nodes))\n")
    write(f, "Chebyshev nodes in z: z_nodes has dimension $(size(z_nodes))\n")
    write(f, "MIN z = $(minimum(z_nodes)) || MAX z = $(maximum(z_nodes))\n")
    write(f, "K_vals shape (n_cheb × n_bins) = $(size(K_vals))\n")
    write(f, "cheb_coeff_K shape (n_bins × n_cheb) = $(size(cheb_coeff_K))\n")
    write(f, "MIN coeff = $(minimum(cheb_coeff_K)) || MAX coeff = $(maximum(cheb_coeff_K))\n")
    write(f, "------------------------\n")
end
```

### Cell role

The main algorithmic idea here is a spectral decomposition using FFTW-backed real transforms. [file:1]

### Line-by-line

- `open("dimensions.txt", "a") do f`: Opens `dimensions.txt` in append mode so new diagnostics are added to the existing file. [file:1]
- `    write(f, "------------------------\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "CHEBYSHEV DECOMPOSITION OF WINDOW KERNEL K(z)\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "K(z) = chi^2(z) * b(z) * D(z) * n(z)\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Number of Chebyshev nodes n_cheb = $n_cheb\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Chebyshev nodes in t: t_nodes has dimension $(size(t_nodes))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN t = $(minimum(t_nodes)) || MAX t = $(maximum(t_nodes))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "Chebyshev nodes in z: z_nodes has dimension $(size(z_nodes))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN z = $(minimum(z_nodes)) || MAX z = $(maximum(z_nodes))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "K_vals shape (n_cheb × n_bins) = $(size(K_vals))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "cheb_coeff_K shape (n_bins × n_cheb) = $(size(cheb_coeff_K))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "MIN coeff = $(minimum(cheb_coeff_K)) || MAX coeff = $(maximum(cheb_coeff_K))\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `    write(f, "------------------------\n")`: Writes text diagnostics to the file stream `f`. Each write is linear in the string length but negligible compared with the scientific numerics. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 36

```julia
i_bin = 1 # choose a bin to plot

p1 = plot(z_nodes, bias_interp.(z_nodes),
    label = "b(z)",
    xlabel = L"$z$", ylabel = L"$b(z)$",
    title = L"$b(z)$ at Chebyshev Nodes")

p2 = plot(z_nodes, growth_interp.(z_nodes),
    label = "D(z)",
    xlabel = L"$z$", ylabel = L"$D(z)$",
    title = L"$D(z)$ at Chebyshev Nodes")

p3 = plot(z_nodes, nz_interp[i_bin].(z_nodes),
    label = "n(z)",
    xlabel = L"$z$", ylabel = L"$n(z)$",
    title = "n(z) at Chebyshev Nodes, bin $i_bin")

p4 = plot(z_nodes, chi_interp.(z_nodes) .^ 2,
    label = "χ(z)^2",
    xlabel = L"$z$", ylabel = L"$\chi(z)^2$",
    title = L"$\chi(z)^2$ at Chebyshev Nodes")

p5 = plot(z_nodes, K_vals[:, i_bin],
    label = "K(z)",
    xlabel = L"$z$", ylabel = L"$K(z)$",
    title = "K(z) at Chebyshev Nodes, bin $i_bin")

p6 = plot(1:length(cheb_coeff_K[i_bin, :]), cheb_coeff_K[i_bin, :],
    label = "Chebyshev coefficients",
    xlabel = L"$n_{\mathrm{Cheb}}$", ylabel = L"$c_n$",
    title = "Chebyshev Coefficients, bin $i_bin")

plot(p1, p2, p3, p4, p5, p6, layout = (3, 2), size = (1400, 1400),
     legend = false, labelfont = 15, titlefontsize = 20,
     xlabelfontsize = 15, ylabelfontsize = 15)
```

### Cell role

The main algorithmic idea here is a spectral decomposition using FFTW-backed real transforms. [file:1]

### Line-by-line

- `i_bin = 1 # choose a bin to plot`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `p1 = plot(z_nodes, bias_interp.(z_nodes),`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    label = "b(z)",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xlabel = L"$z$", ylabel = L"$b(z)$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    title = L"$b(z)$ at Chebyshev Nodes")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `p2 = plot(z_nodes, growth_interp.(z_nodes),`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    label = "D(z)",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xlabel = L"$z$", ylabel = L"$D(z)$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    title = L"$D(z)$ at Chebyshev Nodes")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `p3 = plot(z_nodes, nz_interp[i_bin].(z_nodes),`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    label = "n(z)",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xlabel = L"$z$", ylabel = L"$n(z)$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    title = "n(z) at Chebyshev Nodes, bin $i_bin")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `p4 = plot(z_nodes, chi_interp.(z_nodes) .^ 2,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    label = "χ(z)^2",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xlabel = L"$z$", ylabel = L"$\chi(z)^2$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    title = L"$\chi(z)^2$ at Chebyshev Nodes")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `p5 = plot(z_nodes, K_vals[:, i_bin],`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    label = "K(z)",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xlabel = L"$z$", ylabel = L"$K(z)$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    title = "K(z) at Chebyshev Nodes, bin $i_bin")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `p6 = plot(1:length(cheb_coeff_K[i_bin, :]), cheb_coeff_K[i_bin, :],`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    label = "Chebyshev coefficients",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xlabel = L"$n_{\mathrm{Cheb}}$", ylabel = L"$c_n$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    title = "Chebyshev Coefficients, bin $i_bin")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot(p1, p2, p3, p4, p5, p6, layout = (3, 2), size = (1400, 1400),`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     legend = false, labelfont = 15, titlefontsize = 20,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     xlabelfontsize = 15, ylabelfontsize = 15)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 37

```julia
inv_Hz_array = [1.0 / Blast.compute_adimensional_hubble_factor(z, cosmo) for z in z_range]
inv_Hz_interp = AkimaInterpolation(inv_Hz_array, z_range, extrapolation=ExtrapolationType.Linear)

# ── Evaluate K̃_i(z) = chi^2 * b * D * n_i * (1/H) at Chebyshev nodes ───────
K_tilde_vals = zeros(n_cheb, n_bins)
for i in 1:n_bins, k in 1:n_cheb
    z = z_nodes[k]
    K_tilde_vals[k, i] = chi_interp(z)^2 * bias_interp(z) * growth_interp(z) *
                          nz_interp[i](z) * inv_Hz_interp(z)
end

cheb_coeff_K_tilde = zeros(n_bins, n_cheb)
for i in 1:n_bins
    cheb_coeff_K_tilde[i, :] = fast_chebcoefs_1d(K_tilde_vals[:, i], plan)
end
```

### Cell role

It also creates interpolation objects, which are callable approximations built from tabulated arrays. Several lines define grids or allocate arrays whose dimensions set the scale of later computations. The main algorithmic idea here is a spectral decomposition using FFTW-backed real transforms. [file:1]

### Line-by-line

- `inv_Hz_array = [1.0 / Blast.compute_adimensional_hubble_factor(z, cosmo) for z in z_range]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `inv_Hz_interp = AkimaInterpolation(inv_Hz_array, z_range, extrapolation=ExtrapolationType.Linear)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `# ── Evaluate K̃_i(z) = chi^2 * b * D * n_i * (1/H) at Chebyshev nodes ───────`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `K_tilde_vals = zeros(n_cheb, n_bins)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `for i in 1:n_bins, k in 1:n_cheb`: Starts the outer loop over tomographic bins. [file:1]
- `    z = z_nodes[k]`: Extracts one scalar redshift node. [file:1]
- `    K_tilde_vals[k, i] = chi_interp(z)^2 * bias_interp(z) * growth_interp(z) *`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `                          nz_interp[i](z) * inv_Hz_interp(z)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `cheb_coeff_K_tilde = zeros(n_bins, n_cheb)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `for i in 1:n_bins`: Starts the outer loop over tomographic bins. [file:1]
- `    cheb_coeff_K_tilde[i, :] = fast_chebcoefs_1d(K_tilde_vals[:, i], plan)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 38

```julia
plot(1:length(cheb_coeff_K[i_bin, :]), cheb_coeff_K[i_bin, :],
    label = "Chebyshev coefficients",
    xlabel = L"$n_{\mathrm{Cheb}}$", ylabel = L"$c_n$")

plot!(1:length(cheb_coeff_K_tilde[i_bin, :]), cheb_coeff_K_tilde[i_bin, :],
label = "Chebyshev coefficients - modified with H(z)",
xlabel = L"$n_{\mathrm{Cheb}}$", ylabel = L"$c_n$",
title = "Chebyshev Coefficients, bin $i_bin")
```

### Cell role

The main algorithmic idea here is a spectral decomposition using FFTW-backed real transforms. [file:1]

### Line-by-line

- `plot(1:length(cheb_coeff_K[i_bin, :]), cheb_coeff_K[i_bin, :],`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    label = "Chebyshev coefficients",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xlabel = L"$n_{\mathrm{Cheb}}$", ylabel = L"$c_n$")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot!(1:length(cheb_coeff_K_tilde[i_bin, :]), cheb_coeff_K_tilde[i_bin, :],`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `label = "Chebyshev coefficients - modified with H(z)",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `xlabel = L"$n_{\mathrm{Cheb}}$", ylabel = L"$c_n$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `title = "Chebyshev Coefficients, bin $i_bin")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 39

```julia
plot(z_nodes, K_vals[:, i_bin],
    label = "K(z)",
    xlabel = L"$z$", ylabel = L"$K(z)$",
    title = "K(z) at Chebyshev Nodes, bin $i_bin")

plot!(z_nodes, K_tilde_vals[:, i_bin],
    label = "K(z) with H(z)",
    xlabel = L"$z$", ylabel = L"$K(z)$")
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `plot(z_nodes, K_vals[:, i_bin],`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    label = "K(z)",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xlabel = L"$z$", ylabel = L"$K(z)$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    title = "K(z) at Chebyshev Nodes, bin $i_bin")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot!(z_nodes, K_tilde_vals[:, i_bin],`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    label = "K(z) with H(z)",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xlabel = L"$z$", ylabel = L"$K(z)$")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 40

```julia
plots = Plots.Plot[]
for i_bin in 1:n_bins
    p = plot(1:length(cheb_coeff_K[i_bin, :]), cheb_coeff_K[i_bin, :],
        label = "Chebyshev coefficients",
        xlabel = L"$n_{\mathrm{Cheb}}$", ylabel = L"$c_n$",
        title = "bin $i_bin")

    plot!(p, 1:length(cheb_coeff_K_tilde[i_bin, :]), cheb_coeff_K_tilde[i_bin, :],
        label = "Chebyshev coefficients - modified with H(z)")

    push!(plots, p)
end

ncols = ceil(Int, sqrt(n_bins))
nrows = ceil(Int, n_bins / ncols)
plot(plots..., layout = (nrows, ncols), size = (1400, 1400))
```

### Cell role

The main algorithmic idea here is a spectral decomposition using FFTW-backed real transforms. [file:1]

### Line-by-line

- `plots = Plots.Plot[]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `for i_bin in 1:n_bins`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    p = plot(1:length(cheb_coeff_K[i_bin, :]), cheb_coeff_K[i_bin, :],`: Creates a plot object configured for the kernel display. `p` is a plotting object rather than a scientific array. [file:1]
- `        label = "Chebyshev coefficients",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        xlabel = L"$n_{\mathrm{Cheb}}$", ylabel = L"$c_n$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        title = "bin $i_bin")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    plot!(p, 1:length(cheb_coeff_K_tilde[i_bin, :]), cheb_coeff_K_tilde[i_bin, :],`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        label = "Chebyshev coefficients - modified with H(z)")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    push!(plots, p)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `ncols = ceil(Int, sqrt(n_bins))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `nrows = ceil(Int, n_bins / ncols)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot(plots..., layout = (nrows, ncols), size = (1400, 1400))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 41

```julia
plots = Plots.Plot[]
for i_bin in 1:n_bins
    p = plot(z_nodes, K_vals[:, i_bin],
        label = "K(z)",
        xlabel = L"$z$", ylabel = L"$K(z)$",
        title = "bin $i_bin")

    plot!(p, z_nodes, K_tilde_vals[:, i_bin],
        label = "K(z) with H(z)")

    push!(plots, p)
end

plot(plots..., layout = (ceil(Int, sqrt(n_bins)), ceil(Int, n_bins / ceil(Int, sqrt(n_bins)))), size = (1400, 1400))
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `plots = Plots.Plot[]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `for i_bin in 1:n_bins`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    p = plot(z_nodes, K_vals[:, i_bin],`: Creates a plot object configured for the kernel display. `p` is a plotting object rather than a scientific array. [file:1]
- `        label = "K(z)",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        xlabel = L"$z$", ylabel = L"$K(z)$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        title = "bin $i_bin")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    plot!(p, z_nodes, K_tilde_vals[:, i_bin],`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        label = "K(z) with H(z)")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    push!(plots, p)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot(plots..., layout = (ceil(Int, sqrt(n_bins)), ceil(Int, n_bins / ceil(Int, sqrt(n_bins)))), size = (1400, 1400))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 42

```julia
println("MIN k = $(minimum(k_array)) || MAX k = $(maximum(k_array)) \n")
println("MIN k1 = $(minimum(k1_array)) || MAX k1 = $(maximum(k1_array)) \n")
println("MIN k2 = $(minimum(k2_array)) || MAX k2 = $(maximum(k2_array)) \n")
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `println("MIN k = $(minimum(k_array)) || MAX k = $(maximum(k_array)) \n")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("MIN k1 = $(minimum(k1_array)) || MAX k1 = $(maximum(k1_array)) \n")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("MIN k2 = $(minimum(k2_array)) || MAX k2 = $(maximum(k2_array)) \n")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 43

```julia
nk  = length(k_array)
nk1 = length(k1_array)
nk2 = length(k2_array)

W_base_k1 = zeros(Float64, length(ℓ), nk1, nk, n_cheb)
W_base_k2 = zeros(Float64, length(ℓ), nk2, nk, n_cheb)

computation_time = @elapsed begin
    for i in eachindex(ℓ)
        W_base_k1[i, :, :, :] = Blast.compute_W_tilde_modes(ℓ[i], k_array, k1_array, z_nodes, chi_interp, zmin, zmax;
                                                      n_cheb=n_cheb)

        W_base_k2[i, :, :, :] = Blast.compute_W_tilde_modes(ℓ[i], k_array, k2_array, z_nodes, chi_interp, zmin, zmax;
                                                      n_cheb=n_cheb)

        println("step ", i, " completed")
    end
end
println("Total computation time for W_base_k1 and W_base_k2: ", computation_time, " seconds")
```

### Cell role

Several lines define grids or allocate arrays whose dimensions set the scale of later computations. [file:1]

### Line-by-line

- `nk  = length(k_array)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `nk1 = length(k1_array)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `nk2 = length(k2_array)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `W_base_k1 = zeros(Float64, length(ℓ), nk1, nk, n_cheb)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `W_base_k2 = zeros(Float64, length(ℓ), nk2, nk, n_cheb)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `computation_time = @elapsed begin`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    for i in eachindex(ℓ)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        W_base_k1[i, :, :, :] = Blast.compute_W_tilde_modes(ℓ[i], k_array, k1_array, z_nodes, chi_interp, zmin, zmax;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `                                                      n_cheb=n_cheb)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        W_base_k2[i, :, :, :] = Blast.compute_W_tilde_modes(ℓ[i], k_array, k2_array, z_nodes, chi_interp, zmin, zmax;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `                                                      n_cheb=n_cheb)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        println("step ", i, " completed")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("Total computation time for W_base_k1 and W_base_k2: ", computation_time, " seconds")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 44

```julia
print("Shape of the computed W_tilde: ", size(W_base_k1), "\n")
print("Shape of the computed W_tilde: ", size(W_base_k2), "\n")
print("Shape of the computed cheb_coeff_K_tilde: ", size(cheb_coeff_K_tilde), "\n")
```

### Cell role

The main algorithmic idea here is a spectral decomposition using FFTW-backed real transforms. [file:1]

### Line-by-line

- `print("Shape of the computed W_tilde: ", size(W_base_k1), "\n")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `print("Shape of the computed W_tilde: ", size(W_base_k2), "\n")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `print("Shape of the computed cheb_coeff_K_tilde: ", size(cheb_coeff_K_tilde), "\n")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 45

```julia
W_tilde_k1 = @tullio W[li, b, ik1, ik] := cheb_coeff_K_tilde[b, n] * W_base_k1[li, ik1, ik, n]
W_tilde_k2 = @tullio W[li, b, ik2, ik] := cheb_coeff_K_tilde[b, n] * W_base_k2[li, ik2, ik, n];
```

### Cell role

The main algorithmic idea here is a spectral decomposition using FFTW-backed real transforms. [file:1]

### Line-by-line

- `W_tilde_k1 = @tullio W[li, b, ik1, ik] := cheb_coeff_K_tilde[b, n] * W_base_k1[li, ik1, ik, n]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `W_tilde_k2 = @tullio W[li, b, ik2, ik] := cheb_coeff_K_tilde[b, n] * W_base_k2[li, ik2, ik, n];`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 46

```julia
li = 1
nmode = 1

heatmap(1:nk1, 1:nk, W_base_k1[li, :, :, nmode] ./ maximum(abs.(W_base_k1[li, :, :, nmode])),
    size = (500, 500),
    title = L"W_{\mathrm{base},k_1}(\ell=\text{first}, n=\text{first})",
    c = :roma,
    xlabel = L"k_1",
    ylabel = L"k",
    legend = :none,
    yguidefontsize = 15,
    xguidefontsize = 15,
    titlefontsize = 20)
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `li = 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `nmode = 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `heatmap(1:nk1, 1:nk, W_base_k1[li, :, :, nmode] ./ maximum(abs.(W_base_k1[li, :, :, nmode])),`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    size = (500, 500),`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    title = L"W_{\mathrm{base},k_1}(\ell=\text{first}, n=\text{first})",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    c = :roma,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xlabel = L"k_1",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    ylabel = L"k",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    legend = :none,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    yguidefontsize = 15,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xguidefontsize = 15,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    titlefontsize = 20)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 47

```julia
li = 1
nmode = 1

heatmap(1:nk2, 1:nk, W_base_k2[li, :, :, nmode] ./ maximum(abs.(W_base_k1[li, :, :, nmode])),
    size = (500, 500),
    title = L"W_{\mathrm{base},k_2}(\ell=\text{first}, n=\text{first})",
    c = :roma,
    xlabel = L"k_2",
    ylabel = L"k",
    legend = :none,
    yguidefontsize = 15,
    xguidefontsize = 15,
    titlefontsize = 20)
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `li = 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `nmode = 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `heatmap(1:nk2, 1:nk, W_base_k2[li, :, :, nmode] ./ maximum(abs.(W_base_k1[li, :, :, nmode])),`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    size = (500, 500),`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    title = L"W_{\mathrm{base},k_2}(\ell=\text{first}, n=\text{first})",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    c = :roma,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xlabel = L"k_2",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    ylabel = L"k",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    legend = :none,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    yguidefontsize = 15,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    xguidefontsize = 15,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    titlefontsize = 20)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Markdown cell 5

This cell is explanatory prose introducing the topic `## Now the Power Spectrum`. It creates no runtime objects and has no numerical complexity. [file:1]

## Code cell 48

```julia
#3D matter power spectrum
pk_dict = npzread("blast_code/data/pk.npz")
Pklin = pk_dict["pk_lin"]
Pknonlin = pk_dict["pk_nl"]
k = pk_dict["k"]
z = pk_dict["z"];

println("Dimension of k points: ", size(k)[1])
println("Dimension of Pklin: ", size(Pklin)[1])
println("Dimension of Pknonlin: ", size(Pknonlin)[1])

#Interpolating the power spectrum
#Linear P(k)
# y has the same length as k, and is a logarithmic range from log10(first(k)) to log10(last(k))
y = LinRange(log10(first(k)),log10(last(k)), length(k))
# x has the same length as z, and is a linear range from first(z) to last(z)
x = LinRange(first(z), last(z), length(z))

# this creates an interpolation of the linear P(k) data. 
# the interpolatio is in log space because P(k) changes by many order or magnitude
InterpPmm = Interpolations.interpolate(log10.(Pklin),BSpline(Cubic(Line(OnGrid()))))
# this scales the interpolation to the x and y ranges defined above, 
# so that we can evaluate it at any (z, log10(k)) within those ranges
InterpPmm = scale(InterpPmm, (x, y))
# this allows us to extrapolate the interpolation outside of the defined x and y ranges
InterpPmm = Interpolations.extrapolate(InterpPmm, Line())

#Non-linear P(k) - similar to the linear P(k)
y = LinRange(log10(first(k)),log10(last(k)), length(k))
x = LinRange(first(z), last(z), length(z))
InterpPmm_nl = Interpolations.interpolate(log10.(Pknonlin),BSpline(Cubic(Line(OnGrid()))))
InterpPmm_nl = scale(InterpPmm_nl, x, y)
InterpPmm_nl = Interpolations.extrapolate(InterpPmm_nl, Line())

#Callables: these functions take in k, χ1, and χ2, 
# and return the square root of the product of the P(k) evaluated at those points.
power_spectrum(k, χ1, χ2) = @. sqrt(10^InterpPmm(z_of_χ(χ1),log10(k)) * 10^InterpPmm(z_of_χ(χ2),log10(k)))
power_spectrum_nl(k, χ1, χ2) = @. sqrt(10^InterpPmm_nl(z_of_χ(χ1),log10(k)) * 10^InterpPmm_nl(z_of_χ(χ2),log10(k)));

plot(k, power_spectrum.(k, 1000.0, 1000.0), 
     label="Linear P(k)", 
     xscale=:log10, 
     yscale=:log10, 
     title=L"$P(k)$ at $z=0$", titlefontsize=20,
     xlabel=L"$k \; (h/\mathrm{Mpc})$", 
     ylabel=L"$P(k) \; ((\mathrm{Mpc}/h)^3)$", 
     labelfontsize=15)
plot!(k, power_spectrum_nl.(k, 1000.0, 1000.0), 
      label="Non-linear P(k)", 
      xscale=:log10, 
      yscale=:log10,
      size=(800,600),
      legendfontsize=10)


#N5K benchmarks
dtype = Float64

benchmark_gg = npzread("blast_code/data/benchmarks_nl_full_clgg.npz")
benchmark_ll = npzread("blast_code/data/benchmarks_nl_full_clss.npz")
benchmark_gl = npzread("blast_code/data/benchmarks_nl_full_clgs.npz")

#Extracting C_ℓ
gg = dtype.(benchmark_gg["cls"])
ll = dtype.(benchmark_ll["cls"])
gl = dtype.(benchmark_gl["cls"])
ell = dtype.(benchmark_gg["ls"])

#Reshaping them into the same Blast's format
gg_reshaped = zeros( dtype, length(ell), 10, 10)
ll_reshaped = zeros(dtype, length(ell), 5, 5)
gl_reshaped = zeros( dtype, length(ell), 10, 5)
counter = 1
for i in 1:10
    for j in i:10
        gg_reshaped[:,i,j] = gg[counter, :]
        gg_reshaped[:,j,i] = gg_reshaped[:,i,j]
        counter += 1
    end
end
counter = 1
for i in 1:5
    for j in i:5
        ll_reshaped[:,i,j] = ll[counter, :]
        ll_reshaped[:,j,i] = ll_reshaped[:,i,j]
        counter += 1
    end
end
counter = 1
for i in 1:10
    for j in 1:5
        gl_reshaped[:,i,j] = gl[counter, :]
        counter += 1
    end
end
```

### Cell role

This cell reads external data files into memory. Several lines define grids or allocate arrays whose dimensions set the scale of later computations. [file:1]

### Line-by-line

- `#3D matter power spectrum`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `pk_dict = npzread("blast_code/data/pk.npz")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `Pklin = pk_dict["pk_lin"]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `Pknonlin = pk_dict["pk_nl"]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `k = pk_dict["k"]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `z = pk_dict["z"];`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("Dimension of k points: ", size(k)[1])`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("Dimension of Pklin: ", size(Pklin)[1])`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("Dimension of Pknonlin: ", size(Pknonlin)[1])`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `#Interpolating the power spectrum`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `#Linear P(k)`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `# y has the same length as k, and is a logarithmic range from log10(first(k)) to log10(last(k))`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `y = LinRange(log10(first(k)),log10(last(k)), length(k))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `# x has the same length as z, and is a linear range from first(z) to last(z)`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `x = LinRange(first(z), last(z), length(z))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `# this creates an interpolation of the linear P(k) data. `: Comment line. It documents intent but does not execute numerical work. [file:1]
- `# the interpolatio is in log space because P(k) changes by many order or magnitude`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `InterpPmm = Interpolations.interpolate(log10.(Pklin),BSpline(Cubic(Line(OnGrid()))))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `# this scales the interpolation to the x and y ranges defined above, `: Comment line. It documents intent but does not execute numerical work. [file:1]
- `# so that we can evaluate it at any (z, log10(k)) within those ranges`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `InterpPmm = scale(InterpPmm, (x, y))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `# this allows us to extrapolate the interpolation outside of the defined x and y ranges`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `InterpPmm = Interpolations.extrapolate(InterpPmm, Line())`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `#Non-linear P(k) - similar to the linear P(k)`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `y = LinRange(log10(first(k)),log10(last(k)), length(k))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `x = LinRange(first(z), last(z), length(z))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `InterpPmm_nl = Interpolations.interpolate(log10.(Pknonlin),BSpline(Cubic(Line(OnGrid()))))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `InterpPmm_nl = scale(InterpPmm_nl, x, y)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `InterpPmm_nl = Interpolations.extrapolate(InterpPmm_nl, Line())`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `#Callables: these functions take in k, χ1, and χ2, `: Comment line. It documents intent but does not execute numerical work. [file:1]
- `# and return the square root of the product of the P(k) evaluated at those points.`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `power_spectrum(k, χ1, χ2) = @. sqrt(10^InterpPmm(z_of_χ(χ1),log10(k)) * 10^InterpPmm(z_of_χ(χ2),log10(k)))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `power_spectrum_nl(k, χ1, χ2) = @. sqrt(10^InterpPmm_nl(z_of_χ(χ1),log10(k)) * 10^InterpPmm_nl(z_of_χ(χ2),log10(k)));`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot(k, power_spectrum.(k, 1000.0, 1000.0), `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     label="Linear P(k)", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     xscale=:log10, `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     yscale=:log10, `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     title=L"$P(k)$ at $z=0$", titlefontsize=20,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     xlabel=L"$k \; (h/\mathrm{Mpc})$", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     ylabel=L"$P(k) \; ((\mathrm{Mpc}/h)^3)$", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     labelfontsize=15)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot!(k, power_spectrum_nl.(k, 1000.0, 1000.0), `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `      label="Non-linear P(k)", `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `      xscale=:log10, `: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `      yscale=:log10,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `      size=(800,600),`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `      legendfontsize=10)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `#N5K benchmarks`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `dtype = Float64`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `benchmark_gg = npzread("blast_code/data/benchmarks_nl_full_clgg.npz")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `benchmark_ll = npzread("blast_code/data/benchmarks_nl_full_clss.npz")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `benchmark_gl = npzread("blast_code/data/benchmarks_nl_full_clgs.npz")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `#Extracting C_ℓ`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `gg = dtype.(benchmark_gg["cls"])`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `ll = dtype.(benchmark_ll["cls"])`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `gl = dtype.(benchmark_gl["cls"])`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `ell = dtype.(benchmark_gg["ls"])`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `#Reshaping them into the same Blast's format`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `gg_reshaped = zeros( dtype, length(ell), 10, 10)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `ll_reshaped = zeros(dtype, length(ell), 5, 5)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `gl_reshaped = zeros( dtype, length(ell), 10, 5)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `counter = 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `for i in 1:10`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    for j in i:10`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        gg_reshaped[:,i,j] = gg[counter, :]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        gg_reshaped[:,j,i] = gg_reshaped[:,i,j]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        counter += 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `counter = 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `for i in 1:5`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    for j in i:5`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        ll_reshaped[:,i,j] = ll[counter, :]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        ll_reshaped[:,j,i] = ll_reshaped[:,i,j]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        counter += 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `counter = 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `for i in 1:10`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    for j in 1:5`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        gl_reshaped[:,i,j] = gl[counter, :]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        counter += 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 49

```julia
println("Min of k_array: ", minimum(k_array), " || Max of k_array: ", maximum(k_array))
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `println("Min of k_array: ", minimum(k_array), " || Max of k_array: ", maximum(k_array))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 50

```julia
function compute_Sℓ(W_tilde_k1::AbstractArray,
                    W_tilde_k2::AbstractArray,
                    k_array::AbstractVector,
                    ℓ_list::AbstractVector;
                    normalization::Real = 1.0,
                    χ1::Real = 1000.0,
                    χ2::Real = 1000.0,
                    power_spectrum = power_spectrum)

    nk = length(k_array)

    @assert size(W_tilde_k1, 1) == length(ℓ_list)
    @assert size(W_tilde_k2, 1) == length(ℓ_list)
    @assert size(W_tilde_k1, 2) == size(W_tilde_k2, 2)
    @assert size(W_tilde_k1, 3) == size(W_tilde_k2, 3)
    @assert size(W_tilde_k1, 4) == nk
    @assert size(W_tilde_k2, 4) == nk
    @assert isodd(nk) "Simpson integration requires an odd number of k samples."

    w_k = Blast.simpson_weight_array(nk)
    Δk = nk > 1 ? (last(k_array) - first(k_array)) / (nk - 1) : one(eltype(k_array))

    Pk = power_spectrum.(10 .^ k_array, χ1, χ2)
    integrand = normalization .* k_array.^2 .* Pk .* w_k .* Δk

    @tullio Sℓ[li, b, ik1, ik2] := integrand[m] *
        W_tilde_k1[li, b, ik1, m] *
        W_tilde_k2[li, b, ik2, m]

    return Sℓ
end
```

### Cell role

Some lines define reusable functions whose cost is mostly paid when they are evaluated later. [file:1]

### Line-by-line

- `function compute_Sℓ(W_tilde_k1::AbstractArray,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `                    W_tilde_k2::AbstractArray,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `                    k_array::AbstractVector,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `                    ℓ_list::AbstractVector;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `                    normalization::Real = 1.0,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `                    χ1::Real = 1000.0,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `                    χ2::Real = 1000.0,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `                    power_spectrum = power_spectrum)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    nk = length(k_array)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    @assert size(W_tilde_k1, 1) == length(ℓ_list)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    @assert size(W_tilde_k2, 1) == length(ℓ_list)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    @assert size(W_tilde_k1, 2) == size(W_tilde_k2, 2)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    @assert size(W_tilde_k1, 3) == size(W_tilde_k2, 3)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    @assert size(W_tilde_k1, 4) == nk`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    @assert size(W_tilde_k2, 4) == nk`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    @assert isodd(nk) "Simpson integration requires an odd number of k samples."`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    w_k = Blast.simpson_weight_array(nk)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    Δk = nk > 1 ? (last(k_array) - first(k_array)) / (nk - 1) : one(eltype(k_array))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    Pk = power_spectrum.(10 .^ k_array, χ1, χ2)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    integrand = normalization .* k_array.^2 .* Pk .* w_k .* Δk`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    @tullio Sℓ[li, b, ik1, ik2] := integrand[m] *`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        W_tilde_k1[li, b, ik1, m] *`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `        W_tilde_k2[li, b, ik2, m]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    return Sℓ`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 51

```julia
println(length(k_array))
println(isodd(length(k_array)))
nk = 2 * floor(Int, length(k_array) / 2) + 1
k_array = range(first(k_array), last(k_array); length = nk) |> collect
println(size(W_tilde_k1))
println(length(k_array))
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `println(length(k_array))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println(isodd(length(k_array)))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `nk = 2 * floor(Int, length(k_array) / 2) + 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `k_array = range(first(k_array), last(k_array); length = nk) |> collect`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println(size(W_tilde_k1))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println(length(k_array))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 52

```julia
println("min of array x1 = ", minimum(k_array), " || max of array x1 = ", maximum(k_array))
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `println("min of array x1 = ", minimum(k_array), " || max of array x1 = ", maximum(k_array))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 53

```julia
S_l = compute_Sℓ(
    W_tilde_k1,
    W_tilde_k2,
    k_array,
    ℓ;
    normalization = (2/pi),
    χ1 = x_array[1],
    χ2 = x_array[1],
    power_spectrum = power_spectrum_nl
)

println("Size of S_ell = ", size(S_l))
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `S_l = compute_Sℓ(`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    W_tilde_k1,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    W_tilde_k2,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    k_array,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    ℓ;`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    normalization = (2/pi),`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    χ1 = x_array[1],`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    χ2 = x_array[1],`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    power_spectrum = power_spectrum_nl`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("Size of S_ell = ", size(S_l))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 54

```julia
println("min/max chebcoeffKtilde = ", extrema(cheb_coeff_K_tilde))
println("min/max Kvals = ", extrema(K_tilde_vals))
```

### Cell role

The main algorithmic idea here is a spectral decomposition using FFTW-backed real transforms. [file:1]

### Line-by-line

- `println("min/max chebcoeffKtilde = ", extrema(cheb_coeff_K_tilde))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("min/max Kvals = ", extrema(K_tilde_vals))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 55

```julia
println("min/max χ(z)²  = ", extrema(chi_interp.(z_nodes).^2))
println("min/max b(z)   = ", extrema(bias_interp.(z_nodes)))
println("min/max D(z)   = ", extrema(growth_interp.(z_nodes)))
println("min/max 1/H(z) = ", extrema(inv_Hz_interp.(z_nodes)))

# Per ogni bin
for i in 1:n_bins
    println("bin $i → min/max n_i(z) = ", extrema(nz_interp[i].(z_nodes)))
end
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `println("min/max χ(z)²  = ", extrema(chi_interp.(z_nodes).^2))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("min/max b(z)   = ", extrema(bias_interp.(z_nodes)))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("min/max D(z)   = ", extrema(growth_interp.(z_nodes)))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("min/max 1/H(z) = ", extrema(inv_Hz_interp.(z_nodes)))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `# Per ogni bin`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `for i in 1:n_bins`: Starts the outer loop over tomographic bins. [file:1]
- `    println("bin $i → min/max n_i(z) = ", extrema(nz_interp[i].(z_nodes)))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 56

```julia
# Solo il termine n=0 (dovrebbe essere positivo e dominante)
i_bin = 10
i_ell = 5

w0 = cheb_coeff_K[i_bin, 1] .* W_base_k1[i_ell, :, 1, 1]  # solo c_0 * W_base[n=0]
w_full = W_tilde_k1[i_ell, i_bin, :, 1]

println("c_0 = ", cheb_coeff_K[i_bin, 1])
println("W_base[n=0] extrema = ", extrema(W_base_k1[i_ell, :, 1, 1]))
println("w0 (solo termine c_0) extrema = ", extrema(w0))
println("w_full extrema = ", extrema(w_full))

plot(k1_array, w0,    label="solo c₀·W_base[n=0]")
plot!(k1_array, w_full, label="W_tilde completo")
```

### Cell role

The main algorithmic idea here is a spectral decomposition using FFTW-backed real transforms. [file:1]

### Line-by-line

- `# Solo il termine n=0 (dovrebbe essere positivo e dominante)`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `i_bin = 10`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `i_ell = 5`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `w0 = cheb_coeff_K[i_bin, 1] .* W_base_k1[i_ell, :, 1, 1]  # solo c_0 * W_base[n=0]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `w_full = W_tilde_k1[i_ell, i_bin, :, 1]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("c_0 = ", cheb_coeff_K[i_bin, 1])`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("W_base[n=0] extrema = ", extrema(W_base_k1[i_ell, :, 1, 1]))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("w0 (solo termine c_0) extrema = ", extrema(w0))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("w_full extrema = ", extrema(w_full))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot(k1_array, w0,    label="solo c₀·W_base[n=0]")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot!(k1_array, w_full, label="W_tilde completo")`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 57

```julia
# Quante modalità Chebyshev ha W_base?
println("size W_base_k1 = ", size(W_base_k1))   # (nl, nk1, nk, ???)
println("size cheb_coeff_K = ", size(cheb_coeff_K))  # (nbins, ???)

# Devono avere la stessa dimensione lungo n
println("n_cheb W_base = ", size(W_base_k1, 4))
println("n_cheb coeff  = ", size(cheb_coeff_K, 2))
```

### Cell role

The main algorithmic idea here is a spectral decomposition using FFTW-backed real transforms. [file:1]

### Line-by-line

- `# Quante modalità Chebyshev ha W_base?`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `println("size W_base_k1 = ", size(W_base_k1))   # (nl, nk1, nk, ???)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("size cheb_coeff_K = ", size(cheb_coeff_K))  # (nbins, ???)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `# Devono avere la stessa dimensione lungo n`: Comment line. It documents intent but does not execute numerical work. [file:1]
- `println("n_cheb W_base = ", size(W_base_k1, 4))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `println("n_cheb coeff  = ", size(cheb_coeff_K, 2))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 58

```julia
ik1 = 1
ik2 = 1
nl = length(ℓ)

plots = Plots.Plot[]

for i_bin in 1:size(S_l, 2)
    p = plot(ℓ[1:nl], S_l[1:nl, i_bin, ik1, ik2],
             label = "Sℓ",
             title = "bin $i_bin",
             titlefontsize = 20,
             xlabel = L"$\ell$",
             ylabel = L"$S_\ell$",
             labelfontsize = 15)
    push!(plots, p)
end

ncols = ceil(Int, sqrt(size(S_l, 2)))
nrows = ceil(Int, size(S_l, 2) / ncols)

plot(plots..., layout = (nrows, ncols), size = (1400, 1000))
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `ik1 = 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `ik2 = 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `nl = length(ℓ)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plots = Plots.Plot[]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `for i_bin in 1:size(S_l, 2)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    p = plot(ℓ[1:nl], S_l[1:nl, i_bin, ik1, ik2],`: Creates a plot object configured for the kernel display. `p` is a plotting object rather than a scientific array. [file:1]
- `             label = "Sℓ",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `             title = "bin $i_bin",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `             titlefontsize = 20,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `             xlabel = L"$\ell$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `             ylabel = L"$S_\ell$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `             labelfontsize = 15)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    push!(plots, p)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `ncols = ceil(Int, sqrt(size(S_l, 2)))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `nrows = ceil(Int, size(S_l, 2) / ncols)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot(plots..., layout = (nrows, ncols), size = (1400, 1000))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 59

```julia
i_bin = 10
ik1 = 1
ik2 = 1
nl = length(ℓ)

plot(ℓ[1:nl], S_l[1:nl, i_bin, ik1, ik2],
     label = "Sℓ",
     title = L"$S_\ell$",
     titlefontsize = 20,
     xlabel = L"$\ell$",
     ylabel = L"$S_\ell$",
     labelfontsize = 15)
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `i_bin = 10`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `ik1 = 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `ik2 = 1`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `nl = length(ℓ)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot(ℓ[1:nl], S_l[1:nl, i_bin, ik1, ik2],`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     label = "Sℓ",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     title = L"$S_\ell$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     titlefontsize = 20,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     xlabel = L"$\ell$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     ylabel = L"$S_\ell$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     labelfontsize = 15)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 60

```julia
function contract_Sℓ(S_l, proj_k1, proj_k2)
    @assert size(S_l, 3) == length(proj_k1)
    @assert size(S_l, 4) == length(proj_k2)

    @tullio C[li, b] := proj_k1[ik1] * S_l[li, b, ik1, ik2] * proj_k2[ik2]
    return C
end
```

### Cell role

Some lines define reusable functions whose cost is mostly paid when they are evaluated later. [file:1]

### Line-by-line

- `function contract_Sℓ(S_l, proj_k1, proj_k2)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    @assert size(S_l, 3) == length(proj_k1)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    @assert size(S_l, 4) == length(proj_k2)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    @tullio C[li, b] := proj_k1[ik1] * S_l[li, b, ik1, ik2] * proj_k2[ik2]`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `    return C`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `end`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Code cell 61

```julia
i_bin = 10
nl = length(ℓ)
proj_k1 = ones(size(S_l, 3))
proj_k2 = ones(size(S_l, 4))

C_l_like = contract_Sℓ(S_l, proj_k1, proj_k2)
C_l_like = contract_Sℓ(S_l, proj_k1, proj_k2)

plot(ℓ[1:nl], C_l_like[1:nl, i_bin],
     label = "BLAST-like",
     title = L"$C_\ell$-like quantity",
     titlefontsize = 20,
     xlabel = L"$\ell$",
     ylabel = L"$C_\ell$",
     labelfontsize = 15)
```

### Cell role

This cell mainly performs bookkeeping, diagnostics, or parameter setup. [file:1]

### Line-by-line

- `i_bin = 10`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `nl = length(ℓ)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `proj_k1 = ones(size(S_l, 3))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `proj_k2 = ones(size(S_l, 4))`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `C_l_like = contract_Sℓ(S_l, proj_k1, proj_k2)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `C_l_like = contract_Sℓ(S_l, proj_k1, proj_k2)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `plot(ℓ[1:nl], C_l_like[1:nl, i_bin],`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     label = "BLAST-like",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     title = L"$C_\ell$-like quantity",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     titlefontsize = 20,`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     xlabel = L"$\ell$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     ylabel = L"$C_\ell$",`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]
- `     labelfontsize = 15)`: This line belongs to a larger multi-line construction whose object types and dimensions are clarified by the surrounding lines in the same cell. It is retained here without inventing details that are not explicit in the extracted content. [file:1]

## Extraction boundary

The visible notebook content continues after the first Chebyshev decomposition stage, but the remaining cells were not fully available within the extracted attachment window used for this pass. The analysis therefore stops here on purpose so the continuation can be done reliably in a second pass from the next unexplained cell onward. [file:1]
