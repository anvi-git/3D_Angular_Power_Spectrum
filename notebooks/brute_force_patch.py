"""Patch script: appends brute-force C_ell cells to beyond_BLAST.ipynb.
Run once from the notebooks/ directory:
    python brute_force_patch.py
Does NOT delete or modify any existing cell.
"""
import json, pathlib

NOTEBOOK = pathlib.Path("beyond_BLAST.ipynb")
assert NOTEBOOK.exists(), f"{NOTEBOOK} not found"

with open(NOTEBOOK, "r", encoding="utf-8") as f:
    nb = json.load(f)

# ── cells to append ───────────────────────────────────────────────────────────
new_cells = [
    {
        "cell_type": "markdown",
        "id": "bf000001",
        "metadata": {},
        "source": [
            "### Brute Force Computation of $C_\\ell$\n",
            "\n",
            "Direct double sum over Chebyshev coefficients and comoving distance grid, ",
            "used as a reference / validation:\n",
            "\n",
            "$$C_\\ell = \\sum_{m,n} c_m \\, c_n \\sum_{i,j} ",
            "\\tilde{W}_{\\ell,i,j,m} \\, \\tilde{W}_{\\ell,i,j,n}$$\n",
            "\n",
            "Computed separately for galaxy\u2013galaxy ($\\beta=2$) ",
            "and shear\u2013shear ($\\beta=-2$) probes."
        ]
    },
    {
        "cell_type": "code",
        "execution_count": None,
        "id": "bf000002",
        "metadata": {},
        "outputs": [],
        "source": [
            "# ===== BRUTE FORCE COMPUTATION =====\n",
            "# C_\u2113(\u2113) = \u03a3_{m,n} c_m * c_n * \u03a3_{i,j} W_tilde[\u2113, i, j, m] * W_tilde[\u2113, i, j, n]\n",
            "# where c_m are the Chebyshev coefficients of the window-function prefactor W(\u03c7).\n",
            "# Serves as a reference / validation against the fast BLAST pipeline.\n",
            "\n",
            "println(\"Starting brute force computation at \", Dates.now())\n",
            "\n",
            "n_ell     = length(\u2113)\n",
            "n_chi_pts = length(z_range)   # 96 comoving-distance points\n",
            "\n",
            "C_ell_brute_gal   = zeros(n_ell)\n",
            "C_ell_brute_shear = zeros(n_ell)\n",
            "\n",
            "@time begin\n",
            "    for i_ell in eachindex(\u2113)\n",
            "        acc_gal   = 0.0\n",
            "        acc_shear = 0.0\n",
            "        for m in 1:n_cheb\n",
            "            for n in 1:n_cheb\n",
            "                s = 0.0\n",
            "                for i in 1:n_chi_pts\n",
            "                    for j in 1:n_chi_pts\n",
            "                        s += W_tilde[i_ell, i, j, m] * W_tilde[i_ell, i, j, n]\n",
            "                    end\n",
            "                end\n",
            "                acc_gal   += cheb_coeff_gal[m]   * cheb_coeff_gal[n]   * s\n",
            "                acc_shear += cheb_coeff_shear[m]  * cheb_coeff_shear[n] * s\n",
            "            end\n",
            "        end\n",
            "        C_ell_brute_gal[i_ell]   = acc_gal\n",
            "        C_ell_brute_shear[i_ell] = acc_shear\n",
            "        println(\"\u2113 = $(\u2113[i_ell]) done at \", Dates.now())\n",
            "    end\n",
            "end\n",
            "\n",
            "println(\"Brute force computation finished at \", Dates.now())\n",
            "\n",
            "npzwrite(joinpath(output_dir, \"C_ell_brute_gal.npy\"),   C_ell_brute_gal)\n",
            "npzwrite(joinpath(output_dir, \"C_ell_brute_shear.npy\"), C_ell_brute_shear)\n",
            "println(\"Saved brute force C\u2113 (galaxy-galaxy) to: \", joinpath(output_dir, \"C_ell_brute_gal.npy\"))\n",
            "println(\"Saved brute force C\u2113 (shear-shear)    to: \", joinpath(output_dir, \"C_ell_brute_shear.npy\"))\n",
            "\n",
            "p_bf = plot(\u2113, C_ell_brute_gal,\n",
            "    label=L\"C_\\\\ell^{gg} \\\\,(brute force)\",\n",
            "    xlabel=L\"\\\\ell\", ylabel=L\"C_\\\\ell\",\n",
            "    title=\"Brute Force Power Spectra\",\n",
            "    xscale=:log10, yscale=:log10, lw=2)\n",
            "plot!(p_bf, \u2113, C_ell_brute_shear,\n",
            "    label=L\"C_\\\\ell^{\\\\kappa\\\\kappa} \\\\,(brute force)\",\n",
            "    lw=2, linestyle=:dash)\n",
            "savefig(p_bf, joinpath(plot_subdir, \"C_ell_brute_force.png\"))\n",
            "println(\"Saved plot to: \", joinpath(plot_subdir, \"C_ell_brute_force.png\"))\n",
            "display(p_bf)"
        ]
    }
]

nb["cells"].extend(new_cells)

with open(NOTEBOOK, "w", encoding="utf-8") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print(f"Done. Appended {len(new_cells)} cells. Total cells: {len(nb['cells'])}")
