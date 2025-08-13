import uproot
import numpy as np
from iminuit import Minuit
import matplotlib.pyplot as plt
import pandas as pd
import os
try:
    import ROOT
    has_pyroot = True
except ImportError:
    has_pyroot = False
    print("Warning: PyROOT not available. ROOT file output will be skipped.")

# Check for hist package
try:
    import hist
except ImportError:
    raise ImportError("The 'hist' package is required for uproot histogram conversion. Install it with:\n    pip install hist")

# Check for iminuit
try:
    import iminuit
except ImportError:
    raise ImportError("The 'iminuit' package is required for fitting. Install it with:\n    pip install iminuit")

# Configuration arrays
roads = ["57", "59", "62", "67", "70", "57_70", "5ea", "5eb", "r6", "5_6", "run5"]
targetName = ["all", "LH2", "flask", "LD2", "None", "Fe", "C", "W"]
liveproton = np.zeros((len(roads), len(targetName)))

# roadset 57
liveproton[0][1] = 3.7927e+16
liveproton[0][2] = 4.23106e+15
liveproton[0][3] = 1.89739e+16
liveproton[0][4] = 6.60844e+15
liveproton[0][5] = 4.14517e+15
liveproton[0][6] = 1.28741e+16
liveproton[0][7] = 4.28783e+15

# roadset 59
liveproton[1][1] = 8.43758e+15
liveproton[1][2] = 9.17579e+14
liveproton[1][3] = 3.83527e+15
liveproton[1][4] = 2.00673e+15
liveproton[1][5] = 9.41273e+14
liveproton[1][6] = 2.83152e+15
liveproton[1][7] = 9.43782e+14

# roadset 62
liveproton[2][1] = 5.50924e+16
liveproton[2][2] = 1.12555e+16
liveproton[2][3] = 2.55071e+16
liveproton[2][4] = 1.23643e+16
liveproton[2][5] = 6.62631e+15
liveproton[2][6] = 1.22319e+16
liveproton[2][7] = 6.72392e+15

# roadset 67
liveproton[3][1] = 1.57319e+17
liveproton[3][2] = 3.57904e+16
liveproton[3][3] = 7.51113e+16
liveproton[3][4] = 3.65750e+16
liveproton[3][5] = 1.71730e+16
liveproton[3][6] = 3.50933e+16
liveproton[3][7] = 1.75003e+16

# roadset 70
liveproton[4][1] = 1.76151e+16
liveproton[4][2] = 3.71074e+15
liveproton[4][3] = 8.54473e+15
liveproton[4][4] = 3.79717e+15
liveproton[4][5] = 1.83760e+15
liveproton[4][6] = 3.66525e+15
liveproton[4][7] = 1.82418e+15

# roadset 57_70
liveproton[5][1] = liveproton[0][1] + liveproton[1][1] + liveproton[2][1] + liveproton[3][1] + liveproton[4][1]
liveproton[5][2] = liveproton[0][2] + liveproton[1][2] + liveproton[2][2] + liveproton[3][2] + liveproton[4][2]
liveproton[5][3] = liveproton[0][3] + liveproton[1][3] + liveproton[2][3] + liveproton[3][3] + liveproton[4][3]
liveproton[5][4] = liveproton[0][4] + liveproton[1][4] + liveproton[2][4] + liveproton[3][4] + liveproton[4][4]
liveproton[5][5] = liveproton[0][5] + liveproton[1][5] + liveproton[2][5] + liveproton[3][5] + liveproton[4][5]
liveproton[5][6] = liveproton[0][6] + liveproton[1][6] + liveproton[2][6] + liveproton[3][6] + liveproton[4][6]
liveproton[5][7] = liveproton[0][7] + liveproton[1][7] + liveproton[2][7] + liveproton[3][7] + liveproton[4][7]

component_base = ["dy", "jpsi", "psip", "mix", "flask"]  # Base component names
colors = ["red", "purple", "orange", "green", "black"]  # Matplotlib colors
nMass_default = 35  # Default number of bins for integration

def fitter_new(rInt, target, nMass=nMass_default):
    # File paths and titles
    s1Name = f"step1/{roads[rInt]}_{targetName[target]}_mass.root"
    outName = f"step2/{roads[rInt]}_{targetName[target]}"
    fitTitle = f"{roads[rInt]} {targetName[target]} mass fit;mass (GeV);count/(0.20GeV)"
    pic = ".pdf"
    flaskPos = 2 if target < 4 else 4
    components = component_base.copy()
    if target >= 4:
        components[4] = "none"
    sFlask = liveproton[rInt][target] / liveproton[rInt][flaskPos]

    # Check if ROOT file exists
    if not os.path.exists(s1Name):
        print(f"Error: ROOT file {s1Name} not found!")
        return

    # Read ROOT file
    try:
        with uproot.open(s1Name) as s1File:
            hData = s1File["h_data"].to_hist()
            hMC = [s1File[f"n_{components[i]}"].to_hist() for i in range(3)]
            wMC = [s1File[f"w_{components[i]}"].to_hist() for i in range(3)]
            hMCw = [s1File[f"h_{components[i]}"].to_hist() for i in range(3)]
            hMix = s1File["h_mix"].to_hist()
            hFlask = s1File[f"h_{components[4]}"].to_hist()
    except KeyError as e:
        print(f"Error: Histogram not found in {s1Name}: {e}")
        return

    # Convert histograms to numpy arrays
    data_counts = hData.values()
    data_bins = hData.axes[0].edges
    if len(data_counts) < nMass:
        print(f"Warning: nMass={nMass} exceeds histogram bins ({len(data_counts)}). Adjusting to {len(data_counts)}.")
        nMass = len(data_counts)
    data_counts = data_counts[:nMass]
    data_errors = np.sqrt(hData.variances())[:nMass]  # Use statistical errors only
    data_errors = np.where(data_errors == 0, 1.0, data_errors)  # Avoid division by zero
    bin_centers = (data_bins[:-1] + data_bins[1:]) / 2
    bin_centers = bin_centers[:nMass]

    # Restrict to mass range 2.0-9.0 GeV
    valid_bins = (bin_centers >= 2.0) & (bin_centers <= 9.0)
    if not np.any(valid_bins):
        print("Error: No bins within mass range 2.0-9.0 GeV!")
        return
    data_counts = data_counts[valid_bins]
    data_errors = data_errors[valid_bins]
    bin_centers = bin_centers[valid_bins]
    data_bins = data_bins[np.concatenate([valid_bins, [valid_bins[-1]]])]
    nMass = len(data_counts)

    # Debug: Print histogram integrals and check validity
    print(f"Data integral (2.0-9.0 GeV): {np.sum(data_counts)}")
    mc_counts = [h.values()[:nMass][valid_bins] for h in hMC]
    wmc_counts = [w.values()[:nMass][valid_bins] for w in wMC]
    hmcw_counts = [h.values()[:nMass][valid_bins] for h in hMCw]
    mix_counts = hMix.values()[:nMass][valid_bins]
    flask_counts = hFlask.values()[:nMass][valid_bins]
    print(f"MC integrals: {[np.sum(h) for h in mc_counts]}")
    print(f"Weight integrals: {[np.sum(w) for w in wmc_counts]}")
    print(f"MC weighted integrals: {[np.sum(h) for h in hmcw_counts]}")
    print(f"Mix integral: {np.sum(mix_counts)}")
    print(f"Flask integral: {np.sum(flask_counts)}")
    if np.any(data_counts < 0) or any(np.any(t < 0) for t in mc_counts + [mix_counts, flask_counts]):
        print("Error: Negative values in data or templates!")
        return
    if np.sum(data_counts) == 0 or any(np.sum(t) == 0 for t in mc_counts + [mix_counts, flask_counts]):
        print("Error: Empty histogram detected!")
        return

    # Debug: Check for large weights
    for i, w in enumerate(wmc_counts):
        max_weight = np.max(np.abs(w))
        if max_weight > 100:
            print(f"Warning: Large weights in {components[i]}: max = {max_weight}")
            print(f"Weight values: {w}")

    # Compute integrals and errors
    nData = np.sum(data_counts)
    eData = np.sqrt(np.sum(data_errors**2))
    n0 = [np.sum(hmcw_counts[i]) for i in range(3)] + [np.sum(mix_counts), np.sum(flask_counts)]
    e0 = [np.sqrt(np.sum(np.where(h > 0, h, 1))) for h in hmcw_counts] + [
        np.sqrt(np.sum(np.where(mix_counts > 0, mix_counts, 1))),
        np.sqrt(np.sum(np.where(flask_counts > 0, mix_counts, 1)))
    ]

    # Normalize templates to data integral with safety check
    templates = [mc_counts[i] * wmc_counts[i] for i in range(3)] + [mix_counts, flask_counts]
    templates = [t * nData / np.sum(t) if np.sum(t) > 0 else np.zeros_like(t) for t in templates]

    # Compute fixed flask fraction
    fFlask = sFlask * n0[4] / nData
    print(f"Fixed flask fraction: {fFlask}")

    # Define fitting methods
    def chi2_likelihood_ratio(f0, f1, f2):
        f3 = 1.0 - (f0 + f1 + f2 + fFlask)
        fractions = [f0, f1, f2, f3, fFlask]
        total = sum(f * t for f, t in zip(fractions, templates))
        chi2_terms = np.zeros_like(data_counts)
        valid_bins = (data_counts > 0) & (total > 0)
        chi2_terms[valid_bins] = 2 * (total[valid_bins] - data_counts[valid_bins] + 
                                      data_counts[valid_bins] * np.log(data_counts[valid_bins] / total[valid_bins]))
        chi2_terms[~valid_bins & (data_counts == 0)] = 2 * total[~valid_bins & (data_counts == 0)]
        chi2_terms[~valid_bins & (total <= 0)] = 0
        chi2 = np.sum(chi2_terms)
        if not np.isfinite(chi2):
            print(f"Warning: Non-finite chi2 (LR)! Fractions: {fractions}")
        return chi2

    def chi2_standard(f0, f1, f2):
        f3 = 1.0 - (f0 + f1 + f2 + fFlask)
        fractions = [f0, f1, f2, f3, fFlask]
        
        # Compute Ntotal_i = sum(Ncomponent_i)
        N_components = [
            templates[i] * f for i, f in enumerate(fractions)
        ]
        N_total = sum(N_components)
        
        # Compute error_components_i = sum(MCWeight_i^2) + (MixScale*NMix_i)^2 * (errorMixScale^2/MixScale^2 + 1/NMix_i)
        #                           + (FlaskScale*NFlask_i)^2 * (errorFlaskScale^2/FlaskScale^2 + 1/NFlask_i)
        error_components = np.zeros_like(data_counts)
        for i in range(len(data_counts)):
            # MC components (sum of squares of weighted MC counts)
            mc_error_term = sum((wmc_counts[j][i] * fractions[j])**2 for j in range(3))
            
            # Mix component: (f3 * mix_counts[i])^2 * (f3_error^2 / f3^2 + 1 / mix_counts[i])
            mix_scale = f3
            mix_contrib = (mix_scale * mix_counts[i])**2
            mix_error_term = mix_contrib * (f3_error**2 / (mix_scale**2 + 1e-10) + 1.0 / (mix_counts[i] + 1e-10))
            
            # Flask component: (fFlask * flask_counts[i])^2 * (0^2 / fFlask^2 + 1 / flask_counts[i])
            flask_scale = fFlask
            flask_contrib = (flask_scale * flask_counts[i])**2
            flask_error_term = flask_contrib * (1.0 / (flask_counts[i] + 1e-10))
            
            error_components[i] = mc_error_term + mix_error_term + flask_error_term
        
        # Total delta_i^2 = data_errors^2 + error_components
        delta_i_squared = data_errors**2 + error_components
        delta_i = np.sqrt(np.where(delta_i_squared > 0, delta_i_squared, 1.0))
        
        # Compute d_i = (Ndata_i - Ntotal_i) / delta_i
        d_i = (data_counts - N_total) / delta_i
        
        # Chi-squared = sum(d_i^2)
        chi2 = np.sum(d_i**2)
        if not np.isfinite(chi2):
            print(f"Warning: Non-finite chi2 (Standard)! Fractions: {fractions}")
        return chi2

    def nll(f0, f1, f2):
        f3 = 1.0 - (f0 + f1 + f2 + fFlask)
        fractions = [f0, f1, f2, f3, fFlask]
        total = sum(f * t for f, t in zip(fractions, templates))
        nll_terms = np.zeros_like(data_counts)
        valid_bins = (total > 0)
        nll_terms[valid_bins] = total[valid_bins] - data_counts[valid_bins] * np.log(total[valid_bins])
        nll = np.sum(nll_terms)
        if not np.isfinite(nll):
            print(f"Warning: Non-finite NLL! Fractions: {fractions}")
        return 2 * nll  # Return -2*log(L) for Minuit

    # Initial guess
    initial_guess = [n0[i] / nData for i in range(4)] + [fFlask]
    sum_initial = sum(initial_guess[:4])
    if sum_initial > 0:
        initial_guess[:4] = [f / sum_initial * (1 - fFlask) for f in initial_guess[:4]]
    else:
        initial_guess[:4] = [0.25 * (1 - fFlask)] * 4
    print(f"Initial guess: {initial_guess}, Sum: {sum(initial_guess)}")

    # Perform fits for each method
    fit_methods = [
        ("Likelihood Ratio Chi2", chi2_likelihood_ratio, Minuit.LEAST_SQUARES),
        ("Standard Chi2", chi2_standard, Minuit.LEAST_SQUARES),
        ("Binned Maximum Likelihood", nll, Minuit.LIKELIHOOD)
    ]
    
    results = []
    for method_name, func, errordef in fit_methods:
        m = Minuit(func, f0=initial_guess[0], f1=initial_guess[1], f2=initial_guess[2])
        m.limits = [(0, 1 - fFlask) for _ in range(3)]
        m.errordef = errordef
        m.migrad()
        
        if not m.valid:
            print(f"Fit failed for {method_name}! Reason: {m.fmin}")
            continue

        fractions = list(m.values) + [1.0 - sum(m.values) - fFlask, fFlask]
        hess_inv = m.covariance
        if hess_inv is not None:
            var_f3 = sum(hess_inv[i, i] for i in range(3)) + 2 * sum(hess_inv[i, j] for i in range(3) for j in range(i + 1, 3))
            f3_error = np.sqrt(var_f3) if var_f3 > 0 else 0
        else:
            f3_error = 0
            print(f"Warning: Hessian inverse not available for {method_name}, setting f3 error to 0")
        fraction_errors = list(m.errors) + [f3_error, 0.0]
        
        chi2_value = func(*m.values)
        n_bins_used = np.sum(data_counts > 0)
        n_free_par = len([f for i, f in enumerate(fractions[:4]) if not np.isclose(f, 0, atol=1e-6) and 
                         not np.isclose(f, 1 - fFlask, atol=1e-6)])
        manual_ndf = n_bins_used - n_free_par
        chi2_ndf = chi2_value / manual_ndf if manual_ndf > 0 else 0
        print(f"{method_name} Chi2/NDF: {chi2_value:.1f}/{manual_ndf} = {chi2_ndf:.2f}")
        
        results.append({
            "name": method_name,
            "fractions": fractions,
            "fraction_errors": fraction_errors,
            "chi2": chi2_value,
            "ndf": manual_ndf,
            "chi2_ndf": chi2_ndf
        })

    # Save results for the primary method (Likelihood Ratio Chi2)
    primary_result = results[0] if results else None
    if primary_result:
        n1 = np.array(primary_result["fractions"]) * nData
        en1 = np.array(primary_result["fraction_errors"]) * nData
        s0 = n1 / np.array(n0)
        es0 = en1 / np.array(n0)
        
        with open(f"{outName}.csv", "w") as outFile:
            outFile.write("comp,scale,eScale\n")
            for i in range(5):
                outFile.write(f"{components[i]},{s0[i]:.6f},{es0[i]:.6f}\n")

    # Plotting for each method
    bins = data_bins
    for result in results:
        plt.figure(figsize=(10, 6))
        plt.xlim(2.0, 9.0)
        plt.ylim(5e-1, 60000)
        
        # Data
        plt.errorbar(bin_centers, data_counts, yerr=data_errors, fmt="o", label="data", color="black", markersize=5)
        
        # Total fit
        total_fit = sum(f * t for f, t in zip(result["fractions"], templates))
        plt.step(bins[:-1], total_fit, where="post", linewidth=2, label="Total", color="blue")
        
        # Individual templates
        mc_new = []
        n1 = np.array(result["fractions"]) * nData
        for i in range(5):
            temp = templates[i].copy()
            temp *= n1[i] / np.sum(temp) if np.sum(temp) > 0 else 0
            mc_new.append(temp)
            plt.step(bins[:-1], temp, where="post", linewidth=2, label=components[i], color=colors[i])
        
        plt.title(f"{fitTitle} - {result['name']}")
        plt.xlabel("mass (GeV)")
        plt.ylabel("count/(0.20GeV)")
        plt.legend(loc="upper right", bbox_to_anchor=(1, 1))
        
        text = f"Total = {nData:.0f}\n"
        for i in range(5):
            text += f"{components[i]} = {n1[i]:.0f}±{result['fraction_errors'][i] * nData:.0f}\n"
        text += f"χ²/NDF = {result['chi2']:.1f}/{result['ndf']} = {result['chi2_ndf']:.2f}"
        plt.text(0.8, 0.4, text, transform=plt.gca().transAxes, verticalalignment="bottom", 
                 bbox=dict(boxstyle="square", facecolor="white", edgecolor="black"))
        
        # Save linear and log plots
        plt.savefig(f"{outName}_{result['name'].replace(' ', '_')}{pic}")
        plt.yscale("log")
        plt.xlim(2.0, 9.0)
        plt.savefig(f"{outName}_{result['name'].replace(' ', '_')}_log{pic}")
        plt.close()

    # Save to ROOT file if PyROOT is available (only for primary method)
    if has_pyroot and primary_result:
        with ROOT.TFile(f"{outName}.root", "RECREATE") as outGr:
            h_data = ROOT.TH1D("data", fitTitle, len(data_bins)-1, data_bins)
            for i, (count, error) in enumerate(zip(data_counts, data_errors), 1):
                h_data.SetBinContent(i, count)
                h_data.SetBinError(i, error)
            h_data.Write()
            
            n1 = np.array(primary_result["fractions"]) * nData
            for i in range(5):
                h_mc = ROOT.TH1D(components[i], components[i], len(data_bins)-1, data_bins)
                temp = templates[i].copy()
                temp *= n1[i] / np.sum(temp) if np.sum(temp) > 0 else 0
                for j, count in enumerate(temp, 1):
                    h_mc.SetBinContent(j, count)
                h_mc.Write()

            # Correlation matrix
            if results[0]["name"] == "Likelihood Ratio Chi2" and hess_inv is not None:
                cor = ROOT.TMatrixD(5, 5)
                for i in range(3):
                    for j in range(3):
                        cor[i][j] = hess_inv[i, j] / np.sqrt(hess_inv[i, i] * hess_inv[j, j]) if hess_inv[i, i] * hess_inv[j, j] > 0 else 0
                for j in range(3):
                    cov_f3_fj = -sum(hess_inv[i, j] for i in range(3))
                    cor[3][j] = cor[j][3] = cov_f3_fj / np.sqrt(var_f3 * hess_inv[j, j]) if var_f3 * hess_inv[j, j] > 0 else 0
                cor[3][3] = 1.0
                cor[4][4] = 1.0
                cor.Write("cor_mat")
            else:
                print("Warning: Skipping correlation matrix output due to unavailable Hessian inverse.")

if __name__ == "__main__":
    fitter_new(rInt=3, target=1)