from scipy.interpolate import interp1d
import numpy as np
import pandas as pd
import os


def get_plus_minus(val, x_array):
    """
    Find the nearest lower and upper values in an array for a given value.
    """
    x_array = np.asarray(x_array)

    lower_val = None
    upper_val = None
    lower_idx = None
    upper_idx = None

    for idx, x in enumerate(x_array):
        if x <= val:
            if lower_val is None or x > lower_val:
                lower_val = x
                lower_idx = idx
        if x >= val:
            if upper_val is None or x < upper_val:
                upper_val = x
                upper_idx = idx

    if lower_val is None:
        lower_val = np.nan
        lower_idx = np.nan
    if upper_val is None:
        upper_val = np.nan
        upper_idx = np.nan

    return (lower_val, lower_idx), (upper_val, upper_idx)


def effi_model(d2, e_d2, npz_filepath, model="interpolate"):
    """
    Calculate efficiency and its uncertainty using an interpolation model.
    """
    data = np.load(npz_filepath)

    x = data["x"]
    y = data["y"]
    y_error_low = data["y_error_low"]
    y_error_high = data["y_error_high"]

    f_linear = interp1d(x, y, kind="linear", fill_value="extrapolate")

    d2 = pd.Series(d2)

    if model == "interpolate":
        effi_values = f_linear(d2)
        effi = pd.DataFrame(effi_values, columns=["estimate"])

        effi_bounds = effi["estimate"].apply(lambda val: get_plus_minus(val, y))
        effi_minus = effi_bounds.apply(lambda t: t[0][0])
        idx_effi_minus = effi_bounds.apply(lambda t: t[0][1])
        effi_plus = effi_bounds.apply(lambda t: t[1][0])
        idx_effi_plus = effi_bounds.apply(lambda t: t[1][1])

        d_bounds = d2.apply(lambda val: get_plus_minus(val, x))
        d_minus = d_bounds.apply(lambda t: t[0][0])
        d_plus = d_bounds.apply(lambda t: t[1][0])

        idx_effi_minus = idx_effi_minus.fillna(-1).astype(int)
        idx_effi_plus = idx_effi_plus.fillna(-1).astype(int)

        e_effi_minus = pd.Series(
            [y_error_low[i] if i >= 0 else np.nan for i in idx_effi_minus]
        )
        e_effi_plus = pd.Series(
            [y_error_high[i] if i >= 0 else np.nan for i in idx_effi_plus]
        )

        effi_plus = effi_plus.reset_index(drop=True)
        effi_minus = effi_minus.reset_index(drop=True)
        d2 = d2.reset_index(drop=True)
        e_effi_minus = e_effi_minus.reset_index(drop=True)
        e_effi_plus = e_effi_plus.reset_index(drop=True)
        d_plus = d_plus.reset_index(drop=True)
        d_minus = d_minus.reset_index(drop=True)

        delta_d = d_plus - d_minus

        e_effi_vals = np.zeros_like(d2, dtype=float)
        valid_mask = delta_d != 0

        de_minus = 1 + d2[valid_mask] * (1 / delta_d[valid_mask])
        de_plus = -d2[valid_mask] * (1 / delta_d[valid_mask])
        dd2 = -(effi_plus[valid_mask] - effi_minus[valid_mask]) / delta_d[valid_mask]
        dd_plus = (
            d2[valid_mask]
            * (effi_plus[valid_mask] - effi_minus[valid_mask])
            / (delta_d[valid_mask] ** 2)
        )
        dd_minus = (
            -d2[valid_mask]
            * (effi_plus[valid_mask] - effi_minus[valid_mask])
            / (delta_d[valid_mask] ** 2)
        )

        e_effi_vals[valid_mask] = np.sqrt(
            (de_minus**2) * (e_effi_minus[valid_mask] ** 2)
            + (de_plus**2) * (e_effi_plus[valid_mask] ** 2)
            + (dd2**2) * d2[valid_mask]
            + (dd_plus**2) * d_plus[valid_mask]
            + (dd_minus**2) * d_minus[valid_mask]
        )

        e_effi_vals[~valid_mask] = np.interp(
            d2[~valid_mask], x, (y_error_low + y_error_high) / 2
        )

        e_effi = pd.DataFrame({"e_estimate": e_effi_vals})
        return effi, e_effi


def effi_estimates(df, ana_kin="xF", model="interpolate", var="D2"):
    """
    Calculate average efficiency per bin using efficiency models.
    """
    if ana_kin == "xT":
        xt_bins = np.array([0.1, 0.13, 0.16, 0.195, 0.24, 0.29, 0.35, 0.45])
    elif ana_kin == "xF":
        xt_bins = np.round(np.arange(0.0, 0.9, 0.05), 2)
    else:
        raise ValueError(f"Kinematic variable '{ana_kin}' not defined with bins.")

    estimates = np.zeros((len(xt_bins) - 1))
    e_estimates = np.zeros((len(xt_bins) - 1))
    avg_eff2 = np.zeros((len(xt_bins) - 1))       # <ε²>
    sq_avg_eff = np.zeros((len(xt_bins) - 1))     # (<ε>)²
    error_avg_eff = np.zeros((len(xt_bins) - 1))  # √[(<ε²> - (<ε>)²)/N]
    d2_rejected = np.zeros((len(xt_bins) - 1))
    n_events = np.zeros((len(xt_bins) - 1), dtype=int)  # number of dimuon events per bin

    for i in range(len(xt_bins) - 1):
        lower_bound = xt_bins[i]
        upper_bound = xt_bins[i + 1]
        npz_file_for_bin = f"interpolation_data_{lower_bound}to{upper_bound}.npz"

        print(
            f"--- Processing bin {i+1}/{len(xt_bins)-1}: ({lower_bound}, {upper_bound}] using '{npz_file_for_bin}' ---"
        )

        if not os.path.exists(npz_file_for_bin):
            print(f"Warning: File not found: {npz_file_for_bin}. Skipping this bin.")
            estimates[i] = np.nan
            e_estimates[i] = np.nan
            avg_eff2[i] = np.nan
            sq_avg_eff[i] = np.nan
            error_avg_eff[i] = np.nan
            d2_rejected[i] = np.nan
            n_events[i] = 0
            continue

        interp_data = np.load(npz_file_for_bin)
        lim_d = np.array([interp_data["x"].min(), interp_data["x"].max()])

        df_cands = df[(df[ana_kin] > lower_bound) & (df[ana_kin] <= upper_bound)]
        df_cands = df_cands[~df_cands[["runID", "eventID"]].duplicated(keep="first")]
        df_cands = df_cands[(df_cands[var] >= 8)]

        if model == "interpolate":
            if len(df_cands) != 0:
                d2_reject_perc = (
                    len(df_cands[(df_cands[var] < lim_d[0]) | (df_cands[var] > lim_d[1])])
                    * 100
                    / len(df_cands)
                )
            else:
                d2_reject_perc = 0

            df_cands = df_cands[
                (df_cands[var] >= lim_d[0]) & (df_cands[var] <= lim_d[1])
            ]
            d2_rejected[i] = d2_reject_perc

        n_events[i] = len(df_cands)

        if len(df_cands) == 0:
            print("No events in this bin after cuts. Moving to the next bin.")
            estimates[i] = 0
            e_estimates[i] = 0
            avg_eff2[i] = 0
            sq_avg_eff[i] = 0
            error_avg_eff[i] = 0
            continue

        estimate, e_estimate = effi_model(
            df_cands[var], np.sqrt(df_cands[var]), npz_file_for_bin, model
        )

        eff_vals = estimate.iloc[:, 0].values

        # Average efficiency
        estimates[i] = np.mean(eff_vals)

        # Average of efficiency squared
        avg_eff2[i] = np.mean(eff_vals**2)

        # Square of average efficiency
        sq_avg_eff[i] = estimates[i] ** 2

        # Error of average efficiency (statistical)
        error_avg_eff[i] = np.sqrt((avg_eff2[i] - sq_avg_eff[i]) / len(df_cands))

        # Existing propagated error
        e_estimates[i] = np.sqrt((e_estimate.iloc[:, 0] ** 2).sum()) / len(df_cands)

    with np.errstate(divide="ignore", invalid="ignore"):
        inv_estimates = 1.0 / estimates
        e_inv_estimates = e_estimates / (estimates**2)

    df_estimates = pd.DataFrame(
        {
            "N (dimuon events)": n_events,
            r"$<\epsilon_t>$": estimates,
            r"$\delta <\epsilon_t>$ (propagated)": e_estimates,
            r"$\delta_{stat} <\epsilon_t>$": error_avg_eff,
            r"$1/<\epsilon_t>$": inv_estimates,
            r"$\delta(1/<\epsilon_t>)$": e_inv_estimates,
            r"$<\epsilon_t^2>$": avg_eff2,
            r"$(<\epsilon_t>)^2$": sq_avg_eff,
        }
    )

    df_estimates.replace([np.inf, -np.inf], np.nan, inplace=True)

    if model == "interpolate":
        df_estimates[r"$\%\ of\ events\ removed$"] = d2_rejected

    index_labels = [f"[{xt_bins[j]}, {xt_bins[j+1]})" for j in range(len(xt_bins) - 1)]
    df_estimates.index = index_labels

    return df_estimates


# --- Main execution block ---
var_kin = "xF"

#input_parquet_file = '/root/github/e906-development/ROOTFiles/Hugo/R008_roadset67_0_2111v42_tmp_noPhys_noOcc.parquet'
input_parquet_file = '/root/github/e906-development/ROOTFiles/Hugo/roadset57_70_R008_2111v42_tmp_noPhys.parquet'
output_csv_file = "average_efficiency_xF_bins.csv"

print(f"Reading dimuon events from: {input_parquet_file}")
df_ld2_c_cuts = pd.read_parquet(input_parquet_file)

# Apply targetPos == 1 cut to only select LH2
df_ld2_c_cuts = df_ld2_c_cuts[df_ld2_c_cuts["targetPos"] == 1]

print(f"\nCalculating average efficiency for '{var_kin}' bins...")
effi_ld2 = effi_estimates(df_ld2_c_cuts, var_kin, "interpolate", var="D2")

effi_ld2.to_csv(output_csv_file)

print("\n--- Calculation Complete ---")
print("Results:")
print(effi_ld2)
print(f"\nResults have been successfully saved to: {output_csv_file}")

