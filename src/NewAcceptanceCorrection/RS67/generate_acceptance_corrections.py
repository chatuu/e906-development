import uproot
import awkward as ak
import numpy as np
import os
import ROOT
from array import array

# ==============================================================================
# Helper Function for Loading Nested Tree Branches
# ==============================================================================
def load_and_format_tree(file_path, tree_name):
    """
    Reads the 'event' branch variables from the specified TTree and 
    zips them into a flat awkward array structure for easier processing.
    Includes fallbacks for different ROOT serialization formats (split vs. unsplit).
    """
    print(f"  -> Loading {tree_name} from {file_path} ...")
    tree = uproot.open(file_path)[tree_name]
    keys = tree.keys()
    
    if "mass" in keys and "xF" in keys:
        raw = tree.arrays(["mass", "xF", "pT", "weight"])
        return ak.zip({"mass": raw["mass"], "xF": raw["xF"], "pT": raw["pT"], "weight": raw["weight"]})
        
    elif "event" in keys:
        raw = tree.arrays(["event"])
        return ak.zip({"mass": raw["event"]["mass"], "xF": raw["event"]["xF"], "pT": raw["event"]["pT"], "weight": raw["event"]["weight"]})
        
    elif "event/mass" in keys:
        raw = tree.arrays(["event/mass", "event/xF", "event/pT", "event/weight"])
        return ak.zip({"mass": raw["event/mass"], "xF": raw["event/xF"], "pT": raw["event/pT"], "weight": raw["event/weight"]})
        
    elif "event.mass" in keys:
        raw = tree.arrays(["event.mass", "event.xF", "event.pT", "event.weight"])
        return ak.zip({"mass": raw["event.mass"], "xF": raw["event.xF"], "pT": raw["event.pT"], "weight": raw["event.weight"]})
        
    else:
        raise KeyError(f"Could not locate mass, xF, pT, and weight branches. Keys found in {tree_name}: {keys}")

# ==============================================================================
# Numpy & ROOT Helper Functions
# ==============================================================================
def get_weighted_histogram(data, weights, bins):
    data = np.asarray(data, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if len(data) == 0:
        return np.zeros(len(bins)-1), np.zeros(len(bins)-1)
    hist, _ = np.histogram(data, bins=bins, weights=weights)
    sumw2, _ = np.histogram(data, bins=bins, weights=weights**2)
    return hist, np.sqrt(sumw2)

def calc_binomial_errors(w_sum, p, w_sqr):
    valid = (w_sum > 0) & (p >= 0) & (p <= 1)
    errors = np.zeros_like(w_sum)
    errors[valid] = (w_sqr[valid] / w_sum[valid]) * np.sqrt(p[valid] * (1 - p[valid]))
    return errors

def calc_ratio_errors_independent(val_num, val_den, err_num, err_den):
    valid = (val_den > 0) & (val_num > 0)
    ratio = np.zeros_like(val_num)
    ratio_err = np.zeros_like(err_num)
    ratio[valid] = val_num[valid] / val_den[valid]
    ratio_err[valid] = ratio[valid] * np.sqrt((err_num[valid]/val_num[valid])**2 + (err_den[valid]/val_den[valid])**2)
    return ratio, ratio_err

def make_th1f(name, title, bins_arr, contents, errors, directory=None):
    h = ROOT.TH1F(name, title, len(bins_arr)-1, bins_arr)
    h.Sumw2()
    if directory:
        h.SetDirectory(directory)
    for i in range(len(contents)):
        h.SetBinContent(i+1, float(contents[i]))
        h.SetBinError(i+1, float(errors[i]))
    return h

def format_hist(h, x_title, y_title, color):
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.GetXaxis().SetTitle(x_title)
    h.GetXaxis().CenterTitle(True)
    h.GetYaxis().SetTitle(y_title)
    h.GetYaxis().CenterTitle(True)

def get_hist_max_with_error(h):
    """Calculates the absolute maximum of a TH1 including error bars."""
    max_val = 0.0
    for i in range(1, h.GetNbinsX() + 1):
        val = h.GetBinContent(i) + h.GetBinError(i)
        if val > max_val:
            max_val = val
    return max_val

def get_hist_min_with_error(h):
    """Calculates the absolute minimum of a TH1 including error bars, ignoring 0 bins."""
    min_val = 1e9
    valid = False
    for i in range(1, h.GetNbinsX() + 1):
        if h.GetBinContent(i) > 0: # Ignore empty bins
            val = h.GetBinContent(i) - h.GetBinError(i)
            if val < min_val:
                min_val = val
            valid = True
    return min_val if valid else 0.0

def add_fit_and_band(h_ratio, pad):
    pad.cd()
    h_ratio.SetStats(0) 
    h_ratio.Fit("pol0", "QS")
    pad.Update()
    
    fit_func = h_ratio.GetFunction("pol0")
    if fit_func:
        fit_func.SetLineColor(ROOT.kRed)
        fit_func.SetLineWidth(2)
        p0 = fit_func.GetParameter(0)
        p0_err = fit_func.GetParError(0)
        x_min = h_ratio.GetXaxis().GetXmin()
        x_max = h_ratio.GetXaxis().GetXmax()
        
        error_box = ROOT.TBox(x_min, p0 - p0_err, x_max, p0 + p0_err)
        error_box.SetFillColorAlpha(ROOT.kRed, 0.3) 
        error_box.Draw("SAME")
        pad._error_band = error_box 
        
        fit_func.Draw("SAME")
        h_ratio.Draw("E1 SAME")
        
        x_pos = x_min + (x_max - x_min) * 0.02 
        y_range = h_ratio.GetMaximum() - h_ratio.GetMinimum()
        y_pos = p0 + (y_range * 0.05) 
        
        latex = ROOT.TLatex()
        label_size = h_ratio.GetYaxis().GetLabelSize()
        if label_size == 0: label_size = 0.05 
        latex.SetTextSize(label_size * 0.9)
        latex.SetTextColor(ROOT.kRed)
        latex.SetTextAlign(11) 
        latex.DrawLatex(x_pos, y_pos, f"Fit: {p0:.3f} #pm {p0_err:.3f}")
        
        line = ROOT.TLine(x_min, 1.0, x_max, 1.0)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kGray+2)
        line.Draw("SAME")
        pad._fit_line = line

# ==============================================================================
# Plotting Generation Logic: General Sliced Acceptances
# ==============================================================================
def process_acceptance_sliced(x_var, x_edges, slice_var, slice_edges, t_lh2_th, t_lh2_ac, t_ld2_th, t_ld2_ac, out_file, custom_name=None):
    """Generates acceptance plots for x_var by slicing along slice_var."""
    n_bins = len(slice_edges) - 1
    x_edges_root = array('d', x_edges)
    
    x_titles = {"mass": "Mass [GeV]", "xF": "x_{F}", "pT": "p_{T} [GeV/c]", "pT2": "p_{T}^{2} [(GeV/c)^{2}]"}
    x_title = x_titles.get(x_var, x_var)

    for i in range(n_bins):
        val_low = slice_edges[i]
        val_high = slice_edges[i+1]
        
        if slice_var == "xF": slice_title = f"{val_low:.2f} #leq x_{{F}} < {val_high:.2f}"
        elif slice_var == "pT": slice_title = f"{val_low:.2f} #leq p_{{T}} < {val_high:.2f}"
        elif slice_var == "pT2": slice_title = f"{val_low:.2f} #leq p_{{T}}^{{2}} < {val_high:.2f}"
        else: slice_title = f"{val_low:.2f} #leq {slice_var} < {val_high:.2f}"
            
        title = f"{slice_title}; {x_title}"
        
        if custom_name:
            binName = f"{custom_name}_bin{i}"
        else:
            binName = f"{x_var}_sliced_by_{slice_var}_bin{i}"
            
        print(f"  -> Processing Sliced Bin: {binName}")

        out_dir = out_file.GetDirectory(binName)
        if not out_dir: out_dir = out_file.mkdir(binName)
        out_dir.cd()

        m_lh2_th = (t_lh2_th[slice_var] >= val_low) & (t_lh2_th[slice_var] < val_high)
        m_ld2_th = (t_ld2_th[slice_var] >= val_low) & (t_ld2_th[slice_var] < val_high)
        m_lh2_ac = (t_lh2_ac[slice_var] >= val_low) & (t_lh2_ac[slice_var] < val_high)
        m_ld2_ac = (t_ld2_ac[slice_var] >= val_low) & (t_ld2_ac[slice_var] < val_high)

        lh2_th_h, lh2_th_err = get_weighted_histogram(t_lh2_th[x_var][m_lh2_th], t_lh2_th.weight[m_lh2_th], x_edges)
        lh2_ac_h, lh2_ac_err = get_weighted_histogram(t_lh2_ac[x_var][m_lh2_ac], t_lh2_ac.weight[m_lh2_ac], x_edges)
        lh2_ratio = np.divide(lh2_ac_h, lh2_th_h, out=np.zeros_like(lh2_ac_h), where=lh2_th_h != 0)
        lh2_ratio_err = calc_binomial_errors(lh2_th_h, lh2_ratio, lh2_th_err)

        ld2_th_h, ld2_th_err = get_weighted_histogram(t_ld2_th[x_var][m_ld2_th], t_ld2_th.weight[m_ld2_th], x_edges)
        ld2_ac_h, ld2_ac_err = get_weighted_histogram(t_ld2_ac[x_var][m_ld2_ac], t_ld2_ac.weight[m_ld2_ac], x_edges)
        ld2_ratio = np.divide(ld2_ac_h, ld2_th_h, out=np.zeros_like(ld2_ac_h), where=ld2_th_h != 0)
        ld2_ratio_err = calc_binomial_errors(ld2_th_h, ld2_ratio, ld2_th_err)

        combine_ratio = 0.5 * (lh2_ratio + ld2_ratio)
        combine_ratio_err = 0.5 * np.sqrt(lh2_ratio_err**2 + ld2_ratio_err**2)

        valid_dr = (ld2_ratio > 0) & (lh2_ratio > 0)
        dr = np.zeros_like(lh2_ratio)
        dr_err = np.zeros_like(lh2_ratio)
        dr[valid_dr] = lh2_ratio[valid_dr] / ld2_ratio[valid_dr]
        dr_err[valid_dr] = dr[valid_dr] * np.sqrt((lh2_ratio_err[valid_dr] / lh2_ratio[valid_dr])**2 + (ld2_ratio_err[valid_dr] / ld2_ratio[valid_dr])**2)

        h_LH2_thrown = make_th1f(f"h_LH2_thrown_{binName}", title, x_edges_root, lh2_th_h, lh2_th_err, out_dir)
        h_LH2_accept = make_th1f(f"h_LH2_accept_{binName}", title, x_edges_root, lh2_ac_h, lh2_ac_err, out_dir)
        h_ratio_LH2 = make_th1f(f"h_ratio_LH2_{binName}", title, x_edges_root, lh2_ratio, lh2_ratio_err, out_dir)
        h_LD2_thrown = make_th1f(f"h_LD2_thrown_{binName}", title, x_edges_root, ld2_th_h, ld2_th_err, out_dir)
        h_LD2_accept = make_th1f(f"h_LD2_accept_{binName}", title, x_edges_root, ld2_ac_h, ld2_ac_err, out_dir)
        h_ratio_LD2 = make_th1f(f"h_ratio_LD2_{binName}", title, x_edges_root, ld2_ratio, ld2_ratio_err, out_dir)
        h_ratio_combine = make_th1f(f"h_ratio_combine_{binName}", title, x_edges_root, combine_ratio, combine_ratio_err, out_dir)
        h_ratio_acceptance = make_th1f(f"h_ratio_acceptance_{binName}", f"LH2/LD2 {title}", x_edges_root, dr, dr_err, out_dir)

        format_hist(h_ratio_LH2, x_title, "Acceptance", ROOT.kBlue)
        format_hist(h_ratio_LD2, x_title, "Acceptance", ROOT.kRed)
        format_hist(h_ratio_combine, x_title, "Acceptance", ROOT.kBlack)
        format_hist(h_ratio_acceptance, x_title, "LH2 / LD2 Ratio", ROOT.kBlack)

        # Intelligent Scaling for Acceptances (Top Plots)
        max_acc = max(get_hist_max_with_error(h_ratio_LH2), get_hist_max_with_error(h_ratio_LD2), get_hist_max_with_error(h_ratio_combine))
        if max_acc > 0:
            pad_acc = max_acc * 1.3 # 30% padding for legend
            h_ratio_LH2.SetMaximum(pad_acc)
            h_ratio_LD2.SetMaximum(pad_acc)
            h_ratio_combine.SetMaximum(pad_acc)
        h_ratio_LH2.SetMinimum(0)
        h_ratio_LD2.SetMinimum(0)
        h_ratio_combine.SetMinimum(0)
        
        # Intelligent Scaling for Ratios (Bottom Plots)
        max_ratio = get_hist_max_with_error(h_ratio_acceptance)
        min_ratio = get_hist_min_with_error(h_ratio_acceptance)
        if max_ratio > 0:
            padding = (max_ratio - min_ratio) * 0.2
            if padding == 0: padding = 0.2
            h_ratio_acceptance.SetMaximum(max_ratio + padding)
            h_ratio_acceptance.SetMinimum(max(0.0, min_ratio - padding))

        make_canvas = lambda cname, ctitle: ROOT.TCanvas(cname, ctitle, 800, 600)
        def setup_canvas(c): c.SetTickx(1); c.SetTicky(1)

        c_lh2 = make_canvas(f"c_lh2_{binName}", "LH2 Acceptance"); setup_canvas(c_lh2); h_ratio_LH2.Draw("E1"); c_lh2.Write(); c_lh2.SaveAs(f"acceptance_LH2_{binName}.pdf")
        c_ld2 = make_canvas(f"c_ld2_{binName}", "LD2 Acceptance"); setup_canvas(c_ld2); h_ratio_LD2.Draw("E1"); c_ld2.Write(); c_ld2.SaveAs(f"acceptance_LD2_{binName}.pdf")
        c_comb = make_canvas(f"c_comb_{binName}", "Combined Acceptance"); setup_canvas(c_comb); h_ratio_combine.Draw("E1"); c_comb.Write(); c_comb.SaveAs(f"acceptance_combine_{binName}.pdf")

        c_overlay = make_canvas(f"c_overlay_{binName}", "Acceptances Overlay")
        setup_canvas(c_overlay)
        h_ratio_LH2.Draw("E1")
        h_ratio_LD2.Draw("E1 SAME")
        h_ratio_combine.Draw("E1 SAME")

        leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
        leg.SetBorderSize(0); leg.SetFillColor(ROOT.kWhite)
        leg.AddEntry(h_ratio_LD2, "LD2", "l")
        leg.AddEntry(h_ratio_LH2, "LH2", "l")
        leg.AddEntry(h_ratio_combine, "LH2+LD2", "l")
        leg.Draw()

        c_overlay.Write(); c_overlay.SaveAs(f"acceptance_overlay_{binName}.pdf")

        c_dr = make_canvas(f"c_dr_{binName}", "Double Ratio")
        setup_canvas(c_dr); h_ratio_acceptance.Draw("E1")
        
        add_fit_and_band(h_ratio_acceptance, c_dr)
        c_dr.Write(); c_dr.SaveAs(f"acceptance_ratio_{binName}.pdf")

        h_LH2_thrown.Write(); h_LH2_accept.Write(); h_ratio_LH2.Write()
        h_LD2_thrown.Write(); h_LD2_accept.Write(); h_ratio_LD2.Write()
        h_ratio_combine.Write(); h_ratio_acceptance.Write()

# ==============================================================================
# Fully Integrated Plotting Logic
# ==============================================================================
def process_integrated_1D(var_name, var_edges, t_lh2_th, t_lh2_ac, t_ld2_th, t_ld2_ac, out_file):
    title = f"Integrated Acceptance vs {var_name}"
    binName = f"Integrated_{var_name}"
    print(f"  -> Processing Integrated 1D: {binName}")

    out_dir = out_file.mkdir(binName)
    out_dir.cd()
    var_edges_root = array('d', var_edges)

    x_titles = {"mass": "Mass [GeV]", "xF": "x_{F}", "pT": "p_{T} [GeV/c]", "pT2": "p_{T}^{2} [(GeV/c)^{2}]"}
    x_title = x_titles.get(var_name, var_name)

    lh2_th_h, lh2_th_err = get_weighted_histogram(t_lh2_th[var_name], t_lh2_th.weight, var_edges)
    lh2_ac_h, lh2_ac_err = get_weighted_histogram(t_lh2_ac[var_name], t_lh2_ac.weight, var_edges)
    lh2_ratio = np.divide(lh2_ac_h, lh2_th_h, out=np.zeros_like(lh2_ac_h), where=lh2_th_h != 0)
    lh2_ratio_err = calc_binomial_errors(lh2_th_h, lh2_ratio, lh2_th_err)

    ld2_th_h, ld2_th_err = get_weighted_histogram(t_ld2_th[var_name], t_ld2_th.weight, var_edges)
    ld2_ac_h, ld2_ac_err = get_weighted_histogram(t_ld2_ac[var_name], t_ld2_ac.weight, var_edges)
    ld2_ratio = np.divide(ld2_ac_h, ld2_th_h, out=np.zeros_like(ld2_ac_h), where=ld2_th_h != 0)
    ld2_ratio_err = calc_binomial_errors(ld2_th_h, ld2_ratio, ld2_th_err)

    combine_ratio = 0.5 * (lh2_ratio + ld2_ratio)
    combine_ratio_err = 0.5 * np.sqrt(lh2_ratio_err**2 + ld2_ratio_err**2)

    valid_dr = (ld2_ratio > 0) & (lh2_ratio > 0)
    dr = np.zeros_like(lh2_ratio)
    dr_err = np.zeros_like(lh2_ratio)
    dr[valid_dr] = lh2_ratio[valid_dr] / ld2_ratio[valid_dr]
    dr_err[valid_dr] = dr[valid_dr] * np.sqrt((lh2_ratio_err[valid_dr] / lh2_ratio[valid_dr])**2 + (ld2_ratio_err[valid_dr] / ld2_ratio[valid_dr])**2)

    h_LH2_thrown = make_th1f(f"h_LH2_thrown_{binName}", title, var_edges_root, lh2_th_h, lh2_th_err, out_dir)
    h_LH2_accept = make_th1f(f"h_LH2_accept_{binName}", title, var_edges_root, lh2_ac_h, lh2_ac_err, out_dir)
    h_ratio_LH2 = make_th1f(f"h_ratio_LH2_{binName}", title, var_edges_root, lh2_ratio, lh2_ratio_err, out_dir)
    h_LD2_thrown = make_th1f(f"h_LD2_thrown_{binName}", title, var_edges_root, ld2_th_h, ld2_th_err, out_dir)
    h_LD2_accept = make_th1f(f"h_LD2_accept_{binName}", title, var_edges_root, ld2_ac_h, ld2_ac_err, out_dir)
    h_ratio_LD2 = make_th1f(f"h_ratio_LD2_{binName}", title, var_edges_root, ld2_ratio, ld2_ratio_err, out_dir)
    h_ratio_combine = make_th1f(f"h_ratio_combine_{binName}", title, var_edges_root, combine_ratio, combine_ratio_err, out_dir)
    h_ratio_acceptance = make_th1f(f"h_ratio_acceptance_{binName}", f"LH2/LD2 {title}", var_edges_root, dr, dr_err, out_dir)

    format_hist(h_ratio_LH2, x_title, "Acceptance", ROOT.kBlue)
    format_hist(h_ratio_LD2, x_title, "Acceptance", ROOT.kRed)
    format_hist(h_ratio_combine, x_title, "Acceptance", ROOT.kBlack)
    format_hist(h_ratio_acceptance, x_title, "LH2 / LD2 Ratio", ROOT.kBlack)

    # Intelligent Scaling for Acceptances (Top Plots)
    max_acc = max(get_hist_max_with_error(h_ratio_LH2), get_hist_max_with_error(h_ratio_LD2), get_hist_max_with_error(h_ratio_combine))
    if max_acc > 0:
        pad_acc = max_acc * 1.3 # 30% padding for legend
        h_ratio_LH2.SetMaximum(pad_acc)
        h_ratio_LD2.SetMaximum(pad_acc)
        h_ratio_combine.SetMaximum(pad_acc)
    h_ratio_LH2.SetMinimum(0)
    h_ratio_LD2.SetMinimum(0)
    h_ratio_combine.SetMinimum(0)
    
    # Intelligent Scaling for Ratios (Bottom Plots)
    max_ratio = get_hist_max_with_error(h_ratio_acceptance)
    min_ratio = get_hist_min_with_error(h_ratio_acceptance)
    if max_ratio > 0:
        padding = (max_ratio - min_ratio) * 0.2
        if padding == 0: padding = 0.2
        h_ratio_acceptance.SetMaximum(max_ratio + padding)
        h_ratio_acceptance.SetMinimum(max(0.0, min_ratio - padding))

    make_canvas = lambda cname, ctitle: ROOT.TCanvas(cname, ctitle, 800, 600)
    def setup_canvas(c): c.SetTickx(1); c.SetTicky(1)

    c_lh2 = make_canvas(f"c_lh2_{binName}", "LH2 Acceptance"); setup_canvas(c_lh2); h_ratio_LH2.Draw("E1"); c_lh2.Write(); c_lh2.SaveAs(f"acceptance_LH2_{binName}.pdf")
    c_ld2 = make_canvas(f"c_ld2_{binName}", "LD2 Acceptance"); setup_canvas(c_ld2); h_ratio_LD2.Draw("E1"); c_ld2.Write(); c_ld2.SaveAs(f"acceptance_LD2_{binName}.pdf")
    c_comb = make_canvas(f"c_comb_{binName}", "Combined Acceptance"); setup_canvas(c_comb); h_ratio_combine.Draw("E1"); c_comb.Write(); c_comb.SaveAs(f"acceptance_combine_{binName}.pdf")

    c_overlay = make_canvas(f"c_overlay_{binName}", "Acceptances Overlay")
    setup_canvas(c_overlay)
    h_ratio_LH2.Draw("E1")
    h_ratio_LD2.Draw("E1 SAME")
    h_ratio_combine.Draw("E1 SAME")

    leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
    leg.SetBorderSize(0); leg.SetFillColor(ROOT.kWhite)
    leg.AddEntry(h_ratio_LD2, "LD2", "l")
    leg.AddEntry(h_ratio_LH2, "LH2", "l")
    leg.AddEntry(h_ratio_combine, "LH2+LD2", "l")
    leg.Draw()

    c_overlay.Write()
    c_overlay.SaveAs(f"acceptance_overlay_{binName}.pdf")

    c_dr = make_canvas(f"c_dr_{binName}", "Double Ratio")
    setup_canvas(c_dr)
    h_ratio_acceptance.Draw("E1")
    add_fit_and_band(h_ratio_acceptance, c_dr)

    c_dr.Write()
    c_dr.SaveAs(f"acceptance_ratio_{binName}.pdf")

    h_LH2_thrown.Write(); h_LH2_accept.Write(); h_ratio_LH2.Write()
    h_LD2_thrown.Write(); h_LD2_accept.Write(); h_ratio_LD2.Write()
    h_ratio_combine.Write(); h_ratio_acceptance.Write()


# ==============================================================================
# Yield Ratio Plots (Split Canvas) 
# ==============================================================================
def create_split_ratio_canvas(var_name, var_edges, t_lh2, t_ld2, mask_lh2, mask_ld2, plot_name, title, x_title, out_dir):
    out_dir.cd()
    
    data_lh2 = t_lh2[var_name][mask_lh2]
    weights_lh2 = t_lh2.weight[mask_lh2]
    data_ld2 = t_ld2[var_name][mask_ld2]
    weights_ld2 = t_ld2.weight[mask_ld2]
    
    lh2_h, lh2_err = get_weighted_histogram(data_lh2, weights_lh2, var_edges)
    ld2_h, ld2_err = get_weighted_histogram(data_ld2, weights_ld2, var_edges)
    ratio, ratio_err = calc_ratio_errors_independent(lh2_h, ld2_h, lh2_err, ld2_err)
    
    var_edges_root = array('d', var_edges)
    h_lh2 = make_th1f(f"h_lh2_clean_{plot_name}", title, var_edges_root, lh2_h, lh2_err)
    h_ld2 = make_th1f(f"h_ld2_clean_{plot_name}", title, var_edges_root, ld2_h, ld2_err)
    h_ratio = make_th1f(f"h_ratio_clean_{plot_name}", "", var_edges_root, ratio, ratio_err)
    
    format_hist(h_lh2, x_title, "Clean Yield", ROOT.kBlue)
    format_hist(h_ld2, x_title, "Clean Yield", ROOT.kRed)
    format_hist(h_ratio, x_title, "LH2 / LD2", ROOT.kBlack)
    
    c = ROOT.TCanvas(f"c_{plot_name}", title, 800, 800)
    
    pad1 = ROOT.TPad(f"pad1_{plot_name}", "pad1", 0, 0.35, 1, 1.0)
    pad1.SetBottomMargin(0.15); pad1.SetTickx(1); pad1.SetTicky(1)
    pad1.Draw()
    
    pad2 = ROOT.TPad(f"pad2_{plot_name}", "pad2", 0, 0.0, 1, 0.35)
    pad2.SetTopMargin(0.05); pad2.SetBottomMargin(0.3); pad2.SetTickx(1); pad2.SetTicky(1)
    pad2.Draw()
    
    pad1.cd()
    # Intelligent Scaling for Yields (Top Pad)
    max_y = max(get_hist_max_with_error(h_lh2), get_hist_max_with_error(h_ld2))
    h_lh2.SetMaximum(max_y * 1.3)
    h_lh2.SetMinimum(0)
    
    h_lh2.GetXaxis().SetLabelSize(0.04)
    h_lh2.GetXaxis().SetTitleSize(0.045)
    h_lh2.GetYaxis().SetLabelSize(0.04)
    h_lh2.GetYaxis().SetTitleSize(0.045)
    h_lh2.GetYaxis().SetTitleOffset(1.2)
    h_lh2.Draw("E1")
    h_ld2.Draw("E1 SAME")
    
    leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    leg.SetBorderSize(0); leg.SetFillColor(ROOT.kWhite)
    leg.AddEntry(h_lh2, "LH2 Clean", "lep")
    leg.AddEntry(h_ld2, "LD2 Clean", "lep")
    leg.Draw()
    
    pad2.cd()
    # Intelligent Scaling for Ratios (Bottom Pad)
    valid_mask = ratio > 0
    if np.any(valid_mask):
        max_ratio_val = np.max(ratio[valid_mask] + ratio_err[valid_mask])
        min_ratio_val = np.min(ratio[valid_mask] - ratio_err[valid_mask])
        range_padding = (max_ratio_val - min_ratio_val) * 0.2
        if range_padding == 0: range_padding = 0.2
        
        h_ratio.SetMaximum(max_ratio_val + range_padding)
        h_ratio.SetMinimum(max(0.0, min_ratio_val - range_padding))
    else:
        h_ratio.SetMaximum(2.0); h_ratio.SetMinimum(0.0)
    
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetLabelSize(0.08)
    h_ratio.GetYaxis().SetTitleSize(0.1)
    h_ratio.GetYaxis().SetTitleOffset(0.5)
    
    h_ratio.GetXaxis().SetLabelSize(0.1)
    h_ratio.GetXaxis().SetTitleSize(0.12)
    h_ratio.GetXaxis().SetTitleOffset(1.0)
    h_ratio.Draw("E1")
    
    add_fit_and_band(h_ratio, pad2)
    
    c.Write(); c.SaveAs(f"{plot_name}.pdf")
    h_lh2.Write(); h_ld2.Write(); h_ratio.Write()

# ==============================================================================
# Acceptance Ratio Plots (Split Canvas)
# ==============================================================================
def create_split_acceptance_canvas(var_name, var_edges, t_lh2_th, t_lh2_ac, t_ld2_th, t_ld2_ac,
                                   mask_lh2_th, mask_lh2_ac, mask_ld2_th, mask_ld2_ac, 
                                   plot_name, title, x_title, out_dir):
    out_dir.cd()
    
    lh2_th_h, lh2_th_err = get_weighted_histogram(t_lh2_th[var_name][mask_lh2_th], t_lh2_th.weight[mask_lh2_th], var_edges)
    ld2_th_h, ld2_th_err = get_weighted_histogram(t_ld2_th[var_name][mask_ld2_th], t_ld2_th.weight[mask_ld2_th], var_edges)
    
    lh2_ac_h, lh2_ac_err = get_weighted_histogram(t_lh2_ac[var_name][mask_lh2_ac], t_lh2_ac.weight[mask_lh2_ac], var_edges)
    ld2_ac_h, ld2_ac_err = get_weighted_histogram(t_ld2_ac[var_name][mask_ld2_ac], t_ld2_ac.weight[mask_ld2_ac], var_edges)
    
    lh2_acc = np.divide(lh2_ac_h, lh2_th_h, out=np.zeros_like(lh2_ac_h), where=lh2_th_h != 0)
    lh2_acc_err = calc_binomial_errors(lh2_th_h, lh2_acc, lh2_th_err)
    
    ld2_acc = np.divide(ld2_ac_h, ld2_th_h, out=np.zeros_like(ld2_ac_h), where=ld2_th_h != 0)
    ld2_acc_err = calc_binomial_errors(ld2_th_h, ld2_acc, ld2_th_err)
    
    ratio, ratio_err = calc_ratio_errors_independent(lh2_acc, ld2_acc, lh2_acc_err, ld2_acc_err)
    
    var_edges_root = array('d', var_edges)
    h_lh2 = make_th1f(f"h_lh2_acc_{plot_name}", title, var_edges_root, lh2_acc, lh2_acc_err)
    h_ld2 = make_th1f(f"h_ld2_acc_{plot_name}", title, var_edges_root, ld2_acc, ld2_acc_err)
    h_ratio = make_th1f(f"h_ratio_acc_{plot_name}", "", var_edges_root, ratio, ratio_err)
    
    format_hist(h_lh2, x_title, "Acceptance", ROOT.kBlue)
    format_hist(h_ld2, x_title, "Acceptance", ROOT.kRed)
    format_hist(h_ratio, x_title, "LH2 / LD2 Acc. Ratio", ROOT.kBlack)
    
    c = ROOT.TCanvas(f"c_acc_{plot_name}", title, 800, 800)
    
    pad1 = ROOT.TPad(f"pad1_acc_{plot_name}", "pad1", 0, 0.35, 1, 1.0)
    pad1.SetBottomMargin(0.15); pad1.SetTickx(1); pad1.SetTicky(1)
    pad1.Draw()
    
    pad2 = ROOT.TPad(f"pad2_acc_{plot_name}", "pad2", 0, 0.0, 1, 0.35)
    pad2.SetTopMargin(0.05); pad2.SetBottomMargin(0.3); pad2.SetTickx(1); pad2.SetTicky(1)
    pad2.Draw()
    
    pad1.cd()
    # Intelligent Scaling for Acceptances (Top Pad)
    max_y = max(get_hist_max_with_error(h_lh2), get_hist_max_with_error(h_ld2))
    h_lh2.SetMaximum(max_y * 1.3) 
    h_lh2.SetMinimum(0)
    
    h_lh2.GetXaxis().SetLabelSize(0.04)
    h_lh2.GetXaxis().SetTitleSize(0.045)
    h_lh2.GetYaxis().SetLabelSize(0.04)
    h_lh2.GetYaxis().SetTitleSize(0.045)
    h_lh2.GetYaxis().SetTitleOffset(1.2)
    h_lh2.Draw("E1")
    h_ld2.Draw("E1 SAME")
    
    leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    leg.SetBorderSize(0); leg.SetFillColor(ROOT.kWhite)
    leg.AddEntry(h_lh2, "LH2 Acceptance", "lep")
    leg.AddEntry(h_ld2, "LD2 Acceptance", "lep")
    leg.Draw()
    
    pad2.cd()
    # Intelligent Scaling for Ratios (Bottom Pad)
    valid_mask = ratio > 0
    if np.any(valid_mask):
        max_ratio_val = np.max(ratio[valid_mask] + ratio_err[valid_mask])
        min_ratio_val = np.min(ratio[valid_mask] - ratio_err[valid_mask])
        range_padding = (max_ratio_val - min_ratio_val) * 0.2
        if range_padding == 0: range_padding = 0.2
        
        h_ratio.SetMaximum(max_ratio_val + range_padding)
        h_ratio.SetMinimum(max(0.0, min_ratio_val - range_padding))
    else:
        h_ratio.SetMaximum(2.0); h_ratio.SetMinimum(0.0)
    
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetLabelSize(0.08)
    h_ratio.GetYaxis().SetTitleSize(0.08)
    h_ratio.GetYaxis().SetTitleOffset(0.6)
    
    h_ratio.GetXaxis().SetLabelSize(0.1)
    h_ratio.GetXaxis().SetTitleSize(0.12)
    h_ratio.GetXaxis().SetTitleOffset(1.0)
    h_ratio.Draw("E1")
    
    add_fit_and_band(h_ratio, pad2)
    
    c.Write(); c.SaveAs(f"{plot_name}.pdf")
    h_lh2.Write(); h_ld2.Write(); h_ratio.Write()

# ==============================================================================
# Main Execution
# ==============================================================================
def main():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(1111)
    ROOT.gROOT.SetBatch(True) 

    massEdge = np.array([3.9, 4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.8, 10.0], dtype=float)
    massEdge_root = array('d', massEdge) 
    
    # 18 Bins representing your full fine xF structure
    xFEdge = np.round(np.arange(-0.05, 0.90, 0.05), 2)
    
    pTEdge_fine = np.linspace(0.0, 3.0, 61) 
    pTEdge_user = np.array([0., 0.32, 0.49, 0.63, 0.77, 0.95, 1.18, 1.8], dtype=float)
    pT2Edge_user = np.square(pTEdge_user) 

    lh2_file = "/root/github/e906-development/ROOTFiles/Kenichi/rs67_lh2_acc.root"
    ld2_file = "/root/github/e906-development/ROOTFiles/Kenichi/rs67_ld2_acc.root"

    print("Loading trees from ROOT files...")
    
    t_lh2_thrown = load_and_format_tree(lh2_file, "tree_4pi")
    t_lh2_accept = load_and_format_tree(lh2_file, "tree_acc")
    t_ld2_thrown = load_and_format_tree(ld2_file, "tree_4pi")
    t_ld2_accept = load_and_format_tree(ld2_file, "tree_acc")

    print("Calculating derived pT2...")
    t_lh2_thrown = ak.with_field(t_lh2_thrown, t_lh2_thrown.pT**2, "pT2")
    t_ld2_thrown = ak.with_field(t_ld2_thrown, t_ld2_thrown.pT**2, "pT2")
    t_lh2_accept = ak.with_field(t_lh2_accept, t_lh2_accept.pT**2, "pT2")
    t_ld2_accept = ak.with_field(t_ld2_accept, t_ld2_accept.pT**2, "pT2")

    print("Applying Generator-Level Fiducial Cuts to Thrown Trees (Updated pT up to 3.0)...")
    th_fiducial_lh2 = (
        (t_lh2_thrown.xF >= -0.2) & (t_lh2_thrown.xF <= 1.0) & 
        (t_lh2_thrown.mass >= 3.0) & (t_lh2_thrown.mass <= 12.0) &
        (t_lh2_thrown.pT > 0.0) & (t_lh2_thrown.pT <= 3.0)
    )
    th_fiducial_ld2 = (
        (t_ld2_thrown.xF >= -0.2) & (t_ld2_thrown.xF <= 1.0) & 
        (t_ld2_thrown.mass >= 3.0) & (t_ld2_thrown.mass <= 12.0) &
        (t_ld2_thrown.pT > 0.0) & (t_ld2_thrown.pT <= 3.0)
    )
    t_lh2_thrown = t_lh2_thrown[th_fiducial_lh2]
    t_ld2_thrown = t_ld2_thrown[th_fiducial_ld2]

    out_file = ROOT.TFile("acceptance_mass_xF_unfolding.root", "RECREATE")

    # Applies custom_name="xF" to generate xF_bin0 through xF_bin17
    print("\n--- Starting Mass Acceptances Sliced by xF ---")
    process_acceptance_sliced("mass", massEdge, "xF", xFEdge, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file, custom_name="xF")
    
    # CRITICAL FIX: Applies custom_name="pT" to generate pT_bin0 through pT_bin6
    print("\n--- Starting Mass Acceptances Sliced by pT ---")
    process_acceptance_sliced("mass", massEdge, "pT", pTEdge_user, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file, custom_name="pT")
    
    # Omitted custom_name here on purpose to avoid overwriting the xF_bin directories created above
    print("\n--- Starting pT Acceptances Sliced by xF ---")
    process_acceptance_sliced("pT", pTEdge_user, "xF", xFEdge, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file)

    print("\n--- Starting Fully Integrated 1D Plots ---")
    process_integrated_1D("mass", massEdge, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file)
    process_integrated_1D("xF", xFEdge, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file)
    process_integrated_1D("pT", pTEdge_user, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file)
    process_integrated_1D("pT2", pT2Edge_user, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file)

    # ==================================================================================
    # MASK GENERATION FOR INTEGRATION SCENARIOS
    # ==================================================================================
    m_base_lh2_th = (t_lh2_thrown.mass > -999); m_base_lh2_ac = (t_lh2_accept.mass > -999)
    m_base_ld2_th = (t_ld2_thrown.mass > -999); m_base_ld2_ac = (t_ld2_accept.mass > -999)

    m_xf1_lh2_th = (t_lh2_thrown.xF > 0.0) & (t_lh2_thrown.xF < 0.4); m_xf1_lh2_ac = (t_lh2_accept.xF > 0.0) & (t_lh2_accept.xF < 0.4)
    m_xf1_ld2_th = (t_ld2_thrown.xF > 0.0) & (t_ld2_thrown.xF < 0.4); m_xf1_ld2_ac = (t_ld2_accept.xF > 0.0) & (t_ld2_accept.xF < 0.4)

    m_xf2_lh2_th = (t_lh2_thrown.xF > 0.4) & (t_lh2_thrown.xF < 0.8); m_xf2_lh2_ac = (t_lh2_accept.xF > 0.4) & (t_lh2_accept.xF < 0.8)
    m_xf2_ld2_th = (t_ld2_thrown.xF > 0.4) & (t_ld2_thrown.xF < 0.8); m_xf2_ld2_ac = (t_ld2_accept.xF > 0.4) & (t_ld2_accept.xF < 0.8)

    m_mass1_lh2_th = (t_lh2_thrown.mass > 4.2) & (t_lh2_thrown.mass < 5.5); m_mass1_lh2_ac = (t_lh2_accept.mass > 4.2) & (t_lh2_accept.mass < 5.5)
    m_mass1_ld2_th = (t_ld2_thrown.mass > 4.2) & (t_ld2_thrown.mass < 5.5); m_mass1_ld2_ac = (t_ld2_accept.mass > 4.2) & (t_ld2_accept.mass < 5.5)

    m_mass3_lh2_th = (t_lh2_thrown.mass > 5.5) & (t_lh2_thrown.mass < 8.7); m_mass3_lh2_ac = (t_lh2_accept.mass > 5.5) & (t_lh2_accept.mass < 8.7)
    m_mass3_ld2_th = (t_ld2_thrown.mass > 5.5) & (t_ld2_thrown.mass < 8.7); m_mass3_ld2_ac = (t_ld2_accept.mass > 5.5) & (t_ld2_accept.mass < 8.7)

    # ==================================================================================
    # CLEAN YIELDS (Split Canvas) 
    # ==================================================================================
    print("\n--- Generating Custom Split Canvas YIELD Ratio Plots ---")
    ratio_out_dir = out_file.mkdir("Yield_Ratios_SplitCanvas")
    
    create_split_ratio_canvas("mass", massEdge, t_lh2_accept, t_ld2_accept, m_base_lh2_ac, m_base_ld2_ac, "Split_Mass_All_xF_pT", "Invariant Mass Yields (All x_{F}, p_{T})", "Mass [GeV]", ratio_out_dir)
    create_split_ratio_canvas("xF", xFEdge, t_lh2_accept, t_ld2_accept, m_base_lh2_ac, m_base_ld2_ac, "Split_xF_All_Mass_pT", "x_{F} Yields (All Mass, p_{T})", "x_{F}", ratio_out_dir)
    
    create_split_ratio_canvas("mass", massEdge, t_lh2_accept, t_ld2_accept, m_xf1_lh2_ac, m_xf1_ld2_ac, "Split_Mass_0.0_xF_0.4", "Invariant Mass Yields (0.0 #leq x_{F} < 0.4)", "Mass [GeV]", ratio_out_dir)
    create_split_ratio_canvas("mass", massEdge, t_lh2_accept, t_ld2_accept, m_xf2_lh2_ac, m_xf2_ld2_ac, "Split_Mass_0.4_xF_0.8", "Invariant Mass Yields (0.4 #leq x_{F} < 0.8)", "Mass [GeV]", ratio_out_dir)
    create_split_ratio_canvas("xF", xFEdge, t_lh2_accept, t_ld2_accept, m_mass1_lh2_ac, m_mass1_ld2_ac, "Split_xF_4.2_Mass_5.5", "x_{F} Yields (4.2 < Mass < 5.5)", "x_{F}", ratio_out_dir)
    create_split_ratio_canvas("xF", xFEdge, t_lh2_accept, t_ld2_accept, m_mass3_lh2_ac, m_mass3_ld2_ac, "Split_xF_5.5_Mass_8.7", "x_{F} Yields (5.5 < Mass < 8.7)", "x_{F}", ratio_out_dir)

    create_split_ratio_canvas("pT", pTEdge_fine, t_lh2_accept, t_ld2_accept, m_base_lh2_ac, m_base_ld2_ac, "Split_pT_Fine_All_Mass_xF", "p_{T} Yield (Fine, All Mass, x_{F})", "p_{T} [GeV/c]", ratio_out_dir)
    create_split_ratio_canvas("pT", pTEdge_fine, t_lh2_accept, t_ld2_accept, m_xf1_lh2_ac, m_xf1_ld2_ac, "Split_pT_Fine_0.0_xF_0.4", "p_{T} Yield (Fine, 0.0 #leq x_{F} < 0.4)", "p_{T} [GeV/c]", ratio_out_dir)
    create_split_ratio_canvas("pT", pTEdge_fine, t_lh2_accept, t_ld2_accept, m_xf2_lh2_ac, m_xf2_ld2_ac, "Split_pT_Fine_0.4_xF_0.8", "p_{T} Yield (Fine, 0.4 #leq x_{F} < 0.8)", "p_{T} [GeV/c]", ratio_out_dir)
    create_split_ratio_canvas("pT", pTEdge_fine, t_lh2_accept, t_ld2_accept, m_mass1_lh2_ac, m_mass1_ld2_ac, "Split_pT_Fine_4.2_Mass_5.5", "p_{T} Yield (Fine, 4.2 < Mass < 5.5)", "p_{T} [GeV/c]", ratio_out_dir)
    create_split_ratio_canvas("pT", pTEdge_fine, t_lh2_accept, t_ld2_accept, m_mass3_lh2_ac, m_mass3_ld2_ac, "Split_pT_Fine_5.5_Mass_8.7", "p_{T} Yield (Fine, 5.5 < Mass < 8.7)", "p_{T} [GeV/c]", ratio_out_dir)

    create_split_ratio_canvas("pT", pTEdge_user, t_lh2_accept, t_ld2_accept, m_base_lh2_ac, m_base_ld2_ac, "Split_pT_User_All_Mass_xF", "p_{T} Yield (User Bins, All Mass, x_{F})", "p_{T} [GeV/c]", ratio_out_dir)
    create_split_ratio_canvas("pT", pTEdge_user, t_lh2_accept, t_ld2_accept, m_xf1_lh2_ac, m_xf1_ld2_ac, "Split_pT_User_0.0_xF_0.4", "p_{T} Yield (User Bins, 0.0 #leq x_{F} < 0.4)", "p_{T} [GeV/c]", ratio_out_dir)
    create_split_ratio_canvas("pT", pTEdge_user, t_lh2_accept, t_ld2_accept, m_xf2_lh2_ac, m_xf2_ld2_ac, "Split_pT_User_0.4_xF_0.8", "p_{T} Yield (User Bins, 0.4 #leq x_{F} < 0.8)", "p_{T} [GeV/c]", ratio_out_dir)
    create_split_ratio_canvas("pT", pTEdge_user, t_lh2_accept, t_ld2_accept, m_mass1_lh2_ac, m_mass1_ld2_ac, "Split_pT_User_4.2_Mass_5.5", "p_{T} Yield (User Bins, 4.2 < Mass < 5.5)", "p_{T} [GeV/c]", ratio_out_dir)
    create_split_ratio_canvas("pT", pTEdge_user, t_lh2_accept, t_ld2_accept, m_mass3_lh2_ac, m_mass3_ld2_ac, "Split_pT_User_5.5_Mass_8.7", "p_{T} Yield (User Bins, 5.5 < Mass < 8.7)", "p_{T} [GeV/c]", ratio_out_dir)

    # ==================================================================================
    # ACCEPTANCE CORRECTIONS (Split Canvas) 
    # ==================================================================================
    print("\n--- Generating Custom Split Canvas ACCEPTANCE Ratio Plots ---")
    acc_ratio_out_dir = out_file.mkdir("Acceptance_Ratios_SplitCanvas")
    
    # --- xF Acceptance ---
    create_split_acceptance_canvas("xF", xFEdge, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, m_base_lh2_th, m_base_lh2_ac, m_base_ld2_th, m_base_ld2_ac, "Acceptance_xF_All_Mass_pT", "x_{F} Acceptance (All Mass, p_{T})", "x_{F}", acc_ratio_out_dir)
    create_split_acceptance_canvas("xF", xFEdge, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, m_mass1_lh2_th, m_mass1_lh2_ac, m_mass1_ld2_th, m_mass1_ld2_ac, "Acceptance_xF_4.2_Mass_5.5", "x_{F} Acceptance (4.2 < Mass < 5.5)", "x_{F}", acc_ratio_out_dir)
    create_split_acceptance_canvas("xF", xFEdge, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, m_mass3_lh2_th, m_mass3_lh2_ac, m_mass3_ld2_th, m_mass3_ld2_ac, "Acceptance_xF_5.5_Mass_8.7", "x_{F} Acceptance (5.5 < Mass < 8.7)", "x_{F}", acc_ratio_out_dir)

    # --- pT Acceptance ---
    create_split_acceptance_canvas("pT", pTEdge_user, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, m_base_lh2_th, m_base_lh2_ac, m_base_ld2_th, m_base_ld2_ac, "Acceptance_pT_All_Mass_xF", "p_{T} Acceptance (All Mass, x_{F})", "p_{T} [GeV/c]", acc_ratio_out_dir)
    create_split_acceptance_canvas("pT", pTEdge_user, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, m_xf1_lh2_th, m_xf1_lh2_ac, m_xf1_ld2_th, m_xf1_ld2_ac, "Acceptance_pT_0.0_xF_0.4", "p_{T} Acceptance (0.0 #leq x_{F} < 0.4)", "p_{T} [GeV/c]", acc_ratio_out_dir)
    create_split_acceptance_canvas("pT", pTEdge_user, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, m_xf2_lh2_th, m_xf2_lh2_ac, m_xf2_ld2_th, m_xf2_ld2_ac, "Acceptance_pT_0.4_xF_0.8", "p_{T} Acceptance (0.4 #leq x_{F} < 0.8)", "p_{T} [GeV/c]", acc_ratio_out_dir)
    create_split_acceptance_canvas("pT", pTEdge_user, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, m_mass1_lh2_th, m_mass1_lh2_ac, m_mass1_ld2_th, m_mass1_ld2_ac, "Acceptance_pT_4.2_Mass_5.5", "p_{T} Acceptance (4.2 < Mass < 5.5)", "p_{T} [GeV/c]", acc_ratio_out_dir)
    create_split_acceptance_canvas("pT", pTEdge_user, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, m_mass3_lh2_th, m_mass3_lh2_ac, m_mass3_ld2_th, m_mass3_ld2_ac, "Acceptance_pT_5.5_Mass_8.7", "p_{T} Acceptance (5.5 < Mass < 8.7)", "p_{T} [GeV/c]", acc_ratio_out_dir)

    # --- pT2 Acceptance ---
    create_split_acceptance_canvas("pT2", pT2Edge_user, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, m_base_lh2_th, m_base_lh2_ac, m_base_ld2_th, m_base_ld2_ac, "Acceptance_pT2_All_Mass_xF", "p_{T}^{2} Acceptance (All Mass, x_{F})", "p_{T}^{2} [(GeV/c)^{2}]", acc_ratio_out_dir)

    out_file.Write()
    out_file.Close()
    print("\nDone. Generated all PDFs and saved TH1 data to 'acceptance_mass_xF_unfolding.root'")

if __name__ == "__main__":
    main()