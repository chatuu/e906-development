import uproot
import awkward as ak
import numpy as np
import os
import ROOT
from array import array

# ==============================================================================
# Event Selection (Cuts)
# ==============================================================================
def e906_chuck_cuts(tree: uproot.TTree, cut=4.2, beam_offset: float = 1.6):
    """
    Applies Chuck's event selection criteria. 
    Returns an awkward array of the filtered events.
    """
    events = tree.arrays()
    
    dimuon_cut_2111_v42 = (
        (np.abs(events.dx) < 0.25) & (np.abs(events.dy - beam_offset) < 0.22) &
        (events.dz < -5.) & (events.dz > -280.) &
        (np.abs(events.dpx) < 1.8) & (np.abs(events.dpy) < 2.0) &
        (events.dpx**2 + events.dpy**2 < 5.) &
        (events.dpz < 116.) & (events.dpz > 38.) &
        (events.mass > cut) & (events.mass < 8.8) &
        (events.dx**2 + (events.dy - beam_offset)**2 < 0.06) &
        (events.xF < 0.95) & (events.xF > -0.1) &
        (events.xT > 0.05) & (events.xT <= 0.58) &
        (np.abs(events.costh) < 0.5) &
        (np.abs(events.trackSeparation) < 270.) &
        (events.chisq_dimuon < 18)
    )

    track1_cut_2111_v42 = (
        (events.chisq1_target < 15.) & (events.pz1_st1 > 9.) & (events.pz1_st1 < 75.) &
        (events.nHits1 > 13) &
        (events.x1_t**2 + (events.y1_t - beam_offset)**2 < 320.) &
        (events.x1_d**2 + (events.y1_d - beam_offset)**2 < 1100.) &
        (events.x1_d**2 + (events.y1_d - beam_offset)**2 > 16.) &
        (events.chisq1_target < 1.5 * events.chisq1_upstream) &
        (events.chisq1_target < 1.5 * events.chisq1_dump) &
        (events.z1_v < -5.) & (events.z1_v > -320.) &
        (events.chisq1 / (events.nHits1 - 5) < 12) &
        (events.y1_st1 / events.y1_st3 < 1.) & 
        (np.abs(np.abs(events.px1_st1 - events.px1_st3) - 0.416) < 0.008) &
        (np.abs(events.py1_st1 - events.py1_st3) < 0.008) &
        (np.abs(events.pz1_st1 - events.pz1_st3) < 0.08) &
        (events.y1_st1 * events.y1_st3 > 0.) & 
        (np.abs(events.py1_st1) > 0.02)
    )

    track2_cut_2111_v42 = (
        (events.chisq2_target < 15.) & (events.pz2_st1 > 9.) & (events.pz2_st1 < 75.) &
        (events.nHits2 > 13) &
        (events.x2_t**2 + (events.y2_t - beam_offset)**2 < 320.) &
        (events.x2_d**2 + (events.y2_d - beam_offset)**2 < 1100.) &
        (events.x2_d**2 + (events.y2_d - beam_offset)**2 > 16.) &
        (events.chisq2_target < 1.5 * events.chisq2_upstream) &
        (events.chisq2_target < 1.5 * events.chisq2_dump) &
        (events.z2_v < -5.) & (events.z2_v > -320.) &
        (events.chisq2 / (events.nHits2 - 5) < 12) &
        (events.y2_st1 / events.y2_st3 < 1.) & 
        (np.abs(np.abs(events.px2_st1 - events.px2_st3) - 0.416) < 0.008) &
        (np.abs(events.py2_st1 - events.py2_st3) < 0.008) &
        (np.abs(events.pz2_st1 - events.pz2_st3) < 0.08) &
        (events.y2_st1 * events.y2_st3 > 0.) & 
        (np.abs(events.py2_st1) > 0.02)
    )

    tracks_cut_2111_v42 = (
        (np.abs(events.chisq1_target + events.chisq2_target - events.chisq_dimuon) < 2.) &
        (events.y1_st3 * events.y2_st3 < 0.) & 
        (events.nHits1 + events.nHits2 > 29) &
        (events.nHits1St1 + events.nHits2St1 > 8) &
        (np.abs(events.x1_st1 + events.x2_st1) < 42)
    )

    occ_cut_2111_v42 = (
        (events.D1 < 400) & (events.D2 < 400) & (events.D3 < 400) &
        (events.D1 + events.D2 + events.D3 < 1000)
    )

    mask = track1_cut_2111_v42 & track2_cut_2111_v42 & tracks_cut_2111_v42 & dimuon_cut_2111_v42 & occ_cut_2111_v42
    return events[mask]

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

# ==============================================================================
# Plotting Generation Logic: Mass Acceptances
# ==============================================================================
def process_kinematic_bins(var_name, var_edges, massEdge, massEdge_root, t_lh2_th, t_lh2_ac, t_ld2_th, t_ld2_ac, out_file):
    n_bins = len(var_edges) - 1
    if var_name == "xF":
        acc_min, acc_max = 0.0, 0.1
        ratio_min, ratio_max = 0.0, 1.2
    elif var_name == "pT":
        acc_min, acc_max = 0.0, 0.04
        ratio_min, ratio_max = 0.0, 1.1
    else:
        acc_min, acc_max = 0.0, 0.1
        ratio_min, ratio_max = 0.0, 1.2

    for i in range(n_bins):
        val_low = var_edges[i]
        val_high = var_edges[i+1]
        title = f"{val_low:.2f} <= {var_name} < {val_high:.2f}; Mass [GeV]"
        binName = f"{var_name}_bin{i}"
        print(f"  -> Processing Sliced Bin: {binName}")

        out_dir = out_file.mkdir(binName)
        out_dir.cd()

        m_lh2_th = (t_lh2_th[var_name] >= val_low) & (t_lh2_th[var_name] < val_high)
        m_ld2_th = (t_ld2_th[var_name] >= val_low) & (t_ld2_th[var_name] < val_high)
        m_lh2_ac = (t_lh2_ac[var_name] >= val_low) & (t_lh2_ac[var_name] < val_high)
        m_ld2_ac = (t_ld2_ac[var_name] >= val_low) & (t_ld2_ac[var_name] < val_high)

        lh2_th_h, lh2_th_err = get_weighted_histogram(t_lh2_th.mass[m_lh2_th], t_lh2_th.ReWeight[m_lh2_th], massEdge)
        lh2_ac_h, lh2_ac_err = get_weighted_histogram(t_lh2_ac.mass[m_lh2_ac], t_lh2_ac.ReWeight[m_lh2_ac], massEdge)
        lh2_ratio = np.divide(lh2_ac_h, lh2_th_h, out=np.zeros_like(lh2_ac_h), where=lh2_th_h != 0)
        lh2_ratio_err = calc_binomial_errors(lh2_th_h, lh2_ratio, lh2_th_err)

        ld2_th_h, ld2_th_err = get_weighted_histogram(t_ld2_th.mass[m_ld2_th], t_ld2_th.ReWeight[m_ld2_th], massEdge)
        ld2_ac_h, ld2_ac_err = get_weighted_histogram(t_ld2_ac.mass[m_ld2_ac], t_ld2_ac.ReWeight[m_ld2_ac], massEdge)
        ld2_ratio = np.divide(ld2_ac_h, ld2_th_h, out=np.zeros_like(ld2_ac_h), where=ld2_th_h != 0)
        ld2_ratio_err = calc_binomial_errors(ld2_th_h, ld2_ratio, ld2_th_err)

        combine_ratio = 0.5 * (lh2_ratio + ld2_ratio)
        combine_ratio_err = 0.5 * np.sqrt(lh2_ratio_err**2 + ld2_ratio_err**2)

        valid_dr = (ld2_ratio > 0) & (lh2_ratio > 0)
        dr = np.zeros_like(lh2_ratio)
        dr_err = np.zeros_like(lh2_ratio)
        dr[valid_dr] = lh2_ratio[valid_dr] / ld2_ratio[valid_dr]
        dr_err[valid_dr] = dr[valid_dr] * np.sqrt((lh2_ratio_err[valid_dr] / lh2_ratio[valid_dr])**2 + (ld2_ratio_err[valid_dr] / ld2_ratio[valid_dr])**2)

        h_LH2_thrown = make_th1f(f"h_LH2_thrown_{binName}", title, massEdge_root, lh2_th_h, lh2_th_err, out_dir)
        h_LH2_accept = make_th1f(f"h_LH2_accept_{binName}", title, massEdge_root, lh2_ac_h, lh2_ac_err, out_dir)
        h_ratio_LH2 = make_th1f(f"h_ratio_LH2_{binName}", title, massEdge_root, lh2_ratio, lh2_ratio_err, out_dir)
        h_LD2_thrown = make_th1f(f"h_LD2_thrown_{binName}", title, massEdge_root, ld2_th_h, ld2_th_err, out_dir)
        h_LD2_accept = make_th1f(f"h_LD2_accept_{binName}", title, massEdge_root, ld2_ac_h, ld2_ac_err, out_dir)
        h_ratio_LD2 = make_th1f(f"h_ratio_LD2_{binName}", title, massEdge_root, ld2_ratio, ld2_ratio_err, out_dir)
        h_ratio_combine = make_th1f(f"h_ratio_combine_{binName}", title, massEdge_root, combine_ratio, combine_ratio_err, out_dir)
        h_ratio_acceptance = make_th1f(f"h_ratio_acceptance_{binName}", f"LH2/LD2 {title}", massEdge_root, dr, dr_err, out_dir)

        format_hist(h_ratio_LH2, "Mass [GeV]", "Acceptance", ROOT.kBlue)
        format_hist(h_ratio_LD2, "Mass [GeV]", "Acceptance", ROOT.kRed)
        format_hist(h_ratio_combine, "Mass [GeV]", "Acceptance", ROOT.kBlack)
        format_hist(h_ratio_acceptance, "Mass [GeV]", "LH2 / LD2 Ratio", ROOT.kBlack)

        max_y = max(h_ratio_LH2.GetMaximum(), h_ratio_LD2.GetMaximum())
        if max_y > 0:
            h_ratio_LH2.SetMaximum(max_y * 1.5)
            h_ratio_LD2.SetMaximum(max_y * 1.5)
            h_ratio_combine.SetMaximum(max_y * 1.5)

        h_ratio_LH2.SetMinimum(acc_min)
        h_ratio_LH2.SetMaximum(acc_max)
        h_ratio_LD2.SetMinimum(acc_min)
        h_ratio_LD2.SetMaximum(acc_max)
        h_ratio_combine.SetMinimum(acc_min)
        h_ratio_combine.SetMaximum(acc_max)
        h_ratio_acceptance.SetMinimum(ratio_min)
        h_ratio_acceptance.SetMaximum(ratio_max)

        make_canvas = lambda cname, ctitle: ROOT.TCanvas(cname, ctitle, 800, 600)
        def setup_canvas(c):
            c.SetTickx(1)
            c.SetTicky(1)

        c_lh2 = make_canvas(f"c_lh2_{binName}", "LH2 Acceptance")
        setup_canvas(c_lh2)
        h_ratio_LH2.Draw("E1")
        c_lh2.Write()
        c_lh2.SaveAs(f"acceptance_LH2_{binName}.pdf")

        c_ld2 = make_canvas(f"c_ld2_{binName}", "LD2 Acceptance")
        setup_canvas(c_ld2)
        h_ratio_LD2.Draw("E1")
        c_ld2.Write()
        c_ld2.SaveAs(f"acceptance_LD2_{binName}.pdf")

        c_comb = make_canvas(f"c_comb_{binName}", "Combined Acceptance")
        setup_canvas(c_comb)
        h_ratio_combine.Draw("E1")
        c_comb.Write()
        c_comb.SaveAs(f"acceptance_combine_{binName}.pdf")

        c_overlay = make_canvas(f"c_overlay_{binName}", "Acceptances Overlay")
        setup_canvas(c_overlay)
        h_ratio_LH2.Draw("E1")
        h_ratio_LD2.Draw("E1 SAME")
        h_ratio_combine.Draw("E1 SAME")

        leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillColor(ROOT.kWhite)
        leg.AddEntry(h_ratio_LD2, "LD2", "l")
        leg.AddEntry(h_ratio_LH2, "LH2", "l")
        leg.AddEntry(h_ratio_combine, "LH2+LD2", "l")
        leg.Draw()

        c_overlay.Write()
        c_overlay.SaveAs(f"acceptance_overlay_{binName}.pdf")

        c_dr = make_canvas(f"c_dr_{binName}", "Double Ratio")
        setup_canvas(c_dr)
        h_ratio_acceptance.Draw("E1")
        h_ratio_acceptance.Fit("pol0", "SQ")
        c_dr.Update() 
        
        stats = h_ratio_acceptance.FindObject("stats")
        if stats:
            stats.SetX1NDC(0.65); stats.SetX2NDC(0.88)
            stats.SetY1NDC(0.75); stats.SetY2NDC(0.88)
            stats.SetBorderSize(0); stats.Draw()

        fitter = ROOT.TVirtualFitter.GetFitter()
        if fitter:
            g_error_band = ROOT.TGraphErrors(h_ratio_acceptance)
            fitter.GetConfidenceIntervals(g_error_band)
            g_error_band.SetFillColor(ROOT.kPink - 9)
            g_error_band.SetFillStyle(3004)
            g_error_band.Draw("E3 SAME")

        h_ratio_acceptance.Draw("E1 SAME")
        c_dr.Write()
        c_dr.SaveAs(f"acceptance_ratio_{binName}.pdf")

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

    if var_name == "xF":
        acc_min, acc_max = 0.0, 0.1
        ratio_min, ratio_max = 0.0, 1.2
    elif var_name == "pT":
        acc_min, acc_max = 0.0, 0.04
        ratio_min, ratio_max = 0.0, 1.1
    else: 
        acc_min, acc_max = 0.0, 0.1
        ratio_min, ratio_max = 0.0, 1.2

    x_titles = {"mass": "Mass [GeV]", "xF": "x_{F}", "pT": "p_{T} [GeV/c]"}
    x_title = x_titles.get(var_name, var_name)

    lh2_th_h, lh2_th_err = get_weighted_histogram(t_lh2_th[var_name], t_lh2_th.ReWeight, var_edges)
    lh2_ac_h, lh2_ac_err = get_weighted_histogram(t_lh2_ac[var_name], t_lh2_ac.ReWeight, var_edges)
    lh2_ratio = np.divide(lh2_ac_h, lh2_th_h, out=np.zeros_like(lh2_ac_h), where=lh2_th_h != 0)
    lh2_ratio_err = calc_binomial_errors(lh2_th_h, lh2_ratio, lh2_th_err)

    ld2_th_h, ld2_th_err = get_weighted_histogram(t_ld2_th[var_name], t_ld2_th.ReWeight, var_edges)
    ld2_ac_h, ld2_ac_err = get_weighted_histogram(t_ld2_ac[var_name], t_ld2_ac.ReWeight, var_edges)
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

    h_ratio_LH2.SetMinimum(acc_min); h_ratio_LH2.SetMaximum(acc_max)
    h_ratio_LD2.SetMinimum(acc_min); h_ratio_LD2.SetMaximum(acc_max)
    h_ratio_combine.SetMinimum(acc_min); h_ratio_combine.SetMaximum(acc_max)
    h_ratio_acceptance.SetMinimum(ratio_min); h_ratio_acceptance.SetMaximum(ratio_max)

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
    h_ratio_acceptance.Fit("pol0", "SQ")
    c_dr.Update() 
    
    stats = h_ratio_acceptance.FindObject("stats")
    if stats:
        stats.SetX1NDC(0.65); stats.SetX2NDC(0.88)
        stats.SetY1NDC(0.75); stats.SetY2NDC(0.88)
        stats.SetBorderSize(0); stats.Draw()

    fitter = ROOT.TVirtualFitter.GetFitter()
    if fitter:
        g_error_band = ROOT.TGraphErrors(h_ratio_acceptance)
        fitter.GetConfidenceIntervals(g_error_band)
        g_error_band.SetFillColor(ROOT.kPink - 9)
        g_error_band.SetFillStyle(3004)
        g_error_band.Draw("E3 SAME")

    h_ratio_acceptance.Draw("E1 SAME")
    c_dr.Write()
    c_dr.SaveAs(f"acceptance_ratio_{binName}.pdf")

    h_LH2_thrown.Write(); h_LH2_accept.Write(); h_ratio_LH2.Write()
    h_LD2_thrown.Write(); h_LD2_accept.Write(); h_ratio_LD2.Write()
    h_ratio_combine.Write(); h_ratio_acceptance.Write()


# ==============================================================================
# Separated Kinematic dp_x and dp_y Shape Comparisons
# ==============================================================================
def create_yield_hist(name, title, data, weights, edges, color, style):
    """Creates a TH1F for kinematic yield comparison."""
    h_vals, h_errs = get_weighted_histogram(data, weights, edges)
    edges_root = array('d', edges)
    h = make_th1f(name, title, edges_root, h_vals, h_errs)
    h.SetLineColor(color)
    h.SetLineStyle(style)
    h.SetLineWidth(2)
    return h

def plot_kinematics_separated(var_name, title_suffix, file_prefix, var_edges, t_lh2_th, t_lh2_ac, t_ld2_th, t_ld2_ac, out_dir):
    """Draws Thrown (4pi) and Accepted (Clean) arrays on strictly separate canvases."""
    
    if var_name == "dpx": axis_title = "dp_{x} [GeV/c]"
    elif var_name == "dpy": axis_title = "dp_{y} [GeV/c]"
    elif var_name == "dpx2": axis_title = "dp_{x}^{2} [(GeV/c)^{2}]"
    elif var_name == "dpy2": axis_title = "dp_{y}^{2} [(GeV/c)^{2}]"
    else: axis_title = var_name

    # ----------------------------------------------------
    # 1. Thrown Canvas (4pi)
    # ----------------------------------------------------
    h_lh2_th = create_yield_hist(f"h_{var_name}_lh2_th_{file_prefix}", f"4#pi Yield of {var_name} ({title_suffix})", t_lh2_th[var_name], t_lh2_th.ReWeight, var_edges, ROOT.kBlue, 1)
    h_ld2_th = create_yield_hist(f"h_{var_name}_ld2_th_{file_prefix}", f"4#pi Yield of {var_name} ({title_suffix})", t_ld2_th[var_name], t_ld2_th.ReWeight, var_edges, ROOT.kRed, 1)

    h_lh2_th.GetXaxis().SetTitle(axis_title)
    h_lh2_th.GetYaxis().SetTitle("Weighted Yield") 
    h_lh2_th.GetXaxis().CenterTitle(True)
    h_lh2_th.GetYaxis().CenterTitle(True)

    max_y_th = max(h_lh2_th.GetMaximum(), h_ld2_th.GetMaximum())
    if max_y_th > 0:
        h_lh2_th.SetMaximum(max_y_th * 1.5)

    c_th = ROOT.TCanvas(f"c_{var_name}_th_{file_prefix}", f"Thrown {var_name}", 800, 600)
    c_th.SetTickx(1); c_th.SetTicky(1)
    
    h_lh2_th.Draw("HIST")
    h_ld2_th.Draw("HIST SAME")

    leg_th = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    leg_th.SetBorderSize(0)
    leg_th.SetFillColor(ROOT.kWhite)
    leg_th.AddEntry(h_lh2_th, "LH2 4#pi", "l")
    leg_th.AddEntry(h_ld2_th, "LD2 4#pi", "l")
    leg_th.Draw()

    c_th.Write()
    c_th.SaveAs(f"yield_th_{var_name}_{file_prefix}.pdf")

    # ----------------------------------------------------
    # 2. Accepted Canvas (Clean)
    # ----------------------------------------------------
    h_lh2_ac = create_yield_hist(f"h_{var_name}_lh2_ac_{file_prefix}", f"Clean Yield of {var_name} ({title_suffix})", t_lh2_ac[var_name], t_lh2_ac.ReWeight, var_edges, ROOT.kBlue, 1)
    h_ld2_ac = create_yield_hist(f"h_{var_name}_ld2_ac_{file_prefix}", f"Clean Yield of {var_name} ({title_suffix})", t_ld2_ac[var_name], t_ld2_ac.ReWeight, var_edges, ROOT.kRed, 1)

    h_lh2_ac.GetXaxis().SetTitle(axis_title)
    h_lh2_ac.GetYaxis().SetTitle("Weighted Yield") 
    h_lh2_ac.GetXaxis().CenterTitle(True)
    h_lh2_ac.GetYaxis().CenterTitle(True)

    max_y_ac = max(h_lh2_ac.GetMaximum(), h_ld2_ac.GetMaximum())
    if max_y_ac > 0:
        h_lh2_ac.SetMaximum(max_y_ac * 1.5)

    c_ac = ROOT.TCanvas(f"c_{var_name}_ac_{file_prefix}", f"Accepted {var_name}", 800, 600)
    c_ac.SetTickx(1); c_ac.SetTicky(1)
    
    h_lh2_ac.Draw("HIST")
    h_ld2_ac.Draw("HIST SAME")

    leg_ac = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    leg_ac.SetBorderSize(0)
    leg_ac.SetFillColor(ROOT.kWhite)
    leg_ac.AddEntry(h_lh2_ac, "LH2 Clean", "l")
    leg_ac.AddEntry(h_ld2_ac, "LD2 Clean", "l")
    leg_ac.Draw()

    c_ac.Write()
    c_ac.SaveAs(f"yield_ac_{var_name}_{file_prefix}.pdf")

def process_dp_integrated(dpEdge, dp2Edge, t_lh2_th, t_lh2_ac, t_ld2_th, t_ld2_ac, out_file):
    print("  -> Processing Integrated dp shapes...")
    out_dir = out_file.mkdir("Integrated_dp")
    out_dir.cd()
    
    plot_kinematics_separated("dpx", "Integrated", "integrated", dpEdge, t_lh2_th, t_lh2_ac, t_ld2_th, t_ld2_ac, out_dir)
    plot_kinematics_separated("dpy", "Integrated", "integrated", dpEdge, t_lh2_th, t_lh2_ac, t_ld2_th, t_ld2_ac, out_dir)
    plot_kinematics_separated("dpx2", "Integrated", "integrated", dp2Edge, t_lh2_th, t_lh2_ac, t_ld2_th, t_ld2_ac, out_dir)
    plot_kinematics_separated("dpy2", "Integrated", "integrated", dp2Edge, t_lh2_th, t_lh2_ac, t_ld2_th, t_ld2_ac, out_dir)

def process_dp_binned_xf(dpEdge, dp2Edge, xF_edges, t_lh2_th, t_lh2_ac, t_ld2_th, t_ld2_ac, out_file):
    print("  -> Processing Sliced dp shapes per xF bin...")
    n_bins = len(xF_edges) - 1
    
    for i in range(n_bins):
        val_low = xF_edges[i]
        val_high = xF_edges[i+1]
        binName = f"xF_bin{i}"
        title_suffix = f"{val_low:.2f} <= xF < {val_high:.2f}"

        m_lh2_th = (t_lh2_th["xF"] >= val_low) & (t_lh2_th["xF"] < val_high)
        m_ld2_th = (t_ld2_th["xF"] >= val_low) & (t_ld2_th["xF"] < val_high)
        m_lh2_ac = (t_lh2_ac["xF"] >= val_low) & (t_lh2_ac["xF"] < val_high)
        m_ld2_ac = (t_ld2_ac["xF"] >= val_low) & (t_ld2_ac["xF"] < val_high)

        out_dir = out_file.GetDirectory(binName)
        if not out_dir: out_dir = out_file.mkdir(binName)
        out_dir.cd()

        plot_kinematics_separated("dpx", title_suffix, binName, dpEdge, t_lh2_th[m_lh2_th], t_lh2_ac[m_lh2_ac], t_ld2_th[m_ld2_th], t_ld2_ac[m_ld2_ac], out_dir)
        plot_kinematics_separated("dpy", title_suffix, binName, dpEdge, t_lh2_th[m_lh2_th], t_lh2_ac[m_lh2_ac], t_ld2_th[m_ld2_th], t_ld2_ac[m_ld2_ac], out_dir)
        plot_kinematics_separated("dpx2", title_suffix, binName, dp2Edge, t_lh2_th[m_lh2_th], t_lh2_ac[m_lh2_ac], t_ld2_th[m_ld2_th], t_ld2_ac[m_ld2_ac], out_dir)
        plot_kinematics_separated("dpy2", title_suffix, binName, dp2Edge, t_lh2_th[m_lh2_th], t_lh2_ac[m_lh2_ac], t_ld2_th[m_ld2_th], t_ld2_ac[m_ld2_ac], out_dir)

# ==============================================================================
# Main Execution
# ==============================================================================
def main():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptFit(1111)
    ROOT.gROOT.SetBatch(True) 

    massEdge = np.array([4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7], dtype=float)
    massEdge_root = array('d', massEdge) 
    
    xFEdge = np.round(np.arange(0.0, 0.85, 0.05), 2)
    pTEdge = np.array([0., 0.32, 0.49, 0.63, 0.77, 0.95, 1.18, 1.8], dtype=float)
    
    dpEdge = np.linspace(-2.0, 2.0, 51)
    dp2Edge = np.linspace(0.0, 5.0, 51)

    MCdir = os.path.expanduser("~/github/e906-development/ROOTFiles/Hugo/")
    lh2_files = {
        'thrown': MCdir + "mc_drellyan_LH2_M027_S001_4pi_pTxFweight_v2.root",
        'accept': MCdir + "mc_drellyan_LH2_M027_S001_clean_occ_pTxFweight_v2.root"
    }
    ld2_files = {
        'thrown': MCdir + "mc_drellyan_LD2_M027_S001_4pi_pTxFweight_v2.root",
        'accept': MCdir + "mc_drellyan_LD2_M027_S001_clean_occ_pTxFweight_v2.root"
    }

    print("Loading trees...")
    t_lh2_thrown = uproot.open(lh2_files['thrown'])["Tree"].arrays(["mass", "xF", "dpx", "dpy", "ReWeight"])
    t_ld2_thrown = uproot.open(ld2_files['thrown'])["Tree"].arrays(["mass", "xF", "dpx", "dpy", "ReWeight"])
    t_lh2_accept = e906_chuck_cuts(uproot.open(lh2_files['accept'])["Tree"])
    t_ld2_accept = e906_chuck_cuts(uproot.open(ld2_files['accept'])["Tree"])

    print("Calculating derived kinematics (pT, dpx^2, dpy^2)...")
    t_lh2_thrown = ak.with_field(t_lh2_thrown, np.sqrt(t_lh2_thrown.dpx**2 + t_lh2_thrown.dpy**2), "pT")
    t_ld2_thrown = ak.with_field(t_ld2_thrown, np.sqrt(t_ld2_thrown.dpx**2 + t_ld2_thrown.dpy**2), "pT")
    t_lh2_accept = ak.with_field(t_lh2_accept, np.sqrt(t_lh2_accept.dpx**2 + t_lh2_accept.dpy**2), "pT")
    t_ld2_accept = ak.with_field(t_ld2_accept, np.sqrt(t_ld2_accept.dpx**2 + t_ld2_accept.dpy**2), "pT")

    t_lh2_thrown = ak.with_field(t_lh2_thrown, t_lh2_thrown.dpx**2, "dpx2")
    t_ld2_thrown = ak.with_field(t_ld2_thrown, t_ld2_thrown.dpx**2, "dpx2")
    t_lh2_accept = ak.with_field(t_lh2_accept, t_lh2_accept.dpx**2, "dpx2")
    t_ld2_accept = ak.with_field(t_ld2_accept, t_ld2_accept.dpx**2, "dpx2")

    t_lh2_thrown = ak.with_field(t_lh2_thrown, t_lh2_thrown.dpy**2, "dpy2")
    t_ld2_thrown = ak.with_field(t_ld2_thrown, t_ld2_thrown.dpy**2, "dpy2")
    t_lh2_accept = ak.with_field(t_lh2_accept, t_lh2_accept.dpy**2, "dpy2")
    t_ld2_accept = ak.with_field(t_ld2_accept, t_ld2_accept.dpy**2, "dpy2")

    print("Applying Generator-Level Fiducial Cuts to Thrown Trees...")
    th_fiducial_lh2 = (
        (t_lh2_thrown.xF > -0.1) & (t_lh2_thrown.xF < 0.95) & 
        (t_lh2_thrown.mass > 4.2) & (t_lh2_thrown.mass < 8.8) &
        (t_lh2_thrown.pT > 0.0) & (t_lh2_thrown.pT <= 1.8)
    )
    th_fiducial_ld2 = (
        (t_ld2_thrown.xF > -0.1) & (t_ld2_thrown.xF < 0.95) & 
        (t_ld2_thrown.mass > 4.2) & (t_ld2_thrown.mass < 8.8) &
        (t_ld2_thrown.pT > 0.0) & (t_ld2_thrown.pT <= 1.8)
    )
    t_lh2_thrown = t_lh2_thrown[th_fiducial_lh2]
    t_ld2_thrown = t_ld2_thrown[th_fiducial_ld2]

    out_file = ROOT.TFile("acceptance_mass_xF.root", "RECREATE")

    print("\n--- Starting xF Sliced Bins ---")
    process_kinematic_bins("xF", xFEdge, massEdge, massEdge_root, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file)
    print("\n--- Starting pT Sliced Bins ---")
    process_kinematic_bins("pT", pTEdge, massEdge, massEdge_root, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file)
    print("\n--- Starting Fully Integrated 1D Plots ---")
    process_integrated_1D("mass", massEdge, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file)
    process_integrated_1D("xF", xFEdge, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file)
    process_integrated_1D("pT", pTEdge, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file)
    print("\n--- Starting Kinematic dp_x & dp_y Yield Plots ---")
    process_dp_integrated(dpEdge, dp2Edge, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file)
    process_dp_binned_xf(dpEdge, dp2Edge, xFEdge, t_lh2_thrown, t_lh2_accept, t_ld2_thrown, t_ld2_accept, out_file)

    out_file.Write()
    out_file.Close()
    print("\nDone. Generated PDFs and saved data to 'acceptance_mass_xF.root'")

if __name__ == "__main__":
    main()