import uproot
import numpy as np
import ROOT
import os

# Run ROOT in batch mode to prevent graphical windows from opening
ROOT.gROOT.SetBatch(True)

# --- Consolidated Cut Function ---
def apply_cuts(tree, is_mc=False):
    """Applies standard physics cuts and Messy MC updates to the given TTree arrays in numpy."""
    # Load all branches for the cut calculations
    events = tree.arrays(library="np")
    
    class EventNamespace:
        def __init__(self, data):
            self.__dict__.update(data)
            
    e = EventNamespace(events)

    # Dynamic BO for Data vs Fixed BO for Messy MC
    if is_mc:
        bo = 1.6
    else:
        bo = np.where(e.runID >= 11000, 1.6, 0.4)

    dimuon_cut = (
        (np.abs(e.dx) < 0.25) & (np.abs(e.dy - bo) < 0.22) &
        (e.dz < -5.) & (e.dz > -280.) & (np.abs(e.dpx) < 1.8) & (np.abs(e.dpy) < 2.0) &
        (e.dpx * e.dpx + e.dpy * e.dpy < 5.) & (e.dpz < 116.) & (e.dpz > 38.) &
        (e.dx * e.dx + (e.dy - bo) * (e.dy - bo) < 0.06) &
        (e.xF < 0.95) & (e.xF > -0.1) & (e.xT > 0.05) & (e.xT <= 0.58) &
        (np.abs(e.costh) < 0.5) & (np.abs(e.trackSeparation) < 270.) &
        (e.chisq_dimuon < 18)
    )

    track1_cut = (
        (e.chisq1_target < 15.) & (e.pz1_st1 > 9.) & (e.pz1_st1 < 75.) & (e.nHits1 > 13) &
        (e.x1_t * e.x1_t + (e.y1_t - bo) * (e.y1_t - bo) < 320.) &
        (e.x1_d * e.x1_d + (e.y1_d - bo) * (e.y1_d - bo) < 1100.) &
        (e.x1_d * e.x1_d + (e.y1_d - bo) * (e.y1_d - bo) > 16.) &
        (e.chisq1_target < 1.5 * e.chisq1_upstream) & (e.chisq1_target < 1.5 * e.chisq1_dump) &
        (e.z1_v < -5.) & (e.z1_v > -320.) & (e.chisq1 / (e.nHits1 - 5) < 12) &
        ((e.y1_st1) / (e.y1_st3) < 1.) & (np.abs(np.abs(e.px1_st1 - e.px1_st3) - 0.416) < 0.008) &
        (np.abs(e.py1_st1 - e.py1_st3) < 0.008) & (np.abs(e.pz1_st1 - e.pz1_st3) < 0.08) &
        ((e.y1_st1) * (e.y1_st3) > 0.) & (np.abs(e.py1_st1) > 0.02)
    )

    track2_cut = (
        (e.chisq2_target < 15.) & (e.pz2_st1 > 9.) & (e.pz2_st1 < 75.) & (e.nHits2 > 13) &
        (e.x2_t * e.x2_t + (e.y2_t - bo) * (e.y2_t - bo) < 320.) &
        (e.x2_d * e.x2_d + (e.y2_d - bo) * (e.y2_d - bo) < 1100.) &
        (e.x2_d * e.x2_d + (e.y2_d - bo) * (e.y2_d - bo) > 16.) &
        (e.chisq2_target < 1.5 * e.chisq2_upstream) & (e.chisq2_target < 1.5 * e.chisq2_dump) &
        (e.z2_v < -5.) & (e.z2_v > -320.) & (e.chisq2 / (e.nHits2 - 5) < 12) &
        ((e.y2_st1) / (e.y2_st3) < 1.) & (np.abs(np.abs(e.px2_st1 - e.px2_st3) - 0.416) < 0.008) &
        (np.abs(e.py2_st1 - e.py2_st3) < 0.008) & (np.abs(e.pz2_st1 - e.pz2_st3) < 0.08) &
        ((e.y2_st1) * (e.y2_st3) > 0.) & (np.abs(e.py2_st1) > 0.02)
    )

    tracks_cut = (
        (np.abs(e.chisq1_target + e.chisq2_target - e.chisq_dimuon) < 2.) &
        ((e.y1_st3) * (e.y2_st3) < 0.) & (e.nHits1 + e.nHits2 > 29) &
        (e.nHits1St1 + e.nHits2St1 > 8) & (np.abs(e.x1_st1 + e.x2_st1) < 42)
    )

    occ_cut = (
        (e.D1 < 400) & (e.D2 < 400) & (e.D3 < 400) & (e.D1 + e.D2 + e.D3 < 1000)
    )

    D1_occ_cut = (
        (e.D1 > 20) & (e.D1 < 385)
    )

    xF_cut = (
        (e.xF < 0.80) & (e.xF > 0.0) 
    )

    mass_cut = (
        (e.mass > 4.2) & (e.mass < 8.8)
    )

    total_cut_mask = (track1_cut & track2_cut & tracks_cut & dimuon_cut & occ_cut & D1_occ_cut & xF_cut & mass_cut)
    return total_cut_mask

# --- Utility Functions ---
def extract_kinematics(filepaths, is_mc=False, treename="Tree"):
    """Reads a list of root files, applies cuts, and returns multiple kinematic arrays to save I/O overhead."""
    if isinstance(filepaths, str):
        filepaths = [filepaths]
        
    all_pts = []
    all_xfs = []
    all_masses = []
    all_weights = []
    
    for filepath in filepaths:
        print(f"Processing {filepath} (Tree: {treename})...")
        try:
            tree = uproot.open(f"{filepath}:{treename}")
        except uproot.exceptions.KeyInFileError:
            print(f"WARNING: Tree '{treename}' not found in {filepath}. Attempting fallback...")
            try:
                tree = uproot.open(filepath).values()[0] 
            except Exception as e:
                print(f"Failed to process fallback for {filepath}: {e}")
                continue
        except Exception as e:
            print(f"Failed to process {filepath}: {e}")
            continue
            
        mask = apply_cuts(tree, is_mc=is_mc)
        
        branches_to_load = ["dpx", "dpy", "xF", "mass"]
        
        # Load branches and assign weights
        if is_mc:
            try:
                events = tree.arrays(branches_to_load + ["ReWeight"], library="np")
                weights = events["ReWeight"][mask]
            except uproot.exceptions.KeyInFileError:
                print(f"  -> WARNING: 'ReWeight' branch missing in {filepath}. Defaulting to 1.0.")
                events = tree.arrays(branches_to_load, library="np")
                weights = np.ones(np.count_nonzero(mask), dtype=np.float64)
        else:
            events = tree.arrays(branches_to_load, library="np")
            weights = np.ones(np.count_nonzero(mask), dtype=np.float64)
            
        pt = np.sqrt(events["dpx"][mask]**2 + events["dpy"][mask]**2)
        xf = events["xF"][mask]
        mass = events["mass"][mask]
        
        all_pts.append(pt)
        all_xfs.append(xf)
        all_masses.append(mass)
        all_weights.append(weights)
        
    if len(all_pts) > 0:
        return np.concatenate(all_pts), np.concatenate(all_xfs), np.concatenate(all_masses), np.concatenate(all_weights)
    return np.array([]), np.array([]), np.array([]), np.array([])

def make_normalized_th1(name, title, xtitle, data_array, weight_array, color, bins, xmin, xmax):
    """Creates a TH1F from numpy arrays (value, weight) using ROOT.FillN and area normalizes it."""
    h = ROOT.TH1F(name, title, bins, xmin, xmax)
    h.Sumw2()
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetMarkerStyle(1)
    h.SetLineWidth(2)
    h.SetStats(0)
    
    if len(data_array) > 0:
        x_vals = data_array.astype(np.float64)
        w_vals = weight_array.astype(np.float64)
        h.FillN(len(data_array), x_vals, w_vals)
        
    integral = h.Integral()
    if integral > 0:
        h.Scale(1.0 / integral)
        
    return h

def make_subtracted_normalized_th1(name, title, xtitle, data_arr, w_data, mix_arr, w_mix, color, bins, xmin, xmax):
    """Creates a TH1F for Data, subtracts Mixed background, and area normalizes the result."""
    h_data = ROOT.TH1F(f"{name}_raw", title, bins, xmin, xmax)
    h_data.Sumw2()
    
    h_mix = ROOT.TH1F(f"{name}_mix", title, bins, xmin, xmax)
    h_mix.Sumw2()
    
    if len(data_arr) > 0:
        h_data.FillN(len(data_arr), data_arr.astype(np.float64), w_data.astype(np.float64))
        
    if len(mix_arr) > 0:
        h_mix.FillN(len(mix_arr), mix_arr.astype(np.float64), w_mix.astype(np.float64))
        
    h_final = h_data.Clone(name)
    h_final.Add(h_mix, -1)
    
    h_final.SetLineColor(color)
    h_final.SetMarkerStyle(20)
    h_final.SetMarkerSize(0.8)
    h_final.SetLineWidth(2)
    h_final.SetStats(0)
    
    integral = h_final.Integral()
    if integral > 0:
        h_final.Scale(1.0 / integral)
        
    return h_final

def plot_variable_canvas(target_name, var_name, xtitle, 
                         data_arr, w_data, 
                         mix_arr, w_mix, 
                         dy_arr, w_dy, 
                         bins, xmin, xmax, ymax=None):
    """Generates and saves a split canvas: Upper pad (Distributions) & Lower pad (Ratio)."""
    c_name = f"c_{target_name}_{var_name}"
    c = ROOT.TCanvas(c_name, f"{target_name} {var_name} Distribution", 800, 800)
    
    # --- Pad 1: Main Distributions ---
    pad1 = ROOT.TPad(f"pad1_{c_name}", "pad1", 0.0, 0.3, 1.0, 1.0)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.12)
    pad1.SetTickx(1)
    pad1.SetTicky(1)
    pad1.Draw()
    pad1.cd()
    
    h_data = make_subtracted_normalized_th1(f"h_data_{c_name}", f"{target_name} Area Normalized {var_name}", xtitle, data_arr, w_data, mix_arr, w_mix, ROOT.kBlack, bins, xmin, xmax)
    h_dy   = make_normalized_th1(f"h_dy_{c_name}", "Messy Drell-Yan", xtitle, dy_arr, w_dy, ROOT.kRed, bins, xmin, xmax)

    # Dynamic scaling if ymax is not explicitly provided
    if ymax is None:
        ymax = max(h_data.GetMaximum(), h_dy.GetMaximum()) * 1.3

    h_data.SetMinimum(0.0)
    h_data.SetMaximum(ymax)
    h_data.GetXaxis().SetLabelSize(0)
    h_data.GetXaxis().SetTitleSize(0)
    
    h_data.GetYaxis().SetTitle("Normalized Yield")
    h_data.GetYaxis().CenterTitle(True)
    h_data.GetYaxis().SetTitleFont(43)
    h_data.GetYaxis().SetTitleSize(22)
    h_data.GetYaxis().SetTitleOffset(1.6)
    h_data.GetYaxis().SetLabelFont(43)
    h_data.GetYaxis().SetLabelSize(18)

    h_data.Draw("E1")
    h_dy.Draw("HIST SAME")
    h_dy.Draw("E1 SAME")

    leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_data, f"Data - Mix ({target_name})", "lep")
    leg.AddEntry(h_dy, "Messy DY MC", "le")
    leg.Draw()

    c.cd()
    
    # --- Pad 2: Ratio Plot ---
    pad2 = ROOT.TPad(f"pad2_{c_name}", "pad2", 0.0, 0.0, 1.0, 0.3)
    pad2.SetTopMargin(0.04)
    pad2.SetBottomMargin(0.35)
    pad2.SetLeftMargin(0.12)
    pad2.SetTickx(1)
    pad2.SetTicky(1)
    pad2.Draw()
    pad2.cd()
    
    h_ratio = h_data.Clone(f"h_ratio_{c_name}")
    h_ratio.SetTitle("")
    h_ratio.Divide(h_dy)
    
    h_ratio.GetYaxis().SetTitle("Data / DY")
    h_ratio.GetYaxis().CenterTitle(True)
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetTitleFont(43)
    h_ratio.GetYaxis().SetTitleSize(22)
    h_ratio.GetYaxis().SetTitleOffset(1.6)
    h_ratio.GetYaxis().SetLabelFont(43)
    h_ratio.GetYaxis().SetLabelSize(18)
    
    h_ratio.GetXaxis().SetTitle(xtitle)
    h_ratio.GetXaxis().CenterTitle(True)
    h_ratio.GetXaxis().SetTitleFont(43)
    h_ratio.GetXaxis().SetTitleSize(22)
    h_ratio.GetXaxis().SetTitleOffset(3.2)
    h_ratio.GetXaxis().SetLabelFont(43)
    h_ratio.GetXaxis().SetLabelSize(18)

    h_ratio.SetMinimum(0.0)
    h_ratio.SetMaximum(2.5)
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetMarkerSize(0.8)
    h_ratio.Draw("ep")
    
    line = ROOT.TLine(xmin, 1.0, xmax, 1.0)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw()

    pdf_filename = f"{var_name}_Distribution_{target_name}.pdf"
    c.SaveAs(pdf_filename)
    print(f"Successfully saved: {pdf_filename} (Y-range: [0, {ymax:.3f}])")
    
    return c, pad1, pad2, h_data, h_dy, h_ratio, line, leg

# --- Configuration ---
base_dir = "/root/github/e906-development/ROOTFiles/Hugo"

lh2_files = [
    "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS57/merged_RS57_LH2_recoeff_hodoeff.root",
    "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS59/merged_RS59_LH2_recoeff_hodoeff.root",
    "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/trimmed_RS62_LH2_recoeff_hodoeff.root",
    "/root/github/e906-development/src/HodoEfficiency/RS67/AllTargets/merged_RS67_3089_LH2_recoeff_hodoeff.root",
    "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS70/merged_RS70_LH2_recoeff_hodoeff.root"
]

ld2_files = [
    "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS57/merged_RS57_LD2_recoeff_hodoeff.root",
    "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS59/merged_RS59_LD2_recoeff_hodoeff.root",
    "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS62/trimmed_RS62_LD2_recoeff_hodoeff.root",
    "/root/github/e906-development/src/HodoEfficiency/RS67/AllTargets/merged_RS67_3089_LD2_recoeff_hodoeff.root",
    "/root/github/e906-development/src/HodoEfficiency/RS57-70/RS70/merged_RS70_LD2_recoeff_hodoeff.root"
]

mc_dy_lh2 = os.path.join(base_dir, "mc_drellyan_LH2_M027_S001_messy_occ_pTxFweight_v2.root")
mc_dy_ld2 = os.path.join(base_dir, "mc_drellyan_LD2_M027_S001_messy_occ_pTxFweight_v2.root")

# --- Extract Data Arrays ---
# Notice we now unpack 4 items per function call: pT, xF, Mass, Weights
print("--- Extracting LH2 Data & Mixed Background ---")
pt_data_lh2, xf_data_lh2, mass_data_lh2, w_data_lh2 = extract_kinematics(lh2_files, is_mc=False, treename="result")
pt_mix_lh2,  xf_mix_lh2,  mass_mix_lh2,  w_mix_lh2  = extract_kinematics(lh2_files, is_mc=False, treename="result_mix")

print("\n--- Extracting LD2 Data & Mixed Background ---")
pt_data_ld2, xf_data_ld2, mass_data_ld2, w_data_ld2 = extract_kinematics(ld2_files, is_mc=False, treename="result")
pt_mix_ld2,  xf_mix_ld2,  mass_mix_ld2,  w_mix_ld2  = extract_kinematics(ld2_files, is_mc=False, treename="result_mix")

print("\n--- Extracting LH2 Messy MC ---")
pt_dy_lh2,   xf_dy_lh2,   mass_dy_lh2,   w_dy_lh2   = extract_kinematics(mc_dy_lh2, is_mc=True, treename="Tree")

print("\n--- Extracting LD2 Messy MC ---")
pt_dy_ld2,   xf_dy_ld2,   mass_dy_ld2,   w_dy_ld2   = extract_kinematics(mc_dy_ld2, is_mc=True, treename="Tree")


# --- Generate Plots ---
print("\n--- Generating Canvases (Data - Mixed) and saving to PDFs ---")

plot_refs = []

# --- 1. pT Plots ---
plot_refs.append(plot_variable_canvas(
    "LH2", "pT", "p_{T} [GeV/c]", 
    pt_data_lh2, w_data_lh2, pt_mix_lh2, w_mix_lh2, pt_dy_lh2, w_dy_lh2, 
    bins=50, xmin=0.0, xmax=2.5, ymax=0.07
))
plot_refs.append(plot_variable_canvas(
    "LD2", "pT", "p_{T} [GeV/c]", 
    pt_data_ld2, w_data_ld2, pt_mix_ld2, w_mix_ld2, pt_dy_ld2, w_dy_ld2, 
    bins=50, xmin=0.0, xmax=2.5, ymax=0.07
))

# --- 2. xF Plots ---
# ymax is left out so the function dynamically sets a limit that won't clip the distribution
plot_refs.append(plot_variable_canvas(
    "LH2", "xF", "x_{F}", 
    xf_data_lh2, w_data_lh2, xf_mix_lh2, w_mix_lh2, xf_dy_lh2, w_dy_lh2, 
    bins=40, xmin=0.0, xmax=0.8
))
plot_refs.append(plot_variable_canvas(
    "LD2", "xF", "x_{F}", 
    xf_data_ld2, w_data_ld2, xf_mix_ld2, w_mix_ld2, xf_dy_ld2, w_dy_ld2, 
    bins=40, xmin=0.0, xmax=0.8
))

# --- 3. Mass Plots ---
# ymax is left out here as well for dynamic scaling
plot_refs.append(plot_variable_canvas(
    "LH2", "Mass", "Mass [GeV/c^{2}]", 
    mass_data_lh2, w_data_lh2, mass_mix_lh2, w_mix_lh2, mass_dy_lh2, w_dy_lh2, 
    bins=46, xmin=4.2, xmax=8.8
))
plot_refs.append(plot_variable_canvas(
    "LD2", "Mass", "Mass [GeV/c^{2}]", 
    mass_data_ld2, w_data_ld2, mass_mix_ld2, w_mix_ld2, mass_dy_ld2, w_dy_ld2, 
    bins=46, xmin=4.2, xmax=8.8
))

print("\nDone! Exactly 6 PDFs (LH2 & LD2 distributions for pT, xF, and Mass) with Data-Mixed backgrounds and Ratio plots have been created.")