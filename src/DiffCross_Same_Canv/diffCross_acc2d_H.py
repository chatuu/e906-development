import ROOT
import array
import os
import multiprocessing
import time
from functools import partial

# --- Attempt to import cuts from cuts.py ---
try:
    import cuts
except ImportError:
    print("Error: cuts.py not found. Please create it and define your cuts.")
    print("This script requires the 'cuts.py' file used by the original scripts.")
    exit()

# --- Main Configuration ---
# This script will create and populate these directories.
DIR_ARCHIVED_RESULTS = "result_rootfiles_archived"
DIR_DIFFCROSS_RESULTS = "result_rootfiles_diffcross"
OUTPUT_DIR_PLOTS = "comparison_plots"
SUMMARY_PLOT_FILENAME = "summary_all_bins.pdf"
SUMMARY_ROOT_FILENAME = "summary_all_bins.root" # Added for ROOT file output

# --- Common Binning and Constants ---
XF_BINS = [
    (0.0, 0.05), (0.05, 0.1), (0.1, 0.15), (0.15, 0.2),
    (0.2, 0.25), (0.25, 0.3), (0.3, 0.35), (0.35, 0.4),
    (0.4, 0.45), (0.45, 0.5), (0.5, 0.55), (0.55, 0.6),
    (0.6, 0.65), (0.65, 0.7), (0.7, 0.75), (0.75, 0.8),
    (0.8, 0.85)
]
K_EFF = [
    1.77909, 1.74145, 1.62481, 1.63713, 1.57874, 1.58371, 1.47446,
    1.4607, 1.42335, 1.43693, 1.36176, 1.35986, 1.37244, 1.34141,
    1.33321, 1.30802
]
NBINS_MASS = 11
EDGES_MASS_LIST = [4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7]
EDGES_MASS_ARR = array.array('d', EDGES_MASS_LIST)
NUM_BINS_TO_PROCESS = min(len(XF_BINS), len(K_EFF))

# --- Configurations for each script version ---
CONFIG_ARCHIVED = {
    "name": "Archived",
    "output_dir": DIR_ARCHIVED_RESULTS,
    "file_paths": {
        "dataTree": "/root/github/e906-development/res/ROOTFiles/Shivangi/roadset57_70_R008_2111v42_tmp_noPhys.root",
        "mixTree": "/root/github/e906-development/res/ROOTFiles/Shivangi/mixFPGA4_67_R008_preBin2.root",
        "jpTree": "/root/github/e906-development/res/ROOTFiles/Shivangi/mc_jpsi_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
        "ppTree": "/root/github/e906-development/res/ROOTFiles/Shivangi/mc_psiprime_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
    },
    "tree_names": {"default": "Tree"},
    "norms": {
        "FLASK_NORM": 4.94392, "MIX_NORM": 8.50747, "JP_NORM": 0.0031626, "PP_NORM": 0.00517343
    },
    "analysis_logic": "archived"
}

CONFIG_DIFFCROSS = {
    "name": "DiffCross",
    "output_dir": DIR_DIFFCROSS_RESULTS,
    "file_paths": {
        "dataTree": "/root/github/e906-development/res/ROOTFiles/MixedEvents/merged_RS67_3089LH2.root",
        "mixTree": "/root/github/e906-development/res/ROOTFiles/MixedEvents/merged_RS67_3089LH2.root",
        "emptyFlaskTree": "/root/github/e906-development/res/ROOTFiles/MixedEvents/merged_RS67_3089flask.root",
        "mixFlaskTree": "/root/github/e906-development/res/ROOTFiles/MixedEvents/merged_RS67_3089flask.root",
        "jpTree": "/root/github/e906-development/res/ROOTFiles/Shivangi/mc_jpsi_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
        "ppTree": "/root/github/e906-development/res/ROOTFiles/Shivangi/mc_psiprime_LH2_M027_S001_messy_occ_pTxFweight_v2.root",
    },
    "tree_names": {
        "dataTree": "result", "mixTree": "result_mix", "emptyFlaskTree": "result", "mixFlaskTree": "result_mix", "default": "Tree"
    },
    "norms": {
        "FLASK_NORM": 4.94392, "MIX_NORM": 1.0, "JP_NORM": 2.1878e-03, "PP_NORM": 3.4826e-03
    },
    "analysis_logic": "diffcross"
}

# --- Utility functions for scaling ---
def scale_hist_y(hist, factor):
    """Scales the Y-axis of a histogram including its errors."""
    if not hist: return None
    for i in range(1, hist.GetNbinsX() + 1):
        content = hist.GetBinContent(i)
        error = hist.GetBinError(i)
        hist.SetBinContent(i, content * factor)
        hist.SetBinError(i, error * factor)
    return hist

def scale_graph_y(graph, factor):
    """Scales the Y-axis of a TGraph."""
    if not graph: return None
    n_points = graph.GetN()
    for i in range(n_points):
        # Use array.array for pass-by-reference behavior needed by GetPoint
        x = array.array('d', [0])
        y = array.array('d', [0])
        graph.GetPoint(i, x, y)
        graph.SetPoint(i, x[0], y[0] * factor) # Use the retrieved values
    return graph

# --- Generic Analysis Function ---
def run_analysis_for_bin(binNum, config):
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    binNumStr = str(binNum)
    pid_suffix = f"_{config['name']}_bin{binNumStr}_{os.getpid()}"
    low, high = XF_BINS[binNum]
    xFCut = f"xF >= {low} && xF < {high}"
    print(f"Process {os.getpid()} ({config['name']}): Starting bin {binNumStr} with xF cut: {xFCut}")
    tree_dict = {}
    open_files = []
    for tree_key, path in config["file_paths"].items():
        try:
            file_obj = ROOT.TFile.Open(path, "READ")
            if not file_obj or file_obj.IsZombie():
                print(f"ERROR ({config['name']}): Could not open file {path} for bin {binNumStr}"); return
            tree_name = config["tree_names"].get(tree_key, config["tree_names"]["default"])
            tree = file_obj.Get(tree_name)
            if not tree:
                print(f"ERROR ({config['name']}): Could not get TTree '{tree_name}' from {path} for bin {binNumStr}"); file_obj.Close(); return
            tree_dict[tree_key] = tree
            open_files.append(file_obj)
        except Exception as e:
            print(f"EXCEPTION opening file {path} for bin {binNumStr}: {e}"); return
    h2 = None
    acc_file_path = "/root/github/e906-development/res/ROOTFiles/Shivangi/acceptance_h.root"
    try:
        saveAcc = ROOT.TFile.Open(acc_file_path, "READ")
        h2_name = f"LH2_acc_M_xF_{binNumStr}"
        h2_temp = saveAcc.Get(h2_name)
        h2 = h2_temp.Clone(f"h2_clone{pid_suffix}"); h2.SetDirectory(0); saveAcc.Close()
    except Exception as e:
        print(f"EXCEPTION opening acceptance for bin {binNumStr}: {e}"); return
    hData = ROOT.TH1F("hData" + pid_suffix, "Data", NBINS_MASS, EDGES_MASS_ARR); hData.Sumw2()
    hMC_jp = ROOT.TH1F("hMC_jp" + pid_suffix, "MC J/psi", NBINS_MASS, EDGES_MASS_ARR); hMC_jp.Sumw2()
    hMC_pp = ROOT.TH1F("hMC_pp" + pid_suffix, "MC psi'", NBINS_MASS, EDGES_MASS_ARR); hMC_pp.Sumw2()
    dataTree = tree_dict.get("dataTree"); mixTree = tree_dict.get("mixTree")
    jpTree = tree_dict.get("jpTree"); ppTree = tree_dict.get("ppTree")
    if config["analysis_logic"] == "archived":
        hDataF = ROOT.TH1F("hDataF" + pid_suffix, "Data Flask", NBINS_MASS, EDGES_MASS_ARR); hDataF.Sumw2()
        hDataMix = ROOT.TH1F("hDataMix" + pid_suffix, "Data Mixed Event", NBINS_MASS, EDGES_MASS_ARR); hDataMix.Sumw2()
        cut_data_target = f"mass^3*(({cuts.gmcTemp_charm}) && ({cuts.physics}) && ({cuts.target1}) && ({xFCut}) && ({cuts.occ}) && ({cuts.tInt}))"
        cut_data_flask = f"mass^3*(({cuts.gmcTemp_charm}) && ({cuts.physics}) && ({xFCut}) && ({cuts.flask}) && ({cuts.occ}) && ({cuts.tInt}))"
        cut_data_mix = f"mass^3*(({cuts.gmcTemp_charm}) && ({cuts.physics}) && ({xFCut}) && (({cuts.target1}) || ({cuts.target3})))"
        cut_mc_jp = f"0.99^3*mass^3*sigWeight*(({cuts.gmcTemp_charm}) && ({cuts.mass99}) && ({xFCut}))"
        cut_mc_pp = f"0.99^3*mass^3*sigWeight*(({cuts.gmcTemp_charm}) && ({cuts.mass99}) && ({xFCut}))"
        dataTree.Draw(f"mass>>hData{pid_suffix}", cut_data_target, "goff"); dataTree.Draw(f"mass>>hDataF{pid_suffix}", cut_data_flask, "goff")
        mixTree.Draw(f"mass>>hDataMix{pid_suffix}", cut_data_mix, "goff"); jpTree.Draw(f"0.99*mass>>hMC_jp{pid_suffix}", cut_mc_jp, "goff")
        ppTree.Draw(f"0.99*mass>>hMC_pp{pid_suffix}", cut_mc_pp, "goff"); hDataF.Scale(config["norms"]["FLASK_NORM"])
        hData.Add(hDataF, -1); hDataMix.Scale(config["norms"]["MIX_NORM"]); hData.Add(hDataMix, -1)
    elif config["analysis_logic"] == "diffcross":
        hDataF = ROOT.TH1F("hDataF" + pid_suffix, "Data Flask", NBINS_MASS, EDGES_MASS_ARR); hDataF.Sumw2()
        hDataFMix = ROOT.TH1F("hDataFMix" + pid_suffix, "Data Flask Mix", NBINS_MASS, EDGES_MASS_ARR); hDataFMix.Sumw2()
        hDataMix = ROOT.TH1F("hDataMix" + pid_suffix, "Data Mixed Event", NBINS_MASS, EDGES_MASS_ARR); hDataMix.Sumw2()
        emptyFlaskTree = tree_dict.get("emptyFlaskTree"); mixFlaskTree = tree_dict.get("mixFlaskTree")
        cut_data_target = f"mass^3*(({cuts.gmcTemp_charm}) && ({cuts.physics}) && ({xFCut}) && ({cuts.occ}))"
        cut_data_flask = f"mass^3*(({cuts.gmcTemp_charm}) && ({cuts.physics}) && ({xFCut}) && ({cuts.occ}))"
        cut_data_mix = f"mass^3*(({cuts.gmcTemp_charm}) && ({cuts.physics}) && ({xFCut}) && ({cuts.occ}))"
        cut_mc_jp = f"0.99^3*mass^3*sigWeight*(({cuts.gmcTemp_charm}) && ({cuts.physics}) && ({xFCut}) && ({cuts.occ}) &&  ({cuts.mass99}))"
        cut_mc_pp = f"0.99^3*mass^3*sigWeight*(({cuts.gmcTemp_charm}) && ({cuts.physics}) && ({xFCut}) && ({cuts.occ}) &&  ({cuts.mass99}))"
        dataTree.Draw(f"mass>>hData{pid_suffix}", cut_data_target, "goff"); mixTree.Draw(f"mass>>hDataMix{pid_suffix}", cut_data_mix, "goff")
        emptyFlaskTree.Draw(f"mass>>hDataF{pid_suffix}", cut_data_flask, "goff"); mixFlaskTree.Draw(f"mass>>hDataFMix{pid_suffix}", cut_data_flask, "goff")
        jpTree.Draw(f"0.99*mass>>hMC_jp{pid_suffix}", cut_mc_jp, "goff"); ppTree.Draw(f"0.99*mass>>hMC_pp{pid_suffix}", cut_mc_pp, "goff")
        hDataF.Scale(config["norms"]["FLASK_NORM"]); hDataFMix.Scale(config["norms"]["FLASK_NORM"]); hDataMix.Scale(config["norms"]["MIX_NORM"])
        hData.Add(hDataMix, -1); hData.Add(hDataF, -1); hData.Add(hDataFMix, -1)
    hMC_jp.Scale(config["norms"]["JP_NORM"]); hMC_pp.Scale(config["norms"]["PP_NORM"])
    hData.Add(hMC_jp, -1); hData.Add(hMC_pp, -1)
    hData.Divide(h2); hData.Scale(K_EFF[binNum]); hData.Scale(3.48489e-8); hData.Scale(1.102)
    for k in range(NBINS_MASS):
        bin_width = EDGES_MASS_ARR[k+1] - EDGES_MASS_ARR[k]
        if bin_width > 0:
            hData.SetBinContent(k + 1, hData.GetBinContent(k + 1) / bin_width)
            hData.SetBinError(k + 1, hData.GetBinError(k + 1) / bin_width)
    hData.SetBinContent(1, 0.0)
    hAccCor = hData.Clone("hAccCor")
    output_hist_filename = os.path.join(config["output_dir"], f"LH2_{binNumStr}_updatedPoT.root")
    savehist_file = ROOT.TFile.Open(output_hist_filename, "RECREATE")
    if savehist_file and not savehist_file.IsZombie():
        hAccCor.Write(); savehist_file.Close()
    for f_obj in open_files:
        if f_obj and f_obj.IsOpen(): f_obj.Close()
    print(f"Process {os.getpid()} ({config['name']}): Finished bin {binNumStr}")

# --- Individual Comparison Plotting Function ---
def create_comparison_plot(binNum):
    ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0)
    binNumStr = str(binNum); pid_suffix = f"_bin{binNumStr}_{os.getpid()}"
    print(f"Plotting Process {os.getpid()}: Starting comparison for bin {binNumStr}")
    file_path_archived = os.path.join(DIR_ARCHIVED_RESULTS, f"LH2_{binNumStr}_updatedPoT.root")
    file_path_diffcross = os.path.join(DIR_DIFFCROSS_RESULTS, f"LH2_{binNumStr}_updatedPoT.root")
    
    f_arch = ROOT.TFile.Open(file_path_archived)
    hData_archived = f_arch.Get("hAccCor").Clone(f"hData_archived{pid_suffix}")
    hData_archived.SetDirectory(0)
    f_arch.Close()

    f_diff = ROOT.TFile.Open(file_path_diffcross)
    hData_diffcross = f_diff.Get("hAccCor").Clone(f"hData_diffcross{pid_suffix}")
    hData_diffcross.SetDirectory(0)
    f_diff.Close()
    
    c = ROOT.TCanvas(f"c{pid_suffix}", f"Comparison for xF bin {binNumStr}", 800, 600)
    c.SetLogy(); c.SetTickx(1); c.SetTicky(1)
    hData_archived.SetMarkerColor(ROOT.kRed); hData_archived.SetLineColor(ROOT.kRed); hData_archived.SetMarkerStyle(20)
    hData_archived.GetXaxis().SetRangeUser(4.2, 8.5); hData_archived.GetXaxis().SetTitle("M (GeV)"); hData_archived.GetXaxis().CenterTitle(True)
    hData_archived.GetYaxis().SetTitle("M^{3} d^{2}#sigma/dMdx_{F} (nb GeV^{2})"); hData_archived.GetYaxis().CenterTitle(True)
    hData_archived.SetTitle(f"Differential Cross Section Comparison, xF bin {binNumStr}"); hData_archived.Draw("E1")
    nnpdf4_file_path = "/root/github/e906-development/res/ROOTFiles/Shivangi/NNPDF40_xFnew_p.root"
    ct18_file_path = "/root/github/e906-development/res/ROOTFiles/Shivangi/CT18_xFnew_p.root"; nnpdf4, ct18 = None, None
    try:
        nnpdf4_file_obj = ROOT.TFile.Open(nnpdf4_file_path)
        nnpdf4_graph = nnpdf4_file_obj.Get(f"gr_xFbin{binNumStr}")
        nnpdf4 = nnpdf4_graph.Clone(f"nnpdf4_clone{pid_suffix}")
        nnpdf4.SetLineColor(ROOT.kBlue + 2); nnpdf4.SetFillColorAlpha(38, 0.5); nnpdf4.SetFillStyle(3002)
        nnpdf4.Draw("L3 SAME"); nnpdf4_file_obj.Close()
        ct18_file_obj = ROOT.TFile.Open(ct18_file_path)
        ct18_graph = ct18_file_obj.Get(f"gr_xFbin{binNumStr}")
        ct18 = ct18_graph.Clone(f"ct18_clone{pid_suffix}")
        ct18.SetLineColor(ROOT.kGreen + 2); ct18.SetFillColorAlpha(30, 0.5); ct18.SetFillStyle(3002)
        ct18.Draw("L3 SAME"); ct18_file_obj.Close()
    except Exception as e: print(f"Warning: Could not draw theory curves for bin {binNumStr}. Error: {e}")
    hData_diffcross.SetMarkerColor(ROOT.kBlue); hData_diffcross.SetLineColor(ROOT.kBlue); hData_diffcross.SetMarkerStyle(21); hData_diffcross.Draw("E1 SAME")
    hData_archived.Draw("E1 SAME")
    legend = ROOT.TLegend(0.55, 0.65, 0.88, 0.88); legend.SetBorderSize(0)
    legend.AddEntry(hData_archived, "E906 (LH2)", "lep"); legend.AddEntry(hData_diffcross, "e906 (NMSU)", "lep")
    if nnpdf4: legend.AddEntry(nnpdf4, "NNPDF 4.0", "lf")
    if ct18: legend.AddEntry(ct18, "CT18", "lf")
    legend.Draw()
    plot_filename = os.path.join(OUTPUT_DIR_PLOTS, f"comparison_plot_bin_{binNumStr}.pdf")
    c.SaveAs(plot_filename)
    print(f"Plotting Process {os.getpid()}: Saved comparison plot to {plot_filename}")

# --- Summary Plot Function ---
def create_summary_plot():
    """Creates a single canvas with all bins plotted, scaled vertically."""
    print("Starting summary plot generation...")
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    # This list will hold all ROOT objects that need to stay in memory for drawing.
    drawn_objects = []

    c = ROOT.TCanvas("c_summary", "Summary of all xF Bins", 1200, 1600)
    drawn_objects.append(c)
    c.SetLogy()
    c.SetTicky(1)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.05)
    
    colors = [ROOT.kRed, ROOT.kOrange-3, ROOT.kYellow+1, ROOT.kGreen+1, ROOT.kTeal-5, 
              ROOT.kCyan+1, ROOT.kAzure+1, ROOT.kBlue, ROOT.kBlue+2, ROOT.kViolet-5, 
              ROOT.kMagenta, ROOT.kMagenta+2, ROOT.kPink+1, ROOT.kGray, ROOT.kBlack, ROOT.kRed-10]
    
    # Use user-defined ranges
    frame = c.DrawFrame(2.0, 1e-5, 10.0, 1e30)
    drawn_objects.append(frame)
    frame.SetTitle("Comparison for All x_{F} Bins")
    frame.GetXaxis().SetTitle("M (GeV)")
    frame.GetXaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle("M^{3} d^{2}#sigma/dMdx_{F} (nb GeV^{2})  (arbitrarily scaled)")
    frame.GetYaxis().CenterTitle(True)

    for binNum in range(NUM_BINS_TO_PROCESS):
        scale_factor = 50**binNum # Increased scaling factor
        binNumStr = str(binNum)
        color = colors[binNum]

        f_arch = ROOT.TFile.Open(os.path.join(DIR_ARCHIVED_RESULTS, f"LH2_{binNumStr}_updatedPoT.root"))
        f_ct18 = ROOT.TFile.Open("/root/github/e906-development/res/ROOTFiles/Shivangi/CT18_xFnew_p.root")

        if not f_arch or f_arch.IsZombie() or not f_ct18 or f_ct18.IsZombie():
            print(f"Warning: Skipping bin {binNumStr} due to one or more missing/bad files.")
            if f_arch: f_arch.Close()
            if f_ct18: f_ct18.Close()
            continue
            
        h_exist_orig = f_arch.Get("hAccCor")
        g_ct18_orig = f_ct18.Get(f"gr_xFbin{binNumStr}")

        if not all([h_exist_orig, g_ct18_orig]): 
            print(f"Warning: Skipping bin {binNumStr} due to missing object inside a file.")
            f_arch.Close(); f_ct18.Close()
            continue

        h_exist = h_exist_orig.Clone(f"h_exist_{binNumStr}"); h_exist.SetDirectory(0)
        g_ct18 = g_ct18_orig.Clone(f"g_ct18_{binNumStr}")

        f_arch.Close(); f_ct18.Close()

        scale_hist_y(h_exist, scale_factor)
        h_exist.SetMarkerColor(color); h_exist.SetLineColor(color); h_exist.SetMarkerStyle(20)
        
        scale_graph_y(g_ct18, scale_factor)
        g_ct18.SetLineColor(color); g_ct18.SetLineWidth(2); g_ct18.SetLineStyle(2)
        
        h_exist.Draw("P E1 SAME")
        g_ct18.Draw("L SAME")
        
        low_xf, high_xf = XF_BINS[binNum]
        # Use the Y-value of the first data point for the label position
        label_y = h_exist.GetBinContent(1) if h_exist.GetBinContent(1) > 0 else scale_factor
        label = ROOT.TLatex(2.2, label_y, f"{low_xf:.2f} #leq x_{{F}} < {high_xf:.2f}")
        label.SetTextColor(color); label.SetTextSize(0.02); label.SetTextAlign(12)
        label.Draw()
        
        drawn_objects.extend([h_exist, g_ct18, label])

    # --- Create Legends ---
    leg = ROOT.TLegend(0.65, 0.80, 0.95, 0.94) # Moved to top-right
    leg.SetBorderSize(0); leg.SetFillStyle(0)
    
    dummy_exist = ROOT.TMarker(0,0,20); dummy_exist.SetMarkerColor(ROOT.kBlack)
    dummy_ct18 = ROOT.TLine(0,0,0,0); dummy_ct18.SetLineStyle(2); dummy_ct18.SetLineWidth(2); dummy_ct18.SetLineColor(ROOT.kBlack)

    leg.AddEntry(dummy_exist, "E906 (LH2)", "p")
    leg.AddEntry(dummy_ct18, "CT18", "l")
    leg.Draw()
    
    drawn_objects.extend([leg, dummy_exist, dummy_ct18])
    
    c.SaveAs(SUMMARY_PLOT_FILENAME)
    c.SaveAs(SUMMARY_ROOT_FILENAME) # Save as ROOT file
    print(f"Summary plot saved to {SUMMARY_PLOT_FILENAME} and {SUMMARY_ROOT_FILENAME}")
    
# --- Main Orchestrator ---
if __name__ == '__main__':
    ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0)
    total_start_time = time.time()
    bin_numbers = list(range(NUM_BINS_TO_PROCESS))
    num_processes = multiprocessing.cpu_count()
    
    print("\n" + "="*60 + "\nSTEP 1: Generating results for 'E906 (LH2)' analysis...\n" + "="*60)
    if not os.path.exists(DIR_ARCHIVED_RESULTS): os.makedirs(DIR_ARCHIVED_RESULTS)
    start_time = time.time()
    run_archived_partial = partial(run_analysis_for_bin, config=CONFIG_ARCHIVED)
    with multiprocessing.Pool(processes=num_processes) as pool: pool.map(run_archived_partial, bin_numbers)
    print(f"--- Step 1 finished in {time.time() - start_time:.2f} seconds. ---\n")

    print("="*60 + "\nSTEP 2: Generating results for 'NMSU' analysis...\n" + "="*60)
    if not os.path.exists(DIR_DIFFCROSS_RESULTS): os.makedirs(DIR_DIFFCROSS_RESULTS)
    start_time = time.time()
    run_diffcross_partial = partial(run_analysis_for_bin, config=CONFIG_DIFFCROSS)
    with multiprocessing.Pool(processes=num_processes) as pool: pool.map(run_diffcross_partial, bin_numbers)
    print(f"--- Step 2 finished in {time.time() - start_time:.2f} seconds. ---\n")
    
    print("="*60 + "\nSTEP 3: Generating individual comparison plots...\n" + "="*60)
    if not os.path.exists(OUTPUT_DIR_PLOTS): os.makedirs(OUTPUT_DIR_PLOTS)
    start_time = time.time()
    with multiprocessing.Pool(processes=num_processes) as pool: pool.map(create_comparison_plot, bin_numbers)
    print(f"--- Step 3 finished in {time.time() - start_time:.2f} seconds. ---\n")

    print("="*60 + "\nSTEP 4: Generating final summary plot...\n" + "="*60)
    start_time = time.time()
    create_summary_plot()
    print(f"--- Step 4 finished in {time.time() - start_time:.2f} seconds. ---\n")

    total_end_time = time.time()
    print("="*60)
    print(f"\nAll steps complete. Total execution time: {total_end_time - total_start_time:.2f} seconds.")
    print(f"Final summary plot is saved as '{SUMMARY_PLOT_FILENAME}' and '{SUMMARY_ROOT_FILENAME}'.\n" + "="*60)
