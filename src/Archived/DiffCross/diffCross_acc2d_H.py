import ROOT
import array # For creating TH1F with variable bin edges
import os
import multiprocessing # Import the multiprocessing module
import time # To measure execution time

# --- Import cuts from cuts.py ---
try:
    import cuts
except ImportError:
    print("Error: cuts.py not found. Please create it and define your cuts.")
    exit()
# --- End of Import cuts ---


# --- Global Configuration (Editable by user) ---
# Normalization constants
FLASK_NORM = 4.94392
MIX_NORM = 8.50747
JP_NORM = 0.0031626
PP_NORM = 0.00517343

# File paths
ROOT_FILE_PATH = "/root/github/e906-development/res/ROOTFiles/Shivangi/roadset57_70_R008_2111v42_tmp_noPhys.root"
MIX_FILE_NAME_PATH = "/root/github/e906-development/res/ROOTFiles/Shivangi/mixFPGA4_67_R008_preBin2.root"
JPSI_FILE_NAME_PATH = "/root/github/e906-development/res/ROOTFiles/Shivangi/mc_jpsi_LH2_M027_S001_messy_occ_pTxFweight_v2.root"
PSIP_FILE_NAME_PATH = "/root/github/e906-development/res/ROOTFiles/Shivangi/mc_psiprime_LH2_M027_S001_messy_occ_pTxFweight_v2.root"
ACCEPTANCE_H_FILE_PATH = "/root/github/e906-development/res/ROOTFiles/Shivangi/acceptance_h.root"
NNPDF4_FILE_PATH = "/root/github/e906-development/res/ROOTFiles/Shivangi/NNPDF40_xFnew_p.root"
CT18_FILE_PATH = "/root/github/e906-development/res/ROOTFiles/Shivangi/CT18_xFnew_p.root"
OUTPUT_DIR_RESULTS = "result_rootfiles"
OUTPUT_DIR_PLOTS = "plots" # Set to None if plots are not to be saved by default

# --- Ensure output directories exist ---
if OUTPUT_DIR_RESULTS and not os.path.exists(OUTPUT_DIR_RESULTS):
    os.makedirs(OUTPUT_DIR_RESULTS)
if OUTPUT_DIR_PLOTS and not os.path.exists(OUTPUT_DIR_PLOTS):
    os.makedirs(OUTPUT_DIR_PLOTS)

# --- xF Binning and k-Efficiency Factors ---
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
# Ensure K_EFF has the same number of elements as XF_BINS if they correspond 1-to-1
if len(K_EFF) != len(XF_BINS):
    print(f"Warning: Mismatch between number of xF bins ({len(XF_BINS)}) and kEff values ({len(K_EFF)}).")
    # Decide how to handle this: truncate, error out, or use a default kEff.
    # For now, we'll allow it but it might cause issues if kEff[binNum] is accessed out of bounds.

NBINS_MASS = 11 # Number of mass bins
EDGES_MASS_LIST = [4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7]
EDGES_MASS_ARR = array.array('d', EDGES_MASS_LIST)

def get_xf_cut(bin_num):
    """Returns the xF cut string for a given bin number."""
    if 0 <= bin_num < len(XF_BINS):
        low, high = XF_BINS[bin_num]
        return f"xF >= {low} && xF < {high}"
    else:
        print(f"Error in get_xf_cut: bin_num {bin_num} is out of range for XF_BINS.")
        return "1==0" # Return a cut that selects nothing

def diffCross_acc2d_H(binNum):
    # Set ROOT to batch mode for this process to avoid GUI conflicts
    ROOT.gROOT.SetBatch(True)
    # It's good practice to also set gStyle options here if they are process-specific
    # or ensure they are inherited correctly. Here, gStyle is simple.
    ROOT.gStyle.SetOptStat(0)


    binNumStr = str(binNum)
    xFCut = get_xf_cut(binNum)
    print(f"Process {os.getpid()}: Starting bin {binNumStr} with xF cut: {xFCut}")


    # --- Acceptance Histogram ---
    h2 = None
    try:
        saveAcc = ROOT.TFile.Open(ACCEPTANCE_H_FILE_PATH, "READ")
        if not saveAcc or saveAcc.IsZombie():
            print(f"Process {os.getpid()}: Error - Could not open acceptance file: {ACCEPTANCE_H_FILE_PATH} for bin {binNumStr}")
            return # Essential to return if a critical file is missing
        h2_name = f"LH2_acc_M_xF_{binNumStr}"
        h2_temp = saveAcc.Get(h2_name)
        if not h2_temp:
            print(f"Process {os.getpid()}: Error - Could not retrieve histogram '{h2_name}' from {ACCEPTANCE_H_FILE_PATH} for bin {binNumStr}")
            saveAcc.Close()
            return
        h2 = h2_temp.Clone(f"h2_clone_bin{binNumStr}_{os.getpid()}") # Clone to make it process-local
        h2.SetDirectory(0)
        saveAcc.Close()
    except Exception as e:
        print(f"Process {os.getpid()}: Exception opening/reading acceptance for bin {binNumStr}: {e}")
        return

    # --- Open data and MC files ---
    # It's generally safer to open files within the function for multiprocessing
    # to avoid issues with file descriptors shared across processes.
    tree_dict = {}
    file_paths_to_load = {
        "dataTree": ROOT_FILE_PATH,
        "mixTree": MIX_FILE_NAME_PATH,
        "jpTree": JPSI_FILE_NAME_PATH,
        "ppTree": PSIP_FILE_NAME_PATH
    }
    open_files = [] # Keep track of successfully opened files to close them

    for tree_name, path in file_paths_to_load.items():
        try:
            file_obj = ROOT.TFile.Open(path, "READ")
            if not file_obj or file_obj.IsZombie():
                print(f"Process {os.getpid()}: Error opening {tree_name} file: {path} for bin {binNumStr}")
                for f_obj in open_files: f_obj.Close() # Close already opened files
                return
            tree = file_obj.Get("Tree")
            if not tree:
                print(f"Process {os.getpid()}: Error getting Tree from {tree_name} file: {path} for bin {binNumStr}")
                file_obj.Close()
                for f_obj in open_files: f_obj.Close()
                return
            tree_dict[tree_name] = tree
            # tree.SetDirectory(0) # Not strictly necessary if file is closed after Draw
            open_files.append(file_obj) # Add to list for later closing
        except Exception as e:
            print(f"Process {os.getpid()}: Exception opening/reading {path} for bin {binNumStr}: {e}")
            for f_obj in open_files: f_obj.Close()
            return

    dataTree = tree_dict.get("dataTree")
    mixTree = tree_dict.get("mixTree")
    jpTree = tree_dict.get("jpTree")
    ppTree = tree_dict.get("ppTree")

    if not all([dataTree, mixTree, jpTree, ppTree]):
        print(f"Process {os.getpid()}: One or more TTrees could not be retrieved for bin {binNumStr}.")
        for f_obj in open_files: f_obj.Close()
        return

    # --- Create Histograms (unique names per process and bin) ---
    pid_suffix = f"_bin{binNumStr}_{os.getpid()}" # Unique suffix for histograms
    hData = ROOT.TH1F("hData" + pid_suffix, "Data", NBINS_MASS, EDGES_MASS_ARR)
    hData.Sumw2()
    hDataF = ROOT.TH1F("hDataF" + pid_suffix, "Data Flask", NBINS_MASS, EDGES_MASS_ARR)
    hDataF.Sumw2()
    hDataMix = ROOT.TH1F("hDataMix" + pid_suffix, "Data Mixed Event", NBINS_MASS, EDGES_MASS_ARR)
    hDataMix.Sumw2()
    hMC_jp = ROOT.TH1F("hMC_jp" + pid_suffix, "MC J/psi", NBINS_MASS, EDGES_MASS_ARR)
    hMC_jp.Sumw2()
    hMC_pp = ROOT.TH1F("hMC_pp" + pid_suffix, "MC psi'", NBINS_MASS, EDGES_MASS_ARR)
    hMC_pp.Sumw2()

    # --- Fill Histograms using cuts from cuts.py ---
    cut_data_target = f"mass^3*(({cuts.gmcTemp_charm}) && ({cuts.physics}) && ({cuts.target1}) && ({xFCut}) && ({cuts.occ}) && ({cuts.tInt}))"
    cut_data_flask = f"mass^3*(({cuts.gmcTemp_charm}) && ({cuts.physics}) && ({xFCut}) && ({cuts.flask}) && ({cuts.occ}) && ({cuts.tInt}))"
    cut_data_mix = f"mass^3*(({cuts.gmcTemp_charm}) && ({cuts.physics}) && ({xFCut}) && (({cuts.target1}) || ({cuts.target3})))"
    cut_mc_jp = f"0.99^3*mass^3*sigWeight*(({cuts.gmcTemp_charm}) && ({cuts.mass99}) && ({xFCut}))"
    cut_mc_pp = f"0.99^3*mass^3*sigWeight*(({cuts.gmcTemp_charm}) && ({cuts.mass99}) && ({xFCut}))"

    try:
        dataTree.Draw(f"mass>>hData{pid_suffix}", cut_data_target, "goff") # goff to not draw
        dataTree.Draw(f"mass>>hDataF{pid_suffix}", cut_data_flask, "goff")
        mixTree.Draw(f"mass>>hDataMix{pid_suffix}", cut_data_mix, "goff")
        jpTree.Draw(f"0.99*mass>>hMC_jp{pid_suffix}", cut_mc_jp, "goff")
        ppTree.Draw(f"0.99*mass>>hMC_pp{pid_suffix}", cut_mc_pp, "goff")
    except Exception as e:
        print(f"Process {os.getpid()}: Exception during TTree::Draw for bin {binNumStr}: {e}")
        for f_obj in open_files: f_obj.Close()
        return

    # print(f"Process {os.getpid()}: --- Bin Content Before Scaling (Bin {binNumStr}) ---")
    # for k in range(NBINS_MASS):
    #     print(f"Bin {k+1}: Data={hData.GetBinContent(k+1):.2f}, Flask={hDataF.GetBinContent(k+1):.2f}, Mix={hDataMix.GetBinContent(k+1):.2f}, Jpsi={hMC_jp.GetBinContent(k+1):.2f}, Psip={hMC_pp.GetBinContent(k+1):.2f}, Acc={h2.GetBinContent(k+1):.4f}")

    # --- Scaling and Normalization ---
    hDataF.Scale(FLASK_NORM)
    hData.Add(hDataF, -1)

    hDataMix.Scale(MIX_NORM)
    hData.Add(hDataMix, -1)

    hMC_jp.Scale(JP_NORM)
    hMC_pp.Scale(PP_NORM)
    hData.Add(hMC_jp, -1)
    hData.Add(hMC_pp, -1)

    hData.Divide(h2)
    if 0 <= binNum < len(K_EFF):
        hData.Scale(K_EFF[binNum])
        # print(f"Process {os.getpid()}: Applied kEff[{binNum}]: {K_EFF[binNum]} for bin {binNumStr}")
    else:
        print(f"Process {os.getpid()}: Warning - binNum {binNum} out of range for K_EFF array. No kEff scaling applied for bin {binNumStr}.")

    hData.Scale(3.48489e-8)
    hData.Scale(1.102)

    for k in range(NBINS_MASS):
        bin_width = EDGES_MASS_ARR[k+1] - EDGES_MASS_ARR[k]
        if bin_width > 0:
            hData.SetBinContent(k + 1, hData.GetBinContent(k + 1) / bin_width)
            hData.SetBinError(k + 1, hData.GetBinError(k + 1) / bin_width)
        else:
            hData.SetBinContent(k + 1, 0)
            hData.SetBinError(k + 1, 0)

    for k in range(NBINS_MASS):
        content = hData.GetBinContent(k + 1)
        error = hData.GetBinError(k + 1)
        if content != 0 and abs(error / content) > 0.99 :
            hData.SetBinContent(k + 1, 0.0)
            hData.SetBinError(k + 1, 0.0)

    hData.SetBinContent(1, 0.0)

    # --- Plotting (can be skipped in batch processing or done later) ---
    # For parallel processing, creating canvases and drawing can be tricky
    # unless ROOT.gROOT.SetBatch(True) is strictly enforced.
    # If plots are needed, it's often better to save histograms and plot them in a separate, sequential script.
    if OUTPUT_DIR_PLOTS:
        c = ROOT.TCanvas(f"c{pid_suffix}", f"Canvas for bin {binNumStr}", 800, 600)
        c.SetLogy()
        c.SetTickx(1)
        c.SetTicky(1)
        # ROOT.gStyle.SetOptStat(0) # Already set or inherited

        hData.Draw("E1")
        hData.GetXaxis().SetRangeUser(4.2, 8.5)
        hData.SetMarkerColor(ROOT.kRed)
        hData.SetMarkerStyle(20)
        hData.SetMarkerSize(1)
        hData.SetLineColor(ROOT.kRed)
        hData.GetXaxis().SetTitle("M (GeV)")
        hData.GetXaxis().CenterTitle(True)
        hData.GetYaxis().SetTitle("M^{3} d^{2}#sigma/dMdx_{F} (nb GeV^{2})")
        hData.GetYaxis().CenterTitle(True)
        hData.SetTitle(f"Differential Cross Section, xF bin {binNumStr}")

        # Theory Curves (open files per call)
        nnpdf4 = None
        nnpdf4_file_obj = None
        try:
            nnpdf4_file_obj = ROOT.TFile.Open(NNPDF4_FILE_PATH, "READ")
            if nnpdf4_file_obj and not nnpdf4_file_obj.IsZombie():
                nnpdf4_graph = nnpdf4_file_obj.Get(f"gr_xFbin{binNumStr}")
                if nnpdf4_graph:
                    nnpdf4 = nnpdf4_graph.Clone(f"nnpdf4_clone{pid_suffix}")
                    nnpdf4.SetMarkerSize(0.75)
                    nnpdf4.SetMarkerStyle(21)
                    nnpdf4.SetMarkerColor(ROOT.kBlue + 2)
                    nnpdf4.SetLineColor(ROOT.kBlue + 2)
                    nnpdf4.SetFillColorAlpha(38, 0.5)
                    nnpdf4.SetFillStyle(3002)
                    nnpdf4.Draw("L3 SAME")
                # else:
                    # print(f"Process {os.getpid()}: Warning - Could not get 'gr_xFbin{binNumStr}' from {NNPDF4_FILE_PATH}")
            if nnpdf4_file_obj: nnpdf4_file_obj.Close()
        except Exception as e_plot: print(f"Process {os.getpid()}: Exception plotting NNPDF4 for bin {binNumStr}: {e_plot}")


        ct18 = None
        ct18_file_obj = None
        try:
            ct18_file_obj = ROOT.TFile.Open(CT18_FILE_PATH, "READ")
            if ct18_file_obj and not ct18_file_obj.IsZombie():
                ct18_graph = ct18_file_obj.Get(f"gr_xFbin{binNumStr}")
                if ct18_graph:
                    ct18 = ct18_graph.Clone(f"ct18_clone{pid_suffix}")
                    ct18.SetMarkerSize(0.75)
                    ct18.SetMarkerStyle(21)
                    ct18.SetMarkerColor(ROOT.kGreen + 2)
                    ct18.SetLineColor(ROOT.kGreen + 2)
                    ct18.SetFillColorAlpha(30, 0.5)
                    ct18.SetFillStyle(3002)
                    ct18.Draw("L3 SAME")
                # else:
                #     print(f"Process {os.getpid()}: Warning - Could not get 'gr_xFbin{binNumStr}' from {CT18_FILE_PATH}")
            if ct18_file_obj: ct18_file_obj.Close()
        except Exception as e_plot: print(f"Process {os.getpid()}: Exception plotting CT18 for bin {binNumStr}: {e_plot}")


        hData.Draw("E1 SAME")

        legend = ROOT.TLegend(0.6, 0.7, 0.85, 0.85)
        legend.SetBorderSize(0)
        legend.AddEntry(hData, "E906", "lep")
        if nnpdf4: legend.AddEntry(nnpdf4, "NNPDF 4.0", "lf")
        if ct18: legend.AddEntry(ct18, "CT18", "lf")
        legend.Draw()

        #plot_filename = f"{OUTPUT_DIR_PLOTS}/LH2_{binNumStr}_cs.png" # Changed to png, eps can be slow
        plot_filename = f"{OUTPUT_DIR_PLOTS}/LH2_{binNumStr}_cs.pdf" # Changed to png, eps can be slow
        c.SaveAs(plot_filename)
        # print(f"Process {os.getpid()}: Saved plot to {plot_filename}")
        # ROOT objects like TCanvas, TLegend, etc., will be garbage collected when this function returns.

    # --- Save final histogram ---
    hAccCor = hData.Clone("hAccCor") # Clone for saving, hData might be drawn on a canvas

    output_hist_filename = f"{OUTPUT_DIR_RESULTS}/LH2_{binNumStr}_updatedPoT.root"
    try:
        savehist_file = ROOT.TFile.Open(output_hist_filename, "RECREATE")
        if savehist_file and not savehist_file.IsZombie():
            hAccCor.Write()
            savehist_file.Close()
            print(f"Process {os.getpid()}: Saved histogram to {output_hist_filename} for bin {binNumStr}")
        else:
            print(f"Process {os.getpid()}: Error - Could not create output histogram file {output_hist_filename} for bin {binNumStr}")
    except Exception as e:
        print(f"Process {os.getpid()}: Exception saving output histogram for bin {binNumStr}: {e}")


    # --- Cleanup TFile objects for this process ---
    for f_obj in open_files:
        if f_obj and f_obj.IsOpen():
            f_obj.Close()
    print(f"Process {os.getpid()}: Finished bin {binNumStr}")


if __name__ == '__main__':
    # Enable ROOT's implicit multi-threading if desired (usually good for I/O bound tasks)
    # For CPU-bound tasks with Python's multiprocessing, it might not always provide
    # additional benefits over process-based parallelism due to the GIL,
    # but for ROOT I/O it can be helpful.
    # ROOT.EnableImplicitMT() # Call this once at the beginning if you want ROOT's internal MT

    # Set ROOT to global batch mode for the main process as well,
    # especially if no interactive plots are expected from the main thread.
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0) # Global style setting

    start_time = time.time()

    # Determine the number of bins to process
    # num_bins_to_process = len(XF_BINS) # or len(K_EFF) if that's the limiting factor
    num_bins_to_process = min(len(XF_BINS), len(K_EFF)) # Process only for defined kEff and xF bins
    
    bin_numbers = list(range(num_bins_to_process))

    # Determine number of processes: use number of CPU cores or a fixed number
    num_processes = multiprocessing.cpu_count()
    # num_processes = 4 # Or set a fixed number
    print(f"Starting parallel processing with {num_processes} processes for {len(bin_numbers)} bins.")

    # Create a Pool of worker processes
    # Using 'spawn' or 'forkserver' can be more stable with ROOT than 'fork' on some systems
    # multiprocessing.set_start_method('spawn', force=True) # Optional: try if 'fork' causes issues

    with multiprocessing.Pool(processes=num_processes) as pool:
        # The map function will distribute the bin_numbers to the diffCross_acc2d_H function
        # Each call to diffCross_acc2d_H(bin_num) runs in a separate process.
        pool.map(diffCross_acc2d_H, bin_numbers)

    end_time = time.time()
    print(f"\nParallel processing complete. Total time: {end_time - start_time:.2f} seconds.")

    # If you were not in batch mode and wanted to show something from the main process:
    # if not ROOT.gROOT.IsBatch():
    # print("\nProcessing complete. Main thread can continue or exit.")
    # ROOT.gApplication.Run() # Keeps the application alive if needed