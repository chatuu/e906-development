"""
analyzer.py
Core Object-Oriented Analysis Module for Drell-Yan Cross-Sections.
"""

import os
import sys
import csv
import math
import uproot
import numpy as np
import ROOT
import config

class DYCrossSectionAnalyzer:
    """
    Main Object-Oriented Analyzer for extracting kinematics, performing background 
    subtractions, and calculating absolute Drell-Yan cross-sections.
    """

    def __init__(self, lh2_files, ld2_files, flask_files, out_filename="All_XSec_Objects.root"):
        self._setup_root()
        
        # Convert inputs to lists if single strings are passed
        self.lh2_paths = lh2_files if isinstance(lh2_files, list) else [lh2_files]
        self.ld2_paths = ld2_files if isinstance(ld2_files, list) else [ld2_files]
        self.flask_paths = flask_files if isinstance(flask_files, list) else [flask_files]
        
        self.out_filename = out_filename
        self.out_file = ROOT.TFile(out_filename, "RECREATE")
        
        # Load interpolation map
        try:
            npz_data = np.load(config.INPUT_NPZ_FILE)
            self.x_curve = npz_data['x']
        except Exception as e:
            print(f"Error loading NPZ file for Covariance generation at '{config.INPUT_NPZ_FILE}': {e}")
            sys.exit(1)

        # Output dictionaries to store histograms between stages
        self.hists_lh2 = None
        self.hists_ld2 = None
        self.hists_fl = None
        self.sub_dict_lh2 = None
        self.sub_dict_pd = None

    def _setup_root(self):
        """Configures global ROOT visual options and mutes console output."""
        ROOT.gROOT.SetBatch(True)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kBird)
        ROOT.gStyle.SetEndErrorSize(5) # Globally sets the width of the error bar caps (edges)
        
        # Mute all ROOT terminal output (Info and Warnings) to protect the progress bars
        ROOT.gErrorIgnoreLevel = ROOT.kFatal

    @staticmethod
    def get_or_create_dir(base_dir, name):
        """Retrieves an existing TDirectory or creates a new one."""
        d = base_dir.GetDirectory(name)
        if not d:
            d = base_dir.mkdir(name)
        return d

    @staticmethod
    def apply_cuts(tree, cut=4.2):
        """Applies standard physics cuts to the given TTree arrays in numpy."""
        events = tree.arrays(library="np")
        class EventNamespace:
            def __init__(self, data):
                self.__dict__.update(data)
        e = EventNamespace(events)

        bo = np.where(e.runID >= 11000, 1.6, 0.4)

        dimuon_cut = (
            (np.abs(e.dx) < 0.25) & (np.abs(e.dy - bo) < 0.22) &
            (e.dz < -5.) & (e.dz > -280.) & (np.abs(e.dpx) < 1.8) & (np.abs(e.dpy) < 2.0) &
            (e.dpx * e.dpx + e.dpy * e.dpy < 5.) & (e.dpz < 116.) & (e.dpz > 38.) &
            (e.mass > cut) & (e.mass < 8.8) &
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

        total_cut_mask = (track1_cut & track2_cut & tracks_cut & dimuon_cut & occ_cut)

        filtered_events = {}
        for key, val in events.items():
            filtered_events[key] = val[total_cut_mask]
            
        # Define pT based on standard kinematics: pT = sqrt(px^2 + py^2)
        filtered_events["pT"] = np.sqrt(filtered_events["dpx"]**2 + filtered_events["dpy"]**2)
            
        return filtered_events

    def get_concatenated_events(self, file_paths, tree_name):
        """Iterates through file paths, applies cuts, and merges the numpy arrays in memory."""
        all_events = None
        for fp in file_paths:
            if not os.path.exists(fp):
                print(f"Warning: File {fp} does not exist, skipping.")
                continue
            try:
                with uproot.open(fp) as f:
                    if tree_name not in f:
                        print(f"Warning: Tree '{tree_name}' not found in {fp}, skipping.")
                        continue
                    filtered = self.apply_cuts(f[tree_name])
                    if all_events is None:
                        all_events = {k: [v] for k, v in filtered.items()}
                    else:
                        for k, v in filtered.items():
                            if k in all_events:
                                all_events[k].append(v)
            except Exception as e:
                print(f"Error reading {fp}: {e}")
                
        if all_events is None:
            raise RuntimeError(f"No valid data found for tree '{tree_name}' in provided files.")
        
        return {k: np.concatenate(v) for k, v in all_events.items()}

    @staticmethod
    def create_histograms(file_label):
        hists = {}
        
        def make_th1(name, title, y_title=""):
            h = ROOT.TH1D(f"{name}_{file_label}", f"{title} ({file_label});p_{{T}} [GeV];{y_title}", 
                          len(config.PT_BINS)-1, config.PT_BINS)
            h.Sumw2()
            h.SetStats(0)
            h.GetXaxis().CenterTitle()
            h.GetYaxis().CenterTitle()
            h.GetXaxis().SetTitleOffset(1.3)
            h.GetYaxis().SetTitleOffset(1.7)
            return h

        hists["Y_total"] = make_th1("Y_total", "Total Yield (result)", "Yield")
        hists["Y_mix"] = make_th1("Y_mix", "Mix Yield (result_mix)", "Yield")
        hists["E_total_reco"] = make_th1("E_total_reco", "Avg Reco Eff (Total)", "Efficiency")
        hists["E_mix_reco"] = make_th1("E_mix_reco", "Avg Reco Eff (Mix)", "Efficiency")
        hists["E_total_hodo"] = make_th1("E_total_hodo", "Avg Hodo Eff (Total)", "Efficiency")
        hists["E_mix_hodo"] = make_th1("E_mix_hodo", "Avg Hodo Eff (Mix)", "Efficiency")
        hists["E_total_final"] = make_th1("E_total_final", "Avg Final Eff (Total)", "Efficiency")
        hists["E_mix_final"] = make_th1("E_mix_final", "Avg Final Eff (Mix)", "Efficiency")
        hists["E_final_signal"] = make_th1("E_final_signal", "Avg Signal Efficiency", "Efficiency")
        
        hists["Y_corrected"] = make_th1("Y_corrected", "Corrected Yield (Total Error)", "Yield")
        hists["Y_corrected_stat"] = make_th1("Y_corrected_stat", "Corrected Yield (Stat Error)", "Yield")
        hists["Y_corrected_sys"] = make_th1("Y_corrected_sys", "Corrected Yield (Sys Error)", "Yield")
        hists["Mass_Centroid"] = make_th1("Mass_Centroid", "Data-Driven Mass Centroid", "Mass [GeV]")

        return hists

    @staticmethod
    def get_weighted_mean_and_error(eff_arr, err_arr):
        if len(eff_arr) == 0: return 0.0, 0.0
        mean = np.mean(eff_arr)
        err_on_mean = np.sqrt(np.sum(err_arr**2)) / len(eff_arr)
        return mean, err_on_mean

    @staticmethod
    def get_correlated_mean_and_error(eff_arr, err_arr, loc_arr):
        N = len(eff_arr)
        print(loc_arr)
        if N == 0: return 0.0, 0.0
        mean_eff = np.mean(eff_arr)
        
        diff_matrix = np.abs(loc_arr[:, None] - loc_arr[None, :])
        correl_matrix = np.zeros((N, N))
        correl_matrix[diff_matrix == 0] = 1.0  
        correl_matrix[diff_matrix == 1] = 1.0  
        
        covar_matrix = correl_matrix * np.outer(err_arr, err_arr)
        err_on_mean = np.sqrt(np.sum(covar_matrix)) / N
        return mean_eff, err_on_mean

    @staticmethod
    def add_latex_to_bin(hist, x_center, value, error):
        val_f = float(value)
        err_f = float(error)
        if np.isnan(val_f) or np.isinf(val_f): val_f = 0.0
        if np.isnan(err_f) or np.isinf(err_f): err_f = 0.0
        
        latex_text = f"#splitline{{{val_f:.4f}}}{{#pm {err_f:.4f}}}"
        
        # Determine the top of the point + error, and add an 8% padding
        offset = abs(val_f) * 0.08
        if offset == 0: offset = 0.005
        y_pos = val_f + err_f + offset
        
        l = ROOT.TLatex(x_center, y_pos, latex_text)
        l.SetTextSize(0.022)
        l.SetTextAlign(21) # Bottom center alignment pushes text upward
        l.SetTextColor(ROOT.kBlack)
        
        # Appending it to the functions list instructs ROOT to dynamically draw this on the pad
        hist.GetListOfFunctions().Add(l)

    def extract_target_stats(self, data_tot, data_mix, loc_tot, loc_mix, m_low, m_high, pt_low, pt_high):
        mask_tot = (data_tot["mass"] >= m_low) & (data_tot["mass"] < m_high) & \
                   (data_tot["pT"] >= pt_low) & (data_tot["pT"] < pt_high)
        mask_mix = (data_mix["mass"] >= m_low) & (data_mix["mass"] < m_high) & \
                   (data_mix["pT"] >= pt_low) & (data_mix["pT"] < pt_high)

        N_tot = np.sum(mask_tot)
        N_mix = np.sum(mask_mix)
        sum_mass_tot = np.sum(data_tot["mass"][mask_tot]) if N_tot > 0 else 0.0
        sum_mass_mix = np.sum(data_mix["mass"][mask_mix]) if N_mix > 0 else 0.0

        eff_reco_tot = data_tot["recoeff"][mask_tot]
        err_reco_tot = data_tot["recoeff_error"][mask_tot]
        eff_hodo_tot = data_tot["hodoeff"][mask_tot]
        err_hodo_tot = data_tot["hodoeff_error"][mask_tot]

        eff_reco_mix = data_mix["recoeff"][mask_mix]
        err_reco_mix = data_mix["recoeff_error"][mask_mix]
        eff_hodo_mix = data_mix["hodoeff"][mask_mix]
        err_hodo_mix = data_mix["hodoeff_error"][mask_mix]

        means = {
            "r_tot": self.get_correlated_mean_and_error(eff_reco_tot, err_reco_tot, loc_tot[mask_tot]),
            "h_tot": self.get_weighted_mean_and_error(eff_hodo_tot, err_hodo_tot),
            "r_mix": self.get_correlated_mean_and_error(eff_reco_mix, err_reco_mix, loc_mix[mask_mix]),
            "h_mix": self.get_weighted_mean_and_error(eff_hodo_mix, err_hodo_mix),
        }

        f_tot_mean = means["r_tot"][0] * means["h_tot"][0]
        f_tot_err = np.sqrt((means["h_tot"][0] * means["r_tot"][1])**2 + (means["r_tot"][0] * means["h_tot"][1])**2)
        means["f_tot"] = (f_tot_mean, f_tot_err)

        f_mix_mean = means["r_mix"][0] * means["h_mix"][0]
        f_mix_err = np.sqrt((means["h_mix"][0] * means["r_mix"][1])**2 + (means["r_mix"][0] * means["h_mix"][1])**2)
        means["f_mix"] = (f_mix_mean, f_mix_err)

        val_sig_eff, err_sig_eff = 0.0, 0.0
        diff_yield = N_tot - N_mix
        if diff_yield != 0:
            numerator = (N_tot * means["f_tot"][0]) - (N_mix * means["f_mix"][0])
            val_sig_eff = numerator / diff_yield
            term1 = (N_tot * means["f_tot"][1])**2
            term2 = (N_mix * means["f_mix"][1])**2
            err_sig_eff = (1.0 / diff_yield) * np.sqrt(term1 + term2)

        return N_tot, N_mix, sum_mass_tot, sum_mass_mix, val_sig_eff, err_sig_eff, diff_yield, means

    def fill_histograms(self, hists, root_x, pt_center, N_tot, N_mix, sig_eff, err_sig_eff, diff_yield, means, final_centroid):
        err_N_tot = np.sqrt(N_tot)
        err_N_mix = np.sqrt(N_mix)

        val_corr_yield, err_corr_stat, err_corr_sys, err_corr_total = 0.0, 0.0, 0.0, 0.0
        if sig_eff > 0 and diff_yield > 0:
            val_corr_yield = diff_yield / sig_eff
            err_corr_stat = np.sqrt(N_tot + N_mix) / sig_eff
            err_corr_sys = val_corr_yield * (err_sig_eff / sig_eff)
            err_corr_total = np.sqrt(err_corr_stat**2 + err_corr_sys**2)

        def set_bin(key, val, err):
            hists[key].SetBinContent(root_x, val)
            hists[key].SetBinError(root_x, err)

        set_bin("Y_total", N_tot, err_N_tot)
        set_bin("Y_mix", N_mix, err_N_mix)
        set_bin("E_total_reco", means["r_tot"][0], means["r_tot"][1])
        set_bin("E_mix_reco", means["r_mix"][0], means["r_mix"][1])
        set_bin("E_total_hodo", means["h_tot"][0], means["h_tot"][1])
        set_bin("E_mix_hodo", means["h_mix"][0], means["h_mix"][1])
        set_bin("E_total_final", means["f_tot"][0], means["f_tot"][1])
        set_bin("E_mix_final", means["f_mix"][0], means["f_mix"][1])
        set_bin("E_final_signal", sig_eff, err_sig_eff)
        set_bin("Y_corrected", val_corr_yield, err_corr_total)
        set_bin("Y_corrected_stat", val_corr_yield, err_corr_stat)
        set_bin("Y_corrected_sys", val_corr_yield, err_corr_sys)
        
        hists["Mass_Centroid"].SetBinContent(root_x, final_centroid)

    def save_1d_pdfs(self, hists_dict, label):
        """Saves all generated TH1D histograms as separate PDF files."""
        for name, hist in hists_dict.items():
            if "stat" in name or "sys" in name or "Centroid" in name: continue 
            c = ROOT.TCanvas(f"c_{name}_{label}", "", 1200, 900)
            c.SetRightMargin(0.05)
            c.SetLeftMargin(0.16)
            c.SetBottomMargin(0.14)
            c.SetTickx(1)
            c.SetTicky(1)
            
            # Find true maximums to comfortably fit TLatex splitline text
            max_y = -1e9
            min_y = 1e9
            for i in range(1, hist.GetNbinsX() + 1):
                val = hist.GetBinContent(i)
                err = hist.GetBinError(i)
                if val == 0 and err == 0: continue
                if val + err > max_y: max_y = val + err
                if val - err < min_y: min_y = val - err
                
            if max_y > -1e8:
                hist.SetMaximum(max_y * 1.5) 
                if min_y < 0:
                    hist.SetMinimum(min_y * 1.2)
                else:
                    hist.SetMinimum(0.0)
            else:
                hist.SetMinimum(0.0)
                hist.SetMaximum(1.0)
            
            hist.SetLineColor(ROOT.kBlue)
            hist.SetLineWidth(2)
            hist.Draw("HIST E1")
            
            # Explicitly draw TLatex objects to avoid garbage collection
            latex_draws = []
            offset = (max_y if max_y > -1e8 else 1.0) * 0.08
            for i in range(1, hist.GetNbinsX() + 1):
                val = hist.GetBinContent(i)
                err = hist.GetBinError(i)
                if val == 0 and err == 0: continue
                pt_center = hist.GetBinCenter(i)
                y_pos = val + err + offset
                txt = ROOT.TLatex(pt_center, y_pos, f"#splitline{{{val:.4f}}}{{#pm {err:.4f}}}")
                txt.SetTextSize(0.022)
                txt.SetTextAlign(21)
                txt.SetTextColor(ROOT.kBlack)
                txt.Draw()
                latex_draws.append(txt)

            c.SaveAs(f"{name}_{label}.pdf")
            c.Close()

    def generate_efficiency_csv(self, data_tot, loc_tot, m_low, m_high, pt_low, pt_high, pt_idx, target_label):
        """Filters data for the bin and generates a detailed CSV of event efficiencies."""
        mask = (data_tot["mass"] >= m_low) & (data_tot["mass"] < m_high) & \
               (data_tot["pT"] >= pt_low) & (data_tot["pT"] < pt_high)
        
        reco = data_tot["recoeff"][mask]
        reco_err = data_tot["recoeff_error"][mask]
        hodo = data_tot["hodoeff"][mask]
        hodo_err = data_tot["hodoeff_error"][mask]
        locs = loc_tot[mask]
        
        filename = f"reco_pT_{pt_idx}_{target_label}.csv"
        
        N = len(reco)
        if N == 0:
            avg_reco, reco_err_uncorr, reco_err_corr = 0.0, 0.0, 0.0
            avg_hodo, hodo_err_stat = 0.0, 0.0
            total_eff, total_eff_err = 0.0, 0.0
        else:
            avg_reco, reco_err_uncorr = self.get_weighted_mean_and_error(reco, reco_err)
            _, reco_err_corr = self.get_correlated_mean_and_error(reco, reco_err, locs)
            avg_hodo, hodo_err_stat = self.get_weighted_mean_and_error(hodo, hodo_err)
            
            total_eff = avg_reco * avg_hodo
            total_eff_err = np.sqrt((avg_hodo * reco_err_corr)**2 + (avg_reco * hodo_err_stat)**2)
            
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["recoeff", "recoeff_error", "hodoeff", "hodoeff_error"])
            for i in range(N):
                writer.writerow([reco[i], reco_err[i], hodo[i], hodo_err[i]])
                
            writer.writerow([])
            writer.writerow(["average recoeff", avg_reco])
            writer.writerow(["propagated recoeff_error (without correlations)", reco_err_uncorr])
            writer.writerow(["propagated recoeff_error (with correlations)", reco_err_corr])
            writer.writerow(["average hodoeff", avg_hodo])
            writer.writerow(["propagated hodoeff_error", hodo_err_stat])
            writer.writerow(["total_eff = average recoeff (with correlations) * average hodoeff", total_eff])
            writer.writerow(["total_eff_error = propagated error of total_eff", total_eff_err])

    def process_kinematics(self):
        """Processes uproot data, creates 1D histograms vs pT, and saves PDFs."""
        def get_loc_data(data_dict, is_mix=False):
            if is_mix and 'ptrk_D1' in data_dict and 'ntrk_D1' in data_dict:
                d1_vals = 0.5 * (data_dict['ptrk_D1'] + data_dict['ntrk_D1'])
            else:
                d1_vals = data_dict['D1']
            return np.digitize(d1_vals, self.x_curve) - 1
        
        # Iteratively load, cut, and concatenate all files
        data_lh2_tot = self.get_concatenated_events(self.lh2_paths, "result")
        data_lh2_mix = self.get_concatenated_events(self.lh2_paths, "result_mix")
        
        data_ld2_tot = self.get_concatenated_events(self.ld2_paths, "result")
        data_ld2_mix = self.get_concatenated_events(self.ld2_paths, "result_mix")
        
        data_fl_tot = self.get_concatenated_events(self.flask_paths, "result")
        data_fl_mix = self.get_concatenated_events(self.flask_paths, "result_mix")
        
        loc_lh2_tot = get_loc_data(data_lh2_tot)
        loc_lh2_mix = get_loc_data(data_lh2_mix, True)
        
        loc_ld2_tot = get_loc_data(data_ld2_tot)
        loc_ld2_mix = get_loc_data(data_ld2_mix, True)
        
        loc_fl_tot = get_loc_data(data_fl_tot)
        loc_fl_mix = get_loc_data(data_fl_mix, True)

        dir_kin = self.get_or_create_dir(self.out_file, "Kinematics")

        dir_lh2 = self.get_or_create_dir(dir_kin, "LH2")
        dir_lh2.cd()
        self.hists_lh2 = self.create_histograms("LH2")

        dir_ld2 = self.get_or_create_dir(dir_kin, "LD2")
        dir_ld2.cd()
        self.hists_ld2 = self.create_histograms("LD2")

        dir_fl = self.get_or_create_dir(dir_kin, "Flask")
        dir_fl.cd()
        self.hists_fl = self.create_histograms("Flask")

        dir_csv = self.get_or_create_dir(self.out_file, "CSV_Tables")
        csv_filename_lh2 = "Table_Kinematics_LH2_vs_pT.csv"
        csv_filename_ld2 = "Table_Kinematics_LD2_vs_pT.csv"
        
        csv_rows_lh2 = []
        csv_rows_ld2 = []

        m_low, m_high = config.MASS_BINS[0], config.MASS_BINS[-1]
        m_center_overall = (m_low + m_high) / 2.0

        for i_pt in range(len(config.PT_BINS) - 1):
            pt_low, pt_high = config.PT_BINS[i_pt], config.PT_BINS[i_pt+1]
            pt_center = (pt_low + pt_high) / 2.0
            root_x = i_pt + 1
            
            # Generate the specific event efficiency CSVs requested
            self.generate_efficiency_csv(data_lh2_tot, loc_lh2_tot, m_low, m_high, pt_low, pt_high, i_pt, "LH2")
            self.generate_efficiency_csv(data_ld2_tot, loc_ld2_tot, m_low, m_high, pt_low, pt_high, i_pt, "LD2")
            self.generate_efficiency_csv(data_fl_tot, loc_fl_tot, m_low, m_high, pt_low, pt_high, i_pt, "Flask")

            N_l2_t, N_l2_m, M_l2_t, M_l2_m, eps_ld2, err_eps_ld2, diff_l2, mean_l2 = self.extract_target_stats(data_ld2_tot, data_ld2_mix, loc_ld2_tot, loc_ld2_mix, m_low, m_high, pt_low, pt_high)
            N_lh_t, N_lh_m, M_lh_t, M_lh_m, eps_lh2, err_eps_lh2, diff_lh, mean_lh = self.extract_target_stats(data_lh2_tot, data_lh2_mix, loc_lh2_tot, loc_lh2_mix, m_low, m_high, pt_low, pt_high)
            N_fl_t, N_fl_m, M_fl_t, M_fl_m, eps_fl, err_eps_fl, diff_fl, mean_fl = self.extract_target_stats(data_fl_tot, data_fl_mix, loc_fl_tot, loc_fl_mix, m_low, m_high, pt_low, pt_high)

            if eps_lh2 > 0:
                cY_lh_t, cM_lh_t = N_lh_t / eps_lh2, M_lh_t / eps_lh2
                cY_lh_m, cM_lh_m = N_lh_m / eps_lh2, M_lh_m / eps_lh2
            else:
                cY_lh_t, cM_lh_t, cY_lh_m, cM_lh_m = 0.0, 0.0, 0.0, 0.0

            if eps_ld2 > 0:
                cY_l2_t, cM_l2_t = N_l2_t / eps_ld2, M_l2_t / eps_ld2
                cY_l2_m, cM_l2_m = N_l2_m / eps_ld2, M_l2_m / eps_ld2
            else:
                cY_l2_t, cM_l2_t, cY_l2_m, cM_l2_m = 0.0, 0.0, 0.0, 0.0

            if eps_fl > 0:
                cY_fl_t_lh, cM_fl_t_lh = (config.FLASK_NORM_LH2 * N_fl_t) / eps_fl, (config.FLASK_NORM_LH2 * M_fl_t) / eps_fl
                cY_fl_m_lh, cM_fl_m_lh = (config.FLASK_NORM_LH2 * N_fl_m) / eps_fl, (config.FLASK_NORM_LH2 * M_fl_m) / eps_fl
                cY_fl_t_l2, cM_fl_t_l2 = (config.FLASK_NORM_LD2 * N_fl_t) / eps_fl, (config.FLASK_NORM_LD2 * M_fl_t) / eps_fl
                cY_fl_m_l2, cM_fl_m_l2 = (config.FLASK_NORM_LD2 * N_fl_m) / eps_fl, (config.FLASK_NORM_LD2 * M_fl_m) / eps_fl
            else:
                cY_fl_t_lh, cM_fl_t_lh, cY_fl_m_lh, cM_fl_m_lh = 0.0, 0.0, 0.0, 0.0
                cY_fl_t_l2, cM_fl_t_l2, cY_fl_m_l2, cM_fl_m_l2 = 0.0, 0.0, 0.0, 0.0

            num_LH2 = (cM_lh_t - cM_lh_m) - (cM_fl_t_lh - cM_fl_m_lh)
            den_LH2 = (cY_lh_t - cY_lh_m) - (cY_fl_t_lh - cY_fl_m_lh)
            cent_LH2 = num_LH2 / den_LH2 if den_LH2 != 0 else m_center_overall

            num_LD2 = (cM_l2_t - cM_l2_m) - (cM_fl_t_l2 - cM_fl_m_l2) - (config.LH2_TO_LD2_NORM * num_LH2)
            den_LD2 = (cY_l2_t - cY_l2_m) - (cY_fl_t_l2 - cY_fl_m_l2) - (config.LH2_TO_LD2_NORM * den_LH2)
            cent_LD2 = num_LD2 / den_LD2 if den_LD2 != 0 else m_center_overall

            self.fill_histograms(self.hists_ld2, root_x, pt_center, N_l2_t, N_l2_m, eps_ld2, err_eps_ld2, diff_l2, mean_l2, cent_LD2)
            self.fill_histograms(self.hists_lh2, root_x, pt_center, N_lh_t, N_lh_m, eps_lh2, err_eps_lh2, diff_lh, mean_lh, cent_LH2)
            num_fl = (M_fl_t - M_fl_m) / (eps_fl if eps_fl > 0 else 1e-9)
            den_fl = (N_fl_t - N_fl_m) / (eps_fl if eps_fl > 0 else 1e-9)
            cent_fl = num_fl / den_fl if den_fl != 0 else m_center_overall
            self.fill_histograms(self.hists_fl, root_x, pt_center, N_fl_t, N_fl_m, eps_fl, err_eps_fl, diff_fl, mean_fl, cent_fl)

            csv_rows_lh2.append({
                "pT Bin": f"[{pt_low:.2f}, {pt_high:.2f})",
                "pT Center": pt_center,
                "LH2 Mass Centroid": cent_LH2,
                "N_LH2_total": N_lh_t,
                "N_LH2_mixed": N_lh_m,
                "N_flask_total": N_fl_t,
                "N_flask_mixed": N_fl_m,
                "eps_LH2": eps_lh2,
                "eps_LD2": eps_ld2,
                "eps_flask": eps_fl,
                "Corrected Total Mass LH2 (num)": num_LH2,
                "Corrected Yield LH2 (denom)": den_LH2
            })

            csv_rows_ld2.append({
                "pT Bin": f"[{pt_low:.2f}, {pt_high:.2f})",
                "pT Center": pt_center,
                "LD2 Mass Centroid": cent_LD2,
                "N_LD2_total": N_l2_t,
                "N_LD2_mixed": N_l2_m,
                "N_LH2_total": N_lh_t,
                "N_LH2_mixed": N_lh_m,
                "N_flask_total": N_fl_t,
                "N_flask_mixed": N_fl_m,
                "eps_LH2": eps_lh2,
                "eps_LD2": eps_ld2,
                "eps_flask": eps_fl,
                "Corrected Total Mass LD2 (num)": num_LD2,
                "Corrected Yield LD2 (denom)": den_LD2
            })

        with open(csv_filename_lh2, "w", newline='') as f:
            writer = csv.DictWriter(f, fieldnames=csv_rows_lh2[0].keys())
            writer.writeheader()
            writer.writerows(csv_rows_lh2)
            
        with open(csv_filename_ld2, "w", newline='') as f:
            writer = csv.DictWriter(f, fieldnames=csv_rows_ld2[0].keys())
            writer.writeheader()
            writer.writerows(csv_rows_ld2)
            
        dir_csv.cd()
        macro_lh2 = ROOT.TMacro(csv_filename_lh2)
        macro_lh2.Write(csv_filename_lh2)
        
        macro_ld2 = ROOT.TMacro(csv_filename_ld2)
        macro_ld2.Write(csv_filename_ld2)

        self.save_1d_pdfs(self.hists_lh2, "LH2")
        self.save_1d_pdfs(self.hists_ld2, "LD2")
        self.save_1d_pdfs(self.hists_fl, "Flask")

    def generate_subtracted_plot(self, hists_target, hists_flask, flask_norm, target_label):
        """Generates and saves 1D subtracted yield plot PDFs for LH2."""
        if "Y_corrected_stat" not in hists_target or "Y_corrected_stat" not in hists_flask:
            print("Error: Y_corrected_stat histogram missing.")
            return None

        dir_sub = self.get_or_create_dir(self.out_file, f"Subtracted_Plots_{target_label}")
        dir_sub.cd()
        
        h_target_stat = hists_target["Y_corrected_stat"]
        h_target_sys = hists_target["Y_corrected_sys"]
        h_flask_stat = hists_flask["Y_corrected_stat"]
        h_flask_sys = hists_flask["Y_corrected_sys"]
        
        name = f"Y_corrected_Subtracted_{target_label}"
        title = f"Corrected Yield ({target_label} - Flask)"
        h_sub = ROOT.TH1D(name, f"{title};p_{{T}} [GeV];Yield", 
                          len(config.PT_BINS)-1, config.PT_BINS)
        h_sub.Sumw2()
        h_sub.SetStats(0)
        h_sub.GetXaxis().CenterTitle()
        h_sub.GetYaxis().CenterTitle()
        h_sub.GetXaxis().SetTitleOffset(1.3)
        h_sub.GetYaxis().SetTitleOffset(1.7)

        h_sub_stat = h_sub.Clone(f"{name}_stat")
        h_sub_sys = h_sub.Clone(f"{name}_sys")
        h_sub_centroid = hists_target["Mass_Centroid"].Clone(f"{name}_Mass_Centroid")

        for i_pt in range(len(config.PT_BINS) - 1):
            bin_x = i_pt + 1
            pt_center = (config.PT_BINS[i_pt] + config.PT_BINS[i_pt+1]) / 2.0
            
            y_target = h_target_stat.GetBinContent(bin_x)
            e_target_stat = h_target_stat.GetBinError(bin_x)
            e_target_sys = h_target_sys.GetBinError(bin_x)
            
            y_flask = h_flask_stat.GetBinContent(bin_x)
            e_flask_stat = h_flask_stat.GetBinError(bin_x)
            e_flask_sys = h_flask_sys.GetBinError(bin_x)
            
            val_sub = y_target - (flask_norm * y_flask)
            
            err_sub_stat = np.sqrt(e_target_stat**2 + (flask_norm * e_flask_stat)**2)
            err_sub_sys = np.sqrt(e_target_sys**2 + (flask_norm * e_flask_sys)**2)
            err_sub_total = np.sqrt(err_sub_stat**2 + err_sub_sys**2)
            
            h_sub.SetBinContent(bin_x, val_sub)
            h_sub.SetBinError(bin_x, err_sub_total)

            h_sub_stat.SetBinContent(bin_x, val_sub)
            h_sub_stat.SetBinError(bin_x, err_sub_stat)

            h_sub_sys.SetBinContent(bin_x, val_sub)
            h_sub_sys.SetBinError(bin_x, err_sub_sys)

        c = ROOT.TCanvas(f"c_{name}", f"c_{name}", 1200, 900)
        c.SetRightMargin(0.05); c.SetLeftMargin(0.16); c.SetBottomMargin(0.14)
        c.SetTickx(1); c.SetTicky(1)
        
        max_y = -1e9
        min_y = 1e9
        for ix in range(1, h_sub.GetNbinsX() + 1):
            val = h_sub.GetBinContent(ix)
            err = h_sub.GetBinError(ix)
            if val == 0 and err == 0: continue
            if val + err > max_y: max_y = val + err
            if val - err < min_y: min_y = val - err
        
        if max_y > -1e8:
            h_sub.SetMaximum(max_y * 1.5)
            h_sub.SetMinimum(0.0 if min_y >= 0 else min_y * 1.2)
        else:
            h_sub.SetMinimum(0.0)
            h_sub.SetMaximum(1.0)
            
        h_sub.SetLineColor(ROOT.kBlack)
        h_sub.SetMarkerStyle(20)
        h_sub.Draw("P E1")
        
        latex_draws = []
        offset = (max_y if max_y > -1e8 else 1.0) * 0.08
        for ix in range(1, h_sub.GetNbinsX() + 1):
            val = h_sub.GetBinContent(ix)
            err = h_sub.GetBinError(ix)
            if val == 0 and err == 0: continue
            pt_center = h_sub.GetBinCenter(ix)
            y_pos = val + err + offset
            txt = ROOT.TLatex(pt_center, y_pos, f"#splitline{{{val:.4f}}}{{#pm {err:.4f}}}")
            txt.SetTextSize(0.022)
            txt.SetTextAlign(21)
            txt.SetTextColor(ROOT.kBlack)
            txt.Draw()
            latex_draws.append(txt)

        c.SaveAs(f"{name}.pdf")
        c.Close()
        
        return {"stat": h_sub_stat, "sys": h_sub_sys, "centroid": h_sub_centroid}

    def generate_pd_subtracted_plot(self, hists_ld2, hists_lh2, hists_flask):
        """Generates and saves 1D subtracted yield plot PDFs for LD2 (pd equation)."""
        dir_sub = self.get_or_create_dir(self.out_file, "Subtracted_Plots_LD2_pd")
        dir_sub.cd()

        h_ld2_stat = hists_ld2["Y_corrected_stat"]
        h_ld2_sys = hists_ld2["Y_corrected_sys"]
        
        h_lh2_stat = hists_lh2["Y_corrected_stat"]
        h_lh2_sys = hists_lh2["Y_corrected_sys"]
        
        h_flask_stat = hists_flask["Y_corrected_stat"]
        h_flask_sys = hists_flask["Y_corrected_sys"]

        c_lh2 = config.THD_THH_RATIO * (config.PROTONS_ON_TARGET_LD2 / config.PROTONS_ON_TARGET_LH2)
        c_flask_sub = config.FLASK_NORM_LD2 - (c_lh2 * config.FLASK_NORM_LH2)

        name = "Y_corrected_Subtracted_LD2"
        title = "Corrected Yield (LD2 - LH2 - Flask)"
        h_pd = ROOT.TH1D(name, f"{title};p_{{T}} [GeV];Yield", 
                          len(config.PT_BINS)-1, config.PT_BINS)
        h_pd.Sumw2()
        h_pd.SetStats(0)
        h_pd.GetXaxis().CenterTitle()
        h_pd.GetYaxis().CenterTitle()
        h_pd.GetXaxis().SetTitleOffset(1.3)
        h_pd.GetYaxis().SetTitleOffset(1.7)

        h_pd_stat = h_pd.Clone(f"{name}_stat")
        h_pd_sys = h_pd.Clone(f"{name}_sys")
        h_pd_centroid = hists_ld2["Mass_Centroid"].Clone(f"{name}_Mass_Centroid")

        for i_pt in range(len(config.PT_BINS) - 1):
            bin_x = i_pt + 1
            pt_center = (config.PT_BINS[i_pt] + config.PT_BINS[i_pt+1]) / 2.0
            
            y_ld2 = h_ld2_stat.GetBinContent(bin_x)
            y_lh2 = h_lh2_stat.GetBinContent(bin_x)
            y_flask = h_flask_stat.GetBinContent(bin_x)

            e_ld2_stat = h_ld2_stat.GetBinError(bin_x)
            e_lh2_stat = h_lh2_stat.GetBinError(bin_x)
            e_flask_stat = h_flask_stat.GetBinError(bin_x)
            
            e_ld2_sys = h_ld2_sys.GetBinError(bin_x)
            e_lh2_sys = h_lh2_sys.GetBinError(bin_x)
            e_flask_sys = h_flask_sys.GetBinError(bin_x)
            
            val_pd = y_ld2 - (c_lh2 * y_lh2) - (c_flask_sub * y_flask)
            
            err_pd_stat = np.sqrt(e_ld2_stat**2 + (c_lh2 * e_lh2_stat)**2 + (c_flask_sub * e_flask_stat)**2)
            err_pd_sys = np.sqrt(e_ld2_sys**2 + (c_lh2 * e_lh2_sys)**2 + (c_flask_sub * e_flask_sys)**2)
            err_pd_total = np.sqrt(err_pd_stat**2 + err_pd_sys**2)
            
            h_pd.SetBinContent(bin_x, val_pd)
            h_pd.SetBinError(bin_x, err_pd_total)

            h_pd_stat.SetBinContent(bin_x, val_pd)
            h_pd_stat.SetBinError(bin_x, err_pd_stat)

            h_pd_sys.SetBinContent(bin_x, val_pd)
            h_pd_sys.SetBinError(bin_x, err_pd_sys)

        c = ROOT.TCanvas(f"c_{name}", f"c_{name}", 1200, 900)
        c.SetRightMargin(0.05); c.SetLeftMargin(0.16); c.SetBottomMargin(0.14)
        c.SetTickx(1); c.SetTicky(1)
        
        max_y = -1e9
        min_y = 1e9
        for ix in range(1, h_pd.GetNbinsX() + 1):
            val = h_pd.GetBinContent(ix)
            err = h_pd.GetBinError(ix)
            if val == 0 and err == 0: continue
            if val + err > max_y: max_y = val + err
            if val - err < min_y: min_y = val - err
        
        if max_y > -1e8:
            h_pd.SetMaximum(max_y * 1.5)
            h_pd.SetMinimum(0.0 if min_y >= 0 else min_y * 1.2)
        else:
            h_pd.SetMinimum(0.0)
            h_pd.SetMaximum(1.0)
            
        h_pd.SetLineColor(ROOT.kBlack)
        h_pd.SetMarkerStyle(20)
        h_pd.Draw("P E1")
        
        latex_draws = []
        offset = (max_y if max_y > -1e8 else 1.0) * 0.08
        for ix in range(1, h_pd.GetNbinsX() + 1):
            val = h_pd.GetBinContent(ix)
            err = h_pd.GetBinError(ix)
            if val == 0 and err == 0: continue
            pt_center = h_pd.GetBinCenter(ix)
            y_pos = val + err + offset
            txt = ROOT.TLatex(pt_center, y_pos, f"#splitline{{{val:.4f}}}{{#pm {err:.4f}}}")
            txt.SetTextSize(0.022)
            txt.SetTextAlign(21)
            txt.SetTextColor(ROOT.kBlack)
            txt.Draw()
            latex_draws.append(txt)

        c.SaveAs(f"{name}.pdf")
        c.Close()

        return {"stat": h_pd_stat, "sys": h_pd_sys, "centroid": h_pd_centroid}

    def calculate_and_plot_cross_section(self, h_sub_dict, target_label, global_constant):
        """Builds single differential cross-sections vs pT using 1D Acceptance Corrections."""
        acc_path = "/root/github/e906-development/src/AcceptanceCorrection/acceptance_mass_xF.root"
        psip_path = "All_PsiP_Contaminations_pT.root" 
        
        if not os.path.exists(acc_path):
            raise FileNotFoundError(f"CRITICAL ERROR: Acceptance file '{acc_path}' not found!")

        try:
            acc_file = ROOT.TFile.Open(acc_path)
            f_psip = ROOT.TFile.Open(psip_path) if os.path.exists(psip_path) else None
        except Exception as e:
            print(f"Error opening files: {e}")
            sys.exit(1)

        # Retrieve the new 1D Acceptance TDirectory and TH1D
        acc_dir = acc_file.Get("Integrated_pT")
        if not acc_dir:
            print("CRITICAL ERROR: TDirectory 'Integrated_pT' not found in acceptance file.")
            sys.exit(1)
            
        h_acc_1d = acc_dir.Get(f"h_ratio_{target_label}_Integrated_pT")
        if not h_acc_1d and target_label == "LD2":
            h_acc_1d = acc_dir.Get("h_ratio_LH2_Integrated_pT")
            
        if not h_acc_1d:
            print(f"CRITICAL ERROR: Acceptance histogram for {target_label} not found.")
            sys.exit(1)

        dir_xsec = self.get_or_create_dir(self.out_file, f"CrossSections_{target_label}")
        dir_xsec.cd()
        
        h_sub_stat = h_sub_dict["stat"]
        h_sub_sys = h_sub_dict["sys"]
        h_sub_centroid = h_sub_dict["centroid"]

        n_pt_bins = len(config.PT_BINS) - 1
        
        latex_psip_table_content = r"""\begin{longtable}{|c|c|c|c|c|}
\caption{$\psi'$ Contamination Table for %s} \label{tab:psip_contamination_%s} \\
\hline
\textbf{pT bin} & \textbf{Mass bin (GeV)} & \textbf{$\psi'$ contamination} & \textbf{Contribution to $\sigma$ (nb/GeV)} & \textbf{$\delta\sigma_{\psi'}^{\rm syst.}$ (nb/GeV)} \\
\hline
\endfirsthead

\multicolumn{5}{c}%%
{{\bfseries \tablename\ \thetable{} -- continued from previous page}} \\
\hline
\textbf{pT bin} & \textbf{Mass bin (GeV)} & \textbf{$\psi'$ contamination} & \textbf{Contribution to $\sigma$ (nb/GeV)} & \textbf{$\delta\sigma_{\psi'}^{\rm syst.}$ (nb/GeV)} \\
\hline
\endhead

\hline \multicolumn{5}{|r|}{{Continued on next page}} \\ \hline
\endfoot

\hline
\endlastfoot
""" % (target_label, target_label)

        g_xsec = ROOT.TGraphErrors()
        g_xsec.SetName(f"g_xsec_{target_label}")
        g_sys = ROOT.TGraphErrors()
        g_sys.SetName(f"g_sys_{target_label}")

        h1_xsec = ROOT.TH1D(f"h1_xsec_{target_label}", f"Single Differential Cross Section {target_label};p_{{T}} [GeV];d#sigma/dp_{{T}} [nb/GeV/Nucleus]", n_pt_bins, config.PT_BINS)
        h1_sys = ROOT.TH1D(f"h1_sys_{target_label}", f"Systematic Error {target_label};p_{{T}} [GeV];d#sigma/dp_{{T}} [nb/GeV/Nucleus]", n_pt_bins, config.PT_BINS)

        point_idx = 0
        y_min_data = sys.float_info.max
        y_max_data = -sys.float_info.max

        for i_pt in range(n_pt_bins):
            root_bin_x = i_pt + 1
            pt_min, pt_max = config.PT_BINS[i_pt], config.PT_BINS[i_pt+1]
            pt_width = pt_max - pt_min
            pt_center = (pt_min + pt_max) / 2.0
            
            # Extract 1D acceptance
            acceptance = h_acc_1d.GetBinContent(root_bin_x)
            acceptance_err = h_acc_1d.GetBinError(root_bin_x)
            
            if acceptance <= 0: continue
            
            Y_sub = h_sub_stat.GetBinContent(root_bin_x)
            if Y_sub <= 0: continue
            
            Y_sub_stat_err = h_sub_stat.GetBinError(root_bin_x)
            Y_sub_sys_err = h_sub_sys.GetBinError(root_bin_x)
            
            actual_mass_center = h_sub_centroid.GetBinContent(root_bin_x)
            if actual_mass_center <= 0:
                actual_mass_center = (config.MASS_BINS[0] + config.MASS_BINS[-1]) / 2.0

            # Evaluate PsiP Ratio at the exact mass centroid for this pT bin
            h_ratio_psip = f_psip.Get(f"hRatio_PsiP_DY_pT_{i_pt}") if f_psip else None
            psip_ratio = 0.0
            if h_ratio_psip:
                ratio_bin = h_ratio_psip.FindBin(actual_mass_center)
                psip_ratio = h_ratio_psip.GetBinContent(ratio_bin)
                if psip_ratio > 1.0: psip_ratio = 0.0

            numerator = global_constant * Y_sub
            denominator = pt_width * acceptance
            bin_dsigma_dpt = numerator / denominator
            
            sum_xsec = bin_dsigma_dpt
            
            stat_unc = (Y_sub_stat_err / Y_sub) * bin_dsigma_dpt
            sys_tot_yield = (Y_sub_sys_err / Y_sub) * bin_dsigma_dpt
            sys_psip_cont = psip_ratio * bin_dsigma_dpt
            
            if psip_ratio > 0.0:
                s_pt = f"[{pt_min:.2f}, {pt_max:.2f})"
                s_mass = f"Integrated [{config.MASS_BINS[0]:.2f}, {config.MASS_BINS[-1]:.2f})"
                s_ratio = f"{psip_ratio:.4f}"
                sys_unc_bin = math.sqrt(sys_tot_yield**2 + sys_psip_cont**2)
                s_sigma_psip_col = f"{bin_dsigma_dpt:.4f} $\\pm$ {stat_unc:.4f} $\\pm$ {sys_unc_bin:.4f}"
                latex_psip_table_content += f"{s_pt} & {s_mass} & {s_ratio} & {s_sigma_psip_col} & {sys_psip_cont:.4f} \\\\ \n\\hline\n"
            
            if sum_xsec > 0:
                # The 1D Acceptance error is fully correlated within this pT bin
                sys_acc_total = (acceptance_err / acceptance) * sum_xsec
                
                total_stat_err = stat_unc
                total_sys_err = math.sqrt(sys_tot_yield**2 + sys_psip_cont**2 + sys_acc_total**2)
                
                max_err_for_range = max(total_stat_err, total_sys_err)
                
                if (max_err_for_range / sum_xsec) > 0.99:
                    continue
                
                g_sys.SetPoint(point_idx, pt_center, sum_xsec)
                g_sys.SetPointError(point_idx, pt_width/2.0, total_sys_err)
                g_xsec.SetPoint(point_idx, pt_center, sum_xsec)
                g_xsec.SetPointError(point_idx, 0.0, total_stat_err)
                
                h1_xsec.SetBinContent(root_bin_x, sum_xsec)
                h1_xsec.SetBinError(root_bin_x, total_stat_err)
                h1_sys.SetBinContent(root_bin_x, sum_xsec)
                h1_sys.SetBinError(root_bin_x, total_sys_err)
                
                y_high = sum_xsec + max_err_for_range
                y_low  = sum_xsec - max_err_for_range
                
                if y_low <= 0: y_low = sum_xsec * 0.5 
                if y_high > y_max_data: y_max_data = y_high
                if y_low < y_min_data: y_min_data = y_low
                
                point_idx += 1
                
        if g_xsec.GetN() > 0:
            c_xsec = ROOT.TCanvas(f"c_xsec_{target_label}_vs_pT", "", 800, 600)
            c_xsec.SetLeftMargin(0.16)
            c_xsec.SetBottomMargin(0.14)
            c_xsec.SetTickx(1); c_xsec.SetTicky(1)
            
            mg = ROOT.TMultiGraph()
            mg.SetTitle(f";p_{{T}} [GeV];d#sigma / dp_{{T}} [nb / GeV / Nucleus]")
            
            if y_min_data < y_max_data:
                mg.SetMinimum(0.0)
                mg.SetMaximum(y_max_data * 1.4)
            else:
                mg.SetMinimum(0.0)
                mg.SetMaximum(3.0)
            
            leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
            leg.SetBorderSize(0)
            
            g_sys.SetMarkerSize(0); g_sys.SetLineColor(ROOT.kRed); g_sys.SetFillColorAlpha(ROOT.kPink - 9, 0.5); g_sys.SetFillStyle(1001)
            mg.Add(g_sys, "2"); leg.AddEntry(g_sys, "Systematic Unc.", "f")

            g_xsec.SetMarkerStyle(20); g_xsec.SetMarkerColor(ROOT.kRed); g_xsec.SetLineColor(ROOT.kRed)
            mg.Add(g_xsec, "P"); leg.AddEntry(g_xsec, f"Data ({target_label})", "lep")
            
            mg.Draw("A"); mg.GetXaxis().CenterTitle(); mg.GetYaxis().CenterTitle() 
            mg.GetXaxis().SetLimits(0.0, 2.0)
            mg.GetXaxis().SetTitleOffset(1.3)
            mg.GetYaxis().SetTitleOffset(1.7)
            leg.Draw()
            
            target_prefix = "pp" if target_label == "LH2" else "pd" if target_label == "LD2" else target_label
            internal_title = ROOT.TLatex()
            internal_title.SetNDC(True); internal_title.SetTextFont(42); internal_title.SetTextSize(0.04); internal_title.SetTextAlign(13)
            internal_title.DrawLatex(0.19, 0.85, f"Drell-Yan process in {target_prefix}")

            prelim = ROOT.TLatex()
            prelim.SetNDC(True); prelim.SetTextColor(ROOT.kBlue); prelim.SetTextAlign(33); prelim.SetTextSize(0.05)
            prelim.DrawLatex(0.82, 0.6, "Preliminary")

            plot_y_min = 0.0
            plot_y_max = y_max_data * 1.4 if y_min_data < y_max_data else 3.0
            dynamic_y = plot_y_min + 0.05 * (plot_y_max - plot_y_min)

            lumi_note = ROOT.TLatex()
            lumi_note.SetNDC(False); lumi_note.SetTextFont(42); lumi_note.SetTextColor(ROOT.kBlack)
            lumi_note.SetTextAlign(11); lumi_note.SetTextSize(0.025)
            lumi_note.DrawLatex(0.1, dynamic_y, "10% global uncertainty due to the integrated luminosity is not included in the error bands")
            
            c_xsec.SaveAs(f"CrossSection_{target_label}_vs_pT.pdf")
            
            dir_xsec.cd()
            g_xsec.Write()
            g_sys.Write()
            h1_xsec.Write()
            h1_sys.Write()
            c_xsec.Write()
            c_xsec.Close()

        with open(f"Table_PsiP_Contamination_{target_label}.tex", "w") as f:
            f.write(latex_psip_table_content + r"\end{longtable}" + "\n")

        acc_file.Close()
        if f_psip: f_psip.Close()

    def generate_overlay_plot(self):
        """Combines and saves LH2 and LD2 single differential cross-sections onto a single plot canvas."""
        ROOT.gStyle.SetTitleAlign(23); ROOT.gStyle.SetTitleX(0.5); ROOT.gStyle.SetTitleY(0.99)
        ROOT.gStyle.SetTitleH(0.04); ROOT.gStyle.SetTitleBorderSize(0)

        canvas = ROOT.TCanvas("canvas_overlay_Targets_pT", "Cross-Section Comparison", 800, 600)
        canvas.SetLeftMargin(0.16); canvas.SetBottomMargin(0.14)
        canvas.SetTickx(1); canvas.SetTicky(1)

        legend = ROOT.TLegend(0.65, 0.7, 0.88, 0.88)
        legend.SetBorderSize(0); legend.SetFillStyle(0)
        legend.SetTextFont(43); legend.SetTextSize(18)

        h_frame = canvas.DrawFrame(0.0, 0.0, 2.0, 1.0)
        h_frame.SetTitle("DY Absolute Cross-Section Vs p_{T}")
        h_frame.GetXaxis().SetTitle("p_{T} (GeV)"); h_frame.GetXaxis().CenterTitle()
        h_frame.GetXaxis().SetTitleOffset(1.3); h_frame.GetYaxis().SetTitle("d#sigma/dp_{T} (nb/GeV)")
        h_frame.GetYaxis().CenterTitle(); h_frame.GetYaxis().SetTitleOffset(1.7)  

        targets = [("LH2", ROOT.kBlue), ("LD2", ROOT.kRed)]
        
        y_min = 1e9
        y_max = -1e9

        for target, color in targets:
            dir_xsec = self.out_file.Get(f"CrossSections_{target}")
            if not dir_xsec: continue
            
            g_xsec = dir_xsec.Get(f"g_xsec_{target}")
            g_sys = dir_xsec.Get(f"g_sys_{target}")
            if not g_xsec or not g_sys: continue

            g_sys_clone = g_sys.Clone(f"g_sys_clone_{target}")
            g_sys_clone.SetLineColor(color); g_sys_clone.SetFillColorAlpha(color, 0.35)
            g_sys_clone.SetFillStyle(1001); g_sys_clone.SetMarkerSize(0)

            g_xsec_clone = g_xsec.Clone(f"g_xsec_clone_{target}")
            g_xsec_clone.SetLineColor(color); g_xsec_clone.SetMarkerColor(color)
            g_xsec_clone.SetMarkerStyle(ROOT.kFullCircle); g_xsec_clone.SetMarkerSize(1.0)

            g_sys_clone.Draw("2 SAME"); g_xsec_clone.Draw("P SAME")
            legend.AddEntry(g_xsec_clone, f"{target} Data", "pl")
            
            for i in range(g_xsec.GetN()):
                y = g_xsec.GetY()[i]
                if y > 0:
                    ey = g_sys.GetErrorY(i)
                    if (y - ey) < y_min and (y - ey) > 0: y_min = y - ey
                    if (y + ey) > y_max: y_max = y + ey

        if y_min < y_max:
            h_frame.SetMinimum(0.0)
            h_frame.SetMaximum(y_max * 1.4)

        legend.Draw()
        
        prelim = ROOT.TLatex()
        prelim.SetNDC(True); prelim.SetTextFont(43); prelim.SetTextSize(24)
        prelim.SetTextColor(ROOT.kBlue); prelim.SetTextAlign(33)
        prelim.DrawLatex(0.85, 0.6, "Preliminary")

        lumi_note = ROOT.TLatex()
        lumi_note.SetNDC(True); lumi_note.SetTextFont(43); lumi_note.SetTextSize(18)
        lumi_note.SetTextColor(ROOT.kBlack); lumi_note.SetTextAlign(11)
        lumi_note.DrawLatex(0.18, 0.16, "10% global uncertainty due to the integrated luminosity is not included in the error bands")

        canvas.Update()
        out_pdf = "cross_section_overlay_Targets_vs_pT.pdf"
        canvas.SaveAs(out_pdf)
        
        dir_overlay = self.get_or_create_dir(self.out_file, "Overlays")
        dir_overlay.cd()
        canvas.Write("canvas_overlay_Targets_pT")

    def calculate_cross_sections(self):
        """Runs the subtraction and cross-section logic for all targets."""
        if self.hists_lh2 and self.hists_fl:
            self.sub_dict_lh2 = self.generate_subtracted_plot(self.hists_lh2, self.hists_fl, config.FLASK_NORM_LH2, "LH2")
            if self.sub_dict_lh2:
                self.calculate_and_plot_cross_section(self.sub_dict_lh2, "LH2", config.GLOBAL_CONSTANT_LH2)

        if self.hists_ld2 and self.hists_lh2 and self.hists_fl:
            self.sub_dict_pd = self.generate_pd_subtracted_plot(self.hists_ld2, self.hists_lh2, self.hists_fl)
            if self.sub_dict_pd:
                self.calculate_and_plot_cross_section(self.sub_dict_pd, "LD2", config.GLOBAL_CONSTANT_LD2)

        self.generate_overlay_plot()

    def generate_latex_appendix(self):
        """Generates LaTeX source code for the Mass Centroid derivations."""
        latex_filename = "Appendix_MassCentroids.tex"
        
        with open(latex_filename, "w") as tex_file:
            intro_text = r"""\section{Appendix: Determination of Mass Bin Centroids (pT Binned)}...""" 
            tex_file.write(intro_text)

    def finalize(self):
        """Writes remaining buffers and closes the TFile safely."""
        self.out_file.Write()
        self.out_file.Close()