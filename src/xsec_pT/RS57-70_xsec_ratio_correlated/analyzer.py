"""
analyzer.py
Core Object-Oriented Analysis Module for Drell-Yan Cross-Sections.
Includes dynamic CSV export for error propagation presentation.
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

        # Load dynamic bin-by-bin roadset systematics AND weighted means
        lh2_sys_csv = "/root/github/e906-development/src/xsec_pT/RS57-70_weighted_average/Roadset_Sys_StdDev_LH2_geom.csv"
        ld2_sys_csv = "/root/github/e906-development/src/xsec_pT/RS57-70_weighted_average/Roadset_Sys_StdDev_LD2_geom.csv"
        
        self.roadset_data_lh2 = self._load_systematic_csv(lh2_sys_csv)
        self.roadset_data_ld2 = self._load_systematic_csv(ld2_sys_csv)

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
        ROOT.gStyle.SetEndErrorSize(5) 
        ROOT.gErrorIgnoreLevel = ROOT.kFatal

    def _load_systematic_csv(self, filepath):
        """Reads the Weighted Mean Cross Section and Standard Deviation from the CSVs."""
        sys_dict = {}
        if not os.path.exists(filepath):
            print(f"Warning: Roadset systematic file not found at {filepath}.")
            return sys_dict
        
        with open(filepath, "r") as f:
            reader = csv.reader(f)
            next(reader, None) # Skip metadata header
            next(reader, None) # Skip column headers
            for row in reader:
                if len(row) < 4: continue
                try:
                    bin_idx = int(row[0])
                    mean_val = float(row[1])
                    std_dev = float(row[2])
                    rel_err = float(row[3].replace('%', '').strip()) / 100.0
                    
                    sys_dict[bin_idx] = {
                        "mean_xsec": mean_val,
                        "sys_abs": std_dev,
                        "sys_rel": rel_err
                    }
                except ValueError:
                    continue
        return sys_dict

    @staticmethod
    def get_or_create_dir(base_dir, name):
        d = base_dir.GetDirectory(name)
        if not d:
            d = base_dir.mkdir(name)
        return d

    @staticmethod
    def apply_cuts(tree, cut=4.2):
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

        filtered_events = {}
        for key, val in events.items():
            filtered_events[key] = val[total_cut_mask]
            
        filtered_events["pT"] = np.sqrt(filtered_events["dpx"]**2 + filtered_events["dpy"]**2)
            
        return filtered_events

    def get_concatenated_events(self, file_paths, tree_name):
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
        hists["Pt_Centroid"] = make_th1("Pt_Centroid", "Data-Driven pT Centroid", "p_{T} [GeV]")

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
        if N == 0: return 0.0, 0.0
        mean_eff = np.mean(eff_arr)
        
        diff_matrix = np.abs(loc_arr[:, None] - loc_arr[None, :])
        correl_matrix = np.zeros((N, N))
        correl_matrix[diff_matrix == 0] = 1.0  
        correl_matrix[diff_matrix == 1] = 1.0  
        
        covar_matrix = correl_matrix * np.outer(err_arr, err_arr)
        err_on_mean = np.sqrt(np.sum(covar_matrix)) / N
        return mean_eff, err_on_mean

    def extract_target_stats(self, data_tot, data_mix, loc_tot, loc_mix, m_low, m_high, pt_low, pt_high):
        mask_tot = (data_tot["mass"] >= m_low) & (data_tot["mass"] < m_high) & \
                   (data_tot["pT"] >= pt_low) & (data_tot["pT"] < pt_high)
        
        mask_mix = (data_mix["mass"] >= m_low) & (data_mix["mass"] < m_high) & \
                   (data_mix["pT"] >= pt_low) & (data_mix["pT"] < pt_high)

        N_tot = np.sum(mask_tot)
        N_mix = np.sum(mask_mix)
        
        sum_mass_tot = np.sum(data_tot["mass"][mask_tot]) if N_tot > 0 else 0.0
        sum_mass_mix = np.sum(data_mix["mass"][mask_mix]) if N_mix > 0 else 0.0
        
        sum_pt_tot = np.sum(data_tot["pT"][mask_tot]) if N_tot > 0 else 0.0
        sum_pt_mix = np.sum(data_mix["pT"][mask_mix]) if N_mix > 0 else 0.0

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

        return N_tot, N_mix, sum_mass_tot, sum_mass_mix, sum_pt_tot, sum_pt_mix, val_sig_eff, err_sig_eff, diff_yield, means

    def fill_histograms(self, hists, root_x, pt_center, N_tot, N_mix, sig_eff, err_sig_eff, diff_yield, means, final_centroid, final_pt_centroid):
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
        hists["Pt_Centroid"].SetBinContent(root_x, final_pt_centroid)

    def save_1d_pdfs(self, hists_dict, label):
        """Saves all generated TH1D histograms as separate PDF files."""
        for name, hist in hists_dict.items():
            if "stat" in name or "sys" in name or "Centroid" in name: continue 
            c = ROOT.TCanvas(f"c_{name}_{label}", "", 1200, 900)
            c.SetRightMargin(0.05); c.SetLeftMargin(0.16); c.SetBottomMargin(0.14)
            c.SetTickx(1); c.SetTicky(1)
            
            max_y, min_y = -1e9, 1e9
            for i in range(1, hist.GetNbinsX() + 1):
                val = hist.GetBinContent(i)
                err = hist.GetBinError(i)
                if val == 0 and err == 0: continue
                if val + err > max_y: max_y = val + err
                if val - err < min_y: min_y = val - err
                
            if max_y > -1e8:
                hist.SetMaximum(max_y * 1.5) 
                hist.SetMinimum(min_y * 1.2 if min_y < 0 else 0.0)
            else:
                hist.SetMinimum(0.0); hist.SetMaximum(1.0)
            
            hist.SetLineColor(ROOT.kBlue); hist.SetLineWidth(2)
            hist.Draw("HIST E1")
            
            latex_draws = []
            offset = (max_y if max_y > -1e8 else 1.0) * 0.08
            for i in range(1, hist.GetNbinsX() + 1):
                val = hist.GetBinContent(i)
                err = hist.GetBinError(i)
                if val == 0 and err == 0: continue
                pt_center = hist.GetBinCenter(i)
                y_pos = val + err + offset
                txt = ROOT.TLatex(pt_center, y_pos, f"#splitline{{{val:.4f}}}{{#pm {err:.4f}}}")
                txt.SetTextSize(0.022); txt.SetTextAlign(21); txt.SetTextColor(ROOT.kBlack)
                txt.Draw()
                latex_draws.append(txt)

            c.SaveAs(f"{name}_{label}.pdf")
            c.Close()

    def generate_efficiency_csv(self, data_tot, loc_tot, m_low, m_high, pt_low, pt_high, pt_idx, target_label):
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
        def get_loc_data(data_dict, is_mix=False):
            if is_mix and 'ptrk_D1' in data_dict and 'ntrk_D1' in data_dict:
                d1_vals = 0.5 * (data_dict['ptrk_D1'] + data_dict['ntrk_D1'])
            else:
                d1_vals = data_dict['D1']
            return np.digitize(d1_vals, self.x_curve) - 1
        
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
        dir_lh2.cd(); self.hists_lh2 = self.create_histograms("LH2")

        dir_ld2 = self.get_or_create_dir(dir_kin, "LD2")
        dir_ld2.cd(); self.hists_ld2 = self.create_histograms("LD2")

        dir_fl = self.get_or_create_dir(dir_kin, "Flask")
        dir_fl.cd(); self.hists_fl = self.create_histograms("Flask")

        dir_csv = self.get_or_create_dir(self.out_file, "CSV_Tables")
        csv_filename_lh2 = "Table_Kinematics_LH2_vs_pT.csv"
        csv_filename_ld2 = "Table_Kinematics_LD2_vs_pT.csv"
        
        csv_rows_lh2, csv_rows_ld2 = [], []

        m_low, m_high = config.MASS_BINS[0], config.MASS_BINS[-1]
        m_center_overall = (m_low + m_high) / 2.0

        for i_pt in range(len(config.PT_BINS) - 1):
            pt_low, pt_high = config.PT_BINS[i_pt], config.PT_BINS[i_pt+1]
            pt_center = (pt_low + pt_high) / 2.0
            root_x = i_pt + 1
            
            self.generate_efficiency_csv(data_lh2_tot, loc_lh2_tot, m_low, m_high, pt_low, pt_high, i_pt, "LH2")
            self.generate_efficiency_csv(data_ld2_tot, loc_ld2_tot, m_low, m_high, pt_low, pt_high, i_pt, "LD2")
            self.generate_efficiency_csv(data_fl_tot, loc_fl_tot, m_low, m_high, pt_low, pt_high, i_pt, "Flask")

            N_l2_t, N_l2_m, M_l2_t, M_l2_m, P_l2_t, P_l2_m, eps_ld2, err_eps_ld2, diff_l2, mean_l2 = self.extract_target_stats(data_ld2_tot, data_ld2_mix, loc_ld2_tot, loc_ld2_mix, m_low, m_high, pt_low, pt_high)
            N_lh_t, N_lh_m, M_lh_t, M_lh_m, P_lh_t, P_lh_m, eps_lh2, err_eps_lh2, diff_lh, mean_lh = self.extract_target_stats(data_lh2_tot, data_lh2_mix, loc_lh2_tot, loc_lh2_mix, m_low, m_high, pt_low, pt_high)
            N_fl_t, N_fl_m, M_fl_t, M_fl_m, P_fl_t, P_fl_m, eps_fl, err_eps_fl, diff_fl, mean_fl = self.extract_target_stats(data_fl_tot, data_fl_mix, loc_fl_tot, loc_fl_mix, m_low, m_high, pt_low, pt_high)

            if eps_lh2 > 0:
                cY_lh_t, cM_lh_t, cP_lh_t = N_lh_t / eps_lh2, M_lh_t / eps_lh2, P_lh_t / eps_lh2
                cY_lh_m, cM_lh_m, cP_lh_m = N_lh_m / eps_lh2, M_lh_m / eps_lh2, P_lh_m / eps_lh2
            else:
                cY_lh_t, cM_lh_t, cP_lh_t, cY_lh_m, cM_lh_m, cP_lh_m = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

            if eps_ld2 > 0:
                cY_l2_t, cM_l2_t, cP_l2_t = N_l2_t / eps_ld2, M_l2_t / eps_ld2, P_l2_t / eps_ld2
                cY_l2_m, cM_l2_m, cP_l2_m = N_l2_m / eps_ld2, M_l2_m / eps_ld2, P_l2_m / eps_ld2
            else:
                cY_l2_t, cM_l2_t, cP_l2_t, cY_l2_m, cM_l2_m, cP_l2_m = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

            if eps_fl > 0:
                cY_fl_t_lh, cM_fl_t_lh, cP_fl_t_lh = (config.FLASK_NORM_LH2 * N_fl_t) / eps_fl, (config.FLASK_NORM_LH2 * M_fl_t) / eps_fl, (config.FLASK_NORM_LH2 * P_fl_t) / eps_fl
                cY_fl_m_lh, cM_fl_m_lh, cP_fl_m_lh = (config.FLASK_NORM_LH2 * N_fl_m) / eps_fl, (config.FLASK_NORM_LH2 * M_fl_m) / eps_fl, (config.FLASK_NORM_LH2 * P_fl_m) / eps_fl
                cY_fl_t_l2, cM_fl_t_l2, cP_fl_t_l2 = (config.FLASK_NORM_LD2 * N_fl_t) / eps_fl, (config.FLASK_NORM_LD2 * M_fl_t) / eps_fl, (config.FLASK_NORM_LD2 * P_fl_t) / eps_fl
                cY_fl_m_l2, cM_fl_m_l2, cP_fl_m_l2 = (config.FLASK_NORM_LD2 * N_fl_m) / eps_fl, (config.FLASK_NORM_LD2 * M_fl_m) / eps_fl, (config.FLASK_NORM_LD2 * P_fl_m) / eps_fl
            else:
                cY_fl_t_lh, cM_fl_t_lh, cP_fl_t_lh, cY_fl_m_lh, cM_fl_m_lh, cP_fl_m_lh = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                cY_fl_t_l2, cM_fl_t_l2, cP_fl_t_l2, cY_fl_m_l2, cM_fl_m_l2, cP_fl_m_l2 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

            # LH2 Centroids
            num_LH2 = (cM_lh_t - cM_lh_m) - (cM_fl_t_lh - cM_fl_m_lh)
            num_pt_LH2 = (cP_lh_t - cP_lh_m) - (cP_fl_t_lh - cP_fl_m_lh)
            den_LH2 = (cY_lh_t - cY_lh_m) - (cY_fl_t_lh - cY_fl_m_lh)
            cent_LH2 = num_LH2 / den_LH2 if den_LH2 != 0 else m_center_overall
            cent_pt_LH2 = num_pt_LH2 / den_LH2 if den_LH2 != 0 else pt_center

            # LD2 Centroids
            num_LD2 = (cM_l2_t - cM_l2_m) - (cM_fl_t_l2 - cM_fl_m_l2) - (config.LH2_TO_LD2_NORM * num_LH2)
            num_pt_LD2 = (cP_l2_t - cP_l2_m) - (cP_fl_t_l2 - cP_fl_m_l2) - (config.LH2_TO_LD2_NORM * num_pt_LH2)
            den_LD2 = (cY_l2_t - cY_l2_m) - (cY_fl_t_l2 - cY_fl_m_l2) - (config.LH2_TO_LD2_NORM * den_LH2)
            cent_LD2 = num_LD2 / den_LD2 if den_LD2 != 0 else m_center_overall
            cent_pt_LD2 = num_pt_LD2 / den_LD2 if den_LD2 != 0 else pt_center

            self.fill_histograms(self.hists_ld2, root_x, pt_center, N_l2_t, N_l2_m, eps_ld2, err_eps_ld2, diff_l2, mean_l2, cent_LD2, cent_pt_LD2)
            self.fill_histograms(self.hists_lh2, root_x, pt_center, N_lh_t, N_lh_m, eps_lh2, err_eps_lh2, diff_lh, mean_lh, cent_LH2, cent_pt_LH2)
            
            # Flask Centroids
            num_fl = (M_fl_t - M_fl_m) / (eps_fl if eps_fl > 0 else 1e-9)
            num_pt_fl = (P_fl_t - P_fl_m) / (eps_fl if eps_fl > 0 else 1e-9)
            den_fl = (N_fl_t - N_fl_m) / (eps_fl if eps_fl > 0 else 1e-9)
            cent_fl = num_fl / den_fl if den_fl != 0 else m_center_overall
            cent_pt_fl = num_pt_fl / den_fl if den_fl != 0 else pt_center
            
            self.fill_histograms(self.hists_fl, root_x, pt_center, N_fl_t, N_fl_m, eps_fl, err_eps_fl, diff_fl, mean_fl, cent_fl, cent_pt_fl)

            csv_rows_lh2.append({
                "pT Bin": f"[{pt_low:.2f}, {pt_high:.2f})", "pT Center": pt_center,
                "LH2 Mass Centroid": cent_LH2, "N_LH2_total": N_lh_t, "N_LH2_mixed": N_lh_m,
                "N_flask_total": N_fl_t, "N_flask_mixed": N_fl_m,
                "eps_LH2": eps_lh2, "eps_LD2": eps_ld2, "eps_flask": eps_fl,
                "Corrected Total Mass LH2 (num)": num_LH2, "Corrected Yield LH2 (denom)": den_LH2
            })

            csv_rows_ld2.append({
                "pT Bin": f"[{pt_low:.2f}, {pt_high:.2f})", "pT Center": pt_center,
                "LD2 Mass Centroid": cent_LD2, "N_LD2_total": N_l2_t, "N_LD2_mixed": N_l2_m,
                "N_LH2_total": N_lh_t, "N_LH2_mixed": N_lh_m,
                "N_flask_total": N_fl_t, "N_flask_mixed": N_fl_m,
                "eps_LH2": eps_lh2, "eps_LD2": eps_ld2, "eps_flask": eps_fl,
                "Corrected Total Mass LD2 (num)": num_LD2, "Corrected Yield LD2 (denom)": den_LD2
            })

        with open(csv_filename_lh2, "w", newline='') as f:
            writer = csv.DictWriter(f, fieldnames=csv_rows_lh2[0].keys()); writer.writeheader(); writer.writerows(csv_rows_lh2)
            
        with open(csv_filename_ld2, "w", newline='') as f:
            writer = csv.DictWriter(f, fieldnames=csv_rows_ld2[0].keys()); writer.writeheader(); writer.writerows(csv_rows_ld2)
            
        dir_csv.cd()
        macro_lh2 = ROOT.TMacro(csv_filename_lh2); macro_lh2.Write(csv_filename_lh2)
        macro_ld2 = ROOT.TMacro(csv_filename_ld2); macro_ld2.Write(csv_filename_ld2)

        self.save_1d_pdfs(self.hists_lh2, "LH2")
        self.save_1d_pdfs(self.hists_ld2, "LD2")
        self.save_1d_pdfs(self.hists_fl, "Flask")

    def generate_subtracted_plot(self, hists_target, hists_flask, flask_norm, target_label):
        if "Y_corrected_stat" not in hists_target or "Y_corrected_stat" not in hists_flask: return None

        dir_sub = self.get_or_create_dir(self.out_file, f"Subtracted_Plots_{target_label}")
        dir_sub.cd()
        
        h_target_stat = hists_target["Y_corrected_stat"]
        h_target_sys = hists_target["Y_corrected_sys"]
        h_flask_stat = hists_flask["Y_corrected_stat"]
        h_flask_sys = hists_flask["Y_corrected_sys"]
        
        name = f"Y_corrected_Subtracted_{target_label}"
        title = f"Corrected Yield ({target_label} - Flask)"
        h_sub = ROOT.TH1D(name, f"{title};p_{{T}} [GeV];Yield", len(config.PT_BINS)-1, config.PT_BINS)
        h_sub.Sumw2(); h_sub.SetStats(0)
        h_sub.GetXaxis().CenterTitle(); h_sub.GetYaxis().CenterTitle()
        h_sub.GetXaxis().SetTitleOffset(1.3); h_sub.GetYaxis().SetTitleOffset(1.7)

        h_sub_stat = h_sub.Clone(f"{name}_stat")
        h_sub_sys = h_sub.Clone(f"{name}_sys")
        h_sub_centroid = hists_target["Mass_Centroid"].Clone(f"{name}_Mass_Centroid")
        h_sub_pt_centroid = hists_target["Pt_Centroid"].Clone(f"{name}_Pt_Centroid")

        for i_pt in range(len(config.PT_BINS) - 1):
            bin_x = i_pt + 1
            
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
            
            h_sub.SetBinContent(bin_x, val_sub); h_sub.SetBinError(bin_x, err_sub_total)
            h_sub_stat.SetBinContent(bin_x, val_sub); h_sub_stat.SetBinError(bin_x, err_sub_stat)
            h_sub_sys.SetBinContent(bin_x, val_sub); h_sub_sys.SetBinError(bin_x, err_sub_sys)

        c = ROOT.TCanvas(f"c_{name}", f"c_{name}", 1200, 900)
        c.SetRightMargin(0.05); c.SetLeftMargin(0.16); c.SetBottomMargin(0.14)
        c.SetTickx(1); c.SetTicky(1)
        
        max_y, min_y = -1e9, 1e9
        for ix in range(1, h_sub.GetNbinsX() + 1):
            val = h_sub.GetBinContent(ix); err = h_sub.GetBinError(ix)
            if val == 0 and err == 0: continue
            if val + err > max_y: max_y = val + err
            if val - err < min_y: min_y = val - err
        
        if max_y > -1e8:
            h_sub.SetMaximum(max_y * 1.5)
            h_sub.SetMinimum(0.0 if min_y >= 0 else min_y * 1.2)
        else:
            h_sub.SetMinimum(0.0); h_sub.SetMaximum(1.0)
            
        h_sub.SetLineColor(ROOT.kBlack); h_sub.SetMarkerStyle(20); h_sub.Draw("P E1")
        
        latex_draws = []
        offset = (max_y if max_y > -1e8 else 1.0) * 0.08
        for ix in range(1, h_sub.GetNbinsX() + 1):
            val = h_sub.GetBinContent(ix); err = h_sub.GetBinError(ix)
            if val == 0 and err == 0: continue
            pt_center = h_sub.GetBinCenter(ix)
            y_pos = val + err + offset
            txt = ROOT.TLatex(pt_center, y_pos, f"#splitline{{{val:.4f}}}{{#pm {err:.4f}}}")
            txt.SetTextSize(0.022); txt.SetTextAlign(21); txt.SetTextColor(ROOT.kBlack)
            txt.Draw(); latex_draws.append(txt)

        c.SaveAs(f"{name}.pdf"); c.Close()
        return {"stat": h_sub_stat, "sys": h_sub_sys, "centroid": h_sub_centroid, "pt_centroid": h_sub_pt_centroid}

    def generate_pd_subtracted_plot(self, hists_ld2, hists_lh2, hists_flask):
        dir_sub = self.get_or_create_dir(self.out_file, "Subtracted_Plots_LD2_pd")
        dir_sub.cd()

        h_ld2_stat, h_ld2_sys = hists_ld2["Y_corrected_stat"], hists_ld2["Y_corrected_sys"]
        h_lh2_stat, h_lh2_sys = hists_lh2["Y_corrected_stat"], hists_lh2["Y_corrected_sys"]
        h_flask_stat, h_flask_sys = hists_flask["Y_corrected_stat"], hists_flask["Y_corrected_sys"]

        c_lh2 = config.THD_THH_RATIO * (config.PROTONS_ON_TARGET_LD2 / config.PROTONS_ON_TARGET_LH2)
        c_flask_sub = config.FLASK_NORM_LD2 - (c_lh2 * config.FLASK_NORM_LH2)

        name = "Y_corrected_Subtracted_LD2"
        title = "Corrected Yield (LD2 - LH2 - Flask)"
        h_pd = ROOT.TH1D(name, f"{title};p_{{T}} [GeV];Yield", len(config.PT_BINS)-1, config.PT_BINS)
        h_pd.Sumw2(); h_pd.SetStats(0)
        h_pd.GetXaxis().CenterTitle(); h_pd.GetYaxis().CenterTitle()
        h_pd.GetXaxis().SetTitleOffset(1.3); h_pd.GetYaxis().SetTitleOffset(1.7)

        h_pd_stat = h_pd.Clone(f"{name}_stat")
        h_pd_sys = h_pd.Clone(f"{name}_sys")
        h_pd_centroid = hists_ld2["Mass_Centroid"].Clone(f"{name}_Mass_Centroid")
        h_pd_pt_centroid = hists_ld2["Pt_Centroid"].Clone(f"{name}_Pt_Centroid")

        for i_pt in range(len(config.PT_BINS) - 1):
            bin_x = i_pt + 1
            y_ld2, y_lh2, y_flask = h_ld2_stat.GetBinContent(bin_x), h_lh2_stat.GetBinContent(bin_x), h_flask_stat.GetBinContent(bin_x)
            e_ld2_stat, e_lh2_stat, e_flask_stat = h_ld2_stat.GetBinError(bin_x), h_lh2_stat.GetBinError(bin_x), h_flask_stat.GetBinError(bin_x)
            e_ld2_sys, e_lh2_sys, e_flask_sys = h_ld2_sys.GetBinError(bin_x), h_lh2_sys.GetBinError(bin_x), h_flask_sys.GetBinError(bin_x)
            
            val_pd = y_ld2 - (c_lh2 * y_lh2) - (c_flask_sub * y_flask)
            
            err_pd_stat = np.sqrt(e_ld2_stat**2 + (c_lh2 * e_lh2_stat)**2 + (c_flask_sub * e_flask_stat)**2)
            err_pd_sys = np.sqrt(e_ld2_sys**2 + (c_lh2 * e_lh2_sys)**2 + (c_flask_sub * e_flask_sys)**2)
            err_pd_total = np.sqrt(err_pd_stat**2 + err_pd_sys**2)
            
            h_pd.SetBinContent(bin_x, val_pd); h_pd.SetBinError(bin_x, err_pd_total)
            h_pd_stat.SetBinContent(bin_x, val_pd); h_pd_stat.SetBinError(bin_x, err_pd_stat)
            h_pd_sys.SetBinContent(bin_x, val_pd); h_pd_sys.SetBinError(bin_x, err_pd_sys)

        c = ROOT.TCanvas(f"c_{name}", f"c_{name}", 1200, 900)
        c.SetRightMargin(0.05); c.SetLeftMargin(0.16); c.SetBottomMargin(0.14)
        c.SetTickx(1); c.SetTicky(1)
        
        max_y, min_y = -1e9, 1e9
        for ix in range(1, h_pd.GetNbinsX() + 1):
            val = h_pd.GetBinContent(ix); err = h_pd.GetBinError(ix)
            if val == 0 and err == 0: continue
            if val + err > max_y: max_y = val + err
            if val - err < min_y: min_y = val - err
        
        if max_y > -1e8:
            h_pd.SetMaximum(max_y * 1.5)
            h_pd.SetMinimum(0.0 if min_y >= 0 else min_y * 1.2)
        else:
            h_pd.SetMinimum(0.0); h_pd.SetMaximum(1.0)
            
        h_pd.SetLineColor(ROOT.kBlack); h_pd.SetMarkerStyle(20); h_pd.Draw("P E1")
        
        latex_draws = []
        offset = (max_y if max_y > -1e8 else 1.0) * 0.08
        for ix in range(1, h_pd.GetNbinsX() + 1):
            val = h_pd.GetBinContent(ix); err = h_pd.GetBinError(ix)
            if val == 0 and err == 0: continue
            pt_center = h_pd.GetBinCenter(ix)
            y_pos = val + err + offset
            txt = ROOT.TLatex(pt_center, y_pos, f"#splitline{{{val:.4f}}}{{#pm {err:.4f}}}")
            txt.SetTextSize(0.022); txt.SetTextAlign(21); txt.SetTextColor(ROOT.kBlack)
            txt.Draw(); latex_draws.append(txt)

        c.SaveAs(f"{name}.pdf"); c.Close()
        return {"stat": h_pd_stat, "sys": h_pd_sys, "centroid": h_pd_centroid, "pt_centroid": h_pd_pt_centroid}

    def calculate_and_plot_cross_section(self, h_sub_dict, target_label, global_constant, use_true_pt=False):
        """Builds single differential cross-sections vs pT, incorporating weighted CSV data."""
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
        
        h_sub_stat, h_sub_sys, h_sub_centroid, h_sub_pt_centroid = h_sub_dict["stat"], h_sub_dict["sys"], h_sub_dict["centroid"], h_sub_dict["pt_centroid"]

        n_pt_bins = len(config.PT_BINS) - 1
        suffix = "_true_pt" if use_true_pt else "_geom"
        
        # --- Dynamically pull the correct Roadset Systematic Dictionary for the target ---
        sys_dict = self.roadset_data_lh2 if target_label == "LH2" else (self.roadset_data_ld2 if target_label == "LD2" else {})
        
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

        g_xsec = ROOT.TGraphErrors(); g_xsec.SetName(f"g_xsec_{target_label}{suffix}")
        g_sys = ROOT.TGraphErrors(); g_sys.SetName(f"g_sys_{target_label}{suffix}")

        h1_xsec = ROOT.TH1D(f"h1_xsec_{target_label}{suffix}", f"Single Differential Cross Section {target_label};p_{{T}} [GeV];d#sigma/dp_{{T}} [nb/GeV/Nucleus]", n_pt_bins, config.PT_BINS)
        h1_sys = ROOT.TH1D(f"h1_sys_{target_label}{suffix}", f"Systematic Error {target_label};p_{{T}} [GeV];d#sigma/dp_{{T}} [nb/GeV/Nucleus]", n_pt_bins, config.PT_BINS)

        point_idx = 0
        y_min_data, y_max_data = sys.float_info.max, -sys.float_info.max
        max_sys_err = 0.0

        for i_pt in range(n_pt_bins):
            root_bin_x = i_pt + 1
            pt_min, pt_max = config.PT_BINS[i_pt], config.PT_BINS[i_pt+1]
            pt_width = pt_max - pt_min
            pt_center = (pt_min + pt_max) / 2.0
            
            acceptance = h_acc_1d.GetBinContent(root_bin_x)
            acceptance_err = h_acc_1d.GetBinError(root_bin_x)
            
            if acceptance <= 0: continue
            Y_sub = h_sub_stat.GetBinContent(root_bin_x)
            if Y_sub <= 0: continue
            
            Y_sub_stat_err = h_sub_stat.GetBinError(root_bin_x)
            Y_sub_sys_err = h_sub_sys.GetBinError(root_bin_x)
            
            actual_mass_center = h_sub_centroid.GetBinContent(root_bin_x)
            if actual_mass_center <= 0: actual_mass_center = (config.MASS_BINS[0] + config.MASS_BINS[-1]) / 2.0
                
            actual_pt_center = h_sub_pt_centroid.GetBinContent(root_bin_x)
            if actual_pt_center <= 0: actual_pt_center = pt_center

            h_ratio_psip = f_psip.Get(f"hRatio_PsiP_DY_pT_{i_pt}") if f_psip else None
            psip_ratio = 0.0
            if h_ratio_psip:
                ratio_bin = h_ratio_psip.FindBin(actual_mass_center)
                psip_ratio = h_ratio_psip.GetBinContent(ratio_bin)
                if psip_ratio > 1.0: psip_ratio = 0.0

            # Raw unweighted calculation from current ROOT files
            numerator = global_constant * Y_sub
            denominator = pt_width * acceptance
            raw_dsigma_dpt = numerator / denominator
            
            stat_unc_raw = (Y_sub_stat_err / Y_sub) * raw_dsigma_dpt
            sys_tot_yield_raw = (Y_sub_sys_err / Y_sub) * raw_dsigma_dpt
            sys_psip_cont_raw = psip_ratio * raw_dsigma_dpt
            
            # --- OVERRIDE WITH DYNAMIC CSV DATA ---
            csv_data = sys_dict.get(i_pt, None)
            
            if csv_data is not None and csv_data["mean_xsec"] > 0:
                # Use the inverse-variance weighted mean from CSV
                sum_xsec = csv_data["mean_xsec"]
                # The roadset systematic is the absolute standard deviation from the CSV
                sys_roadset = csv_data["sys_abs"]
                
                # Scale the raw uncertainties to match the new central value
                scale = sum_xsec / raw_dsigma_dpt if raw_dsigma_dpt > 0 else 1.0
                stat_unc = stat_unc_raw * scale
                sys_tot_yield = sys_tot_yield_raw * scale
                sys_psip_cont = sys_psip_cont_raw * scale
            else:
                # Fallback to the raw calculated data
                sum_xsec = raw_dsigma_dpt
                sys_roadset = 0.0
                stat_unc = stat_unc_raw
                sys_tot_yield = sys_tot_yield_raw
                sys_psip_cont = sys_psip_cont_raw
            
            if psip_ratio > 0.0 and not use_true_pt:
                s_pt = f"[{pt_min:.2f}, {pt_max:.2f})"
                s_mass = f"Integrated [{config.MASS_BINS[0]:.2f}, {config.MASS_BINS[-1]:.2f})"
                s_ratio = f"{psip_ratio:.4f}"
                sys_unc_bin = math.sqrt(sys_tot_yield**2 + sys_psip_cont**2 + sys_roadset**2)
                s_sigma_psip_col = f"{sum_xsec:.4f} $\\pm$ {stat_unc:.4f} $\\pm$ {sys_unc_bin:.4f}"
                latex_psip_table_content += f"{s_pt} & {s_mass} & {s_ratio} & {s_sigma_psip_col} & {sys_psip_cont:.4f} \\\\ \n\\hline\n"
            
            if sum_xsec > 0:
                sys_acc_total = (acceptance_err / acceptance) * sum_xsec
                total_stat_err = stat_unc
                sys_lumi = 0.10 * sum_xsec
                
                # Add roadset systematic and lumi in quadrature
                total_sys_err = math.sqrt(sys_tot_yield**2 + sys_psip_cont**2 + sys_acc_total**2 + sys_roadset**2 + sys_lumi**2)
                
                # Keep tracking max sys for baseline visualization
                if total_sys_err > max_sys_err: max_sys_err = total_sys_err
                
                # Used for generic filtering only
                if (max(total_stat_err, total_sys_err) / sum_xsec) > 0.99: continue
                
                plot_x = actual_pt_center if use_true_pt else pt_center
                
                g_sys.SetPoint(point_idx, pt_center, sum_xsec)
                g_sys.SetPointError(point_idx, pt_width/2.0, total_sys_err)
                
                g_xsec.SetPoint(point_idx, plot_x, sum_xsec)
                
                # Adding the horizontal error bar mapping to the TGraphErrors marker
                g_xsec.SetPointError(point_idx, pt_width/2.0, total_stat_err) 
                
                h1_xsec.SetBinContent(root_bin_x, sum_xsec); h1_xsec.SetBinError(root_bin_x, total_stat_err)
                h1_sys.SetBinContent(root_bin_x, sum_xsec); h1_sys.SetBinError(root_bin_x, total_sys_err)
                
                # Dynamic range scaling based exclusively on data points + stat error 
                y_high = sum_xsec + total_stat_err
                y_low  = sum_xsec - total_stat_err
                
                if y_low <= 0: y_low = sum_xsec * 0.5 
                if y_high > y_max_data: y_max_data = y_high
                if y_low < y_min_data: y_min_data = y_low
                
                point_idx += 1
                
        if g_xsec.GetN() > 0:
            c_xsec = ROOT.TCanvas(f"c_xsec_{target_label}{suffix}_vs_pT", "", 800, 600)
            c_xsec.SetLeftMargin(0.16); c_xsec.SetBottomMargin(0.14)
            c_xsec.SetTickx(1); c_xsec.SetTicky(1)
            
            # --- DYNAMIC Y-AXIS & BASELINE CALCULATION ---
            delta_y = y_max_data - y_min_data if y_max_data > y_min_data else (y_max_data * 0.5 if y_max_data > 0 else 1.0)
            gap = delta_y * 0.15  # 15% visual gap between data and sys band
            
            # Ensure the top of the sys band (baseline + max_sys_err) is 'gap' below the lowest data point
            baseline = y_min_data - gap - max_sys_err
            
            # Set the frame minimum slightly below the bottom of the sys band (baseline - max_sys_err)
            plot_y_min = baseline - max_sys_err - (delta_y * 0.05)
            
            # Dynamically set the frame maximum to create empty space at the top.
            # Scaling the data span by 1.6 places the highest point at ~62% of the axis,
            # ensuring it comfortably clears the TLatex text located at NDC y=0.75 to 0.85.
            plot_y_max = plot_y_min + (y_max_data - plot_y_min) * 1.6

            mg = ROOT.TMultiGraph()
            mg.SetTitle(f";p_{{T}} [GeV];d#sigma / dp_{{T}} [nb / GeV / Nucleus]")
            mg.SetMinimum(plot_y_min); mg.SetMaximum(plot_y_max)
            
            leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
            leg.SetBorderSize(0)
            
            # Repurpose cloned g_sys to bottom band level
            g_sys_bottom = g_sys.Clone(f"g_sys_bottom_{target_label}{suffix}")
            for i in range(g_sys_bottom.GetN()):
                g_sys_bottom.SetPoint(i, g_sys_bottom.GetX()[i], baseline)

            g_sys_bottom.SetMarkerSize(0); g_sys_bottom.SetLineColor(ROOT.kRed)
            g_sys_bottom.SetFillColorAlpha(ROOT.kPink - 9, 0.5); g_sys_bottom.SetFillStyle(1001)
            
            mg.Add(g_sys_bottom, "2"); leg.AddEntry(g_sys_bottom, "Syst. Unc. Band", "f")

            g_xsec.SetMarkerStyle(20); g_xsec.SetMarkerColor(ROOT.kRed); g_xsec.SetLineColor(ROOT.kRed)
            mg.Add(g_xsec, "P"); leg.AddEntry(g_xsec, f"Data ({target_label})", "lep")
            
            mg.Draw("A"); mg.GetXaxis().CenterTitle(); mg.GetYaxis().CenterTitle() 
            mg.GetXaxis().SetLimits(0.0, 2.0); mg.GetXaxis().SetTitleOffset(1.3); mg.GetYaxis().SetTitleOffset(1.7)
            
            leg.Draw()
            
            target_prefix = "pp" if target_label == "LH2" else "pd" if target_label == "LD2" else target_label
            internal_title = ROOT.TLatex()
            internal_title.SetNDC(True)
            internal_title.SetTextFont(42)
            internal_title.SetTextSize(0.04)
            internal_title.SetTextAlign(13)

            internal_title.DrawLatex(0.19, 0.85, f"Drell-Yan process in {target_prefix}")
            internal_title.DrawLatex(0.19, 0.80, "0.0 < x_{F} < 0.8")
            internal_title.DrawLatex(0.19, 0.75, "4.2 GeV < M < 8.8 GeV")
            prelim = ROOT.TLatex()
            prelim.SetNDC(True)
            prelim.SetTextColor(ROOT.kBlue)
            prelim.SetTextAlign(33) 

            prelim.SetTextSize(0.05)
            prelim.DrawLatex(0.82, 0.60, "Preliminary")

            prelim.SetTextSize(0.0252) 
            prelim.DrawLatex(0.82, 0.54, "Run Period 2014-2015")

            # -------------------------------------------------------------
            # Updated Luminosity Note
            # lumi_note = ROOT.TLatex()
            # lumi_note.SetNDC(True)
            # lumi_note.SetTextFont(42)
            # lumi_note.SetTextColor(ROOT.kBlack)
            # lumi_note.SetTextAlign(11)
            # lumi_note.SetTextSize(0.025)

            # note_text = "10% global uncertainty due to integrated luminosity is included in the error bands."
            # lumi_note.DrawLatex(0.19, 0.70, note_text)
            # -------------------------------------------------------------
            
            c_xsec.SaveAs(f"CrossSection_{target_label}{suffix}_vs_pT.pdf")
            
            dir_xsec.cd()
            g_xsec.Write(); g_sys.Write(); h1_xsec.Write(); h1_sys.Write(); c_xsec.Write(); c_xsec.Close()

        if not use_true_pt:
            with open(f"Table_PsiP_Contamination_{target_label}.tex", "w") as f:
                f.write(latex_psip_table_content + r"\end{longtable}" + "\n")

        acc_file.Close()
        if f_psip: f_psip.Close()

    def generate_overlay_plot(self, use_true_pt=False):
        suffix = "_true_pt" if use_true_pt else "_geom"
        ROOT.gStyle.SetTitleAlign(23); ROOT.gStyle.SetTitleX(0.5); ROOT.gStyle.SetTitleY(0.99)
        ROOT.gStyle.SetTitleH(0.04); ROOT.gStyle.SetTitleBorderSize(0)

        canvas = ROOT.TCanvas(f"canvas_overlay_Targets_pT{suffix}", "Cross-Section Comparison", 800, 600)
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
        y_min, y_max = 1e9, -1e9
        max_sys_all = 0.0

        for target, color in targets:
            dir_xsec = self.out_file.Get(f"CrossSections_{target}")
            if not dir_xsec: continue
            
            g_xsec = dir_xsec.Get(f"g_xsec_{target}{suffix}")
            g_sys = dir_xsec.Get(f"g_sys_{target}{suffix}")
            if not g_xsec or not g_sys: continue
            
            for i in range(g_xsec.GetN()):
                y = g_xsec.GetY()[i]
                if y > 0:
                    ey = g_xsec.GetErrorY(i)
                    if (y - ey) < y_min and (y - ey) > 0: y_min = y - ey
                    if (y + ey) > y_max: y_max = y + ey
                    
                    sys_ey = g_sys.GetErrorY(i)
                    if sys_ey > max_sys_all: max_sys_all = sys_ey

        plot_y_max = y_max * 1.6 if y_min < y_max else 3.0
        if y_min < y_max: h_frame.SetMinimum(0.0); h_frame.SetMaximum(plot_y_max)

        # Baseline offset for systematic bands
        baseline = plot_y_max * 0.12
        if baseline - max_sys_all < plot_y_max * 0.02:
            baseline = max_sys_all + plot_y_max * 0.02
            
        line_base = ROOT.TLine(0.0, baseline, 2.0, baseline)
        line_base.SetLineStyle(2); line_base.SetLineColor(ROOT.kGray+2)
        line_base.Draw("SAME")

        keepalive = []
        for target, color in targets:
            dir_xsec = self.out_file.Get(f"CrossSections_{target}")
            if not dir_xsec: continue
            g_xsec = dir_xsec.Get(f"g_xsec_{target}{suffix}")
            g_sys = dir_xsec.Get(f"g_sys_{target}{suffix}")
            if not g_xsec or not g_sys: continue

            g_sys_clone = g_sys.Clone(f"g_sys_clone_{target}{suffix}")
            for i in range(g_sys_clone.GetN()):
                g_sys_clone.SetPoint(i, g_sys_clone.GetX()[i], baseline)

            g_sys_clone.SetLineColor(color); g_sys_clone.SetFillColorAlpha(color, 0.35)
            g_sys_clone.SetFillStyle(1001); g_sys_clone.SetMarkerSize(0)

            g_xsec_clone = g_xsec.Clone(f"g_xsec_clone_{target}{suffix}")
            g_xsec_clone.SetLineColor(color); g_xsec_clone.SetMarkerColor(color)
            g_xsec_clone.SetMarkerStyle(ROOT.kFullCircle); g_xsec_clone.SetMarkerSize(1.0)

            g_sys_clone.Draw("2 SAME"); g_xsec_clone.Draw("P SAME")
            legend.AddEntry(g_xsec_clone, f"{target} Data", "pl")
            
            keepalive.extend([g_sys_clone, g_xsec_clone])

        legend.Draw()
        
        prelim = ROOT.TLatex()
        prelim.SetNDC(True)
        prelim.SetTextColor(ROOT.kBlue)
        prelim.SetTextAlign(33) 

        prelim.SetTextSize(0.05)
        prelim.DrawLatex(0.82, 0.60, "Preliminary")

        prelim.SetTextSize(0.0272) 
        prelim.DrawLatex(0.82, 0.54, "Run Period 2014-2015")
        
        # -------------------------------------------------------------
        # Updated Luminosity Note
        # lumi_note = ROOT.TLatex()
        # lumi_note.SetNDC(True); lumi_note.SetTextFont(43); lumi_note.SetTextSize(18)
        # lumi_note.SetTextColor(ROOT.kBlack); lumi_note.SetTextAlign(11)
        # lumi_note.DrawLatex(0.18, 0.16, "10% global uncertainty due to the integrated luminosity is included in the error bands.")
        # -------------------------------------------------------------
        
        canvas.Update(); canvas.SaveAs(f"cross_section_overlay_Targets_vs_pT{suffix}.pdf")
        
        dir_overlay = self.get_or_create_dir(self.out_file, "Overlays")
        dir_overlay.cd(); canvas.Write(f"canvas_overlay_Targets_pT{suffix}")

    def generate_ratio_plot(self, use_true_pt=False):
        suffix = "_true_pt" if use_true_pt else "_geom"
        dir_lh2 = self.out_file.Get("CrossSections_LH2")
        dir_ld2 = self.out_file.Get("CrossSections_LD2")

        if not dir_lh2 or not dir_ld2:
            print("Warning: Cross-section directories not found. Cannot calculate ratio.")
            return

        h1_xsec_lh2 = dir_lh2.Get(f"h1_xsec_LH2{suffix}")
        h1_sys_lh2 = dir_lh2.Get(f"h1_sys_LH2{suffix}")
        h1_xsec_ld2 = dir_ld2.Get(f"h1_xsec_LD2{suffix}")
        h1_sys_ld2 = dir_ld2.Get(f"h1_sys_LD2{suffix}")

        if not h1_xsec_lh2 or not h1_xsec_ld2 or not h1_sys_lh2 or not h1_sys_ld2:
            print(f"Warning: Cross-section histograms not found for suffix '{suffix}'. Cannot calculate ratio.")
            return

        h_pt_lh2 = self.sub_dict_lh2["pt_centroid"] if self.sub_dict_lh2 else None
        h_pt_ld2 = self.sub_dict_pd["pt_centroid"] if self.sub_dict_pd else None

        dir_ratio = self.get_or_create_dir(self.out_file, "Ratio_pd_2pp")
        dir_ratio.cd()

        g_ratio_stat = ROOT.TGraphErrors()
        g_ratio_stat.SetName(f"g_ratio_stat{suffix}")
        
        g_ratio_sys = ROOT.TGraphErrors()
        g_ratio_sys.SetName(f"g_ratio_sys{suffix}")
        
        g_ratio_total = ROOT.TGraphErrors()
        g_ratio_total.SetName(f"g_ratio_total{suffix}")

        pt_idx = 0
        y_max_ratio, y_min_ratio = -1e9, 1e9

        print(f"\n--- Cross-Section Ratio pd/2pp ({'True pT' if use_true_pt else 'Geom pT'}) ---")
        print(f"{'pT Bin [GeV]':<15} | {'Ratio':<10} | {'Stat Err':<10} | {'Sys Err':<10} | {'Total Err':<10}")
        print("-" * 65)

        csv_rows = []

        for i_pt in range(len(config.PT_BINS) - 1):
            bin_x = i_pt + 1
            y_lh2 = h1_xsec_lh2.GetBinContent(bin_x)
            y_ld2 = h1_xsec_ld2.GetBinContent(bin_x)
            
            if y_lh2 <= 0 or y_ld2 <= 0: continue

            err_stat_lh2 = h1_xsec_lh2.GetBinError(bin_x)
            err_sys_lh2 = h1_sys_lh2.GetBinError(bin_x)
            err_stat_ld2 = h1_xsec_ld2.GetBinError(bin_x)
            err_sys_ld2 = h1_sys_ld2.GetBinError(bin_x)

            ratio = y_ld2 / (2.0 * y_lh2)
            
            rel_stat_ld2 = err_stat_ld2 / y_ld2
            rel_stat_lh2 = err_stat_lh2 / y_lh2
            err_ratio_stat = ratio * math.sqrt(rel_stat_ld2**2 + rel_stat_lh2**2)

            rel_sys_ld2 = err_sys_ld2 / y_ld2
            rel_sys_lh2 = err_sys_lh2 / y_lh2
            err_ratio_sys = ratio * abs(rel_sys_ld2 - rel_sys_lh2)

            err_ratio_total = math.sqrt(err_ratio_stat**2 + err_ratio_sys**2)

            pt_min = config.PT_BINS[i_pt]
            pt_max = config.PT_BINS[i_pt+1]
            print(f"[{pt_min:.2f}, {pt_max:.2f})   | {ratio:.4f}     | {err_ratio_stat:.4f}     | {err_ratio_sys:.4f}     | {err_ratio_total:.4f}")

            pt_str = f"[{pt_min:.2f}, {pt_max:.2f})"
            csv_rows.append([pt_str, f"{ratio:.4f}", f"{err_ratio_stat:.4f}", f"{err_ratio_sys:.4f}", f"{err_ratio_total:.4f}"])

            bin_center = (config.PT_BINS[i_pt] + config.PT_BINS[i_pt+1]) / 2.0

            if use_true_pt and h_pt_lh2 and h_pt_ld2:
                x_lh2 = h_pt_lh2.GetBinContent(bin_x)
                x_ld2 = h_pt_ld2.GetBinContent(bin_x)
                plot_x = (x_lh2 + x_ld2) / 2.0
            else:
                plot_x = bin_center

            pt_width = config.PT_BINS[i_pt+1] - config.PT_BINS[i_pt]

            g_ratio_stat.SetPoint(pt_idx, plot_x, ratio)
            g_ratio_stat.SetPointError(pt_idx, 0.0, err_ratio_stat)

            g_ratio_sys.SetPoint(pt_idx, bin_center, ratio)
            g_ratio_sys.SetPointError(pt_idx, pt_width/2.0, err_ratio_sys)
            
            g_ratio_total.SetPoint(pt_idx, plot_x, ratio)
            g_ratio_total.SetPointError(pt_idx, 0.0, err_ratio_total)

            pt_idx += 1

        print("\n")

        if use_true_pt:
            csv_filename = "pd_2pp_errors.csv"
            with open(csv_filename, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["pT Bin [GeV]", "Ratio", "Stat Err", "Sys Err", "Total Err"])
                writer.writerows(csv_rows)
            print(f"[*] Dynamically saved error table to {csv_filename}")

        if g_ratio_total.GetN() > 0:
            c_ratio = ROOT.TCanvas(f"c_ratio_pd_2pp{suffix}", "Ratio pd / 2pp", 800, 600)
            c_ratio.SetLeftMargin(0.16); c_ratio.SetBottomMargin(0.14)
            c_ratio.SetTickx(1); c_ratio.SetTicky(1)

            mg = ROOT.TMultiGraph()
            mg.SetTitle(";p_{T} [GeV];#sigma_{pd} / 2#sigma_{pp}")

            fit_func = ROOT.TF1("fit_ratio", "pol0", 0.0, 2.0)
            g_ratio_total.Fit(fit_func, "Q0")
            fit_val = fit_func.GetParameter(0)
            fit_err = fit_func.GetParError(0)

            g_fit_band = ROOT.TGraphErrors()
            g_fit_band.SetPoint(0, 0.0, fit_val)
            g_fit_band.SetPointError(0, 0.0, fit_err)
            g_fit_band.SetPoint(1, 2.0, fit_val)
            g_fit_band.SetPointError(1, 0.0, fit_err)
            g_fit_band.SetFillColorAlpha(ROOT.kPink, 0.4)
            g_fit_band.SetFillStyle(1001)

            y_high_fit = fit_val + fit_err
            y_low_fit = fit_val - fit_err
            if y_high_fit > y_max_ratio: y_max_ratio = y_high_fit
            if y_low_fit < y_min_ratio: y_min_ratio = y_low_fit
            
            for i in range(g_ratio_total.GetN()):
                y = g_ratio_total.GetY()[i]
                err = g_ratio_total.GetErrorY(i)
                if y + err > y_max_ratio: y_max_ratio = y + err
                if y - err < y_min_ratio: y_min_ratio = y - err

            if y_min_ratio < y_max_ratio:
                y_range = y_max_ratio - y_min_ratio
                mg.SetMinimum(max(0.0, y_min_ratio - 0.5 * y_range))
                mg.SetMaximum(y_max_ratio + 0.6 * y_range)
            else:
                mg.SetMinimum(0.0); mg.SetMaximum(2.0)

            leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
            leg.SetBorderSize(0)

            mg.Add(g_fit_band, "3")

            g_ratio_total.SetMarkerStyle(20)
            g_ratio_total.SetMarkerColor(ROOT.kBlack)
            g_ratio_total.SetLineColor(ROOT.kBlack)
            mg.Add(g_ratio_total, "P")
            leg.AddEntry(g_ratio_total, "Data Ratio (Total Error)", "lep")

            mg.Draw("A"); mg.GetXaxis().CenterTitle(); mg.GetYaxis().CenterTitle()
            mg.GetXaxis().SetLimits(0.0, 2.0)
            mg.GetXaxis().SetTitleOffset(1.3); mg.GetYaxis().SetTitleOffset(1.3)

            line = ROOT.TLine(0.0, 1.0, 2.0, 1.0)
            line.SetLineStyle(2)
            line.SetLineColor(ROOT.kGray+2)
            line.SetLineWidth(2)
            line.Draw("SAME")

            fit_func.SetLineColor(ROOT.kRed)
            fit_func.SetLineWidth(2)
            fit_func.Draw("SAME")

            leg.Draw()

            latex_fit = ROOT.TLatex()
            latex_fit.SetTextFont(42)
            latex_fit.SetTextSize(0.04)
            latex_fit.SetTextColor(ROOT.kRed)
            
            y_text = 0.12 + fit_val + fit_err + ((mg.GetYaxis().GetXmax() - mg.GetYaxis().GetXmin()) * 0.02)
            latex_fit.DrawLatex(0.2, y_text, f"Best Fit: {fit_val:.4f} #pm {fit_err:.4f}")

            prelim = ROOT.TLatex()
            prelim.SetNDC(True)
            prelim.SetTextColor(ROOT.kBlue)
            prelim.SetTextAlign(33)

            prelim.SetTextSize(0.05)
            prelim.DrawLatex(0.82, 0.60, "Preliminary")

            prelim.SetTextSize(0.0252) 
            prelim.DrawLatex(0.82, 0.54, "Run Period 2014-2015")

            c_ratio.SaveAs(f"CrossSection_Ratio_pd_2pp_vs_pT{suffix}.pdf")

            dir_ratio.cd()
            g_ratio_stat.Write()
            g_ratio_sys.Write()
            g_ratio_total.Write()
            c_ratio.Write()
            c_ratio.Close()

    def generate_combined_ratio_overlay_plot(self, use_true_pt=False):
        suffix = "_true_pt" if use_true_pt else "_geom"

        dir_lh2 = self.out_file.Get("CrossSections_LH2")
        dir_ld2 = self.out_file.Get("CrossSections_LD2")
        dir_ratio = self.out_file.Get("Ratio_pd_2pp")

        if not dir_lh2 or not dir_ld2 or not dir_ratio:
            print(f"Warning: Missing directories for combined plot for suffix {suffix}.")
            return

        g_xsec_lh2 = dir_lh2.Get(f"g_xsec_LH2{suffix}")
        g_sys_lh2  = dir_lh2.Get(f"g_sys_LH2{suffix}")
        g_xsec_ld2 = dir_ld2.Get(f"g_xsec_LD2{suffix}")
        g_sys_ld2  = dir_ld2.Get(f"g_sys_LD2{suffix}")

        g_ratio_total = dir_ratio.Get(f"g_ratio_total{suffix}")

        if not all([g_xsec_lh2, g_sys_lh2, g_xsec_ld2, g_sys_ld2, g_ratio_total]):
            print(f"Warning: Missing graphs for combined plot for suffix {suffix}.")
            return

        canvas = ROOT.TCanvas(f"c_combined_{suffix}", "Combined Overlay and Ratio", 800, 800)

        # --- Pad 1: Cross Sections Overlay (Top 65%) ---
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
        pad1.SetBottomMargin(0.15) 
        pad1.SetLeftMargin(0.16)
        pad1.SetRightMargin(0.05)
        pad1.SetTopMargin(0.08)
        pad1.SetTickx(1); pad1.SetTicky(1)
        pad1.Draw()
        pad1.cd()

        mg_top = ROOT.TMultiGraph()
        
        y_max = 0
        max_sys_both = 0
        for i in range(g_xsec_lh2.GetN()):
            if g_xsec_lh2.GetY()[i] + g_xsec_lh2.GetErrorY(i) > y_max: y_max = g_xsec_lh2.GetY()[i] + g_xsec_lh2.GetErrorY(i)
            if g_xsec_ld2.GetY()[i] + g_xsec_ld2.GetErrorY(i) > y_max: y_max = g_xsec_ld2.GetY()[i] + g_xsec_ld2.GetErrorY(i)
            if g_sys_lh2.GetErrorY(i) > max_sys_both: max_sys_both = g_sys_lh2.GetErrorY(i)
            if g_sys_ld2.GetErrorY(i) > max_sys_both: max_sys_both = g_sys_ld2.GetErrorY(i)

        # Calculate the dynamic range and baseline
        delta_y = y_max if y_max > 0 else 1.0
        gap = delta_y * 0.15  # 15% visual gap between data and sys band
        
        # Set the baseline so the top of the sys band sits below the data gap
        # Add a tiny 2% offset to ensure it NEVER touches the true 0.0 line
        baseline = -gap - max_sys_both + (delta_y * 0.02)
        
        # Set the frame minimum slightly below the bottom of the sys band
        plot_y_min = baseline - max_sys_both - (delta_y * 0.05)
        plot_y_max = y_max * 1.6

        mg_top.SetMinimum(plot_y_min)
        mg_top.SetMaximum(plot_y_max)

        g_sys_lh2_clone = g_sys_lh2.Clone()
        for i in range(g_sys_lh2_clone.GetN()): g_sys_lh2_clone.SetPoint(i, g_sys_lh2_clone.GetX()[i], baseline)
        g_sys_lh2_clone.SetLineColor(ROOT.kRed)
        g_sys_lh2_clone.SetFillColorAlpha(ROOT.kPink - 9, 0.5)
        
        g_xsec_lh2_clone = g_xsec_lh2.Clone()
        g_xsec_lh2_clone.SetLineColor(ROOT.kRed)
        g_xsec_lh2_clone.SetMarkerColor(ROOT.kRed)

        g_sys_ld2_clone = g_sys_ld2.Clone()
        for i in range(g_sys_ld2_clone.GetN()): g_sys_ld2_clone.SetPoint(i, g_sys_ld2_clone.GetX()[i], baseline)
        g_sys_ld2_clone.SetLineColor(ROOT.kBlue)
        g_sys_ld2_clone.SetFillColorAlpha(ROOT.kAzure + 1, 0.5) 
        
        g_xsec_ld2_clone = g_xsec_ld2.Clone()
        g_xsec_ld2_clone.SetLineColor(ROOT.kBlue)
        g_xsec_ld2_clone.SetMarkerColor(ROOT.kBlue)

        mg_top.Add(g_sys_lh2_clone, "2")
        mg_top.Add(g_sys_ld2_clone, "2")
        mg_top.Add(g_xsec_lh2_clone, "P")
        mg_top.Add(g_xsec_ld2_clone, "P")

        mg_top.Draw("A")
        mg_top.SetTitle("")
        mg_top.GetXaxis().SetLimits(0.0, 2.0)

        mg_top.GetYaxis().SetTitle("d#sigma/dp_{T} [nb/GeV]")
        mg_top.GetYaxis().CenterTitle()
        mg_top.GetYaxis().SetTitleFont(43); mg_top.GetYaxis().SetTitleSize(22)
        mg_top.GetYaxis().SetTitleOffset(2.0)
        mg_top.GetYaxis().SetLabelFont(43); mg_top.GetYaxis().SetLabelSize(20)

        mg_top.GetXaxis().SetTitle("p_{T} [GeV]")
        mg_top.GetXaxis().CenterTitle()
        mg_top.GetXaxis().SetTitleFont(43); mg_top.GetXaxis().SetTitleSize(22)
        mg_top.GetXaxis().SetTitleOffset(1.2) 
        mg_top.GetXaxis().SetLabelFont(43); mg_top.GetXaxis().SetLabelSize(20)

        line_base = ROOT.TLine(0.0, baseline, 2.0, baseline)
        line_base.SetLineStyle(2); line_base.SetLineColor(ROOT.kGray+2)
        line_base.Draw("SAME")
        
        # Add a zero line to separate the domains cleanly
        line_zero = ROOT.TLine(0.0, 0.0, 2.0, 0.0)
        line_zero.SetLineStyle(1); line_zero.SetLineColor(ROOT.kBlack)
        line_zero.Draw("SAME")

        leg_top = ROOT.TLegend(0.65, 0.65, 0.88, 0.88)
        leg_top.SetBorderSize(0); leg_top.SetFillStyle(0); leg_top.SetTextFont(43); leg_top.SetTextSize(18)
        leg_top.AddEntry(g_xsec_lh2_clone, "LH2 Data", "pl")
        leg_top.AddEntry(g_sys_lh2_clone, "LH2 Sys. Unc.", "f")
        leg_top.AddEntry(g_xsec_ld2_clone, "LD2 Data", "pl")
        leg_top.AddEntry(g_sys_ld2_clone, "LD2 Sys. Unc.", "f")
        leg_top.Draw()

        prelim = ROOT.TLatex()
        prelim.SetNDC(True)
        prelim.SetTextColor(ROOT.kBlue)
        prelim.SetTextAlign(33) 

        prelim.SetTextSize(0.05)
        prelim.DrawLatex(0.33, 0.77, "Preliminary")
        
        prelim.SetTextSize(0.0252) 
        prelim.DrawLatex(0.33, 0.71, "Run Period 2014-2015")

        # --- Pad 2: Ratio (Bottom 35%) ---
        canvas.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.35)
        pad2.SetTopMargin(0.05)
        pad2.SetBottomMargin(0.35) 
        pad2.SetLeftMargin(0.16)
        pad2.SetRightMargin(0.05)
        pad2.SetTickx(1); pad2.SetTicky(1)
        pad2.Draw()
        pad2.cd()

        mg_bottom = ROOT.TMultiGraph()

        fit_func = ROOT.TF1("fit_ratio_comb", "pol0", 0.0, 2.0)
        g_ratio_total.Fit(fit_func, "Q0")
        fit_val = fit_func.GetParameter(0)
        fit_err = fit_func.GetParError(0)

        g_fit_band = ROOT.TGraphErrors()
        g_fit_band.SetPoint(0, 0.0, fit_val); g_fit_band.SetPointError(0, 0.0, fit_err)
        g_fit_band.SetPoint(1, 2.0, fit_val); g_fit_band.SetPointError(1, 0.0, fit_err)
        g_fit_band.SetFillColorAlpha(ROOT.kPink, 0.4); g_fit_band.SetFillStyle(1001)

        g_ratio_total_clone = g_ratio_total.Clone()
        g_ratio_total_clone.SetMarkerStyle(20)
        g_ratio_total_clone.SetMarkerColor(ROOT.kBlack)
        g_ratio_total_clone.SetLineColor(ROOT.kBlack)

        mg_bottom.Add(g_fit_band, "3")
        mg_bottom.Add(g_ratio_total_clone, "P")

        mg_bottom.Draw("A")
        mg_bottom.SetTitle("")
        mg_bottom.GetXaxis().SetLimits(0.0, 2.0)

        mg_bottom.GetYaxis().SetTitle("#sigma_{pd}/2#sigma_{pp}")
        mg_bottom.GetYaxis().CenterTitle()
        mg_bottom.GetYaxis().SetTitleFont(43); mg_bottom.GetYaxis().SetTitleSize(22)
        mg_bottom.GetYaxis().SetTitleOffset(2.0)
        mg_bottom.GetYaxis().SetLabelFont(43); mg_bottom.GetYaxis().SetLabelSize(20)
        mg_bottom.GetYaxis().SetNdivisions(505) 

        mg_bottom.GetXaxis().SetTitle("p_{T} [GeV]")
        mg_bottom.GetXaxis().CenterTitle()
        mg_bottom.GetXaxis().SetTitleFont(43); mg_bottom.GetXaxis().SetTitleSize(22)
        mg_bottom.GetXaxis().SetTitleOffset(1.2) 
        mg_bottom.GetXaxis().SetLabelFont(43); mg_bottom.GetXaxis().SetLabelSize(20)

        y_max_r, y_min_r = -1e9, 1e9
        y_high_fit = fit_val + fit_err; y_low_fit = fit_val - fit_err
        if y_high_fit > y_max_r: y_max_r = y_high_fit
        if y_low_fit < y_min_r: y_min_r = y_low_fit
        
        for i in range(g_ratio_total.GetN()):
            y = g_ratio_total.GetY()[i]
            err = g_ratio_total.GetErrorY(i)
            
            if y + err > y_max_r: y_max_r = y + err
            if y - err < y_min_r: y_min_r = y - err

        if y_min_r < y_max_r:
            y_range = y_max_r - y_min_r
            mg_bottom.SetMinimum(max(0.0, y_min_r - 0.5 * y_range))
            mg_bottom.SetMaximum(y_max_r + 0.6 * y_range)
        else:
            mg_bottom.SetMinimum(0.0); mg_bottom.SetMaximum(2.0)

        line = ROOT.TLine(0.0, 1.0, 2.0, 1.0)
        line.SetLineStyle(2); line.SetLineColor(ROOT.kGray+2); line.SetLineWidth(2); line.Draw("SAME")
        fit_func.SetLineColor(ROOT.kRed); fit_func.SetLineWidth(2); fit_func.Draw("SAME")

        latex_fit = ROOT.TLatex()
        latex_fit.SetTextFont(43); latex_fit.SetTextSize(20); latex_fit.SetTextColor(ROOT.kRed)
        y_text = 0.12 + fit_val + fit_err + ((mg_bottom.GetYaxis().GetXmax() - mg_bottom.GetYaxis().GetXmin()) * 0.05)
        latex_fit.DrawLatex(0.2, y_text, f"Best Fit: {fit_val:.4f} #pm {fit_err:.4f}")

        canvas.SaveAs(f"Combined_XSec_Ratio_vs_pT{suffix}.pdf")

        dir_comb = self.get_or_create_dir(self.out_file, "Combined_Plots")
        dir_comb.cd()
        canvas.Write(f"c_combined_ratio_overlay{suffix}")
        canvas.Close()

    def calculate_cross_sections(self):
        if self.hists_lh2 and self.hists_fl:
            self.sub_dict_lh2 = self.generate_subtracted_plot(self.hists_lh2, self.hists_fl, config.FLASK_NORM_LH2, "LH2")
            
        if self.hists_ld2 and self.hists_lh2 and self.hists_fl:
            self.sub_dict_pd = self.generate_pd_subtracted_plot(self.hists_ld2, self.hists_lh2, self.hists_fl)
            
        for use_true_pt in [False, True]:
            if self.sub_dict_lh2:
                self.calculate_and_plot_cross_section(self.sub_dict_lh2, "LH2", config.GLOBAL_CONSTANT_LH2, use_true_pt)
            if self.sub_dict_pd:
                self.calculate_and_plot_cross_section(self.sub_dict_pd, "LD2", config.GLOBAL_CONSTANT_LD2, use_true_pt)
                
            self.generate_overlay_plot(use_true_pt)
            self.generate_ratio_plot(use_true_pt)
            self.generate_combined_ratio_overlay_plot(use_true_pt) 

    def generate_latex_appendix(self):
        latex_filename = "Appendix_MassCentroids.tex"
        with open(latex_filename, "w") as tex_file:
            intro_text = r"""\section{Appendix: Determination of Mass Bin Centroids (pT Binned)}...""" 
            tex_file.write(intro_text)

    def finalize(self):
        self.out_file.Write()
        self.out_file.Close()