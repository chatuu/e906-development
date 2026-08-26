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
from rich.console import Console

class DYCrossSectionAnalyzer:
    def __init__(self, lh2_files, ld2_files, flask_files, mc_messy_files=None, mc_clean_files=None, out_filename="All_XSec_Objects.root"):
        self.console = Console()
        self._setup_root()
        
        self.lh2_paths = lh2_files if isinstance(lh2_files, list) else [lh2_files]
        self.ld2_paths = ld2_files if isinstance(ld2_files, list) else [ld2_files]
        self.flask_paths = flask_files if isinstance(flask_files, list) else [flask_files]
        
        self.mc_messy_dict = mc_messy_files if mc_messy_files else {}
        self.mc_clean_dict = mc_clean_files if mc_clean_files else {}
        
        self.response_matrices_messy = {}
        self.response_matrices_clean = {}
        self.rm_file = None
        
        self.out_filename = out_filename
        self.out_file = ROOT.TFile(out_filename, "RECREATE")
        
        try:
            npz_data = np.load(config.INPUT_NPZ_FILE)
            self.x_curve = npz_data['x']
        except Exception as e:
            self.console.print(f"[bold red]Error loading NPZ file for Covariance at '{config.INPUT_NPZ_FILE}': {e}[/bold red]")
            sys.exit(1)

        self.hists_lh2 = None
        self.hists_ld2 = None
        self.hists_fl = None
        self.sub_dict_lh2 = None
        self.sub_dict_ld2 = None

    def _setup_root(self):
        ROOT.gROOT.SetBatch(True)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kBird)
        ROOT.gErrorIgnoreLevel = ROOT.kFatal
        
        roounfold_lib = os.path.join(config.ROOUNFOLD_PATH, "build", "libRooUnfold.so")
        if not os.path.exists(roounfold_lib):
            roounfold_lib = os.path.join(config.ROOUNFOLD_PATH, "libRooUnfold.so")
            
        if os.path.exists(roounfold_lib):
            ROOT.gSystem.Load(roounfold_lib)
        else:
            self.console.print(f"[bold red]WARNING: RooUnfold library not found at {roounfold_lib}. Unfolding will fail![/bold red]")

    @staticmethod
    def get_or_create_dir(base_dir, name):
        d = base_dir.GetDirectory(name)
        if not d:
            d = base_dir.mkdir(name)
        return d

    def get_cut_mask(self, events, is_mc=False):
        class EventNamespace:
            def __init__(self, data):
                self.__dict__.update(data)
        e = EventNamespace(events)

        try:
            n_events = len(getattr(e, 'mass', getattr(e, 'mMass', [])))

            if is_mc:
                bo = np.full(n_events, 1.6)
            else:
                if hasattr(e, 'runID'):
                    bo = np.where(e.runID >= 11000, 1.6, 0.4)
                else:
                    bo = np.full(n_events, 0.4)

            dimuon_cut = (
                (np.abs(e.dx) < 0.25) & (np.abs(e.dy - bo) < 0.22) &
                (e.dz < -5.) & (e.dz > -280.) & (np.abs(e.dpx) < 1.8) & (np.abs(e.dpy) < 2.0) &
                (e.dpx * e.dpx + e.dpy * e.dpy < 5.) & (e.dpz < 116.) & (e.dpz > 38.) &
                (e.mass >= 3.0) & (e.mass <= 12.0) & 
                (e.dx * e.dx + (e.dy - bo) * (e.dy - bo) < 0.06) &
                (e.xF >= -0.20) & (e.xF <= 1.0) & 
                (e.xT > 0.05) & (e.xT <= 0.58) &
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

            D1_occ_cut = ((e.D1 > 20) & (e.D1 < 385))

            return (track1_cut & track2_cut & tracks_cut & dimuon_cut & occ_cut & D1_occ_cut)
        except Exception as err:
            self.console.print(f"[yellow]Warning: Falling back to unconstrained skim cuts (Exception: {err})[/yellow]")
            return np.ones(len(getattr(e, 'mass', getattr(e, 'mMass', []))), dtype=bool)

    def apply_cuts(self, tree, is_mc=False):
        events = tree.arrays(library="np")
        total_cut_mask = self.get_cut_mask(events, is_mc=is_mc)
        
        filtered_events = {}
        for key, val in events.items():
            filtered_events[key] = val[total_cut_mask]
            
        return filtered_events

    def get_concatenated_events(self, file_paths, tree_name):
        all_events = None
        for fp in file_paths:
            if not os.path.exists(fp): continue
            try:
                with uproot.open(fp) as f:
                    if tree_name not in f: continue
                    filtered = self.apply_cuts(f[tree_name], is_mc=False)
                    if all_events is None:
                        all_events = {k: [v] for k, v in filtered.items()}
                    else:
                        for k, v in filtered.items():
                            if k in all_events: all_events[k].append(v)
            except Exception as e:
                self.console.print(f"[bold red]Error reading {fp}: {e}[/bold red]")
                
        if all_events is None:
            raise RuntimeError(f"No valid data found for tree '{tree_name}' in provided files.")
        
        return {k: np.concatenate(v) for k, v in all_events.items()}

    @staticmethod
    def create_histograms(file_label):
        hists = {}
        def make_th2(name, title):
            h = ROOT.TH2D(f"{name}_{file_label}", f"{title} ({file_label});Mass [GeV];x_{{F}}", 
                          len(config.MASS_BINS)-1, config.MASS_BINS, 
                          len(config.XF_BINS)-1, config.XF_BINS)
            h.Sumw2(); h.SetStats(0); h.GetXaxis().CenterTitle(); h.GetYaxis().CenterTitle()
            h.GetXaxis().SetTitleOffset(1.2); h.GetYaxis().SetTitleOffset(1.2)
            return h

        hists["Y_total"] = make_th2("Y_total", "Total Yield (result)")
        hists["Y_mix"] = make_th2("Y_mix", "Mix Yield (result_mix)")
        hists["E_total_reco"] = make_th2("E_total_reco", "Avg Reco Eff (Total)")
        hists["E_mix_reco"] = make_th2("E_mix_reco", "Avg Reco Eff (Mix)")
        hists["E_total_hodo"] = make_th2("E_total_hodo", "Avg Hodo Eff (Total)")
        hists["E_mix_hodo"] = make_th2("E_mix_hodo", "Avg Hodo Eff (Mix)")
        hists["E_total_final"] = make_th2("E_total_final", "Avg Final Eff (Total)")
        hists["E_mix_final"] = make_th2("E_mix_final", "Avg Final Eff (Mix)")
        hists["E_final_signal"] = make_th2("E_final_signal", "Avg Signal Efficiency")
        hists["Y_corrected"] = make_th2("Y_corrected", "Corrected Yield (Total Error)")
        hists["Y_corrected_stat"] = make_th2("Y_corrected_stat", "Corrected Yield (Stat Error)")
        hists["Y_corrected_sys"] = make_th2("Y_corrected_sys", "Corrected Yield (Sys Error)")
        hists["Mass_Centroid"] = make_th2("Mass_Centroid", "Data-Driven Mass Centroid")
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
        correl_matrix[diff_matrix == 0] = 1.0; correl_matrix[diff_matrix == 1] = 1.0  
        covar_matrix = correl_matrix * np.outer(err_arr, err_arr)
        err_on_mean = np.sqrt(np.sum(covar_matrix)) / N
        return mean_eff, err_on_mean

    @staticmethod
    def add_latex_to_bin(hist, x_center, y_center, value, error):
        val_f, err_f = float(value), float(error)
        if np.isnan(val_f) or np.isinf(val_f): val_f = 0.0
        if np.isnan(err_f) or np.isinf(err_f): err_f = 0.0
        latex_text = f"#splitline{{{val_f:.3f}}}{{#pm {err_f:.3f}}}"
        l = ROOT.TLatex(x_center, y_center, latex_text)
        l.SetTextSize(0.015); l.SetTextAlign(22); l.SetTextColor(ROOT.kBlack)
        hist.GetListOfFunctions().Add(l)

    def extract_target_stats(self, data_tot, data_mix, loc_tot, loc_mix, m_low, m_high, x_low, x_high):
        mask_tot = (data_tot["mass"] >= m_low) & (data_tot["mass"] < m_high) & (data_tot["xF"] >= x_low) & (data_tot["xF"] < x_high)
        mask_mix = (data_mix["mass"] >= m_low) & (data_mix["mass"] < m_high) & (data_mix["xF"] >= x_low) & (data_mix["xF"] < x_high)

        N_tot, N_mix = np.sum(mask_tot), np.sum(mask_mix)
        sum_mass_tot = np.sum(data_tot["mass"][mask_tot]) if N_tot > 0 else 0.0
        sum_mass_mix = np.sum(data_mix["mass"][mask_mix]) if N_mix > 0 else 0.0

        eff_reco_tot, err_reco_tot = data_tot["recoeff"][mask_tot], data_tot["recoeff_error"][mask_tot]
        eff_hodo_tot, err_hodo_tot = data_tot["hodoeff"][mask_tot], data_tot["hodoeff_error"][mask_tot]
        eff_reco_mix, err_reco_mix = data_mix["recoeff"][mask_mix], data_mix["recoeff_error"][mask_mix]
        eff_hodo_mix, err_hodo_mix = data_mix["hodoeff"][mask_mix], data_mix["hodoeff_error"][mask_mix]

        means = {
            "r_tot": self.get_correlated_mean_and_error(eff_reco_tot, err_reco_tot, loc_tot[mask_tot]),
            "h_tot": self.get_weighted_mean_and_error(eff_hodo_tot, err_hodo_tot),
            "r_mix": self.get_correlated_mean_and_error(eff_reco_mix, err_reco_mix, loc_mix[mask_mix]),
            "h_mix": self.get_weighted_mean_and_error(eff_hodo_mix, err_hodo_mix),
        }

        means["f_tot"] = (means["r_tot"][0] * means["h_tot"][0], np.sqrt((means["h_tot"][0] * means["r_tot"][1])**2 + (means["r_tot"][0] * means["h_tot"][1])**2))
        means["f_mix"] = (means["r_mix"][0] * means["h_mix"][0], np.sqrt((means["h_mix"][0] * means["r_mix"][1])**2 + (means["r_mix"][0] * means["h_mix"][1])**2))

        val_sig_eff, err_sig_eff = 0.0, 0.0
        diff_yield = N_tot - N_mix
        if diff_yield != 0:
            val_sig_eff = ((N_tot * means["f_tot"][0]) - (N_mix * means["f_mix"][0])) / diff_yield
            err_sig_eff = (1.0 / diff_yield) * np.sqrt((N_tot * means["f_tot"][1])**2 + (N_mix * means["f_mix"][1])**2)

        return N_tot, N_mix, sum_mass_tot, sum_mass_mix, val_sig_eff, err_sig_eff, diff_yield, means

    def fill_histograms(self, hists, root_x, root_y, N_tot, N_mix, m_center, x_center, sig_eff, err_sig_eff, diff_yield, means, final_centroid):
        err_N_tot, err_N_mix = np.sqrt(N_tot), np.sqrt(N_mix)

        val_corr_yield, err_corr_stat, err_corr_sys, err_corr_total = 0.0, 0.0, 0.0, 0.0
        if sig_eff > 0 and diff_yield > 0:
            val_corr_yield = diff_yield / sig_eff
            err_corr_stat = np.sqrt(N_tot + N_mix) / sig_eff
            err_corr_sys = val_corr_yield * (err_sig_eff / sig_eff)
            err_corr_total = np.sqrt(err_corr_stat**2 + err_corr_sys**2)

        def set_bin(key, val, err):
            hists[key].SetBinContent(root_x, root_y, val)
            hists[key].SetBinError(root_x, root_y, err)
            if "stat" not in key and "sys" not in key and "Centroid" not in key:
                self.add_latex_to_bin(hists[key], m_center, x_center, val, err)

        set_bin("Y_total", N_tot, err_N_tot); set_bin("Y_mix", N_mix, err_N_mix)
        set_bin("E_total_reco", means["r_tot"][0], means["r_tot"][1]); set_bin("E_mix_reco", means["r_mix"][0], means["r_mix"][1])
        set_bin("E_total_hodo", means["h_tot"][0], means["h_tot"][1]); set_bin("E_mix_hodo", means["h_mix"][0], means["h_mix"][1])
        set_bin("E_total_final", means["f_tot"][0], means["f_tot"][1]); set_bin("E_mix_final", means["f_mix"][0], means["f_mix"][1])
        set_bin("E_final_signal", sig_eff, err_sig_eff)
        set_bin("Y_corrected", val_corr_yield, err_corr_total); set_bin("Y_corrected_stat", val_corr_yield, err_corr_stat)
        set_bin("Y_corrected_sys", val_corr_yield, err_corr_sys)
        hists["Mass_Centroid"].SetBinContent(root_x, root_y, final_centroid)

    def save_2d_pdfs(self, hists_dict, label):
        for name, hist in hists_dict.items():
            if "stat" in name or "sys" in name or "Centroid" in name: continue 
            c = ROOT.TCanvas(f"c_{name}_{label}", "", 1200, 900)
            c.SetRightMargin(0.15); c.SetLeftMargin(0.12); c.SetBottomMargin(0.12); c.SetTickx(1); c.SetTicky(1)
            if hist.GetMinimum() < 0: hist.SetMinimum(0)
            hist.Draw("COLZ"); c.SaveAs(f"{name}_{label}.pdf"); c.Close()

    def build_response_matrix(self):
        """Builds separate 1D Mass Smearing Matrices per xF bin for both Messy and Clean MC."""
        if not self.mc_messy_dict and not self.mc_clean_dict:
            self.console.print("[yellow]WARNING: No MC files provided. Unfolding step will be skipped.[/yellow]")
            return

        self.rm_file = ROOT.TFile("DY_ResponseMatrices.root", "RECREATE")
        n_xf_bins = len(config.XF_BINS) - 1
        n_mass_bins = len(config.MASS_BINS) - 1
        
        def process_mc_dict(mc_dict, label_prefix, store_dict):
            for target_label, paths in mc_dict.items():
                if not paths: continue
                file_list = paths if isinstance(paths, list) else [paths]
                store_dict[target_label] = {}
                
                for i_x in range(n_xf_bins):
                    x_low = config.XF_BINS[i_x]; x_high = config.XF_BINS[i_x+1]
                    h_reco = ROOT.TH1D(f"h_mc_reco_{label_prefix}_{target_label}_xF_{i_x}", f"MC Reco {label_prefix} {target_label} xF [{x_low:.2f}, {x_high:.2f});Mass [GeV]", n_mass_bins, config.MASS_BINS)
                    h_true = ROOT.TH1D(f"h_mc_true_{label_prefix}_{target_label}_xF_{i_x}", f"MC Truth {label_prefix} {target_label} xF [{x_low:.2f}, {x_high:.2f});Mass [GeV]", n_mass_bins, config.MASS_BINS)
                    store_dict[target_label][i_x] = ROOT.RooUnfoldResponse(h_reco, h_true, f"DY_1D_Response_{label_prefix}_{target_label}_xF_{i_x}")

                for mc_file in file_list:
                    if not os.path.exists(mc_file): continue
                    try:
                        with uproot.open(mc_file) as f:
                            if "Tree" not in f: continue
                            mc_tree = f["Tree"] 
                            events = mc_tree.arrays(library="np")
                            mask_reco = self.get_cut_mask(events, is_mc=True)
                            
                            mass_true, xf_true = events["mMass"], events["mxF"]
                            mass_reco, xf_reco = events["mass"], events["xF"]
                            
                            xf_reco_indices = np.digitize(xf_reco, config.XF_BINS) - 1
                            
                            n_filled = 0
                            for i in range(len(mass_true)):
                                if mask_reco[i]:
                                    mr, xr = float(mass_reco[i]), float(xf_reco[i])
                                    mt, xt = float(mass_true[i]), float(xf_true[i])
                                    ix = xf_reco_indices[i]
                                    
                                    if 0 <= ix < n_xf_bins:
                                        if not math.isnan(mr) and not math.isnan(mt) and mr > 0.0 and mt > 0.0:
                                            store_dict[target_label][ix].Fill(mr, mt)
                                            n_filled += 1
                                            
                            self.console.print(f"[cyan][{target_label} - {label_prefix.upper()}] Matrices Trained: {n_filled} events survived cuts out of {len(mass_true)} generated.[/cyan]")
                    except Exception as e:
                        self.console.print(f"[bold red]Error building response matrix from {mc_file}: {e}[/bold red]")
                
                self.rm_file.cd()
                for i_x in range(n_xf_bins):
                    store_dict[target_label][i_x].Write()

        if self.mc_messy_dict: process_mc_dict(self.mc_messy_dict, "messy", self.response_matrices_messy)
        if self.mc_clean_dict: process_mc_dict(self.mc_clean_dict, "clean", self.response_matrices_clean)

    def process_kinematics(self):
        def get_loc_data(data_dict, is_mix=False):
            if is_mix and 'ptrk_D1' in data_dict and 'ntrk_D1' in data_dict: d1_vals = 0.5 * (data_dict['ptrk_D1'] + data_dict['ntrk_D1'])
            else: d1_vals = data_dict['D1']
            return np.digitize(d1_vals, self.x_curve) - 1

        os.makedirs("./MassBinCentroids", exist_ok=True)
        
        data_lh2_tot = self.get_concatenated_events(self.lh2_paths, "result")
        data_lh2_mix = self.get_concatenated_events(self.lh2_paths, "result_mix")
        data_ld2_tot = self.get_concatenated_events(self.ld2_paths, "result")
        data_ld2_mix = self.get_concatenated_events(self.ld2_paths, "result_mix")
        data_fl_tot = self.get_concatenated_events(self.flask_paths, "result")
        data_fl_mix = self.get_concatenated_events(self.flask_paths, "result_mix")
        
        loc_lh2_tot, loc_lh2_mix = get_loc_data(data_lh2_tot), get_loc_data(data_lh2_mix, True)
        loc_ld2_tot, loc_ld2_mix = get_loc_data(data_ld2_tot), get_loc_data(data_ld2_mix, True)
        loc_fl_tot, loc_fl_mix = get_loc_data(data_fl_tot), get_loc_data(data_fl_mix, True)

        dir_kin = self.get_or_create_dir(self.out_file, "Kinematics")

        dir_lh2 = self.get_or_create_dir(dir_kin, "LH2")
        dir_lh2.cd(); self.hists_lh2 = self.create_histograms("LH2")
        dir_ld2 = self.get_or_create_dir(dir_kin, "LD2")
        dir_ld2.cd(); self.hists_ld2 = self.create_histograms("LD2")
        dir_fl = self.get_or_create_dir(dir_kin, "Flask")
        dir_fl.cd(); self.hists_fl = self.create_histograms("Flask")

        colors = {"LD2_tot": ROOT.kRed, "LD2_mix": ROOT.kBlue, "LH2_tot": ROOT.kRed, "LH2_mix": ROOT.kBlue, "Fl_tot": ROOT.kGreen+2, "Fl_mix": ROOT.kOrange+1}

        dir_mass1d = self.get_or_create_dir(self.out_file, "Mass_Distributions_1D")
        dir_csv = self.get_or_create_dir(self.out_file, "CSV_Tables")

        for i_x in range(len(config.XF_BINS) - 1):
            x_low, x_high = config.XF_BINS[i_x], config.XF_BINS[i_x+1]
            x_center = (x_low + x_high) / 2.0
            root_y = i_x + 1

            csv_filename_lh2, csv_filename_ld2 = f"Table_Kinematics_LH2_xF_{x_low:.2f}_{x_high:.2f}.csv", f"Table_Kinematics_LD2_xF_{x_low:.2f}_{x_high:.2f}.csv"
            csv_rows_lh2, csv_rows_ld2 = [], []

            dir_mass1d.cd()
            h1_ld2_tot = ROOT.TH1D(f"h1_ld2_tot_{i_x}", "", len(config.MASS_BINS)-1, config.MASS_BINS)
            h1_ld2_mix = ROOT.TH1D(f"h1_ld2_mix_{i_x}", "", len(config.MASS_BINS)-1, config.MASS_BINS)
            h1_lh2_tot = ROOT.TH1D(f"h1_lh2_tot_{i_x}", "", len(config.MASS_BINS)-1, config.MASS_BINS)
            h1_lh2_mix = ROOT.TH1D(f"h1_lh2_mix_{i_x}", "", len(config.MASS_BINS)-1, config.MASS_BINS)
            h1_fl_tot = ROOT.TH1D(f"h1_fl_tot_{i_x}", "", len(config.MASS_BINS)-1, config.MASS_BINS)
            h1_fl_mix = ROOT.TH1D(f"h1_fl_mix_{i_x}", "", len(config.MASS_BINS)-1, config.MASS_BINS)

            h1s = {"LD2_tot": h1_ld2_tot, "LD2_mix": h1_ld2_mix, "LH2_tot": h1_lh2_tot, "LH2_mix": h1_lh2_mix, "Fl_tot": h1_fl_tot, "Fl_mix": h1_fl_mix}
            for k, h in h1s.items(): h.SetLineColor(colors[k]); h.SetLineWidth(2)
            
            latex_draws_lh2, latex_draws_ld2 = [], []
            max_label_y_lh2, max_label_y_ld2 = 0.5, 0.5

            for i_m in range(len(config.MASS_BINS) - 1):
                m_low, m_high = config.MASS_BINS[i_m], config.MASS_BINS[i_m+1]
                m_center = (m_low + m_high) / 2.0
                root_x = i_m + 1

                N_l2_t, N_l2_m, M_l2_t, M_l2_m, eps_ld2, err_eps_ld2, diff_l2, mean_l2 = self.extract_target_stats(data_ld2_tot, data_ld2_mix, loc_ld2_tot, loc_ld2_mix, m_low, m_high, x_low, x_high)
                N_lh_t, N_lh_m, M_lh_t, M_lh_m, eps_lh2, err_eps_lh2, diff_lh, mean_lh = self.extract_target_stats(data_lh2_tot, data_lh2_mix, loc_lh2_tot, loc_lh2_mix, m_low, m_high, x_low, x_high)
                N_fl_t, N_fl_m, M_fl_t, M_fl_m, eps_fl, err_eps_fl, diff_fl, mean_fl = self.extract_target_stats(data_fl_tot, data_fl_mix, loc_fl_tot, loc_fl_mix, m_low, m_high, x_low, x_high)

                h1_ld2_tot.SetBinContent(root_x, N_l2_t); h1_ld2_mix.SetBinContent(root_x, N_l2_m)
                h1_lh2_tot.SetBinContent(root_x, N_lh_t); h1_lh2_mix.SetBinContent(root_x, N_lh_m)
                h1_fl_tot.SetBinContent(root_x, N_fl_t); h1_fl_mix.SetBinContent(root_x, N_fl_m)

                if eps_lh2 > 0: cY_lh_t, cM_lh_t, cY_lh_m, cM_lh_m = N_lh_t / eps_lh2, M_lh_t / eps_lh2, N_lh_m / eps_lh2, M_lh_m / eps_lh2
                else: cY_lh_t, cM_lh_t, cY_lh_m, cM_lh_m = 0.0, 0.0, 0.0, 0.0

                if eps_ld2 > 0: cY_l2_t, cM_l2_t, cY_l2_m, cM_l2_m = N_l2_t / eps_ld2, M_l2_t / eps_ld2, N_l2_m / eps_ld2, M_l2_m / eps_ld2
                else: cY_l2_t, cM_l2_t, cY_l2_m, cM_l2_m = 0.0, 0.0, 0.0, 0.0

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
                cent_LH2 = num_LH2 / den_LH2 if den_LH2 != 0 else m_center

                num_LD2 = (cM_l2_t - cM_l2_m) - (cM_fl_t_l2 - cM_fl_m_l2) - (config.LH2_TO_LD2_NORM * num_LH2)
                den_LD2 = (cY_l2_t - cY_l2_m) - (cY_fl_t_l2 - cY_fl_m_l2) - (config.LH2_TO_LD2_NORM * den_LH2)
                cent_LD2 = num_LD2 / den_LD2 if den_LD2 != 0 else m_center

                self.fill_histograms(self.hists_ld2, root_x, root_y, N_l2_t, N_l2_m, m_center, x_center, eps_ld2, err_eps_ld2, diff_l2, mean_l2, cent_LD2)
                self.fill_histograms(self.hists_lh2, root_x, root_y, N_lh_t, N_lh_m, m_center, x_center, eps_lh2, err_eps_lh2, diff_lh, mean_lh, cent_LH2)
                
                num_fl = (M_fl_t - M_fl_m) / (eps_fl if eps_fl > 0 else 1e-9)
                den_fl = (N_fl_t - N_fl_m) / (eps_fl if eps_fl > 0 else 1e-9)
                cent_fl = num_fl / den_fl if den_fl != 0 else m_center
                self.fill_histograms(self.hists_fl, root_x, root_y, N_fl_t, N_fl_m, m_center, x_center, eps_fl, err_eps_fl, diff_fl, mean_fl, cent_fl)

                jitter_x = [-0.05, 0.05]
                comp_data_lh2 = [("LH2_tot", cM_lh_t, N_lh_t), ("LH2_mix", cM_lh_m, N_lh_m), ("Fl_tot", cM_fl_t_lh, N_fl_t), ("Fl_mix", cM_fl_m_lh, N_fl_m)]
                valid_lh2 = [(lbl, cM, N) for lbl, cM, N in comp_data_lh2 if N > 0]
                valid_lh2.sort(key=lambda x: x[2])
                
                last_y = 1e-9
                for j, (lbl, cM, N) in enumerate(valid_lh2):
                    desired_y = N * 1.5; actual_y = max(desired_y, last_y * 2.5)
                    x_pos = m_center + jitter_x[j % 2]
                    txt = ROOT.TLatex(x_pos, actual_y, f"#splitline{{{cM:.2f}}}{{({N})}}")
                    txt.SetTextSize(0.025); txt.SetTextColor(colors[lbl]); txt.SetTextAlign(22)
                    latex_draws_lh2.append(txt); last_y = actual_y
                    max_label_y_lh2 = max(max_label_y_lh2, actual_y)

                comp_data_ld2 = [("LD2_tot", cM_l2_t, N_l2_t), ("LD2_mix", cM_l2_m, N_l2_m), ("Fl_tot", cM_fl_t_l2, N_fl_t), ("Fl_mix", cM_fl_m_l2, N_fl_m)]
                valid_ld2 = [(lbl, cM, N) for lbl, cM, N in comp_data_ld2 if N > 0]
                valid_ld2.sort(key=lambda x: x[2])
                
                last_y = 1e-9
                for j, (lbl, cM, N) in enumerate(valid_ld2):
                    desired_y = N * 1.5; actual_y = max(desired_y, last_y * 2.5)
                    x_pos = m_center + jitter_x[j % 2]
                    txt = ROOT.TLatex(x_pos, actual_y, f"#splitline{{{cM:.2f}}}{{({N})}}")
                    txt.SetTextSize(0.025); txt.SetTextColor(colors[lbl]); txt.SetTextAlign(22)
                    latex_draws_ld2.append(txt); last_y = actual_y
                    max_label_y_ld2 = max(max_label_y_ld2, actual_y)

                csv_rows_lh2.append({
                    "Mass Bin": f"[{m_low:.2f}, {m_high:.2f})", "Mass Center": m_center, "LH2 Mass Bin Average": cent_LH2,
                    "N_LH2_total": N_lh_t, "N_LH2_mixed": N_lh_m, "N_flask_total": N_fl_t, "N_flask_mixed": N_fl_m,
                    "eps_LH2": eps_lh2, "eps_LD2": eps_ld2, "eps_flask": eps_fl,
                    "Corrected Total Mass LH2 (num)": num_LH2, "Corrected Yield LH2 (denom)": den_LH2
                })

                csv_rows_ld2.append({
                    "Mass Bin": f"[{m_low:.2f}, {m_high:.2f})", "Mass Center": m_center, "LD2 Mass Bin Average": cent_LD2,
                    "N_LD2_total": N_l2_t, "N_LD2_mixed": N_l2_m, "N_LH2_total": N_lh_t, "N_LH2_mixed": N_lh_m,
                    "N_flask_total": N_fl_t, "N_flask_mixed": N_fl_m,
                    "eps_LH2": eps_lh2, "eps_LD2": eps_ld2, "eps_flask": eps_fl,
                    "Corrected Total Mass LD2 (num)": num_LD2, "Corrected Yield LD2 (denom)": den_LD2
                })

            c_lh2 = ROOT.TCanvas(f"c_1d_mass_lh2_{i_x}", "", 1000, 700)
            c_lh2.SetRightMargin(0.05); c_lh2.SetLogy(); c_lh2.SetTickx(1); c_lh2.SetTicky(1) 
            
            max_y_hist_lh2 = max([h1_lh2_tot.GetMaximum(), h1_lh2_mix.GetMaximum(), h1_fl_tot.GetMaximum(), h1_fl_mix.GetMaximum()])
            overall_max_lh2 = max(max_y_hist_lh2, max_label_y_lh2)
            if overall_max_lh2 <= 0: overall_max_lh2 = 1.0
            h1_lh2_tot.SetMaximum(overall_max_lh2 * 10.0); h1_lh2_tot.SetMinimum(0.5) 
            h1_lh2_tot.SetTitle(f"Mass Distributions (LH2) {x_low:.2f} <= x_{{F}} < {x_high:.2f};Mass [GeV];Counts")
            h1_lh2_tot.GetXaxis().CenterTitle(True); h1_lh2_tot.GetYaxis().CenterTitle(True)
            h1_lh2_tot.Draw("HIST"); h1_lh2_mix.Draw("HIST SAME"); h1_fl_tot.Draw("HIST SAME"); h1_fl_mix.Draw("HIST SAME")
            for l in latex_draws_lh2: l.Draw()
            leg_lh2 = ROOT.TLegend(0.75, 0.72, 0.93, 0.88); leg_lh2.SetBorderSize(0)
            leg_lh2.AddEntry(h1_lh2_tot, "LH2 Total", "l"); leg_lh2.AddEntry(h1_lh2_mix, "LH2 Mix", "l")
            leg_lh2.AddEntry(h1_fl_tot,  "Flask Total", "l"); leg_lh2.AddEntry(h1_fl_mix,  "Flask Mix", "l")
            leg_lh2.Draw()
            c_lh2.SaveAs(f"./MassBinCentroids/MassDist_LH2_xF_{x_low:.2f}_{x_high:.2f}.pdf"); c_lh2.Close()

            c_ld2 = ROOT.TCanvas(f"c_1d_mass_ld2_{i_x}", "", 1000, 700)
            c_ld2.SetRightMargin(0.05); c_ld2.SetLogy(); c_ld2.SetTickx(1); c_ld2.SetTicky(1) 
            
            max_y_hist_ld2 = max([h1_ld2_tot.GetMaximum(), h1_ld2_mix.GetMaximum(), h1_fl_tot.GetMaximum(), h1_fl_mix.GetMaximum()])
            overall_max_ld2 = max(max_y_hist_ld2, max_label_y_ld2)
            if overall_max_ld2 <= 0: overall_max_ld2 = 1.0
            h1_ld2_tot.SetMaximum(overall_max_ld2 * 10.0); h1_ld2_tot.SetMinimum(0.5)
            h1_ld2_tot.SetTitle(f"Mass Distributions (LD2) {x_low:.2f} <= x_{{F}} < {x_high:.2f};Mass [GeV];Counts")
            h1_ld2_tot.GetXaxis().CenterTitle(True); h1_ld2_tot.GetYaxis().CenterTitle(True)
            h1_ld2_tot.Draw("HIST"); h1_ld2_mix.Draw("HIST SAME"); h1_fl_tot.Draw("HIST SAME"); h1_fl_mix.Draw("HIST SAME")
            for l in latex_draws_ld2: l.Draw()
            leg_ld2 = ROOT.TLegend(0.75, 0.72, 0.93, 0.88); leg_ld2.SetBorderSize(0)
            leg_ld2.AddEntry(h1_ld2_tot, "LD2 Total", "l"); leg_ld2.AddEntry(h1_ld2_mix, "LD2 Mix", "l")
            leg_ld2.AddEntry(h1_fl_tot,  "Flask Total", "l"); leg_ld2.AddEntry(h1_fl_mix,  "Flask Mix", "l")
            leg_ld2.Draw()
            c_ld2.SaveAs(f"./MassBinCentroids/MassDist_LD2_xF_{x_low:.2f}_{x_high:.2f}.pdf"); c_ld2.Close()

            with open(csv_filename_lh2, "w", newline='') as f:
                writer = csv.DictWriter(f, fieldnames=csv_rows_lh2[0].keys())
                writer.writeheader(); writer.writerows(csv_rows_lh2)
            with open(csv_filename_ld2, "w", newline='') as f:
                writer = csv.DictWriter(f, fieldnames=csv_rows_ld2[0].keys())
                writer.writeheader(); writer.writerows(csv_rows_ld2)
                
            dir_csv.cd(); ROOT.TMacro(csv_filename_lh2).Write(csv_filename_lh2); ROOT.TMacro(csv_filename_ld2).Write(csv_filename_ld2)

        self.save_2d_pdfs(self.hists_lh2, "LH2"); self.save_2d_pdfs(self.hists_ld2, "LD2"); self.save_2d_pdfs(self.hists_fl, "Flask")

    def generate_subtracted_plot(self, hists_target, hists_flask, flask_norm, target_label):
        dir_sub = self.get_or_create_dir(self.out_file, f"Subtracted_Plots_{target_label}")
        dir_sub.cd()
        
        h_target_stat, h_target_sys = hists_target["Y_corrected_stat"], hists_target["Y_corrected_sys"]
        h_flask_stat, h_flask_sys = hists_flask["Y_corrected_stat"], hists_flask["Y_corrected_sys"]
        
        name = f"Y_corrected_Subtracted_{target_label}"
        h_sub = ROOT.TH2D(name, f"Corrected Yield ({target_label} - Flask);Mass [GeV];x_{{F}}", 
                          len(config.MASS_BINS)-1, config.MASS_BINS, len(config.XF_BINS)-1, config.XF_BINS)
        h_sub.Sumw2(); h_sub.SetStats(0); h_sub.GetXaxis().CenterTitle(); h_sub.GetYaxis().CenterTitle()

        h_sub_stat, h_sub_sys = h_sub.Clone(f"{name}_stat"), h_sub.Clone(f"{name}_sys")
        h_sub_centroid = hists_target["Mass_Centroid"].Clone(f"{name}_Mass_Centroid")

        for i_m in range(len(config.MASS_BINS) - 1):
            m_center = (config.MASS_BINS[i_m] + config.MASS_BINS[i_m+1]) / 2.0
            bin_x = i_m + 1
            for i_x in range(len(config.XF_BINS) - 1):
                x_center = (config.XF_BINS[i_x] + config.XF_BINS[i_x+1]) / 2.0
                bin_y = i_x + 1
                
                y_target, e_target_stat, e_target_sys = h_target_stat.GetBinContent(bin_x, bin_y), h_target_stat.GetBinError(bin_x, bin_y), h_target_sys.GetBinError(bin_x, bin_y)
                y_flask, e_flask_stat, e_flask_sys = h_flask_stat.GetBinContent(bin_x, bin_y), h_flask_stat.GetBinError(bin_x, bin_y), h_flask_sys.GetBinError(bin_x, bin_y)
                
                val_sub = y_target - (flask_norm * y_flask)
                err_sub_stat = np.sqrt(e_target_stat**2 + (flask_norm * e_flask_stat)**2)
                err_sub_sys = np.sqrt(e_target_sys**2 + (flask_norm * e_flask_sys)**2)
                err_sub_total = np.sqrt(err_sub_stat**2 + err_sub_sys**2)
                
                h_sub.SetBinContent(bin_x, bin_y, val_sub); h_sub.SetBinError(bin_x, bin_y, err_sub_total)
                self.add_latex_to_bin(h_sub, m_center, x_center, val_sub, err_sub_total)
                h_sub_stat.SetBinContent(bin_x, bin_y, val_sub); h_sub_stat.SetBinError(bin_x, bin_y, err_sub_stat)
                h_sub_sys.SetBinContent(bin_x, bin_y, val_sub); h_sub_sys.SetBinError(bin_x, bin_y, err_sub_sys)

        c = ROOT.TCanvas(f"c_{name}", f"c_{name}", 1200, 900)
        c.SetRightMargin(0.15); c.SetLeftMargin(0.12); c.SetBottomMargin(0.12); c.SetTickx(1); c.SetTicky(1)
        
        local_min, local_max, has_data = sys.float_info.max, -sys.float_info.max, False
        for ix in range(1, h_sub.GetNbinsX() + 1):
            for iy in range(1, h_sub.GetNbinsY() + 1):
                val = h_sub.GetBinContent(ix, iy)
                if val != 0.0:
                    has_data = True
                    if val < local_min: local_min = val
                    if val > local_max: local_max = val
        if has_data: h_sub.SetMinimum(0.0 if local_min < 0 else local_min); h_sub.SetMaximum(local_max)
        else: h_sub.SetMinimum(0.0); h_sub.SetMaximum(1.0)
            
        h_sub.Draw("COLZ"); c.SaveAs(f"{name}.pdf"); c.Close()
        return {"stat": h_sub_stat, "sys": h_sub_sys, "centroid": h_sub_centroid}

    def calculate_and_plot_cross_section(self, h_sub_dict, target_label, global_constant):
        dir_xsec = self.get_or_create_dir(self.out_file, f"CrossSections_{target_label}")
        dir_xsec.cd()
        
        h_sub_stat = h_sub_dict["stat"]
        h_sub_sys  = h_sub_dict["sys"]
        h_centroid = h_sub_dict["centroid"]

        acc_path = config.ACCEPTANCE_FILE
        psip_path = "All_PsiP_Contaminations.root" 
        
        if target_label == "LD2":
            theory_ct18_path = "CT18_xFnew_d_1sigma.root"
            theory_nnpdf_path = "NNPDF40_xFnew_d.root"
        else:
            theory_ct18_path = "CT18_xFnew_p_1sigma.root"
            theory_nnpdf_path = "NNPDF40_xFnew_p.root"
            
        acc_file, f_ct18, f_nnpdf, f_psip = None, None, None, None

        try:
            if os.path.exists(acc_path): acc_file = ROOT.TFile.Open(acc_path)
            if os.path.exists(theory_ct18_path): f_ct18 = ROOT.TFile.Open(theory_ct18_path)
            if os.path.exists(theory_nnpdf_path): f_nnpdf = ROOT.TFile.Open(theory_nnpdf_path)
            if os.path.exists(psip_path): f_psip = ROOT.TFile.Open(psip_path)
        except Exception as e:
            self.console.print(f"[yellow]Warning: Could not open one or more auxiliary files: {e}[/yellow]")

        n_xf_bins = len(config.XF_BINS) - 1
        n_mass_bins = len(config.MASS_BINS) - 1
        
        latex_xsec_table_content = r"""\begingroup
\renewcommand{\arraystretch}{1.8}
\begin{longtable}{|c|c|c|c|c|c|c|}
\caption{Final TARGET_LABEL Unfolded Cross Section Table} \label{tab:cross_section_final_TARGET_LABEL} \\
\hline
\textbf{xF bin} & \textbf{xF cent} & \textbf{Mass bin} & \textbf{Mass cent} & \textbf{$\sigma_{Raw}$} & \textbf{$\sigma_{Unf (Messy)}$} & \textbf{$\sigma_{Unf (Clean)}$} \\
\hline
\endfirsthead
\multicolumn{7}{c}%
{{\bfseries \tablename\ \thetable{} -- continued from previous page}} \\
\hline
\textbf{xF bin} & \textbf{xF cent} & \textbf{Mass bin} & \textbf{Mass cent} & \textbf{$\sigma_{Raw}$} & \textbf{$\sigma_{Unf (Messy)}$} & \textbf{$\sigma_{Unf (Clean)}$} \\
\hline
\endhead
\hline \multicolumn{7}{|r|}{{Continued on next page}} \\ \hline
\endfoot
\hline
\endlastfoot
""".replace("TARGET_LABEL", target_label)

        def draw_and_save_canvas(plot_type, g_raw_xsec, g_raw_sys, g_unf_m_xsec, g_unf_m_sys, g_unf_c_xsec, g_unf_c_sys, xf_min, xf_max, xf_bin_index, theory_idx, y_min_data, y_max_data, ct18_file, nnpdf_file):
            c_xsec = ROOT.TCanvas(f"c_xsec_{target_label}_{xf_bin_index}_{plot_type}", "", 800, 600)
            c_xsec.SetLogy(); c_xsec.SetTickx(1); c_xsec.SetTicky(1)
            
            mg = ROOT.TMultiGraph()
            mg.SetTitle(f";Mass [GeV];M^{{3}} d^{{2}}\\sigma / dM dx_{{F}} [nb GeV^{{2}}/Nucleus]")
            
            if y_min_data < y_max_data and y_max_data > 0:
                mg.SetMinimum(max(y_min_data * 0.2, 1e-7)) 
                mg.SetMaximum(y_max_data * 15.0)
            else:
                mg.SetMinimum(1e-6); mg.SetMaximum(3.0)
            
            leg = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
            leg.SetBorderSize(0)

            gr_name = f"gr_xFbin{theory_idx}"
            if ct18_file and theory_idx >= 0:
                g_ct18 = ct18_file.Get(gr_name)
                if g_ct18:
                    g_ct18_clone = g_ct18.Clone()
                    g_ct18_clone.SetLineColor(ROOT.kGreen + 2); g_ct18_clone.SetFillColorAlpha(ROOT.kGreen - 5, 0.5); g_ct18_clone.SetFillStyle(3002)
                    mg.Add(g_ct18_clone, "L3"); leg.AddEntry(g_ct18_clone, "CT18 NLO", "lf") 
            
            if nnpdf_file and theory_idx >= 0:
                g_nnpdf = nnpdf_file.Get(gr_name)
                if g_nnpdf:
                    g_nnpdf_clone = g_nnpdf.Clone()
                    g_nnpdf_clone.SetLineColor(ROOT.kMagenta + 2); g_nnpdf_clone.SetFillColorAlpha(ROOT.kMagenta - 9, 0.5); g_nnpdf_clone.SetFillStyle(3002)
                    mg.Add(g_nnpdf_clone, "L3"); leg.AddEntry(g_nnpdf_clone, "NNPDF4.0 NLO", "lf") 

            if g_raw_sys.GetN() > 0:
                g_raw_sys.SetMarkerSize(0); g_raw_sys.SetLineColor(ROOT.kBlue); g_raw_sys.SetFillColorAlpha(ROOT.kAzure - 9, 0.5); g_raw_sys.SetFillStyle(1001)
                mg.Add(g_raw_sys, "2"); leg.AddEntry(g_raw_sys, "Raw Syst. Unc.", "f")

                g_raw_xsec.SetMarkerStyle(24); g_raw_xsec.SetMarkerColor(ROOT.kBlue); g_raw_xsec.SetLineColor(ROOT.kBlue)
                mg.Add(g_raw_xsec, "P"); leg.AddEntry(g_raw_xsec, f"Raw Data ({target_label})", "lep")

            if g_unf_m_sys.GetN() > 0:
                g_unf_m_sys.SetMarkerSize(0); g_unf_m_sys.SetLineColor(ROOT.kRed); g_unf_m_sys.SetFillColorAlpha(ROOT.kPink - 9, 0.5); g_unf_m_sys.SetFillStyle(1001)
                mg.Add(g_unf_m_sys, "2"); leg.AddEntry(g_unf_m_sys, "Unf. (Messy) Syst.", "f")

                g_unf_m_xsec.SetMarkerStyle(20); g_unf_m_xsec.SetMarkerColor(ROOT.kRed); g_unf_m_xsec.SetLineColor(ROOT.kRed)
                mg.Add(g_unf_m_xsec, "P"); leg.AddEntry(g_unf_m_xsec, f"Unfolded Messy", "lep")
                
            if g_unf_c_sys.GetN() > 0:
                g_unf_c_sys.SetMarkerSize(0); g_unf_c_sys.SetLineColor(ROOT.kBlack); g_unf_c_sys.SetFillColorAlpha(ROOT.kGray, 0.5); g_unf_c_sys.SetFillStyle(1001)
                mg.Add(g_unf_c_sys, "2"); leg.AddEntry(g_unf_c_sys, "Unf. (Clean) Syst.", "f")

                g_unf_c_xsec.SetMarkerStyle(21); g_unf_c_xsec.SetMarkerColor(ROOT.kBlack); g_unf_c_xsec.SetLineColor(ROOT.kBlack)
                mg.Add(g_unf_c_xsec, "P"); leg.AddEntry(g_unf_c_xsec, f"Unfolded Clean", "lep")
            
            mg.Draw("A"); mg.GetXaxis().CenterTitle(); mg.GetYaxis().CenterTitle(); mg.GetXaxis().SetLimits(3.9, 10.0)
            leg.Draw()

            target_prefix = "pp" if target_label == "LH2" else "pd" if target_label == "LD2" else target_label
            internal_title = ROOT.TLatex()
            internal_title.SetNDC(True); internal_title.SetTextFont(42); internal_title.SetTextSize(0.04); internal_title.SetTextAlign(13)
            internal_title.DrawLatex(0.14, 0.86, f"Drell-Yan in {target_prefix} at {xf_min:.2f} #leq x_{{F}} < {xf_max:.2f}")

            plot_name = f"CrossSection_{target_label}_xF_{xf_min:.2f}_{xf_max:.2f}_{plot_type}.pdf"
            c_xsec.SaveAs(plot_name)
            
            dir_xsec.cd()
            if g_raw_xsec.GetN() > 0:
                g_raw_xsec.Write(f"g_xsec_{target_label}_{xf_bin_index}_Raw_{plot_type}")
                g_raw_sys.Write(f"g_sys_{target_label}_{xf_bin_index}_Raw_{plot_type}")
            if g_unf_m_xsec.GetN() > 0:
                g_unf_m_xsec.Write(f"g_xsec_{target_label}_{xf_bin_index}_UnfMessy_{plot_type}")
                g_unf_m_sys.Write(f"g_sys_{target_label}_{xf_bin_index}_UnfMessy_{plot_type}")
            if g_unf_c_xsec.GetN() > 0:
                g_unf_c_xsec.Write(f"g_xsec_{target_label}_{xf_bin_index}_UnfClean_{plot_type}")
                g_unf_c_sys.Write(f"g_sys_{target_label}_{xf_bin_index}_UnfClean_{plot_type}")

            c_xsec.Write(f"c_xsec_{target_label}_{xf_bin_index}_{plot_type}")
            c_xsec.Close()

        for i_x in range(n_xf_bins):
            root_bin_y = i_x + 1
            xf_min, xf_max = config.XF_BINS[i_x], config.XF_BINS[i_x+1]
            xf_center = (xf_min + xf_max) / 2.0

            # --- Extract 1D Mass Slice for Unfolding ---
            h_raw_1d = ROOT.TH1D(f"h_raw_1d_{target_label}_xF_{i_x}", "", n_mass_bins, config.MASS_BINS)
            for i_m in range(n_mass_bins):
                val = h_sub_stat.GetBinContent(i_m + 1, root_bin_y)
                err = h_sub_stat.GetBinError(i_m + 1, root_bin_y)
                if val <= 0:
                    h_raw_1d.SetBinContent(i_m + 1, 0.0)
                    h_raw_1d.SetBinError(i_m + 1, max(err, 1e-6))
                else:
                    h_raw_1d.SetBinContent(i_m + 1, val)
                    h_raw_1d.SetBinError(i_m + 1, err)

            # Unfold Messy
            resp_matrix_messy = self.response_matrices_messy.get(target_label, {}).get(i_x)
            if resp_matrix_messy is not None and resp_matrix_messy.Htruth().Integral() > 0:
                unfold_m = ROOT.RooUnfoldBayes(resp_matrix_messy, h_raw_1d, config.UNFOLDING_ITERATIONS)
                h_unf_m = unfold_m.Hunfold()
            else:
                self.console.print(f"[yellow]WARNING: Skipping Messy Unfolding for {target_label} xF bin {i_x}. Matrix empty.[/yellow]")
                h_unf_m = h_raw_1d.Clone(f"h_unf_messy_fallback_{target_label}_{i_x}")
                
            # Unfold Clean
            resp_matrix_clean = self.response_matrices_clean.get(target_label, {}).get(i_x)
            if resp_matrix_clean is not None and resp_matrix_clean.Htruth().Integral() > 0:
                unfold_c = ROOT.RooUnfoldBayes(resp_matrix_clean, h_raw_1d, config.UNFOLDING_ITERATIONS)
                h_unf_c = unfold_c.Hunfold()
            else:
                self.console.print(f"[yellow]WARNING: Skipping Clean Unfolding for {target_label} xF bin {i_x}. Matrix empty.[/yellow]")
                h_unf_c = h_raw_1d.Clone(f"h_unf_clean_fallback_{target_label}_{i_x}")

            acc_folder = f"mass_sliced_by_xF_bin{i_x}"
            acc_hist_name = f"{acc_folder}/h_ratio_{target_label}_{acc_folder}"
            h_acc = acc_file.Get(acc_hist_name) if acc_file else None
            if acc_file and not h_acc: 
                h_acc = acc_file.Get(f"{acc_folder}/h_ratio_LH2_{acc_folder}")
            
            # The older files (Theory, PsiP) still lack the -0.05 bin, so they need a -1 shift
            theory_idx = i_x - 1 

            best_psip_idx = -1
            best_psip_dist = float('inf')
            if f_psip:
                for b_idx in range(25):
                    tmp_p = f_psip.Get(f"hRatio_PsiP_DY_xF_{b_idx}")
                    if tmp_p:
                        old_center = 0.025 + b_idx * 0.05
                        dist = abs(old_center - xf_center)
                        if dist < best_psip_dist:
                            best_psip_dist = dist
                            best_psip_idx = b_idx
            
            h_ratio_psip = None
            if best_psip_idx >= 0:
                h_ratio_psip = f_psip.Get(f"hRatio_PsiP_DY_xF_{best_psip_idx}")

            dir_xsec.cd()
            g_raw_cent_xsec = ROOT.TGraphErrors(); g_raw_cent_xsec.SetName(f"g_xsec_{target_label}_xF_{i_x}_Raw_Centroid")
            g_raw_cent_sys  = ROOT.TGraphErrors(); g_raw_cent_sys.SetName(f"g_sys_{target_label}_xF_{i_x}_Raw_Centroid")
            g_raw_geo_xsec  = ROOT.TGraphErrors(); g_raw_geo_xsec.SetName(f"g_xsec_{target_label}_xF_{i_x}_Raw_GeoCenter")
            g_raw_geo_sys   = ROOT.TGraphErrors(); g_raw_geo_sys.SetName(f"g_sys_{target_label}_xF_{i_x}_Raw_GeoCenter")
            
            g_unf_m_cent_xsec = ROOT.TGraphErrors(); g_unf_m_cent_xsec.SetName(f"g_xsec_{target_label}_xF_{i_x}_UnfMessy_Centroid")
            g_unf_m_cent_sys  = ROOT.TGraphErrors(); g_unf_m_cent_sys.SetName(f"g_sys_{target_label}_xF_{i_x}_UnfMessy_Centroid")
            g_unf_m_geo_xsec  = ROOT.TGraphErrors(); g_unf_m_geo_xsec.SetName(f"g_xsec_{target_label}_xF_{i_x}_UnfMessy_GeoCenter")
            g_unf_m_geo_sys   = ROOT.TGraphErrors(); g_unf_m_geo_sys.SetName(f"g_sys_{target_label}_xF_{i_x}_UnfMessy_GeoCenter")

            g_unf_c_cent_xsec = ROOT.TGraphErrors(); g_unf_c_cent_xsec.SetName(f"g_xsec_{target_label}_xF_{i_x}_UnfClean_Centroid")
            g_unf_c_cent_sys  = ROOT.TGraphErrors(); g_unf_c_cent_sys.SetName(f"g_sys_{target_label}_xF_{i_x}_UnfClean_Centroid")
            g_unf_c_geo_xsec  = ROOT.TGraphErrors(); g_unf_c_geo_xsec.SetName(f"g_xsec_{target_label}_xF_{i_x}_UnfClean_GeoCenter")
            g_unf_c_geo_sys   = ROOT.TGraphErrors(); g_unf_c_geo_sys.SetName(f"g_sys_{target_label}_xF_{i_x}_UnfClean_GeoCenter")

            pt_idx_raw, pt_idx_unfm, pt_idx_unfc = 0, 0, 0
            y_min_raw, y_max_raw = sys.float_info.max, -sys.float_info.max
            y_min_unf, y_max_unf = sys.float_info.max, -sys.float_info.max
            
            for i_m in range(n_mass_bins):
                root_bin_x = i_m + 1
                mass_min, mass_max = config.MASS_BINS[i_m], config.MASS_BINS[i_m+1]
                geometric_center = (mass_min + mass_max) / 2.0
                mass_width = mass_max - mass_min
                
                actual_mass_center = h_centroid.GetBinContent(root_bin_x, root_bin_y)
                if actual_mass_center < mass_min or actual_mass_center > mass_max:
                    actual_mass_center = geometric_center 
                
                Y_raw = h_sub_stat.GetBinContent(root_bin_x, root_bin_y)
                Y_raw_err = h_sub_stat.GetBinError(root_bin_x, root_bin_y)
                Y_sys_err = h_sub_sys.GetBinError(root_bin_x, root_bin_y)

                Y_unf_m = h_unf_m.GetBinContent(root_bin_x)
                Y_unf_m_err = h_unf_m.GetBinError(root_bin_x)

                Y_unf_c = h_unf_c.GetBinContent(root_bin_x)
                Y_unf_c_err = h_unf_c.GetBinError(root_bin_x)
                
                # Dynamic Acceptance via Nearest Neighbor mapping for Mass
                acceptance, acceptance_err = 1.0, 0.0
                if h_acc:
                    acc_bin = h_acc.FindBin(actual_mass_center)
                    acceptance = h_acc.GetBinContent(acc_bin)
                    acceptance_err = h_acc.GetBinError(acc_bin)
                    
                    if acceptance <= 0:
                        best_dist = float('inf')
                        for b in range(1, h_acc.GetNbinsX() + 1):
                            val = h_acc.GetBinContent(b)
                            if val > 0:
                                dist = abs(h_acc.GetBinCenter(b) - actual_mass_center)
                                if dist < best_dist:
                                    best_dist = dist
                                    acceptance = val
                                    acceptance_err = h_acc.GetBinError(b)
                
                psip_ratio = 0.0
                if h_ratio_psip:
                    ratio_bin = h_ratio_psip.FindBin(actual_mass_center)
                    psip_ratio = h_ratio_psip.GetBinContent(ratio_bin)
                    if psip_ratio > 1.0 or psip_ratio <= 0.0:
                        best_dist = float('inf')
                        for b in range(1, h_ratio_psip.GetNbinsX() + 1):
                            val = h_ratio_psip.GetBinContent(b)
                            if 0.0 < val <= 1.0:
                                dist = abs(h_ratio_psip.GetBinCenter(b) - actual_mass_center)
                                if dist < best_dist:
                                    best_dist = dist
                                    psip_ratio = val

                xsec_raw, sys_fraction = 0.0, 0.0
                scaled_xsec_raw, scaled_stat_raw, scaled_sys_raw = 0.0, 0.0, 0.0
                
                # Format variables for LaTeX table
                str_raw, str_unfm, str_unfc = "-", "-", "-"

                if Y_raw > 0 and acceptance > 0:
                    xsec_raw = (global_constant * Y_raw) / (mass_width * acceptance)
                    stat_unc_raw = (Y_raw_err / Y_raw) * xsec_raw
                    sys_acc = (acceptance_err / acceptance) * xsec_raw
                    sys_tot_yield = (Y_sys_err / Y_raw) * xsec_raw
                    sys_psip_cont = psip_ratio * xsec_raw
                    total_sys_raw = np.sqrt(sys_acc**2 + sys_tot_yield**2 + sys_psip_cont**2)
                    
                    sys_fraction = total_sys_raw / xsec_raw
                    
                    scaled_xsec_raw = xsec_raw * (actual_mass_center**3)
                    scaled_stat_raw = stat_unc_raw * (actual_mass_center**3)
                    scaled_sys_raw = total_sys_raw * (actual_mass_center**3)

                    str_raw = f"${scaled_xsec_raw:.4f}^{{ \\pm {scaled_stat_raw:.4f} }}_{{ \\pm {scaled_sys_raw:.4f} }}$"

                    g_raw_cent_xsec.SetPoint(pt_idx_raw, actual_mass_center, scaled_xsec_raw)
                    g_raw_cent_xsec.SetPointError(pt_idx_raw, 0.0, scaled_stat_raw)
                    g_raw_cent_sys.SetPoint(pt_idx_raw, geometric_center, scaled_xsec_raw)
                    g_raw_cent_sys.SetPointError(pt_idx_raw, mass_width/2.0, scaled_sys_raw)

                    g_raw_geo_xsec.SetPoint(pt_idx_raw, geometric_center, scaled_xsec_raw)
                    g_raw_geo_xsec.SetPointError(pt_idx_raw, 0.0, scaled_stat_raw)
                    g_raw_geo_sys.SetPoint(pt_idx_raw, geometric_center, scaled_xsec_raw)
                    g_raw_geo_sys.SetPointError(pt_idx_raw, mass_width/2.0, scaled_sys_raw)

                    y_high_raw = scaled_xsec_raw + max(scaled_stat_raw, scaled_sys_raw)
                    y_low_raw = max(scaled_xsec_raw * 0.5, scaled_xsec_raw - max(scaled_stat_raw, scaled_sys_raw))
                    if y_high_raw > y_max_raw: y_max_raw = y_high_raw
                    if y_low_raw < y_min_raw: y_min_raw = y_low_raw
                    
                    pt_idx_raw += 1
                else:
                    sys_fraction = 0.15 

                # Unfolded Messy Extraction
                if Y_unf_m > 0 and acceptance > 0:
                    xsec_unf_m = (global_constant * Y_unf_m) / (mass_width * acceptance)
                    stat_unc_unf_m = (Y_unf_m_err / Y_unf_m) * xsec_unf_m
                    sys_unf_m = sys_fraction * xsec_unf_m
                    
                    scaled_xsec_unf_m = xsec_unf_m * (actual_mass_center**3)
                    scaled_stat_unf_m = stat_unc_unf_m * (actual_mass_center**3)
                    scaled_sys_unf_m = sys_unf_m * (actual_mass_center**3)
                    
                    str_unfm = f"${scaled_xsec_unf_m:.4f}^{{ \\pm {scaled_stat_unf_m:.4f} }}_{{ \\pm {scaled_sys_unf_m:.4f} }}$"

                    g_unf_m_cent_xsec.SetPoint(pt_idx_unfm, actual_mass_center, scaled_xsec_unf_m)
                    g_unf_m_cent_xsec.SetPointError(pt_idx_unfm, 0.0, scaled_stat_unf_m)
                    g_unf_m_cent_sys.SetPoint(pt_idx_unfm, geometric_center, scaled_xsec_unf_m)
                    g_unf_m_cent_sys.SetPointError(pt_idx_unfm, mass_width/2.0, scaled_sys_unf_m)

                    g_unf_m_geo_xsec.SetPoint(pt_idx_unfm, geometric_center, scaled_xsec_unf_m)
                    g_unf_m_geo_xsec.SetPointError(pt_idx_unfm, 0.0, scaled_stat_unf_m)
                    g_unf_m_geo_sys.SetPoint(pt_idx_unfm, geometric_center, scaled_xsec_unf_m)
                    g_unf_m_geo_sys.SetPointError(pt_idx_unfm, mass_width/2.0, scaled_sys_unf_m)

                    y_high_unf = scaled_xsec_unf_m + max(scaled_stat_unf_m, scaled_sys_unf_m)
                    y_low_unf = max(scaled_xsec_unf_m * 0.5, scaled_xsec_unf_m - max(scaled_stat_unf_m, scaled_sys_unf_m))
                    if y_high_unf > y_max_unf: y_max_unf = y_high_unf
                    if y_low_unf < y_min_unf: y_min_unf = y_low_unf
                    
                    pt_idx_unfm += 1

                # Unfolded Clean Extraction
                if Y_unf_c > 0 and acceptance > 0:
                    xsec_unf_c = (global_constant * Y_unf_c) / (mass_width * acceptance)
                    stat_unc_unf_c = (Y_unf_c_err / Y_unf_c) * xsec_unf_c
                    sys_unf_c = sys_fraction * xsec_unf_c
                    
                    scaled_xsec_unf_c = xsec_unf_c * (actual_mass_center**3)
                    scaled_stat_unf_c = stat_unc_unf_c * (actual_mass_center**3)
                    scaled_sys_unf_c = sys_unf_c * (actual_mass_center**3)
                    
                    str_unfc = f"${scaled_xsec_unf_c:.4f}^{{ \\pm {scaled_stat_unf_c:.4f} }}_{{ \\pm {scaled_sys_unf_c:.4f} }}$"

                    g_unf_c_cent_xsec.SetPoint(pt_idx_unfc, actual_mass_center, scaled_xsec_unf_c)
                    g_unf_c_cent_xsec.SetPointError(pt_idx_unfc, 0.0, scaled_stat_unf_c)
                    g_unf_c_cent_sys.SetPoint(pt_idx_unfc, geometric_center, scaled_xsec_unf_c)
                    g_unf_c_cent_sys.SetPointError(pt_idx_unfc, mass_width/2.0, scaled_sys_unf_c)

                    g_unf_c_geo_xsec.SetPoint(pt_idx_unfc, geometric_center, scaled_xsec_unf_c)
                    g_unf_c_geo_xsec.SetPointError(pt_idx_unfc, 0.0, scaled_stat_unf_c)
                    g_unf_c_geo_sys.SetPoint(pt_idx_unfc, geometric_center, scaled_xsec_unf_c)
                    g_unf_c_geo_sys.SetPointError(pt_idx_unfc, mass_width/2.0, scaled_sys_unf_c)

                    y_high_unf = scaled_xsec_unf_c + max(scaled_stat_unf_c, scaled_sys_unf_c)
                    y_low_unf = max(scaled_xsec_unf_c * 0.5, scaled_xsec_unf_c - max(scaled_stat_unf_c, scaled_sys_unf_c))
                    if y_high_unf > y_max_unf: y_max_unf = y_high_unf
                    if y_low_unf < y_min_unf: y_min_unf = y_low_unf
                    
                    pt_idx_unfc += 1

                if str_raw != "-" or str_unfm != "-" or str_unfc != "-":
                    s_xf_xsec = f"[{xf_min:.2f}, {xf_max:.2f})"
                    s_mass_xsec = f"[{mass_min:.2f}, {mass_max:.2f})"
                    row_xsec = f"{s_xf_xsec} & {xf_center:.2f} & {s_mass_xsec} & {geometric_center:.2f} & {str_raw} & {str_unfm} & {str_unfc} \\\\ \n\\hline\n"
                    latex_xsec_table_content += row_xsec
                
            if g_unf_m_cent_xsec.GetN() > 0 or g_raw_cent_xsec.GetN() > 0:
                overall_min = min(y_min_raw, y_min_unf) if min(y_min_raw, y_min_unf) > 0 else 1e-7
                overall_max = max(y_max_raw, y_max_unf)
                draw_and_save_canvas("Centroid", g_raw_cent_xsec, g_raw_cent_sys, g_unf_m_cent_xsec, g_unf_m_cent_sys, g_unf_c_cent_xsec, g_unf_c_cent_sys, xf_min, xf_max, i_x, theory_idx, overall_min, overall_max, f_ct18, f_nnpdf)
                draw_and_save_canvas("GeoCenter", g_raw_geo_xsec, g_raw_geo_sys, g_unf_m_geo_xsec, g_unf_m_geo_sys, g_unf_c_geo_xsec, g_unf_c_geo_sys, xf_min, xf_max, i_x, theory_idx, overall_min, overall_max, f_ct18, f_nnpdf)
            else:
                self.console.print(f"[dim]Note: Zero surviving events for {target_label} in xF [{xf_min:.2f}, {xf_max:.2f}]. No markers drawn.[/dim]")

        with open(f"Table_CrossSection_Final_Unfolded_{target_label}.tex", "w") as f:
            f.write(latex_xsec_table_content + r"\end{longtable}" + "\n" + r"\endgroup" + "\n")

        if acc_file: acc_file.Close()
        if f_ct18: f_ct18.Close()
        if f_nnpdf: f_nnpdf.Close()
        if f_psip: f_psip.Close()

    def generate_overlay_plot(self, target_label, plot_type, unfold_type):
        def scale_tgrapherrors(g, scale):
            if not g: return
            x_buf = g.GetX(); y_buf = g.GetY()
            for i in range(g.GetN()):
                g.SetPoint(i, x_buf[i], y_buf[i] * scale)
                g.SetPointError(i, g.GetErrorX(i), g.GetErrorY(i) * scale)

        ROOT.gStyle.SetTitleAlign(23); ROOT.gStyle.SetTitleX(0.5); ROOT.gStyle.SetTitleY(0.99)
        ROOT.gStyle.SetTitleH(0.04); ROOT.gStyle.SetTitleBorderSize(0)

        canvas = ROOT.TCanvas(f"canvas_overlay_{target_label}_{unfold_type}_{plot_type}", "Cross-Section Comparison", 1200, 1800)
        canvas.SetLogy(); canvas.SetLeftMargin(0.15); canvas.SetBottomMargin(0.12)
        canvas.SetTickx(1); canvas.SetTicky(1)

        legend = ROOT.TLegend(0.75, 0.45, 0.9, 0.9)
        legend.SetHeader(f"x_{{F}} Bins ({unfold_type})"); legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextFont(43); legend.SetTextSize(18)

        colors = [
            ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kMagenta, ROOT.kCyan,
            ROOT.kOrange + 7, ROOT.kSpring + 5, ROOT.kTeal + 5, ROOT.kAzure + 1,
            ROOT.kGray + 2, ROOT.kViolet - 5, ROOT.kYellow + 2, ROOT.kPink + 1,
            ROOT.kGreen - 9, ROOT.kRed - 9
        ]

        h_frame = canvas.DrawFrame(3.9, 1e-6, 10.0, 1e35)
        h_frame.SetTitle(f"DY Absolute Cross-Section Vs Mass for x_{{F}} bins ({target_label} {unfold_type} {plot_type})")
        h_frame.GetXaxis().SetTitle("Invariant Mass (GeV)"); h_frame.GetXaxis().CenterTitle()
        h_frame.GetXaxis().SetTitleOffset(1.2); h_frame.GetYaxis().SetTitle("M^{3} #frac{d^{2}#sigma}{dMdx_{F}} (nb GeV^{2})")
        h_frame.GetYaxis().CenterTitle(); h_frame.GetYaxis().SetTitleOffset(1.8)  
        
        latex_labels = []
        n_xf_bins = len(config.XF_BINS) - 1
        dir_xsec = self.out_file.Get(f"CrossSections_{target_label}")
        if not dir_xsec: return

        for i in range(n_xf_bins):
            g_xsec = dir_xsec.Get(f"g_xsec_{target_label}_{i}_{unfold_type}_{plot_type}")
            g_sys = dir_xsec.Get(f"g_sys_{target_label}_{i}_{unfold_type}_{plot_type}")
            if not g_xsec or not g_sys: continue
                
            g_xsec_clone = g_xsec.Clone(f"g_xsec_clone_{i}")
            g_sys_clone = g_sys.Clone(f"g_sys_clone_{i}")
            
            scale_factor = 1 * (10**(2*i))
            sf_txt = f"1#times10^{{{2*i}}}"
            scale_tgrapherrors(g_xsec_clone, scale_factor); scale_tgrapherrors(g_sys_clone, scale_factor)

            color = colors[i % len(colors)]
            g_sys_clone.SetLineColor(color); g_sys_clone.SetFillColorAlpha(color, 0.35)
            g_sys_clone.SetFillStyle(1001); g_sys_clone.SetMarkerSize(0)

            g_xsec_clone.SetLineColor(color); g_xsec_clone.SetMarkerColor(color)
            g_xsec_clone.SetMarkerStyle(ROOT.kFullCircle); g_xsec_clone.SetMarkerSize(1.0)

            g_sys_clone.Draw("2 SAME"); g_xsec_clone.Draw("P SAME")

            low_edge, high_edge = config.XF_BINS[i], config.XF_BINS[i+1]
            y_pos = 1.0 * scale_factor
            
            latex = ROOT.TLatex(4.0, y_pos, f"{low_edge:.2f}#leq x_{{F}} < {high_edge:.2f} ({sf_txt})")
            latex.SetTextFont(43); latex.SetTextSize(20); latex.SetTextColor(color)
            latex.Draw()
            latex_labels.append(latex) 
            legend.AddEntry(g_xsec_clone, f"x_{{F}} bin {i}", "pl")

        canvas.Update()
        out_pdf = f"cross_section_overlay_{target_label}_{unfold_type}_{plot_type}.pdf"
        canvas.SaveAs(out_pdf)
        
        dir_overlay = self.get_or_create_dir(self.out_file, "Overlays")
        dir_overlay.cd()
        canvas.Write(f"canvas_overlay_{target_label}_{unfold_type}_{plot_type}")

    def calculate_cross_sections(self):
        if self.hists_lh2 and self.hists_fl:
            self.sub_dict_lh2 = self.generate_subtracted_plot(self.hists_lh2, self.hists_fl, config.FLASK_NORM_LH2, "LH2")
            if self.sub_dict_lh2:
                self.calculate_and_plot_cross_section(self.sub_dict_lh2, "LH2", config.GLOBAL_CONSTANT_LH2)

        if self.hists_ld2 and self.hists_fl:
            self.sub_dict_ld2 = self.generate_subtracted_plot(self.hists_ld2, self.hists_fl, config.FLASK_NORM_LD2, "LD2")
            if self.sub_dict_ld2:
                self.calculate_and_plot_cross_section(self.sub_dict_ld2, "LD2", config.GLOBAL_CONSTANT_LD2)

        for tgt in ["LH2", "LD2"]:
            for p_type in ["Centroid", "GeoCenter"]:
                for u_type in ["Raw", "UnfMessy", "UnfClean"]:
                    self.generate_overlay_plot(tgt, p_type, u_type)

    def generate_latex_appendix(self):
        latex_filename = "Appendix_MassCentroids.tex"
        with open(latex_filename, "w") as tex_file:
            intro_text = r"""\section{Appendix: Determination of Mass Bin Centroids}...""" 
            tex_file.write(intro_text)

    def finalize(self):
        self.out_file.Write()
        self.out_file.Close()
        if self.rm_file:
            self.rm_file.Write()
            self.rm_file.Close()