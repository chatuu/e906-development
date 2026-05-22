"""
analyzer.py
Core Object-Oriented Analysis Module for Drell-Yan Cross-Sections vs pT (Binned by xF).
"""

import os
import sys
import math
import uproot
import numpy as np
import ROOT
import config

class DYCrossSectionAnalyzer:
    def __init__(self, lh2_files, ld2_files, flask_files, out_filename="All_XSec_Objects_Binned_xF.root"):
        self._setup_root()
        
        self.lh2_paths = lh2_files if isinstance(lh2_files, list) else [lh2_files]
        self.ld2_paths = ld2_files if isinstance(ld2_files, list) else [ld2_files]
        self.flask_paths = flask_files if isinstance(flask_files, list) else [flask_files]
        
        self.out_filename = out_filename
        self.out_file = ROOT.TFile(out_filename, "RECREATE")
        
        try:
            npz_data = np.load(config.INPUT_NPZ_FILE)
            self.x_curve = npz_data['x']
        except Exception as e:
            print(f"Error loading NPZ file: {e}")
            sys.exit(1)

        self.n_xf_bins = len(config.XF_BINS) - 1
        
        # Arrays to hold dictionaries of histograms per xF bin
        self.hists_lh2 = [None] * self.n_xf_bins
        self.hists_ld2 = [None] * self.n_xf_bins
        self.hists_fl = [None] * self.n_xf_bins
        self.sub_dict_lh2 = [None] * self.n_xf_bins
        self.sub_dict_pd = [None] * self.n_xf_bins

    def _setup_root(self):
        ROOT.gROOT.SetBatch(True)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetPalette(ROOT.kBird)
        ROOT.gStyle.SetEndErrorSize(5) 
        ROOT.gErrorIgnoreLevel = ROOT.kWarning

    @staticmethod
    def get_or_create_dir(base_dir, name):
        d = base_dir.GetDirectory(name)
        if not d: d = base_dir.mkdir(name)
        return d

    @staticmethod
    def apply_cuts(tree, cut=4.2):
        events = tree.arrays(library="np")
        class EventNamespace:
            def __init__(self, data): self.__dict__.update(data)
        e = EventNamespace(events)

        bo = np.where(e.runID >= 11000, 1.6, 0.4)
        dimuon_cut = (
            (np.abs(e.dx) < 0.25) & (np.abs(e.dy - bo) < 0.22) & (e.dz < -5.) & (e.dz > -280.) & 
            (np.abs(e.dpx) < 1.8) & (np.abs(e.dpy) < 2.0) & (e.dpx**2 + e.dpy**2 < 5.) & 
            (e.dpz < 116.) & (e.dpz > 38.) & (e.mass > cut) & (e.mass < 8.8) &
            (e.dx**2 + (e.dy - bo)**2 < 0.06) & (e.xF < 0.95) & (e.xF > -0.1) & 
            (e.xT > 0.05) & (e.xT <= 0.58) & (np.abs(e.costh) < 0.5) & 
            (np.abs(e.trackSeparation) < 270.) & (e.chisq_dimuon < 18)
        )
        track1_cut = ((e.chisq1_target < 15.) & (e.pz1_st1 > 9.) & (e.pz1_st1 < 75.) & (e.nHits1 > 13) & (e.y1_st1 * e.y1_st3 > 0.) & (np.abs(e.py1_st1) > 0.02))
        track2_cut = ((e.chisq2_target < 15.) & (e.pz2_st1 > 9.) & (e.pz2_st1 < 75.) & (e.nHits2 > 13) & (e.y2_st1 * e.y2_st3 > 0.) & (np.abs(e.py2_st1) > 0.02))
        tracks_cut = ((np.abs(e.chisq1_target + e.chisq2_target - e.chisq_dimuon) < 2.) & ((e.y1_st3) * (e.y2_st3) < 0.) & (e.nHits1 + e.nHits2 > 29))
        occ_cut = ((e.D1 < 400) & (e.D2 < 400) & (e.D3 < 400) & (e.D1 + e.D2 + e.D3 < 1000))

        total_cut_mask = (track1_cut & track2_cut & tracks_cut & dimuon_cut & occ_cut)
        filtered_events = {key: val[total_cut_mask] for key, val in events.items()}
        filtered_events["pT"] = np.sqrt(filtered_events["dpx"]**2 + filtered_events["dpy"]**2)
        return filtered_events

    def get_concatenated_events(self, file_paths, tree_name):
        all_events = None
        for fp in file_paths:
            if not os.path.exists(fp): continue
            try:
                with uproot.open(fp) as f:
                    if tree_name not in f: continue
                    filtered = self.apply_cuts(f[tree_name])
                    if all_events is None:
                        all_events = {k: [v] for k, v in filtered.items()}
                    else:
                        for k, v in filtered.items():
                            if k in all_events: all_events[k].append(v)
            except Exception as e: print(f"Error reading {fp}: {e}")
                
        if all_events is None: raise RuntimeError(f"No valid data found for tree '{tree_name}'.")
        return {k: np.concatenate(v) for k, v in all_events.items()}

    @staticmethod
    def create_histograms(file_label, i_xf):
        hists = {}
        def make_th1(name, title, y_title=""):
            h = ROOT.TH1D(f"{name}_{file_label}_xF{i_xf}", f"{title} ({file_label}, xF Bin {i_xf});p_{{T}} [GeV];{y_title}", 
                          len(config.PT_BINS)-1, config.PT_BINS)
            h.Sumw2(); h.SetStats(0); h.GetXaxis().CenterTitle(); h.GetYaxis().CenterTitle()
            return h

        hists["Y_total"] = make_th1("Y_total", "Total Yield", "Yield")
        hists["Y_mix"] = make_th1("Y_mix", "Mix Yield", "Yield")
        hists["E_total_reco"] = make_th1("E_total_reco", "Avg Reco Eff", "Efficiency")
        hists["E_total_hodo"] = make_th1("E_total_hodo", "Avg Hodo Eff", "Efficiency")
        hists["E_total_final"] = make_th1("E_total_final", "Avg Final Eff", "Efficiency")
        hists["E_final_signal"] = make_th1("E_final_signal", "Avg Signal Efficiency", "Efficiency")
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
        correl_matrix[diff_matrix == 0] = 1.0; correl_matrix[diff_matrix == 1] = 1.0  
        covar_matrix = correl_matrix * np.outer(err_arr, err_arr)
        err_on_mean = np.sqrt(np.sum(covar_matrix)) / N
        return mean_eff, err_on_mean

    def extract_target_stats(self, data_tot, data_mix, loc_tot, loc_mix, m_low, m_high, pt_low, pt_high, xf_low, xf_high):
        mask_tot = (data_tot["mass"] >= m_low) & (data_tot["mass"] < m_high) & (data_tot["pT"] >= pt_low) & (data_tot["pT"] < pt_high) & (data_tot["xF"] >= xf_low) & (data_tot["xF"] < xf_high)
        mask_mix = (data_mix["mass"] >= m_low) & (data_mix["mass"] < m_high) & (data_mix["pT"] >= pt_low) & (data_mix["pT"] < pt_high) & (data_mix["xF"] >= xf_low) & (data_mix["xF"] < xf_high)

        N_tot = np.sum(mask_tot)
        N_mix = np.sum(mask_mix)
        
        sum_mass_tot = np.sum(data_tot["mass"][mask_tot]) if N_tot > 0 else 0.0
        sum_mass_mix = np.sum(data_mix["mass"][mask_mix]) if N_mix > 0 else 0.0
        sum_pt_tot = np.sum(data_tot["pT"][mask_tot]) if N_tot > 0 else 0.0
        sum_pt_mix = np.sum(data_mix["pT"][mask_mix]) if N_mix > 0 else 0.0

        means = {
            "r_tot": self.get_correlated_mean_and_error(data_tot["recoeff"][mask_tot], data_tot["recoeff_error"][mask_tot], loc_tot[mask_tot]),
            "h_tot": self.get_weighted_mean_and_error(data_tot["hodoeff"][mask_tot], data_tot["hodoeff_error"][mask_tot]),
            "r_mix": self.get_correlated_mean_and_error(data_mix["recoeff"][mask_mix], data_mix["recoeff_error"][mask_mix], loc_mix[mask_mix]),
            "h_mix": self.get_weighted_mean_and_error(data_mix["hodoeff"][mask_mix], data_mix["hodoeff_error"][mask_mix]),
        }

        means["f_tot"] = (means["r_tot"][0] * means["h_tot"][0], np.sqrt((means["h_tot"][0] * means["r_tot"][1])**2 + (means["r_tot"][0] * means["h_tot"][1])**2))
        means["f_mix"] = (means["r_mix"][0] * means["h_mix"][0], np.sqrt((means["h_mix"][0] * means["r_mix"][1])**2 + (means["r_mix"][0] * means["h_mix"][1])**2))

        val_sig_eff, err_sig_eff = 0.0, 0.0
        diff_yield = N_tot - N_mix
        if diff_yield != 0:
            val_sig_eff = ((N_tot * means["f_tot"][0]) - (N_mix * means["f_mix"][0])) / diff_yield
            err_sig_eff = (1.0 / diff_yield) * np.sqrt((N_tot * means["f_tot"][1])**2 + (N_mix * means["f_mix"][1])**2)

        return N_tot, N_mix, sum_mass_tot, sum_mass_mix, sum_pt_tot, sum_pt_mix, val_sig_eff, err_sig_eff, diff_yield, means

    def fill_histograms(self, hists, root_x, N_tot, N_mix, sig_eff, err_sig_eff, diff_yield, means, final_centroid, final_pt_centroid):
        val_corr_yield, err_corr_stat, err_corr_sys = 0.0, 0.0, 0.0
        if sig_eff > 0 and diff_yield > 0:
            val_corr_yield = diff_yield / sig_eff
            err_corr_stat = np.sqrt(N_tot + N_mix) / sig_eff
            err_corr_sys = val_corr_yield * (err_sig_eff / sig_eff)

        hists["Y_total"].SetBinContent(root_x, N_tot); hists["Y_total"].SetBinError(root_x, np.sqrt(N_tot))
        hists["Y_mix"].SetBinContent(root_x, N_mix); hists["Y_mix"].SetBinError(root_x, np.sqrt(N_mix))
        hists["E_total_reco"].SetBinContent(root_x, means["r_tot"][0]); hists["E_total_reco"].SetBinError(root_x, means["r_tot"][1])
        hists["E_total_hodo"].SetBinContent(root_x, means["h_tot"][0]); hists["E_total_hodo"].SetBinError(root_x, means["h_tot"][1])
        hists["E_total_final"].SetBinContent(root_x, means["f_tot"][0]); hists["E_total_final"].SetBinError(root_x, means["f_tot"][1])
        hists["E_final_signal"].SetBinContent(root_x, sig_eff); hists["E_final_signal"].SetBinError(root_x, err_sig_eff)
        hists["Y_corrected_stat"].SetBinContent(root_x, val_corr_yield); hists["Y_corrected_stat"].SetBinError(root_x, err_corr_stat)
        hists["Y_corrected_sys"].SetBinContent(root_x, val_corr_yield); hists["Y_corrected_sys"].SetBinError(root_x, err_corr_sys)
        hists["Mass_Centroid"].SetBinContent(root_x, final_centroid)
        hists["Pt_Centroid"].SetBinContent(root_x, final_pt_centroid)

    def process_kinematics(self):
        def get_loc_data(data_dict, is_mix=False):
            d1_vals = 0.5 * (data_dict['ptrk_D1'] + data_dict['ntrk_D1']) if (is_mix and 'ptrk_D1' in data_dict) else data_dict['D1']
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

        m_low, m_high = config.MASS_BINS[0], config.MASS_BINS[-1]
        m_center_overall = (m_low + m_high) / 2.0

        for i_xf in range(self.n_xf_bins):
            xf_low, xf_high = config.XF_BINS[i_xf], config.XF_BINS[i_xf+1]
            
            dir_xf = self.get_or_create_dir(self.out_file, f"xF_Bin_{i_xf}")
            dir_lh2 = self.get_or_create_dir(dir_xf, "LH2"); dir_lh2.cd(); self.hists_lh2[i_xf] = self.create_histograms("LH2", i_xf)
            dir_ld2 = self.get_or_create_dir(dir_xf, "LD2"); dir_ld2.cd(); self.hists_ld2[i_xf] = self.create_histograms("LD2", i_xf)
            dir_fl = self.get_or_create_dir(dir_xf, "Flask"); dir_fl.cd(); self.hists_fl[i_xf] = self.create_histograms("Flask", i_xf)

            for i_pt in range(len(config.PT_BINS) - 1):
                pt_low, pt_high = config.PT_BINS[i_pt], config.PT_BINS[i_pt+1]
                pt_center = (pt_low + pt_high) / 2.0
                root_x = i_pt + 1
                
                N_l2_t, N_l2_m, M_l2_t, M_l2_m, P_l2_t, P_l2_m, eps_ld2, err_eps_ld2, diff_l2, mean_l2 = self.extract_target_stats(data_ld2_tot, data_ld2_mix, loc_ld2_tot, loc_ld2_mix, m_low, m_high, pt_low, pt_high, xf_low, xf_high)
                N_lh_t, N_lh_m, M_lh_t, M_lh_m, P_lh_t, P_lh_m, eps_lh2, err_eps_lh2, diff_lh, mean_lh = self.extract_target_stats(data_lh2_tot, data_lh2_mix, loc_lh2_tot, loc_lh2_mix, m_low, m_high, pt_low, pt_high, xf_low, xf_high)
                N_fl_t, N_fl_m, M_fl_t, M_fl_m, P_fl_t, P_fl_m, eps_fl, err_eps_fl, diff_fl, mean_fl = self.extract_target_stats(data_fl_tot, data_fl_mix, loc_fl_tot, loc_fl_mix, m_low, m_high, pt_low, pt_high, xf_low, xf_high)

                if eps_lh2 > 0:
                    cY_lh_t, cM_lh_t, cP_lh_t = N_lh_t / eps_lh2, M_lh_t / eps_lh2, P_lh_t / eps_lh2
                    cY_lh_m, cM_lh_m, cP_lh_m = N_lh_m / eps_lh2, M_lh_m / eps_lh2, P_lh_m / eps_lh2
                else: cY_lh_t = cM_lh_t = cP_lh_t = cY_lh_m = cM_lh_m = cP_lh_m = 0.0

                if eps_ld2 > 0:
                    cY_l2_t, cM_l2_t, cP_l2_t = N_l2_t / eps_ld2, M_l2_t / eps_ld2, P_l2_t / eps_ld2
                    cY_l2_m, cM_l2_m, cP_l2_m = N_l2_m / eps_ld2, M_l2_m / eps_ld2, P_l2_m / eps_ld2
                else: cY_l2_t = cM_l2_t = cP_l2_t = cY_l2_m = cM_l2_m = cP_l2_m = 0.0

                if eps_fl > 0:
                    cY_fl_t_lh, cM_fl_t_lh, cP_fl_t_lh = (config.FLASK_NORM_LH2 * N_fl_t) / eps_fl, (config.FLASK_NORM_LH2 * M_fl_t) / eps_fl, (config.FLASK_NORM_LH2 * P_fl_t) / eps_fl
                    cY_fl_m_lh, cM_fl_m_lh, cP_fl_m_lh = (config.FLASK_NORM_LH2 * N_fl_m) / eps_fl, (config.FLASK_NORM_LH2 * M_fl_m) / eps_fl, (config.FLASK_NORM_LH2 * P_fl_m) / eps_fl
                    cY_fl_t_l2, cM_fl_t_l2, cP_fl_t_l2 = (config.FLASK_NORM_LD2 * N_fl_t) / eps_fl, (config.FLASK_NORM_LD2 * M_fl_t) / eps_fl, (config.FLASK_NORM_LD2 * P_fl_t) / eps_fl
                    cY_fl_m_l2, cM_fl_m_l2, cP_fl_m_l2 = (config.FLASK_NORM_LD2 * N_fl_m) / eps_fl, (config.FLASK_NORM_LD2 * M_fl_m) / eps_fl, (config.FLASK_NORM_LD2 * P_fl_m) / eps_fl
                else:
                    cY_fl_t_lh = cM_fl_t_lh = cP_fl_t_lh = cY_fl_m_lh = cM_fl_m_lh = cP_fl_m_lh = 0.0
                    cY_fl_t_l2 = cM_fl_t_l2 = cP_fl_t_l2 = cY_fl_m_l2 = cM_fl_m_l2 = cP_fl_m_l2 = 0.0

                num_LH2 = (cM_lh_t - cM_lh_m) - (cM_fl_t_lh - cM_fl_m_lh)
                num_pt_LH2 = (cP_lh_t - cP_lh_m) - (cP_fl_t_lh - cP_fl_m_lh)
                den_LH2 = (cY_lh_t - cY_lh_m) - (cY_fl_t_lh - cY_fl_m_lh)
                cent_LH2 = num_LH2 / den_LH2 if den_LH2 != 0 else m_center_overall
                cent_pt_LH2 = num_pt_LH2 / den_LH2 if den_LH2 != 0 else pt_center

                num_LD2 = (cM_l2_t - cM_l2_m) - (cM_fl_t_l2 - cM_fl_m_l2) - (config.LH2_TO_LD2_NORM * num_LH2)
                num_pt_LD2 = (cP_l2_t - cP_l2_m) - (cP_fl_t_l2 - cP_fl_m_l2) - (config.LH2_TO_LD2_NORM * num_pt_LH2)
                den_LD2 = (cY_l2_t - cY_l2_m) - (cY_fl_t_l2 - cY_fl_m_l2) - (config.LH2_TO_LD2_NORM * den_LH2)
                cent_LD2 = num_LD2 / den_LD2 if den_LD2 != 0 else m_center_overall
                cent_pt_LD2 = num_pt_LD2 / den_LD2 if den_LD2 != 0 else pt_center

                self.fill_histograms(self.hists_ld2[i_xf], root_x, N_l2_t, N_l2_m, eps_ld2, err_eps_ld2, diff_l2, mean_l2, cent_LD2, cent_pt_LD2)
                self.fill_histograms(self.hists_lh2[i_xf], root_x, N_lh_t, N_lh_m, eps_lh2, err_eps_lh2, diff_lh, mean_lh, cent_LH2, cent_pt_LH2)
                
                num_fl = (M_fl_t - M_fl_m) / (eps_fl if eps_fl > 0 else 1e-9)
                num_pt_fl = (P_fl_t - P_fl_m) / (eps_fl if eps_fl > 0 else 1e-9)
                den_fl = (N_fl_t - N_fl_m) / (eps_fl if eps_fl > 0 else 1e-9)
                self.fill_histograms(self.hists_fl[i_xf], root_x, N_fl_t, N_fl_m, eps_fl, err_eps_fl, diff_fl, mean_fl, num_fl/den_fl if den_fl!=0 else m_center_overall, num_pt_fl/den_fl if den_fl!=0 else pt_center)

    def generate_subtracted_plot(self, hists_target, hists_flask, flask_norm, target_label, i_xf):
        dir_sub = self.get_or_create_dir(self.out_file, f"Subtracted_Plots_{target_label}")
        dir_sub.cd()
        
        h_target_stat, h_target_sys = hists_target["Y_corrected_stat"], hists_target["Y_corrected_sys"]
        h_flask_stat, h_flask_sys = hists_flask["Y_corrected_stat"], hists_flask["Y_corrected_sys"]
        
        name = f"Y_corrected_Subtracted_{target_label}_xF{i_xf}"
        h_sub = ROOT.TH1D(name, f"Corrected Yield {target_label} (xF Bin {i_xf});p_{{T}} [GeV];Yield", len(config.PT_BINS)-1, config.PT_BINS)
        h_sub_stat = h_sub.Clone(f"{name}_stat"); h_sub_sys = h_sub.Clone(f"{name}_sys")
        
        for i_pt in range(len(config.PT_BINS) - 1):
            bin_x = i_pt + 1
            y_target, e_target_stat, e_target_sys = h_target_stat.GetBinContent(bin_x), h_target_stat.GetBinError(bin_x), h_target_sys.GetBinError(bin_x)
            y_flask, e_flask_stat, e_flask_sys = h_flask_stat.GetBinContent(bin_x), h_flask_stat.GetBinError(bin_x), h_flask_sys.GetBinError(bin_x)
            
            val_sub = y_target - (flask_norm * y_flask)
            err_sub_stat = np.sqrt(e_target_stat**2 + (flask_norm * e_flask_stat)**2)
            err_sub_sys = np.sqrt(e_target_sys**2 + (flask_norm * e_flask_sys)**2)
            
            h_sub_stat.SetBinContent(bin_x, val_sub); h_sub_stat.SetBinError(bin_x, err_sub_stat)
            h_sub_sys.SetBinContent(bin_x, val_sub); h_sub_sys.SetBinError(bin_x, err_sub_sys)
            h_sub.SetBinContent(bin_x, val_sub); h_sub.SetBinError(bin_x, np.sqrt(err_sub_stat**2 + err_sub_sys**2))

        return {"stat": h_sub_stat, "sys": h_sub_sys, "centroid": hists_target["Mass_Centroid"], "pt_centroid": hists_target["Pt_Centroid"]}

    def generate_pd_subtracted_plot(self, hists_ld2, hists_lh2, hists_flask, i_xf):
        dir_sub = self.get_or_create_dir(self.out_file, "Subtracted_Plots_LD2_pd")
        dir_sub.cd()

        h_ld2_stat, h_ld2_sys = hists_ld2["Y_corrected_stat"], hists_ld2["Y_corrected_sys"]
        h_lh2_stat, h_lh2_sys = hists_lh2["Y_corrected_stat"], hists_lh2["Y_corrected_sys"]
        h_flask_stat, h_flask_sys = hists_flask["Y_corrected_stat"], hists_flask["Y_corrected_sys"]

        c_lh2 = config.THD_THH_RATIO * (config.PROTONS_ON_TARGET_LD2 / config.PROTONS_ON_TARGET_LH2)
        c_flask_sub = config.FLASK_NORM_LD2 - (c_lh2 * config.FLASK_NORM_LH2)

        name = f"Y_corrected_Subtracted_LD2_pd_xF{i_xf}"
        h_pd = ROOT.TH1D(name, f"Corrected Yield LD2 (xF Bin {i_xf});p_{{T}} [GeV];Yield", len(config.PT_BINS)-1, config.PT_BINS)
        h_pd_stat = h_pd.Clone(f"{name}_stat"); h_pd_sys = h_pd.Clone(f"{name}_sys")

        for i_pt in range(len(config.PT_BINS) - 1):
            bin_x = i_pt + 1
            y_ld2, y_lh2, y_flask = h_ld2_stat.GetBinContent(bin_x), h_lh2_stat.GetBinContent(bin_x), h_flask_stat.GetBinContent(bin_x)
            e_ld2_stat, e_lh2_stat, e_flask_stat = h_ld2_stat.GetBinError(bin_x), h_lh2_stat.GetBinError(bin_x), h_flask_stat.GetBinError(bin_x)
            e_ld2_sys, e_lh2_sys, e_flask_sys = h_ld2_sys.GetBinError(bin_x), h_lh2_sys.GetBinError(bin_x), h_flask_sys.GetBinError(bin_x)
            
            val_pd = y_ld2 - (c_lh2 * y_lh2) - (c_flask_sub * y_flask)
            err_pd_stat = np.sqrt(e_ld2_stat**2 + (c_lh2 * e_lh2_stat)**2 + (c_flask_sub * e_flask_stat)**2)
            err_pd_sys = np.sqrt(e_ld2_sys**2 + (c_lh2 * e_lh2_sys)**2 + (c_flask_sub * e_flask_sys)**2)
            
            h_pd_stat.SetBinContent(bin_x, val_pd); h_pd_stat.SetBinError(bin_x, err_pd_stat)
            h_pd_sys.SetBinContent(bin_x, val_pd); h_pd_sys.SetBinError(bin_x, err_pd_sys)
            h_pd.SetBinContent(bin_x, val_pd); h_pd.SetBinError(bin_x, np.sqrt(err_pd_stat**2 + err_pd_sys**2))

        return {"stat": h_pd_stat, "sys": h_pd_sys, "centroid": hists_ld2["Mass_Centroid"], "pt_centroid": hists_ld2["Pt_Centroid"]}

    def calculate_and_plot_cross_section(self, h_sub_dict, target_label, global_constant, i_xf, use_true_pt=False):
        acc_path = "/root/github/e906-development/src/AcceptanceCorrection/acceptance_mass_xF.root"
        psip_path = "/root/github/e906-development/src/psip_contamination/PsiP_Contamination_vs_pT.root" 
        
        acc_file, f_psip = None, None
        if os.path.exists(acc_path): acc_file = ROOT.TFile.Open(acc_path)
        if os.path.exists(psip_path): f_psip = ROOT.TFile.Open(psip_path)

        # =================================================================================================
        # CRITICAL FALLBACK LOGIC: 
        # If the exact directory/histogram names inside your ROOT files don't match the strings below,
        # the code will NO LONGER crash or abort. It will default acceptance to 1.0, print a warning, 
        # and continue generating all plots.
        # =================================================================================================
        h_acc_1d = None
        if acc_file:
            # Try to grab the acceptance histogram for this specific xF bin. 
            # You may need to change "xF_bin_{i_xf}" to match your ROOT file's directory structure!
            h_acc_1d = acc_file.Get(f"xF_bin_{i_xf}/h_ratio_{target_label}_pT")
            if not h_acc_1d and target_label == "LD2": h_acc_1d = acc_file.Get(f"xF_bin_{i_xf}/h_ratio_LH2_pT")
            
        if not h_acc_1d: 
            print(f"[Warning] Acceptance histogram for {target_label} xF Bin {i_xf} not found! Defaulting acceptance to 1.0 to ensure plots generate.")

        h_ratio_psip_1d = None
        if f_psip: h_ratio_psip_1d = f_psip.Get(f"xF_bin_{i_xf}/h_Ratio_PsiP_DY_pT")

        dir_xsec = self.get_or_create_dir(self.out_file, f"CrossSections_{target_label}")
        dir_xsec.cd()
        
        h_sub_stat, h_sub_sys, h_sub_centroid, h_sub_pt_centroid = h_sub_dict["stat"], h_sub_dict["sys"], h_sub_dict["centroid"], h_sub_dict["pt_centroid"]
        n_pt_bins = len(config.PT_BINS) - 1
        suffix = f"_true_pt_xF{i_xf}" if use_true_pt else f"_geom_xF{i_xf}"

        # ----------------------------------------------------------------------
        # PLOT 1 & 2: Individual Differential Cross-Section Plots
        # ----------------------------------------------------------------------
        g_xsec = ROOT.TGraphErrors(); g_xsec.SetName(f"g_xsec_{target_label}{suffix}")
        g_sys = ROOT.TGraphErrors(); g_sys.SetName(f"g_sys_{target_label}{suffix}")
        h1_xsec = ROOT.TH1D(f"h1_xsec_{target_label}{suffix}", f"Cross Section {target_label} xF{i_xf};p_{{T}} [GeV];d#sigma/dp_{{T}}", n_pt_bins, config.PT_BINS)
        h1_sys = ROOT.TH1D(f"h1_sys_{target_label}{suffix}", f"Systematic {target_label} xF{i_xf};p_{{T}} [GeV];d#sigma/dp_{{T}}", n_pt_bins, config.PT_BINS)

        point_idx = 0
        y_max_data = -1e9
        
        for i_pt in range(n_pt_bins):
            root_bin_x = i_pt + 1
            pt_min, pt_max = config.PT_BINS[i_pt], config.PT_BINS[i_pt+1]
            pt_width = pt_max - pt_min
            pt_center = (pt_min + pt_max) / 2.0
            
            # Default to 1.0 if histogram wasn't found
            acceptance = 1.0; acceptance_err = 0.0
            if h_acc_1d:
                acceptance = h_acc_1d.GetBinContent(root_bin_x)
                acceptance_err = h_acc_1d.GetBinError(root_bin_x)

            if acceptance <= 0: continue
            Y_sub = h_sub_stat.GetBinContent(root_bin_x)
            if Y_sub <= 0: continue
            
            actual_pt_center = h_sub_pt_centroid.GetBinContent(root_bin_x)
            if actual_pt_center <= 0: actual_pt_center = pt_center

            psip_ratio = h_ratio_psip_1d.GetBinContent(root_bin_x) if h_ratio_psip_1d else 0.0
            if psip_ratio < 0.0 or psip_ratio > 1.0: psip_ratio = 0.0

            numerator = global_constant * Y_sub
            denominator = pt_width * acceptance
            bin_dsigma_dpt = numerator / denominator
            
            stat_unc = (h_sub_stat.GetBinError(root_bin_x) / Y_sub) * bin_dsigma_dpt
            sys_tot_yield = (h_sub_sys.GetBinError(root_bin_x) / Y_sub) * bin_dsigma_dpt
            sys_psip_cont = psip_ratio * bin_dsigma_dpt
            
            if bin_dsigma_dpt > 0:
                sys_acc_total = (acceptance_err / acceptance) * bin_dsigma_dpt
                total_stat_err = stat_unc
                total_sys_err = math.sqrt(sys_tot_yield**2 + sys_psip_cont**2 + sys_acc_total**2)
                
                plot_x = actual_pt_center if use_true_pt else pt_center
                g_sys.SetPoint(point_idx, pt_center, bin_dsigma_dpt)
                g_sys.SetPointError(point_idx, pt_width/2.0, total_sys_err)
                g_xsec.SetPoint(point_idx, plot_x, bin_dsigma_dpt)
                g_xsec.SetPointError(point_idx, 0.0, total_stat_err)
                
                h1_xsec.SetBinContent(root_bin_x, bin_dsigma_dpt); h1_xsec.SetBinError(root_bin_x, total_stat_err)
                h1_sys.SetBinContent(root_bin_x, bin_dsigma_dpt); h1_sys.SetBinError(root_bin_x, total_sys_err)
                
                if bin_dsigma_dpt + total_sys_err > y_max_data: y_max_data = bin_dsigma_dpt + total_sys_err
                point_idx += 1
                
        # Draw and Save the Canvas
        if g_xsec.GetN() > 0:
            c_xsec = ROOT.TCanvas(f"c_xsec_{target_label}{suffix}", f"Cross-Section {target_label} xF Bin {i_xf}", 800, 600)
            c_xsec.SetLeftMargin(0.16); c_xsec.SetBottomMargin(0.14); c_xsec.SetTickx(1); c_xsec.SetTicky(1)
            
            mg = ROOT.TMultiGraph()
            target_pfx = "pp" if target_label == "LH2" else "pd"
            mg.SetTitle(f"d#sigma_{{{target_pfx}}}/dp_{{T}} (xF Bin {i_xf});p_{{T}} [GeV];d#sigma/dp_{{T}} [nb/GeV/Nucleon]")
            mg.SetMinimum(0.0); mg.SetMaximum(y_max_data * 1.4 if y_max_data > 0 else 3.0)
            
            g_sys.SetMarkerSize(0); g_sys.SetLineColor(ROOT.kRed); g_sys.SetFillColorAlpha(ROOT.kPink - 9, 0.5); g_sys.SetFillStyle(1001)
            mg.Add(g_sys, "2")
            g_xsec.SetMarkerStyle(20); g_xsec.SetMarkerColor(ROOT.kRed); g_xsec.SetLineColor(ROOT.kRed)
            mg.Add(g_xsec, "P")
            
            mg.Draw("A"); mg.GetXaxis().SetLimits(0.0, 2.0)
            
            leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88); leg.SetBorderSize(0)
            leg.AddEntry(g_xsec, f"Data ({target_label})", "lep"); leg.AddEntry(g_sys, "Systematic Unc.", "f")
            leg.Draw()
            
            c_xsec.SaveAs(f"CrossSection_{target_label}_vs_pT{suffix}.pdf")
            c_xsec.Write(); c_xsec.Close()
                
        dir_xsec.cd()
        g_xsec.Write(); g_sys.Write(); h1_xsec.Write(); h1_sys.Write()
        if f_psip: f_psip.Close()
        if acc_file: acc_file.Close()

    def generate_ratio_plot(self, i_xf, use_true_pt=False):
        suffix = f"_true_pt_xF{i_xf}" if use_true_pt else f"_geom_xF{i_xf}"
        dir_lh2 = self.out_file.Get("CrossSections_LH2")
        dir_ld2 = self.out_file.Get("CrossSections_LD2")

        if not dir_lh2 or not dir_ld2: return

        h1_xsec_lh2, h1_sys_lh2 = dir_lh2.Get(f"h1_xsec_LH2{suffix}"), dir_lh2.Get(f"h1_sys_LH2{suffix}")
        h1_xsec_ld2, h1_sys_ld2 = dir_ld2.Get(f"h1_xsec_LD2{suffix}"), dir_ld2.Get(f"h1_sys_LD2{suffix}")

        if not all([h1_xsec_lh2, h1_sys_lh2, h1_xsec_ld2, h1_sys_ld2]): return

        h_pt_lh2 = self.sub_dict_lh2[i_xf]["pt_centroid"] if self.sub_dict_lh2[i_xf] else None
        h_pt_ld2 = self.sub_dict_pd[i_xf]["pt_centroid"] if self.sub_dict_pd[i_xf] else None

        dir_ratio = self.get_or_create_dir(self.out_file, "Ratio_pd_2pp")
        dir_ratio.cd()

        # ----------------------------------------------------------------------
        # PLOT 3: Individual Ratio Plots vs pT
        # ----------------------------------------------------------------------
        g_ratio_stat = ROOT.TGraphErrors(); g_ratio_stat.SetName(f"g_ratio_stat{suffix}")
        g_ratio_sys = ROOT.TGraphErrors(); g_ratio_sys.SetName(f"g_ratio_sys{suffix}")

        pt_idx = 0
        y_max_ratio, y_min_ratio = -1e9, 1e9

        for i_pt in range(len(config.PT_BINS) - 1):
            bin_x = i_pt + 1
            y_lh2, y_ld2 = h1_xsec_lh2.GetBinContent(bin_x), h1_xsec_ld2.GetBinContent(bin_x)
            if y_lh2 <= 0 or y_ld2 <= 0: continue

            ratio = y_ld2 / (2.0 * y_lh2)
            rel_stat_ld2, rel_stat_lh2 = h1_xsec_ld2.GetBinError(bin_x) / y_ld2, h1_xsec_lh2.GetBinError(bin_x) / y_lh2
            err_ratio_stat = ratio * math.sqrt(rel_stat_ld2**2 + rel_stat_lh2**2)

            rel_sys_ld2, rel_sys_lh2 = h1_sys_ld2.GetBinError(bin_x) / y_ld2, h1_sys_lh2.GetBinError(bin_x) / y_lh2
            err_ratio_sys = ratio * math.sqrt(rel_sys_ld2**2 + rel_sys_lh2**2)

            bin_center = (config.PT_BINS[i_pt] + config.PT_BINS[i_pt+1]) / 2.0
            plot_x = (h_pt_lh2.GetBinContent(bin_x) + h_pt_ld2.GetBinContent(bin_x)) / 2.0 if (use_true_pt and h_pt_lh2 and h_pt_ld2) else bin_center
            pt_width = config.PT_BINS[i_pt+1] - config.PT_BINS[i_pt]

            g_ratio_stat.SetPoint(pt_idx, plot_x, ratio); g_ratio_stat.SetPointError(pt_idx, 0.0, err_ratio_stat)
            g_ratio_sys.SetPoint(pt_idx, bin_center, ratio); g_ratio_sys.SetPointError(pt_idx, pt_width/2.0, err_ratio_sys)
            pt_idx += 1

        if g_ratio_stat.GetN() > 0:
            c_ratio = ROOT.TCanvas(f"c_ratio_pd_2pp{suffix}", f"Ratio pd/2pp xF Bin {i_xf}", 800, 600)
            c_ratio.SetLeftMargin(0.16); c_ratio.SetBottomMargin(0.14); c_ratio.SetTickx(1); c_ratio.SetTicky(1)

            mg = ROOT.TMultiGraph()
            mg.SetTitle(f"#sigma_{{pd}} / 2#sigma_{{pp}} Ratio (xF Bin {i_xf});p_{{T}} [GeV];Ratio")

            fit_func = ROOT.TF1("fit_ratio", "pol0", 0.0, 2.0)
            g_ratio_sys.Fit(fit_func, "Q0")
            fit_val, fit_err = fit_func.GetParameter(0), fit_func.GetParError(0)

            g_fit_band = ROOT.TGraphErrors()
            g_fit_band.SetPoint(0, 0.0, fit_val); g_fit_band.SetPointError(0, 0.0, fit_err)
            g_fit_band.SetPoint(1, 2.0, fit_val); g_fit_band.SetPointError(1, 0.0, fit_err)
            g_fit_band.SetFillColorAlpha(ROOT.kPink, 0.4); g_fit_band.SetFillStyle(1001)

            for i in range(g_ratio_stat.GetN()):
                y = g_ratio_stat.GetY()[i]; err = max(g_ratio_stat.GetErrorY(i), g_ratio_sys.GetErrorY(i))
                if y + err > y_max_ratio: y_max_ratio = y + err
                if y - err < y_min_ratio: y_min_ratio = y - err

            if y_min_ratio < y_max_ratio:
                y_range = y_max_ratio - y_min_ratio
                mg.SetMinimum(max(0.0, y_min_ratio - 0.5 * y_range)); mg.SetMaximum(y_max_ratio + 0.6 * y_range)
            else: mg.SetMinimum(0.0); mg.SetMaximum(2.0)

            mg.Add(g_fit_band, "3")
            g_ratio_sys.SetMarkerSize(0); g_ratio_sys.SetLineColor(ROOT.kGreen+2); g_ratio_sys.SetFillColorAlpha(ROOT.kGreen-9, 0.5); g_ratio_sys.SetFillStyle(1001)
            mg.Add(g_ratio_sys, "2")
            g_ratio_stat.SetMarkerStyle(20); g_ratio_stat.SetMarkerColor(ROOT.kBlack); g_ratio_stat.SetLineColor(ROOT.kBlack)
            mg.Add(g_ratio_stat, "P")

            mg.Draw("A"); mg.GetXaxis().SetLimits(0.0, 2.0)
            
            line = ROOT.TLine(0.0, 1.0, 2.0, 1.0); line.SetLineStyle(2); line.SetLineColor(ROOT.kGray+2); line.SetLineWidth(2); line.Draw("SAME")
            fit_func.SetLineColor(ROOT.kRed); fit_func.SetLineWidth(2); fit_func.Draw("SAME")
            
            leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88); leg.SetBorderSize(0)
            leg.AddEntry(g_ratio_sys, "Systematic Unc.", "f"); leg.AddEntry(g_ratio_stat, "Data Ratio", "lep")
            leg.Draw()

            latex_fit = ROOT.TLatex(); latex_fit.SetTextFont(42); latex_fit.SetTextSize(0.04); latex_fit.SetTextColor(ROOT.kRed)
            y_text = fit_val + fit_err + ((mg.GetYaxis().GetXmax() - mg.GetYaxis().GetXmin()) * 0.02)
            latex_fit.DrawLatex(0.2, y_text, f"Best Fit: {fit_val:.4f} #pm {fit_err:.4f}")

            c_ratio.SaveAs(f"CrossSection_Ratio_pd_2pp_vs_pT{suffix}.pdf")
            c_ratio.Write(); c_ratio.Close()

        dir_ratio.cd()
        g_ratio_stat.Write(); g_ratio_sys.Write()

    def generate_combined_ratio_overlay_plot(self, i_xf, use_true_pt=False):
        suffix = f"_true_pt_xF{i_xf}" if use_true_pt else f"_geom_xF{i_xf}"

        dir_lh2 = self.out_file.Get("CrossSections_LH2")
        dir_ld2 = self.out_file.Get("CrossSections_LD2")
        dir_ratio = self.out_file.Get("Ratio_pd_2pp")

        if not all([dir_lh2, dir_ld2, dir_ratio]): return

        g_xsec_lh2 = dir_lh2.Get(f"g_xsec_LH2{suffix}")
        g_sys_lh2  = dir_lh2.Get(f"g_sys_LH2{suffix}")
        g_xsec_ld2 = dir_ld2.Get(f"g_xsec_LD2{suffix}")
        g_sys_ld2  = dir_ld2.Get(f"g_sys_LD2{suffix}")
        g_ratio_stat = dir_ratio.Get(f"g_ratio_stat{suffix}")
        g_ratio_sys  = dir_ratio.Get(f"g_ratio_sys{suffix}")

        if not all([g_xsec_lh2, g_sys_lh2, g_xsec_ld2, g_sys_ld2, g_ratio_stat, g_ratio_sys]): return

        # ----------------------------------------------------------------------
        # PLOT 4: 3-in-1 Combined Canvas (Overlays Top, Ratio Bottom)
        # ----------------------------------------------------------------------
        canvas = ROOT.TCanvas(f"c_combined_{suffix}", f"Combined Overlay xF Bin {i_xf}", 800, 800)

        # Pad 1: Cross Sections Overlay (Upper TFrame)
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.35, 1, 1.0)
        pad1.SetBottomMargin(0.15); pad1.SetLeftMargin(0.16); pad1.SetRightMargin(0.05); pad1.SetTopMargin(0.08)
        pad1.Draw(); pad1.cd()

        mg_top = ROOT.TMultiGraph()
        g_sys_lh2_c = g_sys_lh2.Clone(); g_sys_lh2_c.SetLineColor(ROOT.kRed); g_sys_lh2_c.SetFillColorAlpha(ROOT.kPink - 9, 0.5)
        g_xsec_lh2_c = g_xsec_lh2.Clone(); g_xsec_lh2_c.SetLineColor(ROOT.kRed); g_xsec_lh2_c.SetMarkerColor(ROOT.kRed)

        g_sys_ld2_c = g_sys_ld2.Clone(); g_sys_ld2_c.SetLineColor(ROOT.kBlue); g_sys_ld2_c.SetFillColorAlpha(ROOT.kAzure + 1, 0.5) 
        g_xsec_ld2_c = g_xsec_ld2.Clone(); g_xsec_ld2_c.SetLineColor(ROOT.kBlue); g_xsec_ld2_c.SetMarkerColor(ROOT.kBlue)

        mg_top.Add(g_sys_lh2_c, "2"); mg_top.Add(g_sys_ld2_c, "2")
        mg_top.Add(g_xsec_lh2_c, "P"); mg_top.Add(g_xsec_ld2_c, "P")

        mg_top.Draw("A"); mg_top.GetXaxis().SetLimits(0.0, 2.0)
        mg_top.SetTitle(f"DY Kinematics (xF Bin {i_xf})")
        mg_top.GetYaxis().SetTitle("d#sigma/dp_{T} [nb/GeV]"); mg_top.GetYaxis().CenterTitle(); mg_top.GetYaxis().SetTitleFont(43); mg_top.GetYaxis().SetTitleSize(22); mg_top.GetYaxis().SetTitleOffset(2.0); mg_top.GetYaxis().SetLabelFont(43); mg_top.GetYaxis().SetLabelSize(20)
        mg_top.GetXaxis().SetTitle("p_{T} [GeV]"); mg_top.GetXaxis().CenterTitle(); mg_top.GetXaxis().SetTitleFont(43); mg_top.GetXaxis().SetTitleSize(22); mg_top.GetXaxis().SetTitleOffset(1.2); mg_top.GetXaxis().SetLabelFont(43); mg_top.GetXaxis().SetLabelSize(20)

        y_max = max([g_xsec_lh2.GetY()[i] + g_sys_lh2.GetErrorY(i) for i in range(g_xsec_lh2.GetN())] + [0])
        y_max = max([g_xsec_ld2.GetY()[i] + g_sys_ld2.GetErrorY(i) for i in range(g_xsec_ld2.GetN())] + [y_max])
        mg_top.SetMinimum(0.0); mg_top.SetMaximum(y_max * 1.4 if y_max > 0 else 3.0)

        leg_top = ROOT.TLegend(0.65, 0.65, 0.88, 0.88); leg_top.SetBorderSize(0); leg_top.SetFillStyle(0); leg_top.SetTextFont(43); leg_top.SetTextSize(18)
        leg_top.AddEntry(g_xsec_lh2_c, "LH2 Data", "pl"); leg_top.AddEntry(g_sys_lh2_c, "LH2 Sys. Unc.", "f")
        leg_top.AddEntry(g_xsec_ld2_c, "LD2 Data", "pl"); leg_top.AddEntry(g_sys_ld2_c, "LD2 Sys. Unc.", "f")
        leg_top.Draw()

        # Pad 2: Ratio (Bottom TFrame)
        canvas.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.35)
        pad2.SetTopMargin(0.05); pad2.SetBottomMargin(0.35); pad2.SetLeftMargin(0.16); pad2.SetRightMargin(0.05)
        pad2.Draw(); pad2.cd()

        mg_bottom = ROOT.TMultiGraph()
        fit_func = ROOT.TF1("fit_ratio_comb", "pol0", 0.0, 2.0)
        if g_ratio_sys.GetN() > 0: g_ratio_sys.Fit(fit_func, "Q0")
        fit_val, fit_err = fit_func.GetParameter(0), fit_func.GetParError(0)

        g_fit_band = ROOT.TGraphErrors()
        g_fit_band.SetPoint(0, 0.0, fit_val); g_fit_band.SetPointError(0, 0.0, fit_err)
        g_fit_band.SetPoint(1, 2.0, fit_val); g_fit_band.SetPointError(1, 0.0, fit_err)
        g_fit_band.SetFillColorAlpha(ROOT.kPink, 0.4); g_fit_band.SetFillStyle(1001)

        g_ratio_sys_c = g_ratio_sys.Clone(); g_ratio_sys_c.SetLineColor(ROOT.kGreen+2); g_ratio_sys_c.SetFillColorAlpha(ROOT.kGreen-9, 0.5)
        
        mg_bottom.Add(g_fit_band, "3"); mg_bottom.Add(g_ratio_sys_c, "2"); mg_bottom.Add(g_ratio_stat.Clone(), "P")

        mg_bottom.Draw("A"); mg_bottom.GetXaxis().SetLimits(0.0, 2.0)
        mg_bottom.GetYaxis().SetTitle("#sigma_{pd}/2#sigma_{pp}"); mg_bottom.GetYaxis().CenterTitle(); mg_bottom.GetYaxis().SetTitleFont(43); mg_bottom.GetYaxis().SetTitleSize(22); mg_bottom.GetYaxis().SetTitleOffset(2.0); mg_bottom.GetYaxis().SetLabelFont(43); mg_bottom.GetYaxis().SetLabelSize(20); mg_bottom.GetYaxis().SetNdivisions(505) 
        mg_bottom.GetXaxis().SetTitle("p_{T} [GeV]"); mg_bottom.GetXaxis().CenterTitle(); mg_bottom.GetXaxis().SetTitleFont(43); mg_bottom.GetXaxis().SetTitleSize(22); mg_bottom.GetXaxis().SetTitleOffset(1.2); mg_bottom.GetXaxis().SetLabelFont(43); mg_bottom.GetXaxis().SetLabelSize(20)

        y_max_r, y_min_r = fit_val + fit_err, fit_val - fit_err
        for i in range(g_ratio_stat.GetN()):
            y = g_ratio_stat.GetY()[i]; err = max(g_ratio_stat.GetErrorY(i), g_ratio_sys.GetErrorY(i))
            if y + err > y_max_r: y_max_r = y + err
            if y - err < y_min_r: y_min_r = y - err

        if y_min_r < y_max_r:
            y_range = y_max_r - y_min_r
            mg_bottom.SetMinimum(max(0.0, y_min_r - 0.5 * y_range)); mg_bottom.SetMaximum(y_max_r + 0.6 * y_range)
        else: mg_bottom.SetMinimum(0.0); mg_bottom.SetMaximum(2.0)

        line = ROOT.TLine(0.0, 1.0, 2.0, 1.0); line.SetLineStyle(2); line.SetLineColor(ROOT.kGray+2); line.SetLineWidth(2); line.Draw("SAME")
        fit_func.SetLineColor(ROOT.kRed); fit_func.SetLineWidth(2); fit_func.Draw("SAME")

        latex_fit = ROOT.TLatex(); latex_fit.SetTextFont(43); latex_fit.SetTextSize(20); latex_fit.SetTextColor(ROOT.kRed)
        y_text = fit_val + fit_err + ((mg_bottom.GetYaxis().GetXmax() - mg_bottom.GetYaxis().GetXmin()) * 0.05)
        latex_fit.DrawLatex(0.2, y_text, f"Best Fit: {fit_val:.4f} #pm {fit_err:.4f}")

        canvas.SaveAs(f"Combined_XSec_Ratio_vs_pT{suffix}.pdf")
        dir_comb = self.get_or_create_dir(self.out_file, "Combined_Plots"); dir_comb.cd()
        canvas.Write(f"c_combined_ratio_overlay{suffix}"); canvas.Close()

    def calculate_cross_sections(self):
        for i_xf in range(self.n_xf_bins):
            if self.hists_lh2[i_xf] and self.hists_fl[i_xf]:
                self.sub_dict_lh2[i_xf] = self.generate_subtracted_plot(self.hists_lh2[i_xf], self.hists_fl[i_xf], config.FLASK_NORM_LH2, "LH2", i_xf)
                
            if self.hists_ld2[i_xf] and self.hists_lh2[i_xf] and self.hists_fl[i_xf]:
                self.sub_dict_pd[i_xf] = self.generate_pd_subtracted_plot(self.hists_ld2[i_xf], self.hists_lh2[i_xf], self.hists_fl[i_xf], i_xf)
                
            for use_true_pt in [False, True]:
                if self.sub_dict_lh2[i_xf]: self.calculate_and_plot_cross_section(self.sub_dict_lh2[i_xf], "LH2", config.GLOBAL_CONSTANT_LH2, i_xf, use_true_pt)
                if self.sub_dict_pd[i_xf]:  self.calculate_and_plot_cross_section(self.sub_dict_pd[i_xf], "LD2", config.GLOBAL_CONSTANT_LD2, i_xf, use_true_pt)
                    
                self.generate_ratio_plot(i_xf, use_true_pt)
                self.generate_combined_ratio_overlay_plot(i_xf, use_true_pt)

    def finalize(self):
        self.out_file.Write(); self.out_file.Close()