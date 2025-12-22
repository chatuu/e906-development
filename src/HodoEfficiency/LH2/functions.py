import numpy as np
import pandas as pd
import ROOT
import sys
import array
import os

# Global cache for efficiency map
HODO_EFF_MAP = None
HODO_EFF_FILE = "hodoscope_eff.tsv"

# Define Bins for Histogram axes
MASS_BINS = [4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7]
XF_BINS = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8]

def apply_data_cuts(df: pd.DataFrame, beam_offset: float = 1.604) -> pd.Series:
    """Applies the 2111v42 analysis cuts to a Pandas DataFrame."""
    cut_pos = (
        (df['chisq1_target'] < 15) & (df['pz1_st1'] > 9) & (df['pz1_st1'] < 75) &
        (df['nHits1'] > 13) & ((df['x1_t']**2 + (df['y1_t'] - beam_offset)**2) < 320) &
        ((df['x1_d']**2 + (df['y1_d'] - beam_offset)**2) < 1100) &
        ((df['x1_d']**2 + (df['y1_d'] - beam_offset)**2) > 8) &
        (df['chisq1_target'] < 1.5 * df['chisq1_upstream']) &
        (df['chisq1_target'] < 1.5 * df['chisq1_dump']) & (df['z1_v'] < -5) & 
        (df['z1_v'] > -320) & ((df['chisq1'] / (df['nHits1'] - 5)) < 12) &
        ((df['y1_st1'] / df['y1_st3']) < 1) &
        (np.abs(np.abs(df['px1_st1'] - df['px1_st3']) - 0.416) < 0.008) &
        (np.abs(df['py1_st1'] - df['py1_st3']) < 0.008) &
        (np.abs(df['pz1_st1'] - df['pz1_st3']) < 0.08) &
        ((df['y1_st1'] * df['y1_st3']) > 0) & (np.abs(df['py1_st1']) > 0.02)
    )
    cut_neg = (
        (df['chisq2_target'] < 15) & (df['pz2_st1'] > 9) & (df['pz2_st1'] < 75) &
        (df['nHits2'] > 13) & ((df['x2_t']**2 + (df['y2_t'] - beam_offset)**2) < 320) &
        ((df['x2_d']**2 + (df['y2_d'] - beam_offset)**2) < 1100) &
        ((df['x2_d']**2 + (df['y2_d'] - beam_offset)**2) > 8) &
        (df['chisq2_target'] < 1.5 * df['chisq2_upstream']) &
        (df['chisq2_target'] < 1.5 * df['chisq2_dump']) & (df['z2_v'] < -5) & 
        (df['z2_v'] > -320) & ((df['chisq2'] / (df['nHits2'] - 5)) < 12) &
        ((df['y2_st1'] / df['y2_st3']) < 1) &
        (np.abs(np.abs(df['px2_st1'] - df['px2_st3']) - 0.416) < 0.008) &
        (np.abs(df['py2_st1'] - df['py2_st3']) < 0.008) &
        (np.abs(df['pz2_st1'] - df['pz2_st3']) < 0.08) &
        ((df['y2_st1'] * df['y2_st3']) > 0) & (np.abs(df['py2_st1']) > 0.02)
    )
    cut_dimuon = (
        (np.abs(df['dx']) < 0.25) & (np.abs(df['dy'] - beam_offset) < 0.22) &
        (df['dz'] > -280) & (df['dz'] < -5) & (np.abs(df['dpx']) < 1.8) &
        (np.abs(df['dpy']) < 2) & ((df['dpx']**2 + df['dpy']**2) < 5) &
        (df['dpz'] > 38) & (df['dpz'] < 116) &
        ((df['dx']**2 + (df['dy'] - beam_offset)**2) < 0.06) &
        (np.abs(df['trackSeparation']) < 270) & (df['chisq_dimuon'] < 18) &
        (np.abs(df['chisq1_target'] + df['chisq2_target'] - df['chisq_dimuon']) < 2) &
        ((df['y1_st3'] * df['y2_st3']) < 0) & ((df['nHits1'] + df['nHits2']) > 29) &
        ((df['nHits1St1'] + df['nHits2St1']) > 8) & (np.abs(df['x1_st1'] + df['x2_st1']) < 42)
    )
    cut_physics = (
        (df['mass'] > 4.2) & (df['mass'] < 8.8) & (df['xF'] < 0.95) &
        (df['xF'] > -0.1) & (df['xT'] > 0.05) & (df['xT'] < 0.55) & (np.abs(df['costh']) < 0.5)
    )
    cut_occ = (
        (df['D1'] < 400) & (df['D2'] < 400) & (df['D3'] < 400) &
        ((df['D1'] + df['D2'] + df['D3']) < 1000)
    )
    return cut_pos & cut_neg & cut_dimuon & cut_physics & cut_occ & (df['mass'] < 8.8)

def get_cut_columns():
    """Returns list of columns needed for the cut logic."""
    return [
        'chisq1_target', 'pz1_st1', 'nHits1', 'x1_t', 'y1_t', 'x1_d', 'y1_d',
        'chisq1_upstream', 'chisq1_dump', 'z1_v', 'chisq1', 'y1_st1', 'y1_st3',
        'px1_st1', 'px1_st3', 'py1_st1', 'py1_st3', 'pz1_st3',
        'chisq2_target', 'pz2_st1', 'nHits2', 'x2_t', 'y2_t', 'x2_d', 'y2_d',
        'chisq2_upstream', 'chisq2_dump', 'z2_v', 'chisq2', 'y2_st1', 'y2_st3',
        'px2_st1', 'px2_st3', 'py2_st1', 'py2_st3', 'pz2_st3',
        'dx', 'dy', 'dz', 'dpx', 'dpy', 'dpz', 'trackSeparation', 'chisq_dimuon',
        'nHits1St1', 'nHits2St1', 'x1_st1', 'x2_st1',
        'mass', 'xF', 'xT', 'costh',
        'D1', 'D2', 'D3'
    ]

def load_data_from_root(filepath, treename, columns, filter_str=""):
    """
    Loads specific columns from a ROOT tree into a Pandas DataFrame.
    """
    try:
        rdf = ROOT.RDataFrame(treename, filepath)
        if filter_str:
            rdf = rdf.Filter(filter_str)
        data_dict = rdf.AsNumpy(columns=columns)
        return pd.DataFrame(data_dict)
    except Exception as e:
        print(f"Error loading {treename} from {filepath}: {e}", file=sys.stderr)
        return pd.DataFrame()

def decode_hodo_ids(road_id):
    rid = int(road_id)
    abs_rid = abs(rid)
    h1_id = int((abs_rid - 1) / (16**3) + 1)
    h2_id = int(((abs_rid - 1) / (16**2)) % 16 + 1)
    h3_id = int(((abs_rid - 1) / 16) % 16 + 1)
    h4_id = int((abs_rid - 1) % 16 + 1)
    if rid < 0: labels = {"H1": "H1B", "H2": "H2B", "H3": "H3B", "H4": "H4B"}
    else: labels = {"H1": "H1T", "H2": "H2T", "H3": "H3T", "H4": "H4T"}
    ids = {"H1": h1_id, "H2": h2_id, "H3": h3_id, "H4": h4_id}
    return ids, labels

def load_hodo_eff_map():
    global HODO_EFF_MAP
    if HODO_EFF_MAP is not None: return HODO_EFF_MAP
    mapping = {}
    try:
        df = pd.read_csv(HODO_EFF_FILE, sep=r'\s+', header=None)
        for _, row in df.iterrows():
            mapping[(str(row.iloc[0]).strip(), int(row.iloc[1]))] = (float(row.iloc[2]), float(row.iloc[3]), float(row.iloc[4]))
        HODO_EFF_MAP = mapping
    except: HODO_EFF_MAP = {}
    return HODO_EFF_MAP

def process_single_bin(df_merged, mass_low, mass_high, xf_low, xf_high):
    """Filters for a specific bin, prints results, and fills ROOT histograms."""
    
    # Load Efficiency Map (Lazy load)
    eff_map = load_hodo_eff_map()

    # Filter for Mass and xF Bin
    df_mass_bin = df_merged[(df_merged['mass'] >= mass_low) & (df_merged['mass'] < mass_high)]
    df_final = df_mass_bin[(df_mass_bin['xF'] >= xf_low) & (df_mass_bin['xF'] < xf_high)]
    
    if not df_final.empty:
        total_entries = len(df_final)
        print(f"BIN_RESULT: Mass[{mass_low}-{mass_high}] xF[{xf_low}-{xf_high}] Count={total_entries}")
        
        # Determine xF Bin Index for filename
        xf_index = -1
        for i, val in enumerate(XF_BINS):
            if abs(val - xf_low) < 1e-5:
                xf_index = i
                break
        
        root_filename = f"HodoEff_xF_{xf_index}.root"
        root_file = ROOT.TFile(root_filename, "UPDATE")
        
        # Prepare arrays for variable binning
        mass_bins_arr = array.array('d', MASS_BINS)
        xf_bins_arr = array.array('d', XF_BINS)

        # -------------------------------------------------------------
        # 1. TH2I for Number of Entries (Total Counts)
        # -------------------------------------------------------------
        hist_counts_name = "h_bin_counts"
        h_counts = root_file.Get(hist_counts_name)
        if not h_counts:
            h_counts = ROOT.TH2I(hist_counts_name, "Number of Entries per Bin; Mass (GeV); xF", 
                                 len(MASS_BINS)-1, mass_bins_arr, 
                                 len(XF_BINS)-1, xf_bins_arr)
        
        # Find global bin for the current Mass/xF center
        mass_center = (mass_low + mass_high) / 2.0
        xf_center = (xf_low + xf_high) / 2.0
        global_bin = h_counts.FindBin(mass_center, xf_center)
        
        # Fill the total count for this bin
        h_counts.SetBinContent(global_bin, total_entries)
        
        # -------------------------------------------------------------
        # 1D Histogram Setup
        # -------------------------------------------------------------
        h_name_suffix = f"Mass{mass_low}_{mass_high}_xF{xf_low}_{xf_high}".replace('.', 'p')
        hist_name = f"h_{h_name_suffix}"
        
        h_eff = root_file.Get(hist_name)
        if not h_eff:
            h_eff = ROOT.TH1D(hist_name, "", 110, 0.0, 1.1)
            h_eff.Sumw2()
            
        # =========================================================================
        # NEW: Element ID Histograms Setup
        # =========================================================================
        station_labels = ['H1T', 'H2T', 'H3T', 'H4T', 'H1B', 'H2B', 'H3B', 'H4B']
        h_elem_counts = {}

        for label in station_labels:
            # Consistent naming scheme for the element ID histograms
            h_elem_name = f"h_ElementID_{label}_{h_name_suffix}"
            h_temp = root_file.Get(h_elem_name)
            
            if not h_temp:
                # Define binning 0 to 50 (should cover typical element IDs)
                h_temp = ROOT.TH1I(h_elem_name, f"Element ID Counts {label}; Element ID; Counts", 50, 0.5, 50.5)
            
            h_elem_counts[label] = h_temp
        # =========================================================================

        zero_eff_count = 0 # Counter for zero efficiency events

        # Loop over Events
        for index, row in df_final.iterrows():
            # Basic Event Info
            run, spill, event = int(row['runID']), int(row['spillID']), int(row['eventID'])
            
            # Decode Pos and Neg Roads
            pos_ids, pos_lbls = decode_hodo_ids(row['posRoad'])
            neg_ids, neg_lbls = decode_hodo_ids(row['negRoad'])
            
            # =========================================================================
            # NEW: Fill Element ID Histograms
            # =========================================================================
            # Loop through the decoded stations (H1..H4) for both tracks
            for station_key in ['H1', 'H2', 'H3', 'H4']:
                # Positive Track: Get Label (e.g., 'H1T') and ID (e.g., 5)
                p_lbl = pos_lbls[station_key]
                p_id = pos_ids[station_key]
                if p_lbl in h_elem_counts:
                    h_elem_counts[p_lbl].Fill(p_id)

                # Negative Track: Get Label (e.g., 'H1B') and ID (e.g., 12)
                n_lbl = neg_lbls[station_key]
                n_id = neg_ids[station_key]
                if n_lbl in h_elem_counts:
                    h_elem_counts[n_lbl].Fill(n_id)
            # =========================================================================

            # --- Helper to calculate efficiency ---
            def calculate_eff(road_val, ids, lbls):
                eff_product = 1.0
                for h in ['H1', 'H2', 'H3', 'H4']:
                    key = (lbls[h], ids[h])
                    if key in eff_map:
                        eff_product *= eff_map[key][0]
                return eff_product

            #print(f"Positive: labels={pos_lbls}, ids={pos_ids} | Negative: labels={neg_lbls}, ids={neg_ids}")
            # if (((pos_lbls['H1'] == 'H1T' and pos_ids['H1']==5) or (pos_lbls['H1'] == 'H1T' and pos_ids['H1']==19)) or ((pos_lbls['H1'] == 'H1B' and pos_ids['H1']==5) or (pos_lbls['H1'] == 'H1B' and pos_ids['H1']==7))):
            #     continue
            # if (((neg_lbls['H1'] == 'H1T' and neg_ids['H1']==5) or (neg_lbls['H1'] == 'H1T' and neg_ids['H1']==19)) or ((neg_lbls['H1'] == 'H1B' and neg_ids['H1']==5) or (neg_lbls['H1'] == 'H1B' and neg_ids['H1']==7))):
            #     continue

            pos_eff = calculate_eff(row['posRoad'], pos_ids, pos_lbls)
            neg_eff = calculate_eff(row['negRoad'], neg_ids, neg_lbls)

            # --- Check for Zero Efficiency ---
            if pos_eff == 0.0 or neg_eff == 0.0:
                print(f"Positive labels={pos_lbls}, ids={pos_ids} | Negative labels={neg_lbls}, ids={neg_ids}")
                zero_eff_count += 1
                
                # CSV Export Logic
                csv_filename = f"EffZeroEvents_xF_{xf_index}.csv"
                zero_data = {
                    'run': [run], 'spill': [spill], 'event': [event],
                    'posRoad': [row['posRoad']], 'negRoad': [row['negRoad']],
                    'pos_eff': [pos_eff], 'neg_eff': [neg_eff]
                }
                df_zero = pd.DataFrame(zero_data)
                header_needed = not os.path.exists(csv_filename)
                df_zero.to_csv(csv_filename, mode='a', header=header_needed, index=False)
                
            total_eff = pos_eff * neg_eff
            h_eff.Fill(total_eff)
        
        # -------------------------------------------------------------
        # 2. TH2I for Number of Zero Efficiency Entries (NEW)
        # -------------------------------------------------------------
        hist_zero_counts_name = "h_zero_eff_counts"
        h_zero_counts = root_file.Get(hist_zero_counts_name)
        if not h_zero_counts:
            h_zero_counts = ROOT.TH2I(hist_zero_counts_name, "Count of Zero Eff Events; Mass (GeV); xF", 
                                      len(MASS_BINS)-1, mass_bins_arr, 
                                      len(XF_BINS)-1, xf_bins_arr)
        
        # Fill the zero count for this bin
        h_zero_counts.SetBinContent(global_bin, zero_eff_count)

        # -------------------------------------------------------------
        # 3. TH2D for Percentage of Zero Efficiency Events
        # -------------------------------------------------------------
        hist_zero_name = "h_zero_eff_percentage"
        h_zero_perc = root_file.Get(hist_zero_name)
        if not h_zero_perc:
            h_zero_perc = ROOT.TH2D(hist_zero_name, "Percentage of Zero Eff Events; Mass (GeV); xF", 
                                    len(MASS_BINS)-1, mass_bins_arr, 
                                    len(XF_BINS)-1, xf_bins_arr)
        
        # Calculate Percentage
        perc_zero = (zero_eff_count / total_entries * 100.0) if total_entries > 0 else 0.0
        h_zero_perc.SetBinContent(global_bin, perc_zero)

        # -------------------------------------------------------------
        # Update 1D Histogram Stats and Write All
        # -------------------------------------------------------------
        mean_eff = h_eff.GetMean()
        std_dev = h_eff.GetStdDev()
        h_eff.SetTitle(f"Mean Efficiency {mean_eff:.4f} #pm {std_dev:.4f}")
        
        h_counts.Write("", ROOT.TObject.kOverwrite)
        h_zero_counts.Write("", ROOT.TObject.kOverwrite)
        h_zero_perc.Write("", ROOT.TObject.kOverwrite)
        h_eff.Write("", ROOT.TObject.kOverwrite)

        # =========================================================================
        # NEW: Write Element ID Histograms
        # =========================================================================
        for label, h_obj in h_elem_counts.items():
            h_obj.Write("", ROOT.TObject.kOverwrite)
        # =========================================================================

        # Handle 2D Mean Efficiency Map (Existing logic)
        hist_2d_name = "hodoeffcorr"
        h_2d = root_file.Get(hist_2d_name)
        if not h_2d:
            h_2d = ROOT.TH2D(hist_2d_name, "Hodoscope Efficiency Correction; Mass (GeV); xF", 
                             len(MASS_BINS)-1, mass_bins_arr, 
                             len(XF_BINS)-1, xf_bins_arr)

        h_2d.SetBinContent(global_bin, mean_eff)
        h_2d.SetBinError(global_bin, std_dev)
        h_2d.Write("", ROOT.TObject.kOverwrite)
        
        root_file.Close()