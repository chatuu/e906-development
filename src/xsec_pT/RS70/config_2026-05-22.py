"""
config.py
Configuration and constants for Drell-Yan Cross-Section Analysis.
"""

import math
import numpy as np
from rich.console import Console
from rich.table import Table

# ==========================================
# Physics Constants
# ==========================================
INPUT_NPZ_FILE = "/root/github/e906-development/src/kTrackerEfficiency/RS67/GlobalEfficiencyCurve/interpolation_data_d1.npz"

PROTONS_ON_TARGET_LH2 = (3.533324e+16 + 9.365986e+15 + 5.654075e+16 + 1.611435e+17 + 1.785745e+16)
PROTONS_ON_TARGET_LD2 = (1.768358e+16 + 4.319952e+15 + 2.517904e+16 + 7.694541e+16 + 8.752588e+15)
PROTONS_ON_TARGET_FLASK = (3.918550e+15 + 1.010350e+15 + 1.176106e+16 + 3.662417e+16 + 3.84128e+15)

LH2_TARGET_DENSITY_MOL_CM2 = 3.5966
LD2_TARGET_DENSITY_MOL_CM2 = 8.0431

LH2_TARGET_LENGTH_CM = 50.8
LD2_TARGET_LENGTH_CM = 50.8

LH2_TARGET_DENSITY_MOL_CM3 = 0.0708
LD2_TARGET_DENSITY_MOL_CM3 = 0.163

AVOGADRO_CONSTANT = 6.022e23
NUCLEONS_PER_NUCLEUS_LH2 = 1.008
NUCLEONS_PER_NUCLEUS_LD2 = 2.014

NUCLEAR_INTERACTION_LENGTH_LH2_GPERCM2 = 52.0
NUCLEAR_INTERACTION_LENGTH_LD2_GPERCM2 = 54.7

TARGET_THICKNESS_THD = 0.1084
TARGET_THICKNESS_THH = 3.5966
THD_THH_RATIO = TARGET_THICKNESS_THD / TARGET_THICKNESS_THH

# ==========================================
# Derived Normalizations
# ==========================================
val_exp_LH2 = -(LH2_TARGET_LENGTH_CM * LH2_TARGET_DENSITY_MOL_CM3) / NUCLEAR_INTERACTION_LENGTH_LH2_GPERCM2
BEAM_ATTENUATION_LH2 = (NUCLEAR_INTERACTION_LENGTH_LH2_GPERCM2 / (LH2_TARGET_DENSITY_MOL_CM3 * LH2_TARGET_LENGTH_CM)) * (1.0 - math.exp(val_exp_LH2))

val_exp_LD2 = -(LD2_TARGET_LENGTH_CM * LD2_TARGET_DENSITY_MOL_CM3) / NUCLEAR_INTERACTION_LENGTH_LD2_GPERCM2
BEAM_ATTENUATION_LD2 = (NUCLEAR_INTERACTION_LENGTH_LD2_GPERCM2 / (LD2_TARGET_DENSITY_MOL_CM3 * LD2_TARGET_LENGTH_CM)) * (1.0 - math.exp(val_exp_LD2))

GLOBAL_CONSTANT_LH2 = (NUCLEONS_PER_NUCLEUS_LH2 * 1e33) / (
    LH2_TARGET_DENSITY_MOL_CM2 * AVOGADRO_CONSTANT * PROTONS_ON_TARGET_LH2 * BEAM_ATTENUATION_LH2
)

GLOBAL_CONSTANT_LD2 = (NUCLEONS_PER_NUCLEUS_LD2 * 1e33) / (
    LD2_TARGET_DENSITY_MOL_CM2 * AVOGADRO_CONSTANT * PROTONS_ON_TARGET_LD2 * BEAM_ATTENUATION_LD2
)

FLASK_NORM_LH2 = PROTONS_ON_TARGET_LH2 / PROTONS_ON_TARGET_FLASK
FLASK_NORM_LD2 = PROTONS_ON_TARGET_LD2 / PROTONS_ON_TARGET_FLASK
LH2_TO_LD2_NORM = THD_THH_RATIO * (PROTONS_ON_TARGET_LD2 / PROTONS_ON_TARGET_LH2)

# ==========================================
# Kinematic Bins
# ==========================================
MASS_BINS = np.array([4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.7], dtype=float)
PT_BINS = np.array([0., 0.32, 0.49, 0.63, 0.77, 0.95, 1.18, 1.8], dtype=float)
#XF_BINS = np.array([0.1, 0.16, 0.22, 0.28, 0.34, 0.40, 0.46, 0.58, 0.70, 0.95], dtype=float)
XF_BINS = np.round(np.arange(0.0, 0.85, 0.05), 2)

def print_physics_constants():
    console = Console()
    table = Table(title="Drell-Yan Physics Constants & Normalizations", header_style="bold magenta")
    table.add_column("Parameter", style="cyan", justify="right")
    table.add_column("Value", style="green", justify="left")

    table.add_row("PoT LH2", f"{PROTONS_ON_TARGET_LH2:.4e}")
    table.add_row("PoT LD2", f"{PROTONS_ON_TARGET_LD2:.4e}")
    table.add_row("PoT Flask", f"{PROTONS_ON_TARGET_FLASK:.4e}")
    table.add_row("pT Bins Definition", str(PT_BINS))
    table.add_row("xF Bins Definition", str(XF_BINS))
    table.add_row("Global Constant LH2", f"{GLOBAL_CONSTANT_LH2:.4e}")
    table.add_row("Global Constant LD2", f"{GLOBAL_CONSTANT_LD2:.4e}")

    console.print(table)