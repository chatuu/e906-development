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

PROTONS_ON_TARGET_LH2 = 5.654075e+16
PROTONS_ON_TARGET_LD2 = 2.517904e+16
PROTONS_ON_TARGET_FLASK = 1.176106e+16

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

def print_physics_constants():
    """
    Renders a colorful terminal table displaying the initialized physics constants.
    """
    console = Console()
    table = Table(title="Drell-Yan Physics Constants & Normalizations", header_style="bold magenta")
    
    table.add_column("Parameter", style="cyan", justify="right")
    table.add_column("Value", style="green", justify="left")

    # Raw Constants
    table.add_row("PoT LH2", f"{PROTONS_ON_TARGET_LH2:.4e}")
    table.add_row("PoT LD2", f"{PROTONS_ON_TARGET_LD2:.4e}")
    table.add_row("PoT Flask", f"{PROTONS_ON_TARGET_FLASK:.4e}")
    table.add_row("LH2 Target Density (mol/cm²)", f"{LH2_TARGET_DENSITY_MOL_CM2:.4f}")
    table.add_row("LD2 Target Density (mol/cm²)", f"{LD2_TARGET_DENSITY_MOL_CM2:.4f}")
    table.add_row("LH2 Target Length (cm)", f"{LH2_TARGET_LENGTH_CM:.1f}")
    table.add_row("LD2 Target Length (cm)", f"{LD2_TARGET_LENGTH_CM:.1f}")
    table.add_row("LH2 Target Density (mol/cm³)", f"{LH2_TARGET_DENSITY_MOL_CM3:.4f}")
    table.add_row("LD2 Target Density (mol/cm³)", f"{LD2_TARGET_DENSITY_MOL_CM3:.4f}")
    table.add_row("Avogadro Constant", f"{AVOGADRO_CONSTANT:.3e}")
    table.add_row("Nucleons/Nucleus LH2", f"{NUCLEONS_PER_NUCLEUS_LH2:.3f}")
    table.add_row("Nucleons/Nucleus LD2", f"{NUCLEONS_PER_NUCLEUS_LD2:.3f}")
    table.add_row("Nuclear Int. Length LH2 (g/cm²)", f"{NUCLEAR_INTERACTION_LENGTH_LH2_GPERCM2:.1f}")
    table.add_row("Nuclear Int. Length LD2 (g/cm²)", f"{NUCLEAR_INTERACTION_LENGTH_LD2_GPERCM2:.1f}")
    table.add_row("Target Thickness THD", f"{TARGET_THICKNESS_THD:.4f}")
    table.add_row("Target Thickness THH", f"{TARGET_THICKNESS_THH:.4f}")
    table.add_row("pT Bins Definition", str(PT_BINS))
    
    # Derived Normalizations
    table.add_row("Global Constant LH2 (w/o pT width)", f"{GLOBAL_CONSTANT_LH2:.4e}")
    table.add_row("Global Constant LD2 (w/o pT width)", f"{GLOBAL_CONSTANT_LD2:.4e}")
    table.add_row("Flask Norm LH2", f"{FLASK_NORM_LH2:.4f}")
    table.add_row("Flask Norm LD2", f"{FLASK_NORM_LD2:.4f}")
    table.add_row("LH2 to LD2 Norm", f"{LH2_TO_LD2_NORM:.4f}")

    console.print(table)