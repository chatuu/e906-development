"""
config.py
Configuration and constants for Drell-Yan Cross-Section Analysis.
"""

import math
import numpy as np
from rich.console import Console
from rich.table import Table

# ==========================================
# Individual Roadset POTs
# ==========================================
POT_LH2_MAP = {
    "RS57": 3.533324e+16,
    "RS59": 9.365986e+15,
    "RS62": 5.281659e+16,
    "RS67": 1.611435e+17,
    "RS70": 1.785745e+16
}

POT_LD2_MAP = {
    "RS57": 1.768358e+16,
    "RS59": 4.319952e+15,
    "RS62": 2.382899e+16,
    "RS67": 7.694541e+16,
    "RS70": 8.752588e+15
}

POT_FLASK_MAP = {
    "RS57": 3.918550e+15,
    "RS59": 1.010350e+15,
    "RS62": 1.097456e+16,
    "RS67": 3.662417e+16,
    "RS70": 3.841280e+15
}

# Default to Total (Sum) if not explicitly overridden
PROTONS_ON_TARGET_LH2 = sum(POT_LH2_MAP.values())
PROTONS_ON_TARGET_LD2 = sum(POT_LD2_MAP.values())
PROTONS_ON_TARGET_FLASK = sum(POT_FLASK_MAP.values())

# ==========================================
# Physics Constants
# ==========================================
INPUT_NPZ_FILE = "/root/github/e906-development/src/kTrackerEfficiency/RS67/GlobalEfficiencyCurve/interpolation_data_d1.npz"

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
# Base Attenuations
# ==========================================
val_exp_LH2 = -(LH2_TARGET_LENGTH_CM * LH2_TARGET_DENSITY_MOL_CM3) / NUCLEAR_INTERACTION_LENGTH_LH2_GPERCM2
BEAM_ATTENUATION_LH2 = (NUCLEAR_INTERACTION_LENGTH_LH2_GPERCM2 / (LH2_TARGET_DENSITY_MOL_CM3 * LH2_TARGET_LENGTH_CM)) * (1.0 - math.exp(val_exp_LH2))

val_exp_LD2 = -(LD2_TARGET_LENGTH_CM * LD2_TARGET_DENSITY_MOL_CM3) / NUCLEAR_INTERACTION_LENGTH_LD2_GPERCM2
BEAM_ATTENUATION_LD2 = (NUCLEAR_INTERACTION_LENGTH_LD2_GPERCM2 / (LD2_TARGET_DENSITY_MOL_CM3 * LD2_TARGET_LENGTH_CM)) * (1.0 - math.exp(val_exp_LD2))

# ==========================================
# Derived Normalizations (Defaults)
# ==========================================
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
MASS_BINS = np.array([4.2, 4.5, 4.8, 5.1, 5.4, 5.7, 6.0, 6.3, 6.6, 6.9, 7.5, 8.8], dtype=float)
PT_BINS = np.array([0., 0.32, 0.49, 0.63, 0.77, 0.95, 1.18, 1.8], dtype=float)


def set_roadset(rs):
    """Overrides global POTs and recalculates all derived normalizations for a specific roadset."""
    global PROTONS_ON_TARGET_LH2, PROTONS_ON_TARGET_LD2, PROTONS_ON_TARGET_FLASK
    global GLOBAL_CONSTANT_LH2, GLOBAL_CONSTANT_LD2, FLASK_NORM_LH2, FLASK_NORM_LD2, LH2_TO_LD2_NORM
    
    if rs == "Combined" or rs == "All":
        PROTONS_ON_TARGET_LH2 = sum(POT_LH2_MAP.values())
        PROTONS_ON_TARGET_LD2 = sum(POT_LD2_MAP.values())
        PROTONS_ON_TARGET_FLASK = sum(POT_FLASK_MAP.values())
    elif rs in POT_LH2_MAP:
        PROTONS_ON_TARGET_LH2 = POT_LH2_MAP[rs]
        PROTONS_ON_TARGET_LD2 = POT_LD2_MAP[rs]
        PROTONS_ON_TARGET_FLASK = POT_FLASK_MAP[rs]
    else:
        return

    GLOBAL_CONSTANT_LH2 = (NUCLEONS_PER_NUCLEUS_LH2 * 1e33) / (
        LH2_TARGET_DENSITY_MOL_CM2 * AVOGADRO_CONSTANT * PROTONS_ON_TARGET_LH2 * BEAM_ATTENUATION_LH2
    )
    GLOBAL_CONSTANT_LD2 = (NUCLEONS_PER_NUCLEUS_LD2 * 1e33) / (
        LD2_TARGET_DENSITY_MOL_CM2 * AVOGADRO_CONSTANT * PROTONS_ON_TARGET_LD2 * BEAM_ATTENUATION_LD2
    )
    
    FLASK_NORM_LH2 = PROTONS_ON_TARGET_LH2 / PROTONS_ON_TARGET_FLASK
    FLASK_NORM_LD2 = PROTONS_ON_TARGET_LD2 / PROTONS_ON_TARGET_FLASK
    LH2_TO_LD2_NORM = THD_THH_RATIO * (PROTONS_ON_TARGET_LD2 / PROTONS_ON_TARGET_LH2)


def print_physics_constants():
    """
    Renders a colorful terminal table displaying the initialized physics constants.
    """
    console = Console()
    table = Table(title="Drell-Yan Physics Constants & Normalizations", header_style="bold magenta")
    
    table.add_column("Parameter", style="cyan", justify="right")
    table.add_column("Value", style="green", justify="left")

    table.add_row("PoT LH2", f"{PROTONS_ON_TARGET_LH2:.4e}")
    table.add_row("PoT LD2", f"{PROTONS_ON_TARGET_LD2:.4e}")
    table.add_row("PoT Flask", f"{PROTONS_ON_TARGET_FLASK:.4e}")
    table.add_row("LH2 Target Density (mol/cm²)", f"{LH2_TARGET_DENSITY_MOL_CM2:.4f}")
    table.add_row("LD2 Target Density (mol/cm²)", f"{LD2_TARGET_DENSITY_MOL_CM2:.4f}")
    table.add_row("LH2 Target Length (cm)", f"{LH2_TARGET_LENGTH_CM:.1f}")
    table.add_row("LD2 Target Length (cm)", f"{LD2_TARGET_LENGTH_CM:.1f}")
    table.add_row("Nuclear Int. Length LH2 (g/cm²)", f"{NUCLEAR_INTERACTION_LENGTH_LH2_GPERCM2:.1f}")
    table.add_row("Nuclear Int. Length LD2 (g/cm²)", f"{NUCLEAR_INTERACTION_LENGTH_LD2_GPERCM2:.1f}")
    table.add_row("pT Bins Definition", str(PT_BINS))
    table.add_row("Global Constant LH2 (w/o pT width)", f"{GLOBAL_CONSTANT_LH2:.4e}")
    table.add_row("Global Constant LD2 (w/o pT width)", f"{GLOBAL_CONSTANT_LD2:.4e}")
    table.add_row("Flask Norm LH2", f"{FLASK_NORM_LH2:.4f}")
    table.add_row("Flask Norm LD2", f"{FLASK_NORM_LD2:.4f}")
    table.add_row("LH2 to LD2 Norm", f"{LH2_TO_LD2_NORM:.4f}")

    console.print(table)