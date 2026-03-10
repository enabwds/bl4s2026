"""
add_x0_z.py
-----------
Adds X0_cm and Z columns to shower_output.csv based on the material column.
Also adds absorber_thickness_cm for convenience.

Usage:
    python3 add_x0_z.py                              # reads/writes shower_output.csv
    python3 add_x0_z.py --input myfile.csv           # custom input
    python3 add_x0_z.py --input a.csv --output b.csv # custom input and output
"""

import argparse
import pandas as pd

# ── Material properties (from PDG / your experiment spec) ──────────────────
X0_map = {
    'G4_Al': 8.897,   # cm
    'G4_Fe': 1.757,   # cm
    'G4_Cu': 1.436,   # cm
    'G4_Pb': 0.5612,  # cm
    'Brass': 1.54,    # cm (70% Cu, 30% Zn) — add after unblinding run_blind.mac
}

Z_map = {
    'G4_Al': 13,
    'G4_Fe': 26,
    'G4_Cu': 29,
    'G4_Pb': 82,
    'Brass': 27.4,    # effective Z (70% Cu @ Z=29, 30% Zn @ Z=30)
}

density_map = {          # g/cm³ — useful for later analysis
    'G4_Al': 2.70,
    'G4_Fe': 7.87,
    'G4_Cu': 8.96,
    'G4_Pb': 11.35,
    'Brass': 8.5,     # typical brass density
}

critical_energy_map = {  # MeV — electron critical energy (PDG 2023)
    'G4_Al': 42.3,
    'G4_Fe': 21.2,
    'G4_Cu': 19.0,
    'G4_Pb':  7.43,  # corrected from 6.71 (PDG 2023 electron value)
    'Brass': 19.2,    # interpolated between Cu and Zn
}

moliere_radius_map = {   # cm — sets lateral shower scale
    'G4_Al': 4.4,
    'G4_Fe': 1.72,
    'G4_Cu': 1.57,
    'G4_Pb': 1.60,
    'Brass': 1.59,    # very similar to Cu (high-Z, similar shower structure)
}

# ── Argument parsing ───────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument('--input',  default='shower_output.csv')
parser.add_argument('--output', default=None)   # defaults to overwriting input
args = parser.parse_args()

output_file = args.output if args.output else args.input

# ── Load ───────────────────────────────────────────────────────────────────
print(f"Loading {args.input} ...")
df = pd.read_csv(args.input)
print(f"  {len(df):,} rows, {len(df.columns)} columns")
print(f"  Materials found: {sorted(df['material'].unique())}")

# ── Check for unknown materials ────────────────────────────────────────────
unknown = set(df['material'].unique()) - set(X0_map.keys())
if unknown:
    print(f"\nWARNING: Unknown materials (will get NaN): {unknown}")
    print("  Add them to the maps at the top of this script if needed.")

# ── Add columns ────────────────────────────────────────────────────────────
df['X0_cm']              = df['material'].map(X0_map)
df['Z']                  = df['material'].map(Z_map)
df['density_g_cm3']      = df['material'].map(density_map)
df['Ec_MeV']             = df['material'].map(critical_energy_map)
df['moliere_radius_cm']  = df['material'].map(moliere_radius_map)
df['thickness_cm']       = df['absorber_thickness_mm'] / 10.0

# Sanity check: recompute absorber_x0 and compare to GEANT4's value.
# NOTE: absorber_x0 in the CSV uses GEANT4's internal Tsai-formula X0,
# which can differ slightly from the PDG table values in X0_map above
# (e.g. G4 NIST X0 for Fe = 1.7588 cm vs PDG 1.757 cm, ~0.1% difference).
# The tolerance of 0.01 accommodates this systematic offset.
# For custom materials (LeadGlass, Brass), absorber_x0 will not appear
# in this CSV unless those materials are used as absorbers — the check
# will trivially pass (no rows) and should not be read as validation.
df['X0_check'] = df['thickness_cm'] / df['X0_cm']
discrepancy = (df['X0_check'] - df['absorber_x0']).abs().max()
print(f"\nSanity check — max discrepancy between recomputed and G4 absorber_x0: {discrepancy:.6f}")
if discrepancy < 0.01:
    print("  OK — within tolerance (G4 Tsai X0 vs PDG table values differ by ~0.1%)")
else:
    print("  FAIL — large discrepancy; check X0 values in the map above")

# Drop the check column (optional — remove this line to keep it)
df.drop(columns=['X0_check'], inplace=True)

# ── Save ───────────────────────────────────────────────────────────────────
df.to_csv(output_file, index=False)
print(f"\nSaved to {output_file}")
print(f"  New columns added: X0_cm, Z, density_g_cm3, Ec_MeV, moliere_radius_cm, thickness_cm")
print(f"\nColumn list:")
for col in df.columns:
    print(f"  {col}")
