#!/usr/bin/env python3
"""
analyse_shower.py — BL4S shower analysis
Usage: python analyse_shower.py shower_output.csv
       python analyse_shower.py shower_output.csv --outdir plots/

Produces 7 physics validation plots:
  1. Energy resolution vs beam energy (per material)
  2. Shower attenuation curves — Etotal vs absorber X0 (per material & energy)
  3. Core fraction vs absorber X0 (per material & energy)
  4. Shower width (2D radial RMS) vs absorber X0
  5. Block energy map — mean energy per calorimeter block
  6. Lateral shower profile — energy per column, per material
  7. Correlation matrix of all observables
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# ── Material properties ───────────────────────────────────────────────────────
X0_map  = {
    'G4_Al': 8.897, 'G4_Ti': 3.560, 'G4_Fe': 1.757,
    'G4_Cu': 1.436, 'G4_Sn': 1.210, 'G4_W':  0.350, 'G4_Pb': 0.5612
}
Z_map   = {
    'G4_Al': 13, 'G4_Ti': 22, 'G4_Fe': 26,
    'G4_Cu': 29, 'G4_Sn': 50, 'G4_W':  74, 'G4_Pb': 82
}
COLOURS = {
    'G4_Al': '#4C9BE8', 'G4_Ti': '#E84C4C', 'G4_Fe': '#E8834C',
    'G4_Cu': '#B05CE8', 'G4_Sn': '#4CE8B0', 'G4_W':  '#E8D84C', 'G4_Pb': '#4CE87A'
}
LABELS  = {
    'G4_Al': 'Aluminium (Z=13)', 'G4_Ti': 'Titanium (Z=22)',
    'G4_Fe': 'Iron (Z=26)',      'G4_Cu': 'Copper (Z=29)',
    'G4_Sn': 'Tin (Z=50)',       'G4_W':  'Tungsten (Z=74)',
    'G4_Pb': 'Lead (Z=82)'
}
BLOCK_COLS = [f'E_block_{i}_GeV' for i in range(16)]

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument('csv', nargs='?', default='shower_output.csv')
parser.add_argument('--outdir', default='plots')
args = parser.parse_args()
os.makedirs(args.outdir, exist_ok=True)

# ── Load ──────────────────────────────────────────────────────────────────────
print(f"Loading {args.csv} ...")
df = pd.read_csv(args.csv)
print(f"  {len(df):,} events")
print(f"  Materials : {sorted(df['material'].unique())}")
print(f"  Energies  : {sorted(df['beam_energy_GeV'].round(1).unique())} GeV")

if 'X0_cm' not in df.columns:
    df['X0_cm'] = df['material'].map(X0_map)
if 'Z' not in df.columns:
    df['Z'] = df['material'].map(Z_map)

materials = sorted(df['material'].unique(), key=lambda m: Z_map.get(m, 0))
energies  = sorted(df['beam_energy_GeV'].round(1).unique())

def save(fig, name):
    path = os.path.join(args.outdir, name)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved {path}")

# ═════════════════════════════════════════════════════════════════════════════
# Plot 1 — Energy resolution vs beam energy
# ═════════════════════════════════════════════════════════════════════════════
print("\nPlot 1: Energy resolution ...")
fig, ax = plt.subplots(figsize=(8, 5))
df_thin = df[df['absorber_x0'] < 0.75]
for mat in materials:
    sub = df_thin[df_thin['material'] == mat]
    res, evals = [], []
    for E in energies:
        ev = sub[sub['beam_energy_GeV'].round(1) == round(E, 1)]['Etotal_GeV']
        if len(ev) < 50:
            continue
        mu = ev.mean()
        res.append(ev.std() / mu * 100 if mu > 0 else np.nan)
        evals.append(E)
    if evals:
        ax.plot(evals, res, 'o-', color=COLOURS[mat], label=LABELS[mat],
                linewidth=2, markersize=7)
E_ref = np.linspace(0.8, 4.5, 100)
ax.plot(E_ref, 6.3 / np.sqrt(E_ref), 'k--', linewidth=1.5,
        label='BL4S spec: 6.3%/√E', alpha=0.7)
ax.set_xlabel('Beam Energy (GeV)', fontsize=12)
ax.set_ylabel('Energy Resolution σ/E (%)', fontsize=12)
ax.set_title('Energy Resolution vs Beam Energy (absorber < 0.75 X0)', fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(0.5, 5)
ax.set_ylim(0)
save(fig, '1_energy_resolution.png')

# ═════════════════════════════════════════════════════════════════════════════
# Plot 2 — Shower attenuation: Etotal vs absorber X0
# ═════════════════════════════════════════════════════════════════════════════
print("Plot 2: Shower attenuation ...")
fig, axes = plt.subplots(1, len(energies), figsize=(5*len(energies), 5), sharey=True)
if len(energies) == 1: axes = [axes]
for ax, E in zip(axes, energies):
    for mat in materials:
        sub = df[(df['material'] == mat) &
                 (df['beam_energy_GeV'].round(1) == round(E, 1))]
        grp = sub.groupby('absorber_x0')['Etotal_GeV'].agg(['mean','std']).reset_index()
        grp = grp.sort_values('absorber_x0')
        ax.errorbar(grp['absorber_x0'], grp['mean'], yerr=grp['std'],
                    fmt='o-', color=COLOURS[mat], label=LABELS[mat],
                    linewidth=2, markersize=6, capsize=3)
    ax.set_xlabel('Absorber Thickness (X₀)', fontsize=11)
    ax.set_title(f'{E:.0f} GeV beam', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
axes[0].set_ylabel('Mean Total Energy in Calorimeter (GeV)', fontsize=11)
fig.suptitle('Shower Attenuation: Etotal vs Absorber Thickness', fontsize=13, y=1.02)
fig.tight_layout()
save(fig, '2_shower_attenuation.png')

# ═════════════════════════════════════════════════════════════════════════════
# Plot 3 — Core fraction vs absorber X0
# ═════════════════════════════════════════════════════════════════════════════
print("Plot 3: Core fraction ...")
fig, axes = plt.subplots(1, len(energies), figsize=(5*len(energies), 5), sharey=True)
if len(energies) == 1: axes = [axes]
for ax, E in zip(axes, energies):
    for mat in materials:
        sub = df[(df['material'] == mat) &
                 (df['beam_energy_GeV'].round(1) == round(E, 1))]
        grp = sub.groupby('absorber_x0')['CoreFraction'].agg(['mean','std']).reset_index()
        grp = grp.sort_values('absorber_x0')
        ax.errorbar(grp['absorber_x0'], grp['mean'], yerr=grp['std'],
                    fmt='s-', color=COLOURS[mat], label=LABELS[mat],
                    linewidth=2, markersize=6, capsize=3)
    ax.set_xlabel('Absorber Thickness (X₀)', fontsize=11)
    ax.set_title(f'{E:.0f} GeV beam', fontsize=12)
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
axes[0].set_ylabel('Mean Core Fraction (central 2×2 / total)', fontsize=11)
fig.suptitle('Core Fraction vs Absorber Thickness\n(higher Z → more compact shower)',
             fontsize=13, y=1.02)
fig.tight_layout()
save(fig, '3_core_fraction.png')

# ═════════════════════════════════════════════════════════════════════════════
# Plot 4 — Shower width vs absorber X0
# ═════════════════════════════════════════════════════════════════════════════
print("Plot 4: Shower width ...")
fig, axes = plt.subplots(1, len(energies), figsize=(5*len(energies), 5), sharey=True)
if len(energies) == 1: axes = [axes]
for ax, E in zip(axes, energies):
    for mat in materials:
        sub = df[(df['material'] == mat) &
                 (df['beam_energy_GeV'].round(1) == round(E, 1))]
        grp = sub.groupby('absorber_x0')['ShowerWidth_cm'].agg(['mean','std']).reset_index()
        grp = grp.sort_values('absorber_x0')
        ax.errorbar(grp['absorber_x0'], grp['mean'], yerr=grp['std'],
                    fmt='^-', color=COLOURS[mat], label=LABELS[mat],
                    linewidth=2, markersize=6, capsize=3)
    ax.set_xlabel('Absorber Thickness (X₀)', fontsize=11)
    ax.set_title(f'{E:.0f} GeV beam', fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
axes[0].set_ylabel('Mean Shower Width — 2D Radial RMS (cm)', fontsize=11)
fig.suptitle('Lateral Shower Width vs Absorber Thickness\n(proxy for Molière radius)',
             fontsize=13, y=1.02)
fig.tight_layout()
save(fig, '4_shower_width.png')

# ═════════════════════════════════════════════════════════════════════════════
# Plot 5 — Block energy heatmap (4×4)
# ═════════════════════════════════════════════════════════════════════════════
print("Plot 5: Block energy maps ...")
ncols = 2
nrows = (len(materials) + 1) // 2
fig, axes = plt.subplots(nrows, ncols, figsize=(10, 5*nrows))
axes = axes.flatten()
tick_labels = ['-15cm', '-5cm', '+5cm', '+15cm']
for ax, mat in zip(axes, materials):
    sub = df[(df['material'] == mat) & (df["absorber_x0"] < df["absorber_x0"].quantile(0.25))]
    if len(sub) == 0:
        sub = df[df['material'] == mat]
    mean_E = sub[BLOCK_COLS].mean().values.reshape(4, 4)
    im = ax.imshow(mean_E, cmap='hot', origin='lower', aspect='equal')
    plt.colorbar(im, ax=ax, label='Mean Energy (GeV)')
    ax.set_title(f'{LABELS[mat]}\n(absorber < 25th percentile X0)', fontsize=11)
    ax.set_xlabel('Column (X)', fontsize=10)
    ax.set_ylabel('Row (Y)', fontsize=10)
    ax.set_xticks([0,1,2,3]); ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_yticks([0,1,2,3]); ax.set_yticklabels(tick_labels, fontsize=8)
    vmax = mean_E.max()
    for r in range(4):
        for c in range(4):
            colour = 'white' if mean_E[r,c] > 0.5*vmax else 'black'
            ax.text(c, r, f'{mean_E[r,c]:.3f}', ha='center', va='center',
                    fontsize=7, color=colour)
for ax in axes[len(materials):]:
    ax.set_visible(False)
fig.suptitle('Mean Energy per Calorimeter Block', fontsize=13)
fig.tight_layout()
save(fig, '5_block_energy_map.png')

# ═════════════════════════════════════════════════════════════════════════════
# Plot 6 — Lateral profile per material
# ═════════════════════════════════════════════════════════════════════════════
print("Plot 6: Lateral profiles ...")
fig, axes = plt.subplots(1, len(energies), figsize=(5*len(energies), 5), sharey=True)
if len(energies) == 1: axes = [axes]
col_centres = np.array([-15., -5., 5., 15.])
for ax, E in zip(axes, energies):
    for mat in materials:
        sub = df[(df['material'] == mat) &
                 (df['beam_energy_GeV'].round(1) == round(E, 1)) &
                 (df["absorber_x0"] < df["absorber_x0"].quantile(0.25))]
        if len(sub) == 0:
            continue
        col_E = np.array([
            sub[[f'E_block_{row*4+col}_GeV' for row in range(4)]].values.sum()
            for col in range(4)
        ])
        # Note: normalisation is by column sums (4-column total), not 16-block total.
        # For well-contained showers these are identical; edge events may differ slightly.
        total = col_E.sum()
        if total > 0:
            col_E /= total
        ax.plot(col_centres, col_E, 'o-', color=COLOURS[mat],
                label=LABELS[mat], linewidth=2, markersize=7)
    ax.set_xlabel('Block Column Centre (cm)', fontsize=11)
    ax.set_title(f'{E:.0f} GeV beam', fontsize=12)
    ax.set_xticks(col_centres)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
axes[0].set_ylabel('Normalised Energy Fraction per Column', fontsize=11)
fig.suptitle('Lateral Shower Profile per Material\n(absorber < 25th percentile X0 — narrower = higher Z / smaller Moliere radius)',
             fontsize=13, y=1.02)
fig.tight_layout()
save(fig, '6_lateral_profile.png')

# ═════════════════════════════════════════════════════════════════════════════
# Plot 7 — Correlation matrix
# ═════════════════════════════════════════════════════════════════════════════
print("Plot 7: Correlation matrix ...")
obs_cols = BLOCK_COLS + ['Etotal_GeV', 'CoreFraction', 'ShowerWidth_cm',
                          'beam_energy_GeV', 'absorber_x0']
sample = df[obs_cols].sample(min(50000, len(df)), random_state=42)
corr = sample.corr()
rename = {f'E_block_{i}_GeV': f'B{i}' for i in range(16)}
rename.update({'Etotal_GeV': 'Etotal', 'CoreFraction': 'Core',
               'ShowerWidth_cm': 'Width', 'beam_energy_GeV': 'BeamE',
               'absorber_x0': 'AbsX0'})
corr = corr.rename(index=rename, columns=rename)
fig, ax = plt.subplots(figsize=(14, 12))
mask = np.triu(np.ones_like(corr, dtype=bool))
sns.heatmap(corr, mask=mask, ax=ax, cmap='RdBu_r', center=0,
            vmin=-1, vmax=1, square=True, linewidths=0.3,
            annot=False, cbar_kws={'label': 'Pearson r'})
ax.set_title('Observable Correlation Matrix\n(B0–B15 = calorimeter blocks)',
             fontsize=13)
fig.tight_layout()
save(fig, '7_correlation_matrix.png')

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'─'*55}")
print(f"  All 7 plots saved to: {os.path.abspath(args.outdir)}/")
print(f"{'─'*55}")
print("\nQuick data summary:")
print(df.groupby('material')[['Etotal_GeV','CoreFraction','ShowerWidth_cm']].mean().round(3))
