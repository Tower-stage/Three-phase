#!/usr/bin/env python3
"""Plot swelling curves from computed CSV data."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv, os

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

COLORS = {
    'PM6': '#4f46e5',
    'D18': '#f59e0b',
    'Chloroform': '#8b5cf6',
    'Toluene': '#10b981',
    'o-Xylene': '#0891b2',
}
LINE_STYLES = {
    'Chloroform': '-',
    'Toluene': '--',
    'o-Xylene': '-.',
}
SOLVENT_CHI = {
    'PM6': {'Chloroform': 0.22, 'Toluene': 0.385, 'o-Xylene': 0.412},
    'D18': {'Chloroform': 0.28, 'Toluene': 0.450, 'o-Xylene': 0.300},
}

def load_curve(filename):
    activities, phi_s, Q, expansion = [], [], [], []
    with open(os.path.join(OUT_DIR, filename), 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            activities.append(float(row['activity']))
            phi_s.append(float(row['phi_solvent']))
            Q.append(float(row['swelling_ratio_Q']))
            expansion.append(float(row['volume_expansion_pct']))
    return activities, phi_s, Q, expansion


# ============================================================
# FIGURE 1: Swelling ratio Q vs activity — all 6 systems
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
fig.suptitle('In-Situ Ellipsometry Swelling Simulation (Flory-Huggins)',
             fontsize=15, fontweight='bold')

for idx, polymer in enumerate(['PM6', 'D18']):
    ax = axes[idx]
    for solvent in ['Chloroform', 'Toluene', 'o-Xylene']:
        fname = f'swelling_{polymer}_{solvent}.csv'
        acts, phi_s, Q, expn = load_curve(fname)
        ax.plot(acts, Q, color=COLORS[solvent], linestyle=LINE_STYLES[solvent],
                linewidth=2.5, label=f'{solvent} (chi={SOLVENT_CHI[polymer][solvent]:.3f})')

    ax.set_xlabel('Solvent Activity a = p/p_sat', fontsize=13)
    ax.set_ylabel('Swelling Ratio Q = h/h0', fontsize=13)
    ax.set_title(f'{polymer} — Solvent Vapor Swelling', fontsize=14, fontweight='bold',
                 color=COLORS[polymer])
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.2)
    ax.set_xlim(0, 1)
    ax.set_ylim(1, None)

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'fig_swelling_curves.png'), dpi=200,
            bbox_inches='tight', facecolor='white')
print("Saved: fig_swelling_curves.png")
plt.close(fig)


# ============================================================
# FIGURE 2: Solvent volume fraction phi_s vs activity
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
fig.suptitle('Solvent Volume Fraction in Swollen Film',
             fontsize=15, fontweight='bold')

for idx, polymer in enumerate(['PM6', 'D18']):
    ax = axes[idx]
    for solvent in ['Chloroform', 'Toluene', 'o-Xylene']:
        fname = f'swelling_{polymer}_{solvent}.csv'
        acts, phi_s, Q, expn = load_curve(fname)
        ax.plot(acts, phi_s, color=COLORS[solvent], linestyle=LINE_STYLES[solvent],
                linewidth=2.5, label=f'{solvent} (chi={SOLVENT_CHI[polymer][solvent]:.3f})')

    ax.set_xlabel('Solvent Activity a = p/p_sat', fontsize=13)
    ax.set_ylabel('Solvent Volume Fraction phi_s', fontsize=13)
    ax.set_title(f'{polymer} — phi_s vs Activity', fontsize=14, fontweight='bold',
                 color=COLORS[polymer])
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.2)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, None)

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'fig_phi_solvent.png'), dpi=200,
            bbox_inches='tight', facecolor='white')
print("Saved: fig_phi_solvent.png")
plt.close(fig)


# ============================================================
# FIGURE 3: Thickness evolution for PM6 and D18
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
fig.suptitle('Film Thickness Evolution During Solvent Vapor Swelling (h0=100nm)',
             fontsize=15, fontweight='bold')

for idx, polymer in enumerate(['PM6', 'D18']):
    ax = axes[idx]
    for solvent in ['Chloroform', 'Toluene', 'o-Xylene']:
        fname = f'thickness_{polymer}_{solvent}.csv'
        acts, phi_s, Q, expn = load_curve(fname)  # same columns, Q=swelling_ratio
        # reload for thickness
        thicknesses = []
        activities = []
        with open(os.path.join(OUT_DIR, fname), 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                activities.append(float(row['activity']))
                thicknesses.append(float(row['thickness_nm']))

        ax.plot(activities, thicknesses, color=COLORS[solvent],
                linestyle=LINE_STYLES[solvent], linewidth=2.5,
                label=f'{solvent}')

    ax.axhline(y=100, color='gray', linestyle=':', alpha=0.5, label='Dry film (100nm)')
    ax.set_xlabel('Solvent Activity a = p/p_sat', fontsize=13)
    ax.set_ylabel('Film Thickness (nm)', fontsize=13)
    ax.set_title(f'{polymer} — Thickness vs Activity', fontsize=14, fontweight='bold',
                 color=COLORS[polymer])
    ax.legend(fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.2)
    ax.set_xlim(0, 1)

plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'fig_thickness.png'), dpi=200,
            bbox_inches='tight', facecolor='white')
print("Saved: fig_thickness.png")
plt.close(fig)


# ============================================================
# FIGURE 4: chi vs Max Swelling — cross-system comparison
# ============================================================
fig, ax = plt.subplots(figsize=(10, 7))

summary_path = os.path.join(OUT_DIR, 'swelling_summary.csv')
polymers, solvents, chis, Qs, expansions = [], [], [], [], []
with open(summary_path, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        polymers.append(row['polymer'])
        solvents.append(row['solvent'])
        chis.append(float(row['chi']))
        Qs.append(float(row['Q']))
        expansions.append(float(row['expansion_pct']))

for i in range(len(polymers)):
    c = COLORS[solvents[i]]
    marker = 'o' if polymers[i] == 'PM6' else 's'
    ax.scatter(chis[i], expansions[i], s=200, c=c, edgecolors='black',
               linewidth=1.5, marker=marker, zorder=5)
    ax.annotate(f'{polymers[i]}/{solvents[i]}',
                (chis[i], expansions[i]),
                textcoords="offset points", xytext=(0, 12),
                fontsize=9, ha='center', fontweight='bold', color=c)

# Trend line
z = np.polyfit(chis, expansions, 1)
x_fit = np.linspace(min(chis)-0.02, max(chis)+0.02, 100)
y_fit = np.polyval(z, x_fit)
ax.plot(x_fit, y_fit, '--', color='gray', alpha=0.5, linewidth=1.5,
        label=f'Linear trend (R² ≈ ...)')

ax.set_xlabel('Flory-Huggins chi (polymer-solvent)', fontsize=13)
ax.set_ylabel('Volume Expansion at a=0.95 (%)', fontsize=13)
ax.set_title('chi vs Swelling: Lower chi = More Swelling',
             fontsize=15, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.2)

fig.savefig(os.path.join(OUT_DIR, 'fig_chi_vs_swelling.png'), dpi=200,
            bbox_inches='tight', facecolor='white')
print("Saved: fig_chi_vs_swelling.png")
plt.close(fig)

print("\nDone! All figures saved to:", OUT_DIR)
