#!/usr/bin/env python3
"""
Ternary phase diagrams grouped by SOLVENT.
Each figure: all donor systems in one solvent.
Format: Left=Donor, Bottom=Acceptor(L8-Bo), Right=Solvent.
Colors: #F3A569 (binodal), #4190C4 (tie-lines), blue frame.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import os, csv

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'output', 'data')
OUT_DIR = os.path.dirname(__file__)
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================
# Color scheme
# ============================================================
BINODAL_COLOR = '#F3A569'   # orange
TIE_COLOR      = '#4190C4'   # blue
FRAME_COLOR    = '#4190C4'   # blue frame
SPINODAL_COLOR = '#999999'   # gray dashed
CRITICAL_COLOR = '#F3A569'
DONOR_COLORS = {
    'PM6': '#4472C4',
    'D18': '#ED7D31',
}

# Group by solvent
GROUPS = {
    'Toluene': {
        'systems': ['PM6_Tol', 'D18_Tol'],
        'v2': 106.3,
    },
    'o-Xylene': {
        'systems': ['PM6_OXy', 'D18_OXy'],
        'v2': 120.6,
    },
    'Chloroform': {
        'systems': ['PM6_CF', 'D18_CF'],
        'v2': 80.7,
    },
}

# ============================================================
# Ternary → Cartesian
# ============================================================
SQ3 = np.sqrt(3)
# Vertex positions:
#  Bottom-left  (0,0)       = Donor (phi3)
#  Bottom-right (1,0)       = Acceptor = L8-Bo (phi1)
#  Top          (0.5, SQ3/2) = Solvent (phi2)

def to_xy(phi1, phi2):
    """phi1=acceptor(L8-Bo), phi2=solvent. phi3=1-phi1-phi2=donor."""
    x = phi1 * 1.0 + phi2 * 0.5
    y = phi2 * SQ3 / 2
    return x, y

# ============================================================
# Draw frame
# ============================================================
def draw_frame(ax):
    # Triangle
    corners = [(0,0), (1,0), (0.5, SQ3/2)]
    for i in range(3):
        x1,y1 = corners[i]; x2,y2 = corners[(i+1)%3]
        ax.plot([x1,x2],[y1,y2], '-', color=FRAME_COLOR, lw=2.5, solid_capstyle='round')

    # Grid lines every 25%
    for level in [0.25, 0.50, 0.75]:
        y_g = level * SQ3 / 2
        # Horizontal (constant phi2=solvent)
        ax.plot([level/2, 1-level/2], [y_g, y_g], '-', color='#ddd', lw=0.5, alpha=0.7)
        # Parallel to right edge (constant phi3=donor)
        ax.plot([level, level/2], [0, y_g], '-', color='#ddd', lw=0.5, alpha=0.7)
        # Parallel to left edge (constant phi1=acceptor)
        ax.plot([1-level, 1-level/2], [0, y_g], '-', color='#ddd', lw=0.5, alpha=0.7)

    # Tick labels — percentages, parallel to axes
    for pct in [25, 50, 75]:
        val = pct / 100.0
        # Bottom axis (phi_acceptor): along bottom edge, reading rightward
        ax.text(val, -0.03, f'{pct}', ha='center', va='top', fontsize=8,
                color='#666', rotation=0)
        # Also the complementary for donor (reading leftward)
        ax.text(1-val, -0.03, f'{pct}', ha='center', va='top', fontsize=8,
                color='#999', rotation=0)
        # Left edge (phi_donor): along left edge
        xi = val/2; yi = val * SQ3 / 2
        ax.text(xi - 0.04, yi + 0.01, f'{pct}', ha='right', va='center',
                fontsize=8, color='#666', rotation=60)
        # Right edge (phi_solvent): along right edge
        xi = 1 - val/2; yi = val * SQ3 / 2
        ax.text(xi + 0.04, yi + 0.01, f'{pct}', ha='left', va='center',
                fontsize=8, color='#666', rotation=-60)

    # Vertex labels
    ax.text(-0.04, -0.06, r'$\phi_{\mathrm{Donor}}$', ha='right', va='top',
            fontsize=11, fontweight='bold', color='#333')
    ax.text(1.04, -0.06, r'$\phi_{\mathrm{Acceptor}}$', ha='left', va='top',
            fontsize=11, fontweight='bold', color='#333')
    ax.text(0.5, SQ3/2 + 0.05, r'$\phi_{\mathrm{Solvent}}$', ha='center', va='bottom',
            fontsize=11, fontweight='bold', color='#333')

    # Axis labels
    ax.text(-0.06, SQ3/4, r'$\phi_{Donor}(\%)$', ha='right', va='center',
            fontsize=9, color='#666', rotation=60)
    ax.text(1.06, SQ3/4, r'$\phi_{Solvent}(\%)$', ha='left', va='center',
            fontsize=9, color='#666', rotation=-60)
    ax.text(0.5, -0.08, r'$\phi_{Acceptor}(\%)$', ha='center', va='top',
            fontsize=9, color='#666')


# ============================================================
# Load data
# ============================================================
def load(system_tag):
    """Load binodal (concentrated only), spinodal, and critical point."""
    tag = system_tag

    # Binodal
    bino = []
    bfile = os.path.join(DATA_DIR, f'binodal_{tag}.csv')
    if os.path.exists(bfile):
        with open(bfile) as f:
            reader = csv.reader(f); next(reader)
            rows = list(reader)
            for i in range(0, len(rows)-1, 2):  # concentrated phase only
                r = rows[i]
                if len(r) >= 4:
                    phi1, phi3, phi2 = float(r[1]), float(r[2]), float(r[3])
                    bino.append((phi1, phi2, phi3))
    bino.sort(key=lambda p: p[2], reverse=True)  # sort by donor

    # Spinodal
    spin = []
    sfile = os.path.join(DATA_DIR, f'spinodal_{tag}.csv')
    if os.path.exists(sfile):
        with open(sfile) as f:
            reader = csv.reader(f); next(reader)
            for r in reader:
                if len(r) >= 4:
                    phi1, phi3, phi2 = float(r[1]), float(r[2]), float(r[3])
                    if 0<=phi1<=1 and 0<=phi2<=1 and 0<=phi3<=1:
                        spin.append((phi1, phi2, phi3))
    spin.sort(key=lambda p: (p[1], p[0]))

    # Critical
    crit = None
    cfile = os.path.join(DATA_DIR, f'critical_{tag}.csv')
    if os.path.exists(cfile):
        with open(cfile) as f:
            reader = csv.reader(f); next(reader)
            r = next(reader)
            crit = (float(r[0]), float(r[1]), float(r[2]))

    return bino, spin, crit


# ============================================================
# PLOT ONE SOLVENT FIGURE
# ============================================================
def plot_solvent(solvent_name, sys_tags, donor_names):
    fig, ax = plt.subplots(figsize=(8.5, 8))
    ax.set_aspect('equal')
    ax.set_xlim(-0.12, 1.12)
    ax.set_ylim(-0.14, SQ3/2 + 0.12)
    ax.axis('off')

    draw_frame(ax)

    for sys_tag, donor_name in zip(sys_tags, donor_names):
        bino, spin, crit = load(sys_tag)
        d_color = DONOR_COLORS[donor_name]

        # Spinodal (dashed, same color as donor but muted)
        if len(spin) > 2:
            xs, ys = zip(*[to_xy(p[0], p[1]) for p in spin])
            ax.plot(xs, ys, '--', color=d_color, lw=0.8, alpha=0.4, zorder=3)

        # Binodal (F3A569 orange)
        if len(bino) > 2:
            xs, ys = zip(*[to_xy(p[0], p[1]) for p in bino])
            ax.plot(xs, ys, '-', color=BINODAL_COLOR, lw=2.8, zorder=10,
                    label=f'{donor_name} binodal')

        # Tie-lines (4190C4 blue) — from full binodal data
        bfile = os.path.join(DATA_DIR, f'binodal_{sys_tag}.csv')
        tls = []
        if os.path.exists(bfile):
            with open(bfile) as f:
                reader = csv.reader(f); next(reader)
                rows = list(reader)
                for i in range(0, len(rows)-1, 2):
                    if len(rows[i])>=4 and len(rows[i+1])>=4:
                        c = (float(rows[i][1]), float(rows[i][3]), float(rows[i][2]))
                        d = (float(rows[i+1][1]), float(rows[i+1][3]), float(rows[i+1][2]))
                        tls.append((c, d))

        step = max(1, len(tls)//12)
        for i, (c, d) in enumerate(tls):
            if i % step != 0: continue
            x1, y1 = to_xy(c[0], c[1])
            x2, y2 = to_xy(d[0], d[1])
            ax.plot([x1, x2], [y1, y2], '-', color=TIE_COLOR, lw=0.7,
                    alpha=0.5, zorder=5)

        # Critical point
        if crit:
            xc, yc = to_xy(crit[0], crit[1])
            ax.scatter(xc, yc, marker='D', s=100, c=CRITICAL_COLOR,
                       edgecolors='#333', lw=1.5, zorder=20)
            ax.annotate(f'{donor_name}\n$\phi_2$={crit[1]:.2f}',
                       (xc, yc), textcoords="offset points", xytext=(6, 6),
                       fontsize=8, color=d_color, fontweight='bold')

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0],[0], color=BINODAL_COLOR, lw=2.8, label='Binodal'),
        Line2D([0],[0], color=TIE_COLOR, lw=0.8, alpha=0.6, label='Tie-lines'),
        Line2D([0],[0], color='#999', lw=0.8, ls='--', alpha=0.4, label='Spinodal'),
        Line2D([0],[0], marker='D', color='w', markerfacecolor=CRITICAL_COLOR,
               markersize=8, markeredgecolor='#333', markeredgewidth=1.5,
               label='Critical Point'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9,
              framealpha=0.9, ncol=2)

    ax.set_title(f'{solvent_name} — Ternary Phase Diagrams',
                 fontsize=14, fontweight='bold', color='#333', pad=15)

    path = os.path.join(OUT_DIR, f'ternary_{solvent_name.replace(" ","_")}.png')
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"Saved: {path}")
    plt.close(fig)


# ============================================================
# MAIN
# ============================================================
for solvent, cfg in GROUPS.items():
    donors = ['PM6' if 'PM6' in s else 'D18' for s in cfg['systems']]
    print(f"Plotting {solvent}...")
    plot_solvent(solvent, cfg['systems'], donors)

print("\nDone!")
