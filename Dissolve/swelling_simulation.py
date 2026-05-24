#!/usr/bin/env python3
"""
In-Situ Ellipsometry Swelling Simulation
=========================================
Back-calculate solvent vapor swelling curves for PM6 and D18 in CF, Toluene, o-Xylene.

Physics:
  Binary Flory-Huggins swelling equilibrium:
    ln(a) = ln(phi_s) + (1 - phi_s) + chi*(1 - phi_s)^2

  where:
    a = p/p_sat (solvent activity)
    phi_s = solvent volume fraction in swollen film
    chi = Flory-Huggins polymer-solvent interaction parameter (= g23 from phase diagram)

  Swelling ratio: Q = h/h0 = 1/(1 - phi_s)
  Volume expansion: (Q - 1) * 100%

Data stored in Three-phase/Dissolve/
"""

import numpy as np
from scipy.optimize import fsolve
import os, csv

# ============================================================
# OUTPUT DIR
# ============================================================
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ============================================================
# POLYMER-SOLVENT CHI PARAMETERS (from ternary phase diagram g23)
# ============================================================
# chi = g23 (binary polymer-solvent interaction)
SYSTEMS = {
    'PM6': {
        'Chloroform': 0.22,    # estimated from HSP
        'Toluene':    0.385,   # from phase diagram optimization
        'o-Xylene':   0.412,   # from phase diagram optimization
        'color': '#4f46e5',
    },
    'D18': {
        'Chloroform': 0.28,    # estimated from HSP
        'Toluene':    0.450,   # from phase diagram optimization
        'o-Xylene':   0.300,   # from phase diagram optimization
        'color': '#f59e0b',
    },
}

# Solvent properties
SOLVENT_INFO = {
    'Chloroform': {'bp': 61.2, 'V_molar': 80.7,  'p_sat_25C': 26.3},   # kPa at 25C
    'Toluene':    {'bp': 110.6,'V_molar': 106.3, 'p_sat_25C': 3.8},
    'o-Xylene':   {'bp': 144.4,'V_molar': 120.6, 'p_sat_25C': 0.9},
}


def swelling_eq(phi_s, a, chi):
    """FH swelling equilibrium: ln(a) = ln(phi_s) + (1-phi_s) + chi*(1-phi_s)^2
    Return residual (should be 0 at equilibrium)"""
    if phi_s <= 1e-10 or phi_s >= 1.0:
        return 1e10
    lhs = np.log(max(a, 1e-10))
    rhs = np.log(phi_s) + (1 - phi_s) + chi * (1 - phi_s)**2
    return lhs - rhs


def solve_swelling(a, chi):
    """Solve for equilibrium phi_s at given activity."""
    # Bracket: phi_s is between 0 and 1
    # For poor solvents (chi > 0.5), swelling is limited
    # For good solvents (chi < 0.5), swelling can be large

    # Initial guess based on chi
    if chi < 0.3:
        phi0 = 0.5  # good solvent, significant swelling
    elif chi < 0.5:
        phi0 = 0.3  # moderate solvent
    else:
        phi0 = 0.15  # poor solvent, limited swelling

    try:
        sol = fsolve(lambda phi: swelling_eq(phi, a, chi), phi0, maxfev=1000)
        phi_s = float(sol[0])
        if 0 < phi_s < 0.999:
            return phi_s
    except:
        pass

    # Try different initial guesses
    for guess in [0.05, 0.1, 0.2, 0.4, 0.6, 0.8]:
        try:
            sol = fsolve(lambda phi: swelling_eq(phi, a, chi), guess, maxfev=1000)
            phi_s = float(sol[0])
            if 0 < phi_s < 0.999:
                return phi_s
        except:
            continue

    return None


def compute_swelling_curve(chi, n_pts=50):
    """Compute swelling curve: activity -> phi_s, Q, expansion%"""
    # Activities from very dilute vapor to near saturation
    activities = np.logspace(-3, -0.02, n_pts)  # 0.001 to ~0.95

    results = []
    for a in activities:
        phi_s = solve_swelling(a, chi)
        if phi_s is not None:
            Q = 1.0 / (1.0 - phi_s)  # swelling ratio h/h0
            expansion = (Q - 1.0) * 100  # volume expansion %
            results.append({
                'activity': a,
                'p_p_sat': a,
                'phi_solvent': phi_s,
                'phi_polymer': 1.0 - phi_s,
                'swelling_ratio_Q': Q,
                'volume_expansion_pct': expansion,
            })

    return results


def compute_thickness_evolution(chi, initial_thickness_nm=100, n_pts=50):
    """Model film thickness as a function of solvent vapor activity."""
    curve = compute_swelling_curve(chi, n_pts)
    for pt in curve:
        pt['thickness_nm'] = initial_thickness_nm * pt['swelling_ratio_Q']
        pt['delta_thickness_nm'] = pt['thickness_nm'] - initial_thickness_nm
    return curve


# ============================================================
# MAIN COMPUTATION
# ============================================================
def main():
    print("=" * 70)
    print("  In-Situ Ellipsometry Swelling Simulation")
    print("  Flory-Huggins binary swelling model")
    print("=" * 70)

    for polymer, solvents in SYSTEMS.items():
        for solvent, chi in solvents.items():
            if solvent == 'color':
                continue

            print(f"\n  {polymer} / {solvent}: chi = {chi:.3f}")

            # Full swelling curve
            curve = compute_swelling_curve(chi, n_pts=80)
            if len(curve) < 5:
                print(f"    WARNING: Only {len(curve)} points converged")
                continue

            # Key metrics
            # Activity at which phi_s = 0.1 (10% solvent)
            phi10_act = None
            phi30_act = None
            max_swelling = curve[-1]['swelling_ratio_Q']
            max_phi = curve[-1]['phi_solvent']

            for pt in curve:
                if phi10_act is None and pt['phi_solvent'] >= 0.10:
                    phi10_act = pt['activity']
                if phi30_act is None and pt['phi_solvent'] >= 0.30:
                    phi30_act = pt['activity']

            print(f"    Max swelling (a=0.95): Q={max_swelling:.3f}, phi_s={max_phi:.3f}")
            if phi10_act:
                print(f"    Activity for 10% solvent uptake: a={phi10_act:.4f}")
            if phi30_act:
                print(f"    Activity for 30% solvent uptake: a={phi30_act:.4f}")

            # Save CSV
            csv_path = os.path.join(OUT_DIR, f'swelling_{polymer}_{solvent.replace(" ","_")}.csv')
            with open(csv_path, 'w', newline='') as f:
                w = csv.DictWriter(f, fieldnames=[
                    'activity', 'p_p_sat', 'phi_solvent', 'phi_polymer',
                    'swelling_ratio_Q', 'volume_expansion_pct'
                ])
                w.writeheader()
                w.writerows(curve)
            print(f"    Saved: {csv_path}")

            # Thickness evolution
            thickness = compute_thickness_evolution(chi, initial_thickness_nm=100, n_pts=80)
            csv_path_t = os.path.join(OUT_DIR, f'thickness_{polymer}_{solvent.replace(" ","_")}.csv')
            with open(csv_path_t, 'w', newline='') as f:
                w = csv.DictWriter(f, fieldnames=[
                    'activity', 'p_p_sat', 'phi_solvent', 'phi_polymer',
                    'swelling_ratio_Q', 'volume_expansion_pct',
                    'thickness_nm', 'delta_thickness_nm'
                ])
                w.writeheader()
                w.writerows(thickness)
            print(f"    Saved: {csv_path_t}")

    # ============================================================
    # SUMMARY TABLE
    # ============================================================
    print(f"\n{'='*70}")
    print(f"  SUMMARY: Swelling at a=0.95 (near saturation)")
    print(f"{'='*70}")
    print(f"  {'System':<25s} {'chi':>8s} {'phi_s':>8s} {'Q':>8s} {'Expansion%':>10s}")
    print(f"  {'-'*25} {'-'*8} {'-'*8} {'-'*8} {'-'*10}")

    summary_data = []
    for polymer in ['PM6', 'D18']:
        for solvent in ['Chloroform', 'Toluene', 'o-Xylene']:
            chi = SYSTEMS[polymer][solvent]
            curve = compute_swelling_curve(chi, n_pts=80)
            if len(curve) >= 5:
                pt = curve[-1]
                print(f"  {polymer}/{solvent:<18s} {chi:>8.3f} {pt['phi_solvent']:>8.3f} "
                      f"{pt['swelling_ratio_Q']:>8.3f} {pt['volume_expansion_pct']:>10.1f}%")
                summary_data.append({
                    'polymer': polymer, 'solvent': solvent, 'chi': chi,
                    'phi_s': pt['phi_solvent'], 'Q': pt['swelling_ratio_Q'],
                    'expansion_pct': pt['volume_expansion_pct'],
                })

    # Save summary
    summary_path = os.path.join(OUT_DIR, 'swelling_summary.csv')
    with open(summary_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['polymer','solvent','chi','phi_s','Q','expansion_pct'])
        w.writeheader()
        w.writerows(summary_data)
    print(f"\n  Summary saved: {summary_path}")

    return True


if __name__ == '__main__':
    main()
