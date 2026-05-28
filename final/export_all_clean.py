#!/usr/bin/env python3
"""
Export ALL ternary data with concentrated/dilute binodal SEPARATED.
Also regenerate Python ternary plots by solvent.
"""
import os, csv, shutil

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'output', 'data')
OUT_DIR = os.path.join(os.path.dirname(__file__), 'clean_data')
os.makedirs(OUT_DIR, exist_ok=True)

SYSTEMS = [
    'PM6_Tol', 'PM6_OXy', 'D18_Tol', 'D18_OXy', 'PM6_CF', 'D18_CF'
]

for tag in SYSTEMS:
    # ---- Binodal: split concentrated / dilute ----
    bin_file = os.path.join(DATA_DIR, f'binodal_{tag}.csv')
    conc_rows = []
    dil_rows = []
    if os.path.exists(bin_file):
        with open(bin_file, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            rows = list(reader)
            for i in range(0, len(rows)-1, 2):
                conc_rows.append(rows[i])
                dil_rows.append(rows[i+1])

        # Save concentrated phase
        with open(os.path.join(OUT_DIR, f'binodal_{tag}_concentrated.csv'), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['residual','phi1_L8Bo','phi2_Solvent','phi3_Donor'])
            for r in conc_rows:
                if len(r) >= 4:
                    w.writerow([r[0], r[1], r[3], r[2]])

        # Save dilute phase
        with open(os.path.join(OUT_DIR, f'binodal_{tag}_dilute.csv'), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['residual','phi1_L8Bo','phi2_Solvent','phi3_Donor'])
            for r in dil_rows:
                if len(r) >= 4:
                    w.writerow([r[0], r[1], r[3], r[2]])

    # ---- Spinodal: just copy, add solvent column header ----
    spin_file = os.path.join(DATA_DIR, f'spinodal_{tag}.csv')
    if os.path.exists(spin_file):
        with open(spin_file, 'r') as f:
            content = f.read()
        with open(os.path.join(OUT_DIR, f'spinodal_{tag}.csv'), 'w') as f:
            f.write(content)

    # ---- Critical: just copy ----
    crit_file = os.path.join(DATA_DIR, f'critical_{tag}.csv')
    if os.path.exists(crit_file):
        shutil.copy(crit_file, os.path.join(OUT_DIR, f'critical_{tag}.csv'))

    n_conc = len(conc_rows)
    n_dil = len(dil_rows)
    print(f"  {tag}: conc={n_conc}pts  dil={n_dil}pts  exported")

# ---- Summary index ----
with open(os.path.join(OUT_DIR, '_FILE_INDEX.txt'), 'w') as f:
    f.write("Concentrated / Dilute Binodal Data — Separated\n")
    f.write("================================================\n\n")
    f.write("For each system, THREE files:\n")
    f.write("  binodal_{system}_concentrated.csv  — concentrated phase (donor-rich)\n")
    f.write("  binodal_{system}_dilute.csv        — dilute phase (acceptor-rich)\n")
    f.write("  spinodal_{system}.csv              — spinodal curve\n")
    f.write("  critical_{system}.csv              — critical point\n\n")
    f.write("Columns in binodal files:\n")
    f.write("  phi1_L8Bo = Acceptor volume fraction\n")
    f.write("  phi2_Solvent = Solvent volume fraction\n")
    f.write("  phi3_Donor = Donor volume fraction\n\n")
    f.write("In Origin ternary plot: X=phi1, Y=phi3, Z=phi2\n")

print(f"\nAll clean data: {OUT_DIR}")
