#!/usr/bin/env python3
"""Compute final CF ternary diagrams with optimized parameters."""
import numpy as np
from scipy.optimize import fsolve
import os, csv

from optimize_cf import chempot_eq, compute_spinodal_cp, compute_binodal

v1 = 1132.1; v2_cf = 80.7
OUT = os.path.join(os.path.dirname(__file__), '..', 'output', 'data')
os.makedirs(OUT, exist_ok=True)

# Optimized parameters (from sweep)
CF_PARAMS = {
    'PM6_CF':  {'v3':100000, 'X13':0.90, 'g23':0.18, 'p':[0,0,0,-0.12,0.40]},
    'D18_CF':  {'v3': 60000, 'X13':0.90, 'g23':0.18, 'p':[0,0,0,-0.12,0.40]},
}

for tag, cfg in CF_PARAMS.items():
    donor = 'PM6' if 'PM6' in tag else 'D18'
    v3=cfg['v3']; X13=cfg['X13']; g23=cfg['g23']; p=cfg['p']
    print(f"\n{'='*60}\n  {donor}/L8-Bo/CF  X13={X13} g23={g23} v3={v3} p={p}\n{'='*60}")

    # Full-res spinodal + CP
    segs,x1c,x3c = compute_spinodal_cp(v1,v2_cf,v3,X13,g23,p,n_g=1200)
    x2c=1-x1c-x3c
    print(f"  CP: phi1={x1c:.4f} phi2={x2c:.4f} phi3={x3c:.4f}")

    # Save spinodal
    all_pts=[]
    for seg in segs:
        for pt in seg:
            x1s,x3s=pt[0],pt[1]; x2s=max(1-x1s-x3s,0)
            all_pts.append([0.0,x1s,x3s,x2s])
    with open(os.path.join(OUT,f'spinodal_{tag}.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['residual','phi1','phi3','phi2']); w.writerows(all_pts)
    print(f"  Spinodal: {len(all_pts)} pts")

    # Full binodal (120 points)
    bino=compute_binodal(v1,v2_cf,v3,X13,g23,p,x3c,n_loop=120)
    np_f=len(bino)//2
    print(f"  Binodal: {np_f} pairs")

    with open(os.path.join(OUT,f'binodal_{tag}.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['residual','phi1_L8Bo','phi3_Donor','phi2_Chloroform']); w.writerows(bino)

    with open(os.path.join(OUT,f'critical_{tag}.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['phi1','phi2','phi3']); w.writerow([x1c,x2c,x3c])

    # Quality check
    if np_f>=3:
        c3=[bino[2*k][2] for k in range(np_f)]
        d3=[bino[2*k+1][2] for k in range(np_f)]
        print(f"  Quality: conc_mono={all(c3[i]>c3[i+1] for i in range(len(c3)-1))} "
              f"dil_mono={all(d3[i]<d3[i+1] for i in range(len(d3)-1))}")
        print(f"  Conc phi3: {c3[0]:.4f}->{c3[-1]:.4f}  Dil phi3: {d3[0]:.6f}->{d3[-1]:.4f}")

print("\nDone!")
