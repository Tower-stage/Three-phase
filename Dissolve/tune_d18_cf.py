"""Quick D18-CF parameter tuning."""
from optimize_cf import chempot_eq, compute_spinodal_cp, compute_binodal
import numpy as np

v1=1132.1; v2=80.7; v3=60000

tests = [
    # X13, g23, p4, p5
    (0.90, 0.20, -0.12, 0.40),
    (0.90, 0.22, -0.12, 0.40),
    (0.90, 0.25, -0.12, 0.40),
    (0.90, 0.28, -0.12, 0.40),
    (0.95, 0.20, -0.12, 0.40),
    (0.95, 0.22, -0.15, 0.45),
    (0.95, 0.28, -0.15, 0.45),
    (0.95, 0.25, -0.15, 0.45),
    (1.00, 0.22, -0.12, 0.40),
    (1.00, 0.25, -0.12, 0.40),
    (0.90, 0.20, -0.15, 0.50),
    (0.90, 0.18, -0.12, 0.40),
]

for X13,g23,p4,p5 in tests:
    p=[0,0,0,p4,p5]
    segs,x1c,x3c=compute_spinodal_cp(v1,v2,v3,X13,g23,p,n_g=600)
    if x3c is None: print(f"  X13={X13} g23={g23} p4={p4} p5={p5} -> NO CP"); continue
    x2c=1-x1c-x3c
    bino=compute_binodal(v1,v2,v3,X13,g23,p,x3c,n_loop=60)
    np_f=len(bino)//2
    marker="***" if np_f>=50 else ("**" if np_f>=30 else ("*" if np_f>=20 else ""))
    c3_max=bino[0][2] if np_f>=1 else 0
    print(f"  X13={X13:.2f} g23={g23:.3f} p4={p4:.2f} p5={p5:.2f} -> {np_f:3d}pairs CP_phi2={x2c:.3f} phi3max={c3_max:.3f} {marker}")
