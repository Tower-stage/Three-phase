#!/usr/bin/env python3
"""
CF Parameter Sweep — find optimal (X13, g23, v3, p) for chloroform systems.
Equivalent to param_sweep_v2.m but in Python.

Sweeps: X13 × g23 × v3 × (g12 slope)
Scores: binodal success rate (n_pairs / 120)
"""

import numpy as np
from scipy.optimize import fsolve
import os, csv

# ============================================================
# FH KERNEL (self-contained, no cross-file imports)
# ============================================================
v1 = 1132.1       # L8-Bo molar volume
v2_cf = 80.7      # Chloroform molar volume

def chempot_eq(x, x3d, v1, v2, v3, s, r, X13, g23, p):
    x1c = max(x[0], 1e-10); x2c = max(x[1], 1e-10)
    x3c = max(1 - x1c - x2c, 1e-10)
    x1d = max(x[2], 1e-10)
    x2d = max(1 - x[2] - x3d, 1e-10)
    x3dv = max(x3d, 1e-10)
    u1c = x1c/(x1c+x2c); u2c = x2c/(x1c+x2c)
    g12c = p[0]*u2c**4+p[1]*u2c**3+p[2]*u2c**2+p[3]*u2c+p[4]
    dgc = 4*p[0]*u2c**3+3*p[1]*u2c**2+2*p[2]*u2c+p[3]
    u1d = x1d/(x1d+x2d); u2d = x2d/(x1d+x2d)
    g12d = p[0]*u2d**4+p[1]*u2d**3+p[2]*u2d**2+p[3]*u2d+p[4]
    dgd = 4*p[0]*u2d**3+3*p[1]*u2d**2+2*p[2]*u2d+p[3]
    F1 = (np.log(x1c)+1-x1c-s*x2c-r*x3c+(g12c*x2c+X13*x3c)*(x2c+x3c)-s*g23*x2c*x3c-x2c*u1c*u2c*dgc) \
       - (np.log(x1d)+1-x1d-s*x2d-r*x3dv+(g12d*x2d+X13*x3dv)*(x2d+x3dv)-s*g23*x2d*x3dv-x2d*u1d*u2d*dgd)
    F2 = (s*np.log(x2c)+s-x1c-s*x2c-r*x3c+(g12c*x1c+g23*s*x3c)*(x1c+x3c)-X13*x1c*x3c+x1c*u1c*u2c*dgc) \
       - (s*np.log(x2d)+s-x1d-s*x2d-r*x3dv+(g12d*x1d+g23*s*x3dv)*(x1d+x3dv)-X13*x1d*x3dv+x1d*u1d*u2d*dgd)
    F3 = (r*np.log(x3c)+r-x1c-s*x2c-r*x3c+(X13*x1c+s*g23*x2c)*(x1c+x2c)-g12c*x1c*x2c) \
       - (r*np.log(x3dv)+r-x1d-s*x2d-r*x3dv+(X13*x1d+s*g23*x2d)*(x1d+x2d)-g12d*x1d*x2d)
    return [F1, F2, F3]


def compute_spinodal_cp(v1, v2, v3, X13, g23, p, n_g=400):
    s = v1/v2; r = v1/v3
    x1v = np.linspace(1e-8, 0.999, n_g); x3v = np.linspace(1e-8, 0.999, n_g)
    X1, X3 = np.meshgrid(x1v, x3v); X2 = 1-X1-X3
    bad = (X1<=1e-8)|(X2<=1e-8)|(X3<=1e-8)
    u1 = X1/(X1+X2+1e-300); u2 = X2/(X1+X2+1e-300)
    g12=p[0]*u2**4+p[1]*u2**3+p[2]*u2**2+p[3]*u2+p[4]
    dg=4*p[0]*u2**3+3*p[1]*u2**2+2*p[2]*u2+p[3]
    d2g=12*p[0]*u2**2+6*p[1]*u2+2*p[2]
    G22=1/X1+s/X2-2*g12+2*(1-2*u2)*dg+u1*u2*d2g
    G23=1/X1-(g12+X13)+s*g23+u2*(1-3*u2)*dg+u1*u2**2*d2g
    G33=1/X1+r/X3-2*X13-2*u2**3*dg+u2**3*u1*d2g
    detG=G22*G33-G23**2; detG[bad]=np.nan

    # extract contour points
    valid = ~np.isnan(detG)
    if valid.sum()<100: return [], None, None
    spin_pts = []
    for i in range(n_g-1):
        for j in range(n_g-1):
            v00=detG[i,j]; v10=detG[i+1,j]; v01=detG[i,j+1]; v11=detG[i+1,j+1]
            if np.isnan(v00) or np.isnan(v10) or np.isnan(v01) or np.isnan(v11): continue
            signs=[np.sign(v00),np.sign(v10),np.sign(v01),np.sign(v11)]
            if all(s==0 for s in signs): continue
            if min(signs)<0 and max(signs)>=0:
                if v00*v10<=0:
                    t=abs(v00)/(abs(v00)+abs(v10)+1e-300)
                    spin_pts.append((x1v[j], x3v[i]+t*(x3v[i+1]-x3v[i])))
                if v00*v01<=0:
                    t=abs(v00)/(abs(v00)+abs(v01)+1e-300)
                    spin_pts.append((x1v[j]+t*(x1v[j+1]-x1v[j]), x3v[i]))
    if len(spin_pts)<10: return [], None, None
    pts=np.array(spin_pts)
    sort_idx=np.argsort(pts[:,1]); pts=pts[sort_idx]

    # Critical point
    x1_crit,x3_crit=None,None
    xs1=pts[:,0]; xs3=pts[:,1]; xs2=1-xs1-xs3
    np_pts=len(xs1); cv=np.full(np_pts, np.nan)
    for i in range(np_pts):
        x1i=xs1[i];x3i=xs3[i];x2i=xs2[i]
        if x1i<=1e-8 or x2i<=1e-8 or x3i<=1e-8: continue
        u1i=x1i/(x1i+x2i);u2i=x2i/(x1i+x2i)
        gi=p[0]*u2i**4+p[1]*u2i**3+p[2]*u2i**2+p[3]*u2i+p[4]
        dgi=4*p[0]*u2i**3+3*p[1]*u2i**2+2*p[2]*u2i+p[3]
        d2gi=12*p[0]*u2i**2+6*p[1]*u2i+2*p[2]
        G22i=1/x1i+s/x2i-2*gi+2*(1-2*u2i)*dgi+u1i*u2i*d2gi
        G23i=1/x1i-(gi+X13)+s*g23+u2i*(1-3*u2i)*dgi+u1i*u2i**2*d2gi
        G33i=1/x1i+r/x3i-2*X13-2*u2i**3*dgi+u2i**3*u1i*d2gi
        G222=1/x1i**2-s/x2i**2-6*u2i/x2i*dgi+(3-6*u2i)*u2i/x2i*d2gi
        G223=1/x1i**2-6*u2i**2/x2i*dgi+3*(1-2*u2i)*u2i**2/x2i*d2gi
        G233=1/x1i**2-6*u2i**3/x2i*dgi+(3*u2i-6*u2i**2)*u2i**2/x2i*d2gi
        G333=1/x1i**2-r/x3i**2+6*u2i**4/x2i*dgi+(3*u2i**2-2*u2i**3)*u2i**2/x2i*d2gi
        cv[i]=G222*G33i**2-3*G223*G23i*G33i+3*G233*G23i**2-G22i*G23i*G333
    ok=~np.isnan(cv)
    if ok.sum()>=2:
        cv_ok=cv[ok]
        for k in range(len(cv_ok)-1):
            if cv_ok[k]*cv_ok[k+1]<0:
                x1o=xs1[ok];x3o=xs3[ok]
                t_=abs(cv_ok[k])/(abs(cv_ok[k])+abs(cv_ok[k+1]))
                x1_crit=x1o[k]+t_*(x1o[k+1]-x1o[k])
                x3_crit=x3o[k]+t_*(x3o[k+1]-x3o[k])
                break
    if x3_crit is None and len(pts)>0:
        x3_crit=float(np.max(pts[:,1])); x1_crit=float(pts[np.argmax(pts[:,1]),0])
    return [(pts.tolist())], x1_crit, x3_crit


def compute_binodal(v1, v2, v3, X13, g23, p, x3_crit, n_loop=40):
    s=v1/v2; r=v1/v3
    if x3_crit is None or np.isnan(x3_crit): x3_crit=0.01
    x3d_max=max(0.005, x3_crit*0.90)
    x3d_vals=np.logspace(-5, np.log10(x3d_max), n_loop)
    guesses=[[0.08,0.85,0.12],[0.10,0.80,0.15],[0.12,0.75,0.18],
             [0.15,0.70,0.20],[0.18,0.65,0.22],[0.20,0.60,0.25],
             [0.25,0.50,0.30],[0.15,0.72,0.19]]
    x_sol=[]; prev=None
    for x3d in x3d_vals:
        best_f=np.inf; best=None
        for x0 in guesses:
            if 1-x0[0]-x0[1]<=0: continue
            try:
                xt,infodict,ier,msg=fsolve(
                    lambda x: chempot_eq(x,x3d,v1,v2,v3,s,r,X13,g23,p),
                    x0, full_output=True, xtol=1e-12, maxfev=5000)
                fv=chempot_eq(xt,x3d,v1,v2,v3,s,r,X13,g23,p)
                fv2=sum(np.array(fv)**2); x3c=1-xt[0]-xt[1]
                if fv2<1e-10 and x3c>x3d+0.003 and abs(xt[0]-xt[2])>0.005:
                    if fv2<best_f: best_f=fv2; best=xt
            except: pass
        if best is None and prev is not None:
            try:
                xt,infodict,ier,msg=fsolve(
                    lambda x: chempot_eq(x,x3d,v1,v2,v3,s,r,X13,g23,p),
                    prev, full_output=True, xtol=1e-12, maxfev=5000)
                fv=chempot_eq(xt,x3d,v1,v2,v3,s,r,X13,g23,p)
                fv2=sum(np.array(fv)**2); x3c=1-xt[0]-xt[1]
                if fv2<1e-10 and x3c>x3d+0.003: best_f=fv2; best=xt
            except: pass
        if best is not None:
            x3c=1-best[0]-best[1]; x2d=1-best[2]-x3d
            x_sol.append([best_f,best[0],x3c,best[1]])
            x_sol.append([best_f,best[2],x3d,x2d]); prev=best
    if len(x_sol)<6: return []
    x_sol=np.array(x_sol)
    mask=~np.isnan(x_sol[:,0]); x_sol=x_sol[mask]
    mask=x_sol[:,0]<1e-8; x_sol=x_sol[mask]
    np_=len(x_sol)//2
    if np_<3: return x_sol.tolist()
    dilute_mask=np.arange(1,2*np_,2)
    sort_idx=np.argsort(x_sol[dilute_mask,2])
    tmp=x_sol.copy()
    for k in range(np_): x_sol[2*k:2*k+2]=tmp[2*sort_idx[k]:2*sort_idx[k]+2]
    keep=np.ones(np_,dtype=bool); pc,pd=np.inf,-np.inf
    for k in range(np_):
        c3=x_sol[2*k,2]; d3=x_sol[2*k+1,2]
        if c3>=pc or c3<=0 or d3<=pd: keep[k]=False
        else: pc,pd=c3,d3
    filtered=[]
    for k in range(np_):
        if keep[k]: filtered.extend([x_sol[2*k].tolist(),x_sol[2*k+1].tolist()])
    return filtered


def score_system(X13, g23, v3, p, n_loop=60):
    """Quick evaluation: binodal success count."""
    spin_segs, x1_crit, x3_crit = compute_spinodal_cp(v1, v2_cf, v3, X13, g23, p, n_g=400)
    if x3_crit is None or np.isnan(x3_crit):
        return 0, None, None
    x2_crit = 1 - x1_crit - x3_crit

    binodal = compute_binodal(v1, v2_cf, v3, X13, g23, p, x3_crit, n_loop=n_loop)
    np_final = len(binodal) // 2
    return np_final, x2_crit, x3_crit


def sweep_pm6():
    """Sweep PM6-CF parameters."""
    print("="*60)
    print("  Sweeping PM6 / L8-Bo / Chloroform")
    print("="*60)

    X13_vals = [0.85, 0.90, 0.95, 1.00]
    g23_vals = [0.15, 0.20, 0.25, 0.30, 0.35]
    v3_vals  = [80000, 100000, 120000, 150000]
    p4_vals  = [-0.10, -0.15, -0.20]
    p5_vals  = [0.40, 0.45, 0.50, 0.55]

    best_score = 0
    best_params = None
    results = []

    total = len(X13_vals)*len(g23_vals)*len(v3_vals)*len(p4_vals)*len(p5_vals)
    count = 0

    for X13 in X13_vals:
        for g23 in g23_vals:
            for v3 in v3_vals:
                for p4 in p4_vals:
                    for p5 in p5_vals:
                        count += 1
                        p = [0, 0, 0, p4, p5]
                        np_final, x2_crit, x3_crit = score_system(X13, g23, v3, p, n_loop=40)
                        results.append((np_final, x2_crit, X13, g23, v3, p4, p5))

                        if np_final > best_score:
                            best_score = np_final
                            best_params = (X13, g23, v3, p4, p5, x2_crit)

                        if count % 20 == 0:
                            print(f"  [{count}/{total}] best={best_score} pairs | "
                                  f"X13={best_params[0]:.2f} g23={best_params[1]:.3f} "
                                  f"v3={best_params[2]} p4={best_params[3]:.2f} p5={best_params[4]:.2f} "
                                  f"CP_phi2={best_params[5]:.3f}")

    # Sort and show top 10
    results.sort(key=lambda x: x[0], reverse=True)
    print(f"\n  Top 10 PM6-CF parameter sets:")
    print(f"  {'Rank':<5s} {'Pairs':<8s} {'CP phi2':<10s} {'X13':<8s} {'g23':<8s} {'v3':<10s} {'p4':<8s} {'p5':<8s}")
    for i, (np_f, cp2, X13, g23, v3, p4, p5) in enumerate(results[:10]):
        print(f"  {i+1:<5d} {np_f:<8d} {cp2:<10.4f} {X13:<8.2f} {g23:<8.3f} {v3:<10d} {p4:<8.2f} {p5:<8.2f}")

    return results


def sweep_d18():
    """Sweep D18-CF parameters."""
    print("\n" + "="*60)
    print("  Sweeping D18 / L8-Bo / Chloroform")
    print("="*60)

    X13_vals = [0.85, 0.90, 0.95]
    g23_vals = [0.20, 0.25, 0.30, 0.35, 0.40]
    v3_vals  = [50000, 60000, 80000, 100000]
    p4_vals  = [-0.10, -0.15, -0.20]
    p5_vals  = [0.40, 0.45, 0.50, 0.55]

    best_score = 0
    best_params = None
    results = []

    total = len(X13_vals)*len(g23_vals)*len(v3_vals)*len(p4_vals)*len(p5_vals)
    count = 0

    for X13 in X13_vals:
        for g23 in g23_vals:
            for v3 in v3_vals:
                for p4 in p4_vals:
                    for p5 in p5_vals:
                        count += 1
                        p = [0, 0, 0, p4, p5]
                        np_final, x2_crit, x3_crit = score_system(X13, g23, v3, p, n_loop=40)
                        results.append((np_final, x2_crit, X13, g23, v3, p4, p5))

                        if np_final > best_score:
                            best_score = np_final
                            best_params = (X13, g23, v3, p4, p5, x2_crit)

                        if count % 15 == 0:
                            print(f"  [{count}/{total}] best={best_score} pairs | "
                                  f"X13={best_params[0]:.2f} g23={best_params[1]:.3f} "
                                  f"v3={best_params[2]} p4={best_params[3]:.2f} p5={best_params[4]:.2f} "
                                  f"CP_phi2={best_params[5]:.3f}")

    results.sort(key=lambda x: x[0], reverse=True)
    print(f"\n  Top 10 D18-CF parameter sets:")
    print(f"  {'Rank':<5s} {'Pairs':<8s} {'CP phi2':<10s} {'X13':<8s} {'g23':<8s} {'v3':<10s} {'p4':<8s} {'p5':<8s}")
    for i, (np_f, cp2, X13, g23, v3, p4, p5) in enumerate(results[:10]):
        print(f"  {i+1:<5d} {np_f:<8d} {cp2:<10.4f} {X13:<8.2f} {g23:<8.3f} {v3:<10d} {p4:<8.2f} {p5:<8.2f}")

    return results


if __name__ == '__main__':
    pm6_results = sweep_pm6()
    d18_results = sweep_d18()

    # Save sweep results
    sweep_path = os.path.join(DATA_OUT, '..', 'Dissolve', 'cf_sweep_results.csv')
    # Actually save to Dissolve folder
    sweep_dir = os.path.join(os.path.dirname(__file__))
    sweep_path = os.path.join(sweep_dir, 'cf_sweep_results.csv')
    with open(sweep_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['system','rank','pairs','CP_phi2','X13','g23','v3','p4','p5'])
        for i, (np_f, cp2, X13, g23, v3, p4, p5) in enumerate(pm6_results[:20]):
            w.writerow(['PM6_CF', i+1, np_f, f'{cp2:.4f}', X13, g23, v3, p4, p5])
        for i, (np_f, cp2, X13, g23, v3, p4, p5) in enumerate(d18_results[:20]):
            w.writerow(['D18_CF', i+1, np_f, f'{cp2:.4f}', X13, g23, v3, p4, p5])
    print(f"\nSweep results saved: {sweep_path}")
