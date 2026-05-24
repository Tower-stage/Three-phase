#!/usr/bin/env python3
"""Focused CF parameter optimization — v3 fixed, tune X13, g23, p4, p5."""
import numpy as np
from scipy.optimize import fsolve
import os, csv, sys, time

# ============================================================
# FH KERNEL
# ============================================================
v1 = 1132.1; v2_cf = 80.7

def chempot_eq(x, x3d, v1, v2, v3, s, r, X13, g23, p):
    x1c=max(x[0],1e-10); x2c=max(x[1],1e-10); x3c=max(1-x1c-x2c,1e-10)
    x1d=max(x[2],1e-10); x2d=max(1-x[2]-x3d,1e-10); x3dv=max(x3d,1e-10)
    u1c=x1c/(x1c+x2c); u2c=x2c/(x1c+x2c)
    g12c=p[0]*u2c**4+p[1]*u2c**3+p[2]*u2c**2+p[3]*u2c+p[4]
    dgc=4*p[0]*u2c**3+3*p[1]*u2c**2+2*p[2]*u2c+p[3]
    u1d=x1d/(x1d+x2d); u2d=x2d/(x1d+x2d)
    g12d=p[0]*u2d**4+p[1]*u2d**3+p[2]*u2d**2+p[3]*u2d+p[4]
    dgd=4*p[0]*u2d**3+3*p[1]*u2d**2+2*p[2]*u2d+p[3]
    F1=(np.log(x1c)+1-x1c-s*x2c-r*x3c+(g12c*x2c+X13*x3c)*(x2c+x3c)-s*g23*x2c*x3c-x2c*u1c*u2c*dgc)-(np.log(x1d)+1-x1d-s*x2d-r*x3dv+(g12d*x2d+X13*x3dv)*(x2d+x3dv)-s*g23*x2d*x3dv-x2d*u1d*u2d*dgd)
    F2=(s*np.log(x2c)+s-x1c-s*x2c-r*x3c+(g12c*x1c+g23*s*x3c)*(x1c+x3c)-X13*x1c*x3c+x1c*u1c*u2c*dgc)-(s*np.log(x2d)+s-x1d-s*x2d-r*x3dv+(g12d*x1d+g23*s*x3dv)*(x1d+x3dv)-X13*x1d*x3dv+x1d*u1d*u2d*dgd)
    F3=(r*np.log(x3c)+r-x1c-s*x2c-r*x3c+(X13*x1c+s*g23*x2c)*(x1c+x2c)-g12c*x1c*x2c)-(r*np.log(x3dv)+r-x1d-s*x2d-r*x3dv+(X13*x1d+s*g23*x2d)*(x1d+x2d)-g12d*x1d*x2d)
    return [F1,F2,F3]

def compute_spinodal_cp(v1,v2,v3,X13,g23,p,n_g=600):
    s=v1/v2; r=v1/v3
    x1v=np.linspace(1e-8,0.999,n_g); x3v=np.linspace(1e-8,0.999,n_g)
    X1,X3=np.meshgrid(x1v,x3v); X2=1-X1-X3
    bad=(X1<=1e-8)|(X2<=1e-8)|(X3<=1e-8)
    u1=X1/(X1+X2+1e-300); u2=X2/(X1+X2+1e-300)
    g12=p[0]*u2**4+p[1]*u2**3+p[2]*u2**2+p[3]*u2+p[4]
    dg=4*p[0]*u2**3+3*p[1]*u2**2+2*p[2]*u2+p[3]
    d2g=12*p[0]*u2**2+6*p[1]*u2+2*p[2]
    G22=1/X1+s/X2-2*g12+2*(1-2*u2)*dg+u1*u2*d2g
    G23=1/X1-(g12+X13)+s*g23+u2*(1-3*u2)*dg+u1*u2**2*d2g
    G33=1/X1+r/X3-2*X13-2*u2**3*dg+u2**3*u1*d2g
    detG=G22*G33-G23**2; detG[bad]=np.nan
    valid=~np.isnan(detG)
    if valid.sum()<100: return [],None,None
    spin_pts=[]
    for i in range(n_g-1):
        for j in range(n_g-1):
            v00=detG[i,j];v10=detG[i+1,j];v01=detG[i,j+1];v11=detG[i+1,j+1]
            if np.isnan(v00)|np.isnan(v10)|np.isnan(v01)|np.isnan(v11): continue
            signs=[np.sign(v00),np.sign(v10),np.sign(v01),np.sign(v11)]
            if all(s==0 for s in signs): continue
            if min(signs)<0 and max(signs)>=0:
                if v00*v10<=0:
                    t=abs(v00)/(abs(v00)+abs(v10)+1e-300)
                    spin_pts.append((x1v[j],x3v[i]+t*(x3v[i+1]-x3v[i])))
                if v00*v01<=0:
                    t=abs(v00)/(abs(v00)+abs(v01)+1e-300)
                    spin_pts.append((x1v[j]+t*(x1v[j+1]-x1v[j]),x3v[i]))
    if len(spin_pts)<10: return [],None,None
    pts=np.array(spin_pts); sort_idx=np.argsort(pts[:,1]); pts=pts[sort_idx]
    xs1=pts[:,0];xs3=pts[:,1];xs2=1-xs1-xs3; np_pts=len(xs1)
    cv=np.full(np_pts,np.nan)
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
    ok=~np.isnan(cv); x1_crit=None; x3_crit=None
    if ok.sum()>=2:
        cv_ok=cv[ok]; x1o=xs1[ok];x3o=xs3[ok]
        for k in range(len(cv_ok)-1):
            if cv_ok[k]*cv_ok[k+1]<0:
                t_=abs(cv_ok[k])/(abs(cv_ok[k])+abs(cv_ok[k+1]))
                x1_crit=x1o[k]+t_*(x1o[k+1]-x1o[k])
                x3_crit=x3o[k]+t_*(x3o[k+1]-x3o[k]); break
    if x3_crit is None and len(pts)>0:
        x3_crit=float(np.max(pts[:,1])); x1_crit=float(pts[np.argmax(pts[:,1]),0])
    return [(pts.tolist())], x1_crit, x3_crit

def compute_binodal(v1,v2,v3,X13,g23,p,x3_crit,n_loop=60):
    s=v1/v2; r=v1/v3
    if x3_crit is None or np.isnan(x3_crit): x3_crit=0.01
    x3d_max=max(0.005,x3_crit*0.90)
    x3d_vals=np.logspace(-5,np.log10(x3d_max),n_loop)
    guesses=[[0.08,0.85,0.12],[0.10,0.80,0.15],[0.12,0.75,0.18],
             [0.15,0.70,0.20],[0.18,0.65,0.22],[0.20,0.60,0.25],
             [0.25,0.50,0.30],[0.15,0.72,0.19]]
    x_sol=[]; prev=None
    for x3d in x3d_vals:
        best_f=np.inf; best=None
        for x0 in guesses:
            if 1-x0[0]-x0[1]<=0: continue
            try:
                xt,infodict,ier,msg=fsolve(lambda x:chempot_eq(x,x3d,v1,v2,v3,s,r,X13,g23,p),x0,full_output=True,xtol=1e-12,maxfev=5000)
                fv=chempot_eq(xt,x3d,v1,v2,v3,s,r,X13,g23,p); fv2=sum(np.array(fv)**2); x3c=1-xt[0]-xt[1]
                if fv2<1e-10 and x3c>x3d+0.003 and abs(xt[0]-xt[2])>0.005:
                    if fv2<best_f: best_f=fv2; best=xt
            except: pass
        if best is None and prev is not None:
            try:
                xt,infodict,ier,msg=fsolve(lambda x:chempot_eq(x,x3d,v1,v2,v3,s,r,X13,g23,p),prev,full_output=True,xtol=1e-12,maxfev=5000)
                fv=chempot_eq(xt,x3d,v1,v2,v3,s,r,X13,g23,p); fv2=sum(np.array(fv)**2); x3c=1-xt[0]-xt[1]
                if fv2<1e-10 and x3c>x3d+0.003: best_f=fv2; best=xt
            except: pass
        if best is not None:
            x3c=1-best[0]-best[1]; x2d=1-best[2]-x3d
            x_sol.append([best_f,best[0],x3c,best[1]]); x_sol.append([best_f,best[2],x3d,x2d]); prev=best
    if len(x_sol)<6: return []
    x_sol=np.array(x_sol)
    mask=~np.isnan(x_sol[:,0]); x_sol=x_sol[mask]; mask=x_sol[:,0]<1e-8; x_sol=x_sol[mask]
    np_=len(x_sol)//2
    if np_<3: return x_sol.tolist()
    dilute_mask=np.arange(1,2*np_,2); sort_idx=np.argsort(x_sol[dilute_mask,2])
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

# ============================================================
# FOCUSED SWEEP
# ============================================================
def test_one(name, X13, g23, v3, p4, p5, n_g=600, n_loop=60):
    p=[0,0,0,p4,p5]
    segs,x1c,x3c=compute_spinodal_cp(v1,v2_cf,v3,X13,g23,p,n_g=n_g)
    if x3c is None or np.isnan(x3c): return 0,None
    x2c=1-x1c-x3c
    bino=compute_binodal(v1,v2_cf,v3,X13,g23,p,x3c,n_loop=n_loop)
    np_f=len(bino)//2
    return np_f,x2c

def run():
    sys.stdout.write("="*60+"\n  CF Parameter Optimization (focused)\n"+"="*60+"\n")
    sys.stdout.flush()

    # PM6-CF: v3 fixed at 100k, sweep X13 × g23 × (p4,p5)
    v3_pm6=100000
    combos_pm6=[]
    for X13 in [0.90,0.95,1.00]:
      for g23 in [0.18,0.22,0.25,0.28,0.32]:
        for p4 in [-0.12,-0.15,-0.18,-0.20]:
          for p5 in [0.40,0.45,0.50,0.55]:
            combos_pm6.append((X13,g23,p4,p5))

    # D18-CF: v3 fixed at 60k
    v3_d18=60000
    combos_d18=[]
    for X13 in [0.90,0.95,1.00]:
      for g23 in [0.22,0.26,0.30,0.34,0.38]:
        for p4 in [-0.12,-0.15,-0.18,-0.20]:
          for p5 in [0.40,0.45,0.50,0.55]:
            combos_d18.append((X13,g23,p4,p5))

    total=len(combos_pm6)+len(combos_d18)
    sys.stdout.write(f"Total combos: PM6={len(combos_pm6)} D18={len(combos_d18)} = {total}\n")
    sys.stdout.write(f"Estimated time: ~{total*2//60} min\n\n")
    sys.stdout.flush()

    results_pm6=[]
    t0=time.time()
    for k,(X13,g23,p4,p5) in enumerate(combos_pm6):
        np_f,cp2=test_one('PM6',X13,g23,v3_pm6,p4,p5)
        results_pm6.append((np_f,cp2,X13,g23,p4,p5))
        if (k+1)%20==0:
            best=max(results_pm6,key=lambda x:x[0])
            elapsed=time.time()-t0
            sys.stdout.write(f"  PM6 [{k+1}/{len(combos_pm6)}] best={best[0]}pairs CP={best[1]:.3f} X13={best[2]} g23={best[3]} p4={best[4]} p5={best[5]} ({elapsed:.0f}s)\n")
            sys.stdout.flush()

    results_d18=[]
    for k,(X13,g23,p4,p5) in enumerate(combos_d18):
        np_f,cp2=test_one('D18',X13,g23,v3_d18,p4,p5)
        results_d18.append((np_f,cp2,X13,g23,p4,p5))
        if (k+1)%20==0:
            best=max(results_d18,key=lambda x:x[0])
            elapsed=time.time()-t0
            sys.stdout.write(f"  D18 [{k+1}/{len(combos_d18)}] best={best[0]}pairs CP={best[1]:.3f} X13={best[2]} g23={best[3]} p4={best[4]} p5={best[5]} ({elapsed:.0f}s)\n")
            sys.stdout.flush()

    # Top results
    results_pm6.sort(key=lambda x:x[0],reverse=True)
    results_d18.sort(key=lambda x:x[0],reverse=True)

    sys.stdout.write(f"\n{'='*60}\n  TOP PM6-CF RESULTS\n{'='*60}\n")
    for i,(np_f,cp2,X13,g23,p4,p5) in enumerate(results_pm6[:8]):
        sys.stdout.write(f"  #{i+1}: {np_f}pairs CP_phi2={cp2:.3f} X13={X13:.2f} g23={g23:.3f} v3={v3_pm6} p=[0,0,0,{p4},{p5}]\n")

    sys.stdout.write(f"\n{'='*60}\n  TOP D18-CF RESULTS\n{'='*60}\n")
    for i,(np_f,cp2,X13,g23,p4,p5) in enumerate(results_d18[:8]):
        sys.stdout.write(f"  #{i+1}: {np_f}pairs CP_phi2={cp2:.3f} X13={X13:.2f} g23={g23:.3f} v3={v3_d18} p=[0,0,0,{p4},{p5}]\n")

    # Save
    out_dir=os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(out_dir,'cf_optimized_params.csv'),'w',newline='') as f:
        w=csv.writer(f)
        w.writerow(['system','rank','pairs','CP_phi2','X13','g23','v3','p4','p5'])
        for i,(np_f,cp2,X13,g23,p4,p5) in enumerate(results_pm6[:10]):
            w.writerow(['PM6_CF',i+1,np_f,f'{cp2:.3f}',X13,g23,v3_pm6,p4,p5])
        for i,(np_f,cp2,X13,g23,p4,p5) in enumerate(results_d18[:10]):
            w.writerow(['D18_CF',i+1,np_f,f'{cp2:.3f}',X13,g23,v3_d18,p4,p5])

    elapsed=time.time()-t0
    sys.stdout.write(f"\nDone in {elapsed:.0f}s ({elapsed/60:.1f} min)\n")
    sys.stdout.flush()

if __name__=='__main__':
    run()
