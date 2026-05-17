import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

v1, v2, v3 = 1132.1, 106.3, 1743900
s, r = v1/v2, v1/v3
g23 = 0.3852
p = [0, 0, 0, -0.2000, 0.6000]
X13 = 0.62

def fun(x, ix3d):
    x3d = ix3d
    x1c, x2c, x1d = x
    x3c = 1 - x1c - x2c
    x2d = 1 - x1d - x3d
    eps = 1e-10
    x1c = max(x1c, eps); x2c = max(x2c, eps); x3c = max(x3c, eps)
    x1d = max(x1d, eps); x2d = max(x2d, eps); x3d = max(x3d, eps)
    def g12(u2):
        return p[0]*u2**4 + p[1]*u2**3 + p[2]*u2**2 + p[3]*u2 + p[4]
    def dg(u2):
        return 4*p[0]*u2**3 + 3*p[1]*u2**2 + 2*p[2]*u2 + p[3]
    u1c = x1c/(x1c+x2c); u2c = x2c/(x1c+x2c)
    u1d = x1d/(x1d+x2d); u2d = x2d/(x1d+x2d)
    g12c, dg12c = g12(u2c), dg(u2c)
    g12d, dg12d = g12(u2d), dg(u2d)
    mu1c = np.log(x1c)+1-x1c-s*x2c-r*x3c+(g12c*x2c+X13*x3c)*(x2c+x3c)-s*g23*x2c*x3c-x2c*u1c*u2c*dg12c
    mu1d = np.log(x1d)+1-x1d-s*x2d-r*x3d+(g12d*x2d+X13*x3d)*(x2d+x3d)-s*g23*x2d*x3d-x2d*u1d*u2d*dg12d
    mu2c = s*np.log(x2c)+s-x1c-s*x2c-r*x3c+(g12c*x1c+g23*s*x3c)*(x1c+x3c)-X13*x1c*x3c+x1c*u1c*u2c*dg12c
    mu2d = s*np.log(x2d)+s-x1d-s*x2d-r*x3d+(g12d*x1d+g23*s*x3d)*(x1d+x3d)-X13*x1d*x3d+x1d*u1d*u2d*dg12d
    mu3c = r*np.log(x3c)+r-x1c-s*x2c-r*x3c+(X13*x1c+s*g23*x2c)*(x1c+x2c)-g12c*x1c*x2c
    mu3d = r*np.log(x3d)+r-x1d-s*x2d-r*x3d+(X13*x1d+s*g23*x2d)*(x1d+x2d)-g12d*x1d*x2d
    return (mu1c-mu1d)**2+(mu2c-mu2d)**2+(mu3c-mu3d)**2

def solve(use_warmstart=True, nloop=80):
    ix3d_values = np.logspace(-5, np.log10(0.3), nloop)
    x0_default = np.array([0.2, 0.5, 0.5])
    x0 = x0_default.copy()
    results = []
    success_count = 0
    for i, ix3d in enumerate(ix3d_values):
        cons = [
            {'type': 'ineq', 'fun': lambda x, ix=ix3d: 1-1e-6-x[0]-x[1]},
            {'type': 'ineq', 'fun': lambda x, ix=ix3d: 1-ix-1e-6-x[2]},
            {'type': 'ineq', 'fun': lambda x: x[0]-1e-6},
            {'type': 'ineq', 'fun': lambda x: x[1]-1e-6},
            {'type': 'ineq', 'fun': lambda x: x[2]-1e-6},
        ]
        res = minimize(lambda x: fun(x, ix3d), x0, method='SLSQP', constraints=cons, options={'ftol':1e-10,'maxiter':10000,'disp':False})
        x = res.x
        f = res.fun
        x2d = 1-x[2]-ix3d
        pd = abs(x[0]-x[2])+abs(x[1]-x2d)
        cd = abs((1-x[0]-x[1])-ix3d)
        ok = f < 1e-7 and pd > 1e-3 and cd > 0.01
        if ok:
            success_count += 1
            if use_warmstart:
                x0 = x.copy()
            else:
                x0 = x0_default.copy()
            results.append({'x1c':x[0],'x3c':1-x[0]-x[1],'x1d':x[2],'x3d':ix3d})
        else:
            x0 = x0_default.copy()
    return results, success_count/nloop

fig, axes = plt.subplots(1,2, figsize=(12,5))
for idx, ws in enumerate([False, True]):
    res, rate = solve(use_warmstart=ws, nloop=80)
    ax = axes[idx]
    if len(res) > 3:
        ax.plot([r['x1c'] for r in res], [r['x3c'] for r in res], 'bo-', markersize=3, label='Conc')
        ax.plot([r['x1d'] for r in res], [r['x3d'] for r in res], 'ro-', markersize=3, label='Dilute')
        for item in res:
            ax.plot([item['x1c'], item['x1d']], [item['x3c'], item['x3d']], 'g--', alpha=0.3, lw=0.6)
    ax.set_title(f'Warmstart={ws} rate={rate:.1%}')
    ax.set_xlim(0,1); ax.set_ylim(0,0.5)
    ax.set_xlabel('x1'); ax.set_ylabel('x3')
    ax.legend(); ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('smooth_compare_Tol.png', dpi=200)
print('Saved smooth_compare_Tol.png')
plt.show()
