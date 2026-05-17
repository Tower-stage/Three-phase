"""
生成最终高质量的 Binodal 相图（Tol 和 O-Xy）
参数：Tol X13=0.62, O-Xy X13=0.65
策略：无热启动（曲线更光滑）
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def solve(system, X13, nloop=100, x3max=0.35):
    if system == 'Tol':
        v1, v2, v3 = 1132.1, 106.3, 1743900
        g23 = 0.3852
        p = [0, 0, 0, -0.2000, 0.6000]
    else:
        v1, v2, v3 = 1132.1, 120.6, 1743900
        g23 = 0.4120
        p = [0, 0, 0, -0.2000, 0.7500]
    s, r = v1/v2, v1/v3

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

    ix3d_values = np.logspace(-5, np.log10(x3max), nloop)
    x0_default = np.array([0.2, 0.5, 0.5])
    results = []
    for ix3d in ix3d_values:
        cons = [
            {'type': 'ineq', 'fun': lambda x: 1-1e-6-x[0]-x[1]},
            {'type': 'ineq', 'fun': lambda x, ix=ix3d: 1-ix-1e-6-x[2]},
            {'type': 'ineq', 'fun': lambda x: x[0]-1e-6},
            {'type': 'ineq', 'fun': lambda x: x[1]-1e-6},
            {'type': 'ineq', 'fun': lambda x: x[2]-1e-6},
        ]
        res = minimize(lambda x: fun(x, ix3d), x0_default, method='SLSQP', constraints=cons, options={'ftol':1e-10,'maxiter':10000,'disp':False})
        x = res.x; f = res.fun
        x2d = 1-x[2]-ix3d
        pd = abs(x[0]-x[2])+abs(x[1]-x2d)
        cd = abs((1-x[0]-x[1])-ix3d)
        if f < 1e-7 and pd > 1e-3 and cd > 0.01:
            results.append({'x1c':x[0],'x2c':x[1],'x3c':1-x[0]-x[1],'x1d':x[2],'x2d':x2d,'x3d':ix3d})
    return results

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for idx, (system, X13) in enumerate([('Tol', 0.62), ('O-Xy', 0.65)]):
    res = solve(system, X13, nloop=120, x3max=0.35)
    ax = axes[idx]

    x1c = [r['x1c'] for r in res]
    x3c = [r['x3c'] for r in res]
    x1d = [r['x1d'] for r in res]
    x3d = [r['x3d'] for r in res]

    # 画 tie line（每隔几个画一条，避免太密）
    for i, item in enumerate(res):
        if i % 5 == 0:
            ax.plot([item['x1c'], item['x1d']], [item['x3c'], item['x3d']], 'g--', alpha=0.4, linewidth=0.8)

    ax.plot(x1c, x3c, 'b-', linewidth=2, label='Concentrated phase arm')
    ax.plot(x1d, x3d, 'r-', linewidth=2, label='Dilute phase arm')
    ax.scatter(x1c, x3c, c='blue', s=15, zorder=5)
    ax.scatter(x1d, x3d, c='red', s=15, zorder=5)

    ax.set_xlabel('Volume fraction of L8-Bo ($\\phi_1$)', fontsize=12)
    ax.set_ylabel('Volume fraction of PM6 ($\\phi_3$)', fontsize=12)
    ax.set_title(f'{system} system  (X13={X13})', fontsize=13)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.35)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # 标注关键点
    if len(res) > 0:
        # 近似临界点（x3 最大处）
        max_idx = np.argmax(x3c)
        ax.annotate('Critical region', xy=(x1c[max_idx], x3c[max_idx]),
                    xytext=(x1c[max_idx]+0.15, x3c[max_idx]+0.05),
                    arrowprops=dict(arrowstyle='->', color='gray', alpha=0.7),
                    fontsize=10, color='gray')

plt.tight_layout()
plt.savefig('final_binodal.png', dpi=300, bbox_inches='tight')
print('最终相图已保存: final_binodal.png')
plt.show()
