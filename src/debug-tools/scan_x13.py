"""
批量扫描 X13，找到能产生正常 binodal 形状的参数
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import os

# =================== 参数设置（Tol 体系）=====================
v1 = 1132.1
v2 = 106.3           # Toluene
v3 = 1743900
s = v1 / v2
r = v1 / v3

g23 = 0.3852
p1, p2, p3, p4, p5 = 0, 0, 0, -0.2000, 0.6000

def fun(x, ix3d, iX13):
    x3d = ix3d
    x1c, x2c, x1d = x
    x3c = 1 - x1c - x2c
    x2d = 1 - x1d - x3d
    eps = 1e-10
    x1c = max(x1c, eps); x2c = max(x2c, eps); x3c = max(x3c, eps)
    x1d = max(x1d, eps); x2d = max(x2d, eps); x3d = max(x3d, eps)

    def calc_g12(u2):
        g = p1*u2**4 + p2*u2**3 + p3*u2**2 + p4*u2 + p5
        dg = 4*p1*u2**3 + 3*p2*u2**2 + 2*p3*u2 + p4
        return g, dg

    u1c = x1c / (x1c + x2c)
    u2c = x2c / (x1c + x2c)
    u1d = x1d / (x1d + x2d)
    u2d = x2d / (x1d + x2d)

    g12c, dgdu2c = calc_g12(u2c)
    g12d, dgdu2d = calc_g12(u2d)

    X13c = iX13; X13d = iX13
    dX13c = 0; dX13d = 0
    g23c = g23; g23d = g23
    dg23c = 0; dg23d = 0

    mu1c = (np.log(x1c) + 1 - x1c - s*x2c - r*x3c +
            (g12c*x2c + X13c*x3c)*(x2c+x3c) -
            s*g23c*x2c*x3c - x2c*u1c*u2c*dgdu2c -
            x1c*x3c**2*dX13c - s*x2c*x3c**2*dg23c)
    mu1d = (np.log(x1d) + 1 - x1d - s*x2d - r*x3d +
            (g12d*x2d + X13d*x3d)*(x2d+x3d) -
            s*g23d*x2d*x3d - x2d*u1d*u2d*dgdu2d -
            x1d*x3d**2*dX13d - s*x2d*x3d**2*dg23d)
    mu2c = (s*np.log(x2c) + s - x1c - s*x2c - r*x3c +
            (g12c*x1c + g23c*s*x3c)*(x1c+x3c) -
            X13c*x1c*x3c + x1c*u1c*u2c*dgdu2c -
            x1c*x3c**2*dX13c - s*x2c*x3c**2*dg23c)
    mu2d = (s*np.log(x2d) + s - x1d - s*x2d - r*x3d +
            (g12d*x1d + g23d*s*x3d)*(x1d+x3d) -
            X13d*x1d*x3d + x1d*u1d*u2d*dgdu2d -
            x1d*x3d**2*dX13d - s*x2d*x3d**2*dg23d)
    mu3c = (r*np.log(x3c) + r - x1c - s*x2c - r*x3c +
            (X13c*x1c + s*g23c*x2c)*(x1c + x2c) -
            g12c*x1c*x2c + (x1c*dX13c + s*x2c*dg23c)*x3c*(x1c+x2c))
    mu3d = (r*np.log(x3d) + r - x1d - s*x2d - r*x3d +
            (X13d*x1d + s*g23d*x2d)*(x1d + x2d) -
            g12d*x1d*x2d + (x1d*dX13d + s*x2d*dg23d)*x3d*(x1d+x2d))
    return (mu1c - mu1d)**2 + (mu2c - mu2d)**2 + (mu3c - mu3d)**2


def solve_binodal(X13, nloop=80):
    ix3d_values = np.logspace(-5, np.log10(0.3), nloop)
    x0_default = np.array([0.2, 0.5, 0.5])
    x0 = x0_default.copy()
    results = []
    success_count = 0
    consecutive_fail = 0

    for i, ix3d in enumerate(ix3d_values):
        cons = [
            {'type': 'ineq', 'fun': lambda x: 1 - 1e-6 - x[0] - x[1]},
            {'type': 'ineq', 'fun': lambda x: 1 - ix3d - 1e-6 - x[2]},
            {'type': 'ineq', 'fun': lambda x: x[0] - 1e-6},
            {'type': 'ineq', 'fun': lambda x: x[1] - 1e-6},
            {'type': 'ineq', 'fun': lambda x: x[2] - 1e-6},
        ]
        res = minimize(lambda x: fun(x, ix3d, X13), x0, method='SLSQP',
                       constraints=cons, options={'ftol': 1e-10, 'maxiter': 10000, 'disp': False})
        x = res.x
        f = res.fun
        x2d_calc = 1 - x[2] - ix3d
        phase_diff = abs(x[0] - x[2]) + abs(x[1] - x2d_calc)
        comp_diff = abs((1 - x[0] - x[1]) - ix3d)
        is_success = f < 1e-7 and phase_diff > 1e-3 and comp_diff > 0.01
        if is_success:
            success_count += 1
            x0 = x.copy()
            consecutive_fail = 0
            results.append({
                'ix3d': ix3d, 'x1c': x[0], 'x2c': x[1], 'x3c': 1-x[0]-x[1],
                'x1d': x[2], 'x2d': x2d_calc, 'x3d': ix3d
            })
        else:
            consecutive_fail += 1
            if consecutive_fail >= 5:
                x0 = x0_default.copy()
                consecutive_fail = 0
    return results, success_count / nloop


# =================== 批量扫描 X13 ==========================
X13_candidates = [0.45, 0.50, 0.55, 0.58, 0.60, 0.62, 0.65]
fig, axes = plt.subplots(2, 4, figsize=(18, 10))
axes = axes.flatten()

best_X13 = None
best_score = -1
best_results = None

for idx, X13 in enumerate(X13_candidates):
    results, rate = solve_binodal(X13, nloop=80)
    ax = axes[idx]

    if len(results) > 3:
        x1c_list = [r['x1c'] for r in results]
        x3c_list = [r['x3c'] for r in results]
        x1d_list = [r['x1d'] for r in results]
        x3d_list = [r['x3d'] for r in results]

        ax.plot(x1c_list, x3c_list, 'bo-', markersize=3, label='Concentrated')
        ax.plot(x1d_list, x3d_list, 'ro-', markersize=3, label='Dilute')
        for res_item in results:
            ax.plot([res_item['x1c'], res_item['x1d']], [res_item['x3c'], res_item['x3d']], 'g--', alpha=0.3, linewidth=0.6)

        # 评分：成功率 + 曲线光滑度（用相邻点距离方差倒数近似）
        score = rate
        if len(x1c_list) > 5:
            dists_c = np.diff(np.array(x1c_list))**2 + np.diff(np.array(x3c_list))**2
            dists_d = np.diff(np.array(x1d_list))**2 + np.diff(np.array(x3d_list))**2
            smoothness = 1.0 / (1.0 + np.var(dists_c) + np.var(dists_d))
            score = rate * 0.7 + smoothness * 0.3

        if score > best_score:
            best_score = score
            best_X13 = X13
            best_results = results
    else:
        ax.text(0.5, 0.5, f'No data\n(rate={rate:.1%})', ha='center', va='center', transform=ax.transAxes)

    ax.set_title(f'X13={X13}  rate={rate:.1%}')
    ax.set_xlabel('x1 (L8-Bo)')
    ax.set_ylabel('x3 (PM6)')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.5)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)

# 最后一个子图画最佳结果的放大图
ax = axes[-1]
if best_results:
    x1c_list = [r['x1c'] for r in best_results]
    x3c_list = [r['x3c'] for r in best_results]
    x1d_list = [r['x1d'] for r in best_results]
    x3d_list = [r['x3d'] for r in best_results]
    ax.plot(x1c_list, x3c_list, 'bo-', markersize=4)
    ax.plot(x1d_list, x3d_list, 'ro-', markersize=4)
    for res_item in best_results:
        ax.plot([res_item['x1c'], res_item['x1d']], [res_item['x3c'], res_item['x3d']], 'g--', alpha=0.4, linewidth=0.8)
    ax.set_title(f'BEST: X13={best_X13}')
    ax.set_xlabel('x1 (L8-Bo)')
    ax.set_ylabel('x3 (PM6)')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.5)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('scan_x13_Tol.png', dpi=200)
print(f'最佳 X13 (Tol) = {best_X13}, 评分={best_score:.3f}')
print('图像已保存为 scan_x13_Tol.png')
plt.show()
