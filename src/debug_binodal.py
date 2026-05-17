"""
Binodal 计算调试脚本（Python 验证版）
用途：在不依赖 MATLAB 的情况下，验证参数合理性、测试收敛性、绘制相图轮廓
运行方式：python debug_binodal.py
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import csv

# =================== 参数设置（O-Xy 体系，可切换为 Tol）=====================
v1 = 1132.1          # L8-Bo
v2 = 120.6           # o-Xylene
v3 = 1743900         # PM6 (注意：已修正为 1743900)
s = v1 / v2
r = v1 / v3

X13 = 0.55           # L8-Bo/PM6
g23 = 0.4120         # o-Xy/PM6

# L8-Bo / o-Xylene 多项式系数
p1, p2, p3, p4, p5 = 16.6738, -40.8682, 39.6207, -19.3621, -1.5993

# ========================= 核心计算函数 ==========================
def fun(x, ix3d, iX13):
    """
    目标函数：两相化学势差值的平方和
    x = [x1c, x2c, x1d]
    """
    x3d = ix3d
    x1c, x2c, x1d = x
    x3c = 1 - x1c - x2c
    x2d = 1 - x1d - x3d

    # 防爆保护
    eps = 1e-10
    x1c = max(x1c, eps); x2c = max(x2c, eps); x3c = max(x3c, eps)
    x1d = max(x1d, eps); x2d = max(x2d, eps); x3d = max(x3d, eps)

    # g12 及其导数
    def calc_g12(u2):
        g = p1*u2**4 + p2*u2**3 + p3*u2**2 + p4*u2 + p5
        dg = 4*p1*u2**3 + 3*p2*u2**2 + 2*p3*u2 + p4
        d2g = 12*p1*u2**2 + 6*p2*u2 + 2*p3
        return g, dg, d2g

    u1c = x1c / (x1c + x2c)
    u2c = x2c / (x1c + x2c)
    u1d = x1d / (x1d + x2d)
    u2d = x2d / (x1d + x2d)

    g12c, dgdu2c, _ = calc_g12(u2c)
    g12d, dgdu2d, _ = calc_g12(u2d)

    g23c = g23
    g23d = g23
    dg23c = 0
    dg23d = 0
    dX13c = 0
    dX13d = 0

    mu1c = (np.log(x1c) + 1 - x1c - s*x2c - r*x3c +
            (g12c*x2c + iX13*x3c)*(x2c+x3c) -
            s*g23c*x2c*x3c - x2c*u1c*u2c*dgdu2c -
            x1c*x3c**2*dX13c - s*x2c*x3c**2*dg23c)

    mu1d = (np.log(x1d) + 1 - x1d - s*x2d - r*x3d +
            (g12d*x2d + iX13*x3d)*(x2d+x3d) -
            s*g23d*x2d*x3d - x2d*u1d*u2d*dgdu2d -
            x1d*x3d**2*dX13d - s*x2d*x3d**2*dg23d)

    mu2c = (s*np.log(x2c) + s - x1c - s*x2c - r*x3c +
            (g12c*x1c + g23c*s*x3c)*(x1c+x3c) -
            iX13*x1c*x3c + x1c*u1c*u2c*dgdu2c -
            x1c*x3c**2*dX13c - s*x2c*x3c**2*dg23c)

    mu2d = (s*np.log(x2d) + s - x1d - s*x2d - r*x3d +
            (g12d*x1d + g23d*s*x3d)*(x1d+x3d) -
            iX13*x1d*x3d + x1d*u1d*u2d*dgdu2d -
            x1d*x3d**2*dX13d - s*x2d*x3d**2*dg23d)

    mu3c = (r*np.log(x3c) + r - x1c - s*x2c - r*x3c +
            (iX13*x1c + s*g23c*x2c)*(x1c + x2c) -
            g12c*x1c*x2c + (x1c*dX13c + s*x2c*dg23c)*x3c*(x1c+x2c))

    mu3d = (r*np.log(x3d) + r - x1d - s*x2d - r*x3d +
            (iX13*x1d + s*g23d*x2d)*(x1d + x2d) -
            g12d*x1d*x2d + (x1d*dX13d + s*x2d*dg23d)*x3d*(x1d+x2d))

    return (mu1c - mu1d)**2 + (mu2c - mu2d)**2 + (mu3c - mu3d)**2


# ========================= 诊断 1：检查 g12 的取值范围 ==========================
print("=" * 60)
print("诊断 1：g12(u2) 的取值范围")
print("=" * 60)
u2_test = np.linspace(0, 1, 11)
for u2 in u2_test:
    g = p1*u2**4 + p2*u2**3 + p3*u2**2 + p4*u2 + p5
    print(f"  u2 = {u2:.2f}  ->  g12 = {g:+.4f}")

print("\n注意：g12 在所有组成下均为负值且绝对值很大（-1.6 ~ -5.8）。")
print("这可能导致 Gibbs 自由能表面异常平坦，使数值优化困难。\n")


# ========================= 诊断 2：逐点求解并观察成功率 ==========================
print("=" * 60)
print("诊断 2：Binodal 逐点求解")
print("=" * 60)

nloop = 60
ix3d_values = np.logspace(-5, np.log10(0.3), nloop)
x0_default = np.array([0.2, 0.5, 0.5])
x0 = x0_default.copy()

results = []
success_count = 0

for i, ix3d in enumerate(ix3d_values):
    # 约束：x1c + x2c <= 1 - 1e-6, x1d <= 1 - ix3d - 1e-6
    cons = [
        {'type': 'ineq', 'fun': lambda x: 1 - 1e-6 - x[0] - x[1]},  # x3c > 0
        {'type': 'ineq', 'fun': lambda x: 1 - ix3d - 1e-6 - x[2]},   # x2d > 0
        {'type': 'ineq', 'fun': lambda x: x[0] - 1e-6},
        {'type': 'ineq', 'fun': lambda x: x[1] - 1e-6},
        {'type': 'ineq', 'fun': lambda x: x[2] - 1e-6},
        {'type': 'ineq', 'fun': lambda x: 1 - 1e-6 - x[0]},
        {'type': 'ineq', 'fun': lambda x: 1 - 1e-6 - x[1]},
        {'type': 'ineq', 'fun': lambda x: 1 - 1e-6 - x[2]},
    ]

    res = minimize(lambda x: fun(x, ix3d, X13), x0, method='SLSQP',
                   constraints=cons, options={'ftol': 1e-10, 'maxiter': 10000, 'disp': False})

    x = res.x
    f = res.fun
    x2d_calc = 1 - x[2] - ix3d
    phase_diff = abs(x[0] - x[2]) + abs(x[1] - x2d_calc)

    is_success = f < 1e-7 and phase_diff > 1e-3
    if is_success:
        success_count += 1
        x0 = x.copy()  # 热启动
        results.append({
            'ix3d': ix3d, 'f': f, 'x1c': x[0], 'x2c': x[1], 'x3c': 1-x[0]-x[1],
            'x1d': x[2], 'x2d': x2d_calc, 'x3d': ix3d
        })
    else:
        x0 = x0_default.copy()

    status = "SUCCESS" if is_success else f"FAIL(f={f:.2e},pd={phase_diff:.2e})"
    print(f"  i={i+1:2d}, ix3d={ix3d:.4e} -> {status}")

print(f"\n成功率: {success_count}/{nloop} = {success_count/nloop*100:.1f}%")


# ========================= 诊断 3：绘制相图 ==========================
if len(results) > 0:
    print("\n" + "=" * 60)
    print("诊断 3：绘制相图")
    print("=" * 60)

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))

    # 浓相 arm
    x1c_list = [r['x1c'] for r in results]
    x3c_list = [r['x3c'] for r in results]
    # 稀相 arm
    x1d_list = [r['x1d'] for r in results]
    x3d_list = [r['x3d'] for r in results]

    ax.plot(x1c_list, x3c_list, 'bo-', label='Concentrated phase arm', markersize=4)
    ax.plot(x1d_list, x3d_list, 'ro-', label='Dilute phase arm', markersize=4)

    # 画 tie lines
    for r in results:
        ax.plot([r['x1c'], r['x1d']], [r['x3c'], r['x3d']], 'g--', alpha=0.4, linewidth=0.8)

    ax.set_xlabel('Volume fraction of L8-Bo (x1)')
    ax.set_ylabel('Volume fraction of PM6 (x3)')
    ax.set_title(f'Binodal Curve (success {success_count}/{nloop})')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.5)

    plt.tight_layout()
    plt.savefig('debug_binodal_OXy.png', dpi=200)
    print("相图已保存为 debug_binodal_OXy.png")
    plt.show()
else:
    print("\n警告：未找到任何成功点，无法绘制相图。")


# ========================= 诊断 4：参数敏感性测试 ==========================
print("\n" + "=" * 60)
print("诊断 4：参数敏感性 —— 测试不同 X13 下的首个 ix3d 收敛情况")
print("=" * 60)

ix3d_test = 1e-5
x0_test = np.array([0.2, 0.5, 0.5])

for X13_test in [0.40, 0.45, 0.447, 0.50, 0.55, 0.60, 0.70, 1.0, 1.7]:
    cons = [
        {'type': 'ineq', 'fun': lambda x: 1 - 1e-6 - x[0] - x[1]},
        {'type': 'ineq', 'fun': lambda x: 1 - ix3d_test - 1e-6 - x[2]},
        {'type': 'ineq', 'fun': lambda x: x[0] - 1e-6},
        {'type': 'ineq', 'fun': lambda x: x[1] - 1e-6},
        {'type': 'ineq', 'fun': lambda x: x[2] - 1e-6},
        {'type': 'ineq', 'fun': lambda x: 1 - 1e-6 - x[0]},
        {'type': 'ineq', 'fun': lambda x: 1 - 1e-6 - x[1]},
        {'type': 'ineq', 'fun': lambda x: 1 - 1e-6 - x[2]},
    ]
    res = minimize(lambda x: fun(x, ix3d_test, X13_test), x0_test, method='SLSQP',
                   constraints=cons, options={'ftol': 1e-10, 'maxiter': 10000, 'disp': False})
    f = res.fun
    x = res.x
    x2d = 1 - x[2] - ix3d_test
    pd = abs(x[0] - x[2]) + abs(x[1] - x2d)
    status = "OK" if f < 1e-7 and pd > 1e-3 else "FAIL"
    print(f"  X13={X13_test:.3f}: f={f:.2e}, pd={pd:.2e}, x=[{x[0]:.3f},{x[1]:.3f},{x[2]:.3f}] -> {status}")

print("\n建议：若 X13 过大导致相分离区过大，或过小导致无法分离，可据此调整。")
