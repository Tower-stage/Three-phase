function f = fun(x, ix3d, iX13)
% PM6 - L8-Bo - Toluene 三元体系 Binodal 计算函数
% 基于 Flory-Huggins 理论，改编自 Gemini 重构版

x3d = ix3d;
x1c = x(1);
x2c = x(2);
x3c = 1 - x1c - x2c;
x1d = x(3);
x2d = 1 - x1d - x3d;

% 【数值防爆保护】防止 log(负数) 导致复数报错
x1c = max(x1c, 1e-10); x2c = max(x2c, 1e-10); x3c = max(x3c, 1e-10);
x1d = max(x1d, 1e-10); x2d = max(x2d, 1e-10); x3d = max(x3d, 1e-10);

% ======== 物理参数区（PM6-L8Bo-Tol 体系）========
v1 = 1132.1;      % L8-Bo 摩尔体积 [cm^3/mol] (MW~1478, rho~1.3)
v2 = 106.3;       % Toluene 摩尔体积
v3 = 1743900;     % PM6 摩尔体积 (1500 * 1162.6，对应聚合度~1500)
s = v1/v2;
r = v1/v3;

% --- L8-Bo(1) - PM6(3) 相互作用参数 ---
X13c = iX13;
X13d = iX13;
dX13c = 0;
dX13d = 0;

% --- Toluene(2) - PM6(3) 相互作用参数 (g23) ---
kk = 0;
bb = 0.3852;      % PM6/Toluene, 源自 UNIFAC 估算
g23c = kk*x3c + bb;
g23d = kk*x3d + bb;
dg23c = kk;
dg23d = kk;

% --- L8-Bo(1) - Toluene(2) 相互作用参数 (g12, p1-p5 多项式) ---
% 基于 HSP 估算重新计算：chi ~ 0.47，考虑溶解度限制后上调至平均 0.50
% 假设 g12 随溶剂浓度增加线性降低（从 0.60 到 0.40）
p1 = 0;
p2 = 0;
p3 = 0;
p4 = -0.2000;
p5 = 0.6000;

u1c = x1c/(x1c + x2c);
u2c = x2c/(x1c + x2c);
u1d = x1d/(x1d + x2d);
u2d = x2d/(x1d + x2d);

g12c = p1*u2c^4 + p2*u2c^3 + p3*u2c^2 + p4*u2c + p5;
g12d = p1*u2d^4 + p2*u2d^3 + p3*u2d^2 + p4*u2d + p5;
dgdu2c = 4*p1*u2c^3 + 3*p2*u2c^2 + 2*p3*u2c + p4;
dgdu2d = 4*p1*u2d^3 + 3*p2*u2d^2 + 2*p3*u2d + p4;

% ======== 化学势相等方程 ========
mu1c = log(x1c) + 1 - x1c - s*x2c - r*x3c + ...
    (g12c*x2c + X13c*x3c)*(x2c+x3c) - ...
    s*g23c*x2c*x3c - x2c*u1c*u2c*dgdu2c - ...
    x1c*x3c^2*dX13c - s*x2c*x3c^2*dg23c;

mu1d = log(x1d) + 1 - x1d - s*x2d - r*x3d + ...
    (g12d*x2d + X13d*x3d)*(x2d+x3d) - ...
    s*g23d*x2d*x3d - x2d*u1d*u2d*dgdu2d - ...
    x1d*x3d^2*dX13d - s*x2d*x3d^2*dg23d;

mu2c = s*log(x2c) + s - x1c - s*x2c - r*x3c ...
    + (g12c*x1c + g23c*s*x3c)*(x1c+x3c) - ...
    X13c*x1c*x3c + x1c*u1c*u2c*dgdu2c - ...
    x1c*x3c^2*dX13c - s*x2c*x3c^2*dg23c;

mu2d = s*log(x2d) + s - x1d - s*x2d - r*x3d ...
    + (g12d*x1d + g23d*s*x3d)*(x1d+x3d) - ...
    X13d*x1d*x3d + x1d*u1d*u2d*dgdu2d - ...
    x1d*x3d^2*dX13d - s*x2d*x3d^2*dg23d;

mu3c = r*log(x3c) + r - x1c - s*x2c - r*x3c ...
    + (X13c*x1c + s*g23c*x2c)*(x1c + x2c) - ...
    g12c*x1c*x2c + (x1c*dX13c + s*x2c*dg23c)*x3c*(x1c+x2c);

mu3d = r*log(x3d) + r - x1d - s*x2d - r*x3d ...
    + (X13d*x1d + s*g23d*x2d)*(x1d + x2d) - ...
    g12d*x1d*x2d + (x1d*dX13d + s*x2d*dg23d)*x3d*(x1d+x2d);

% 目标函数：化学势差值平方和
f = (mu1c - mu1d)^2 + (mu2c - mu2d)^2 + (mu3c - mu3d)^2;
end
