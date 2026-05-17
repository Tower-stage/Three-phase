clear; clc;

%% 临界点计算脚本（PM6 - L8-Bo - o-Xylene 体系）
% 策略：先计算 Spinodal 曲线，再在 Spinodal 上寻找三阶条件变号点

%% 1. 物理参数（与 binodal/spinodal 保持一致）
v1 = 1132.1;       % L8-Bo 摩尔体积 [cm^3/mol]
v2 = 120.6;        % o-Xylene 摩尔体积
v3 = 1743900;      % PM6 摩尔体积 (1500 * 1162.6)
s = v1/v2;
r = v1/v3;

X13 = 0.55;        % L8-Bo/PM6 相互作用参数
                   % 调试提示：若相图异常，可尝试 0.447 (spinodal-gemini 中使用的值)
g23 = 0.4120;      % o-Xylene/PM6

% L8-Bo / o-Xylene 多项式系数 (g12)
p1 = 16.6738; p2 = -40.8682; p3 = 39.6207; p4 = -19.3621; p5 = -1.5993;

%% 2. 生成网格并计算 Spinodal 判别式 det(G'') = 0
x1_val = linspace(1e-8, 0.999, 800);
x3_val = linspace(1e-8, 0.999, 800);
[X1, X3] = meshgrid(x1_val, x3_val);
X2 = 1 - X1 - X3;

% 屏蔽非法区域
invalid = (X1 <= 1e-8) | (X2 <= 1e-8) | (X3 <= 1e-8);
X1(invalid) = 0.33; X2(invalid) = 0.33; X3(invalid) = 0.34;

% 计算 g12 及其导数
u1 = X1 ./ (X1 + X2);
u2 = X2 ./ (X1 + X2);
g12 = p1.*u2.^4 + p2.*u2.^3 + p3.*u2.^2 + p4.*u2 + p5;
dg12du = 4*p1.*u2.^3 + 3*p2.*u2.^2 + 2*p3.*u2 + p4;
d2g12du = 12*p1.*u2.^2 + 6*p2.*u2 + 2*p3;

% Gibbs 自由能二阶导数矩阵元素
G22 = 1./X1 + s./X2 - 2.*g12 + 2.*(1-2.*u2).*dg12du + u1.*u2.*d2g12du;
G23 = 1./X1 - (g12+X13) + s.*g23 + u2.*(1-3.*u2).*dg12du + u1.*u2.^2.*d2g12du;
G33 = 1./X1 + r./X3 - 2.*X13 - 2.*u2.^3.*dg12du + u2.^3.*u1.*d2g12du;

detG = G22.*G33 - G23.^2;
detG(invalid) = NaN;

%% 3. 提取 Spinodal 曲线 (det = 0)
fig = figure('Visible','off');
[C_spin, ~] = contour(X1, X3, detG, [0 0], 'b-', 'LineWidth', 1);
close(fig);

num_points = C_spin(2,1);
x1_spin = C_spin(1, 2:1+num_points);
x3_spin = C_spin(2, 2:1+num_points);
x2_spin = 1 - x1_spin - x3_spin;

%% 4. 在 Spinodal 曲线上计算三阶条件
crit_val = zeros(1, num_points);

for i = 1:num_points
    x1 = x1_spin(i); x3 = x3_spin(i); x2 = x2_spin(i);
    if x1 <= 1e-8 || x2 <= 1e-8 || x3 <= 1e-8
        crit_val(i) = NaN;
        continue;
    end

    u1 = x1/(x1+x2); u2 = x2/(x1+x2);
    g12 = p1*u2^4 + p2*u2^3 + p3*u2^2 + p4*u2 + p5;
    dg12du = 4*p1*u2^3 + 3*p2*u2^2 + 2*p3*u2 + p4;
    d2g12du = 12*p1*u2^2 + 6*p2*u2 + 2*p3;

    G22 = 1/x1 + s/x2 - 2*g12 + 2*(1-2*u2)*dg12du + u1*u2*d2g12du;
    G23 = 1/x1 - (g12+X13) + s*g23 + u2*(1-3*u2)*dg12du + u1*u2^2*d2g12du;
    G33 = 1/x1 + r/x3 - 2*X13 - 2*u2^3*dg12du + u2^3*u1*d2g12du;

    G222 = 1/x1^2 - s/x2^2 - 6*u2/x2*dg12du + (3-6*u2)*u2/x2*d2g12du;
    G223 = 1/x1^2 - 6*u2^2/x2*dg12du + 3*(1-2*u2)*u2^2/x2*d2g12du;
    G233 = 1/x1^2 - 6*u2^3/x2*dg12du + (3*u2-6*u2^2)*u2^2/x2*d2g12du;
    G333 = 1/x1^2 - r/x3^2 + 6*u2^4/x2*dg12du + (3*u2^2-2*u2^3)*u2^2/x2*d2g12du;

    crit_val(i) = G222*G33^2 - 3*G223*G23*G33 + 3*G233*G23^2 - G22*G23*G333;
end

%% 5. 寻找 crit_val 变号点并插值得到临界点
valid_idx = ~isnan(crit_val);
x1_v = x1_spin(valid_idx);
x3_v = x3_spin(valid_idx);
crit_v = crit_val(valid_idx);

sign_changes = find(diff(sign(crit_v)) ~= 0);

if isempty(sign_changes)
    warning('未在 Spinodal 曲线上找到临界点（crit 未变号）。请检查 X13 是否合理，或增大网格分辨率。');
    fprintf('crit_val 范围: %.4e ~ %.4e\n', min(crit_v), max(crit_v));
    return;
end

idx = sign_changes(1);
t = abs(crit_v(idx)) / (abs(crit_v(idx)) + abs(crit_v(idx+1)));
x1_crit = x1_v(idx) + t*(x1_v(idx+1) - x1_v(idx));
x3_crit = x3_v(idx) + t*(x3_v(idx+1) - x3_v(idx));
x2_crit = 1 - x1_crit - x3_crit;

fprintf('========================================\n');
fprintf('PM6-L8Bo-o-Xylene 体系临界点计算完成\n');
fprintf('========================================\n');
fprintf('  phi1 (L8-Bo)    = %.6f\n', x1_crit);
fprintf('  phi2 (o-Xylene) = %.6f\n', x2_crit);
fprintf('  phi3 (PM6)      = %.6f\n', x3_crit);
fprintf('========================================\n');

%% 6. 保存结果
crit_result = table(x1_crit, x2_crit, x3_crit, ...
    'VariableNames', {'phi1_L8Bo','phi2_oXy','phi3_PM6'});
writetable(crit_result, 'critical_point_OXy.csv');
fprintf('结果已保存至 critical_point_OXy.csv\n');
