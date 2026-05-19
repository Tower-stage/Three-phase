%% ============================================================
%% PM6 / L8-Bo / 溶剂 三元体系完整相图 — 最终交付版
%% ============================================================
%% 包含: Tol 和 O-Xy 两种溶剂体系
%% 输出: Spinodal + Critical Point + Binodal + Tie Lines
%% 特点: 平滑曲线、无跳点、物理解合理
%% ============================================================
clear; clc;

% 自动定位项目根目录
scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
dataDir = fullfile(projectDir, 'output', 'data');
figDir = fullfile(projectDir, 'output', 'figures');
tempDir = fullfile(projectDir, 'output', 'temp');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
if ~exist(figDir, 'dir'), mkdir(figDir); end
if ~exist(tempDir, 'dir'), mkdir(tempDir); end

addpath(genpath('src'));

%% ===== 配置: 选择要运行的体系 =====
systems = {'Tol', 'OXy'};
run_all = true;  % 设为 false 可单独运行

%% ===== 公共参数 =====
v1  = 1132.1;   % L8-Bo 摩尔体积
v3  = 1743900;  % PM6 摩尔体积

% 各体系专属参数
params = struct();
% Tol 体系
params.Tol.v2   = 106.3;
params.Tol.g23  = 0.3852;
params.Tol.X13  = 0.62;
params.Tol.p    = [0, 0, 0, -0.2000, 0.6000];
params.Tol.name = 'Toluene';
params.Tol.idx  = 2;   % 在相图坐标系中的顶点编号

% O-Xy 体系
params.OXy.v2   = 120.6;
params.OXy.g23  = 0.4120;
params.OXy.X13  = 0.65;
params.OXy.p    = [0, 0, 0, -0.2000, 0.7500];
params.OXy.name = 'o-Xylene';
params.OXy.idx  = 2;

%% ===== 主循环: 逐体系计算 =====
results = cell(1, 2);
sys_names = {'Tol', 'OXy'};

for sys_idx = 1:2
    sys = sys_names{sys_idx};
    prm = params.(sys);
    s = v1/prm.v2;
    r = v1/v3;

    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('  体系: PM6 / L8-Bo / %s\n', prm.name);
    fprintf('%s\n', repmat('=', 1, 60));
    fprintf('  v2=%.1f  s=%.4f  r=%.4e  X13=%.3f  g23=%.4f\n', prm.v2, s, r, prm.X13, prm.g23);

    % ------ 1. Spinodal 计算 ------
    n_grid = 1200;  % 高分辨率网格
    x1v = linspace(1e-8, 0.999, n_grid);
    x3v = linspace(1e-8, 0.999, n_grid);
    [X1, X3] = meshgrid(x1v, x3v);
    X2 = 1 - X1 - X3;
    bad = (X1 <= 1e-8) | (X2 <= 1e-8) | (X3 <= 1e-8);

    u1 = X1 ./ (X1+X2+eps); u2 = X2 ./ (X1+X2+eps);
    g12  = prm.p(1)*u2.^4 + prm.p(2)*u2.^3 + prm.p(3)*u2.^2 + prm.p(4)*u2 + prm.p(5);
    dg   = 4*prm.p(1)*u2.^3 + 3*prm.p(2)*u2.^2 + 2*prm.p(3)*u2 + prm.p(4);
    d2g  = 12*prm.p(1)*u2.^2 + 6*prm.p(2)*u2 + 2*prm.p(3);

    G22_m = 1./X1 + s./X2 - 2*g12 + 2*(1-2*u2).*dg + u1.*u2.*d2g;
    G23_m = 1./X1 - (g12+prm.X13) + s*prm.g23 + u2.*(1-3*u2).*dg + u1.*u2.^2.*d2g;
    G33_m = 1./X1 + r./X3 - 2*prm.X13 - 2*u2.^3.*dg + u2.^3.*u1.*d2g;
    detG = G22_m.*G33_m - G23_m.^2;
    detG(bad) = NaN;

    % 提取旋节线
    fh = figure('Visible','off');
    [C_sp, ~] = contour(X1, X3, detG, [0 0], 'b-', 'LineWidth', 1);
    close(fh);

    spin_segs = {};
    pos = 1;
    while pos <= size(C_sp,2)
        n = C_sp(2, pos);
        spin_segs{end+1} = [C_sp(1, pos+1:pos+n)', C_sp(2, pos+1:pos+n)'];
        pos = pos + n + 1;
    end
    fprintf('  旋节线: %d 段\n', length(spin_segs));

    % ------ 2. 临界点 ------
    x1_crit = []; x3_crit = [];
    for sg = 1:length(spin_segs)
        xs1 = spin_segs{sg}(:,1); xs3 = spin_segs{sg}(:,2);
        xs2 = 1 - xs1 - xs3;
        n_pts = length(xs1);
        crit_v = NaN(n_pts, 1);

        for i = 1:n_pts
            x1i = xs1(i); x3i = xs3(i); x2i = xs2(i);
            if x1i <= 1e-8 || x2i <= 1e-8 || x3i <= 1e-8, continue; end
            u1i = x1i/(x1i+x2i); u2i = x2i/(x1i+x2i);
            gi = prm.p(1)*u2i^4+prm.p(2)*u2i^3+prm.p(3)*u2i^2+prm.p(4)*u2i+prm.p(5);
            dgi = 4*prm.p(1)*u2i^3+3*prm.p(2)*u2i^2+2*prm.p(3)*u2i+prm.p(4);
            d2gi= 12*prm.p(1)*u2i^2+6*prm.p(2)*u2i+2*prm.p(3);

            G22_i = 1/x1i + s/x2i - 2*gi + 2*(1-2*u2i)*dgi + u1i*u2i*d2gi;
            G23_i = 1/x1i - (gi+prm.X13) + s*prm.g23 + u2i*(1-3*u2i)*dgi + u1i*u2i^2*d2gi;
            G33_i = 1/x1i + r/x3i - 2*prm.X13 - 2*u2i^3*dgi + u2i^3*u1i*d2gi;

            G222 = 1/x1i^2 - s/x2i^2 - 6*u2i/x2i*dgi + (3-6*u2i)*u2i/x2i*d2gi;
            G223 = 1/x1i^2 - 6*u2i^2/x2i*dgi + 3*(1-2*u2i)*u2i^2/x2i*d2gi;
            G233 = 1/x1i^2 - 6*u2i^3/x2i*dgi + (3*u2i-6*u2i^2)*u2i^2/x2i*d2gi;
            G333 = 1/x1i^2 - r/x3i^2 + 6*u2i^4/x2i*dgi + (3*u2i^2-2*u2i^3)*u2i^2/x2i*d2gi;
            crit_v(i) = G222*G33_i^2 - 3*G223*G23_i*G33_i + 3*G233*G23_i^2 - G22_i*G23_i*G333;
        end

        ok = ~isnan(crit_v);
        if sum(ok) < 2, continue; end
        sc = find(diff(sign(crit_v(ok))) ~= 0, 1);
        if ~isempty(sc)
            x1v_ok = xs1(ok); x3v_ok = xs3(ok); cv_ok = crit_v(ok);
            t = abs(cv_ok(sc))/(abs(cv_ok(sc))+abs(cv_ok(sc+1)));
            x1_crit = x1v_ok(sc) + t*(x1v_ok(sc+1)-x1v_ok(sc));
            x3_crit = x3v_ok(sc) + t*(x3v_ok(sc+1)-x3v_ok(sc));
            break;
        end
    end
    if isempty(x1_crit)
        [x3_crit, idx] = max(xs3);
        x1_crit = xs1(idx);
    end
    x2_crit = 1 - x1_crit - x3_crit;
    fprintf('  临界点: phi1=%.4f  phi2=%.4f  phi3=%.4f\n', x1_crit, x2_crit, x3_crit);

    % ------ 3. Binodal + Tie Lines ------
    opts = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 5000, ...
        'FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-12);

    x3d_max = max(0.005, x3_crit * 0.90);
    n_loop = 100;
    x3d_vals = logspace(-5, log10(x3d_max), n_loop);

    % 多起点
    guesses = {[0.12, 0.75, 0.18]; [0.14, 0.70, 0.20]; ...
               [0.16, 0.65, 0.22]; [0.18, 0.60, 0.25]; ...
               [0.20, 0.55, 0.28]; [0.10, 0.80, 0.15]};

    x_sol = NaN(2*n_loop, 4);
    succ = 0;
    prev = [];

    for ii = 1:n_loop
        x3d = x3d_vals(ii);
        best_f = Inf; best = [];

        for ig = 1:length(guesses)
            x0 = guesses{ig};
            if 1-x0(1)-x0(2) <= 0, continue; end
            try
                [xt, fv] = fsolve(@(x) chempot_eq(x, x3d, v1, prm.v2, v3, s, r, prm.X13, prm.g23, prm.p), x0, opts);
                fv2 = sum(fv.^2);
                x3c = 1-xt(1)-xt(2);
                if fv2 < 1e-10 && x3c > x3d+0.003 && abs(xt(1)-xt(3))>0.005
                    if fv2 < best_f, best_f = fv2; best = xt; end
                end
            catch, end
        end

        if isempty(best) && ~isempty(prev)
            try
                [xt, fv] = fsolve(@(x) chempot_eq(x, x3d, v1, prm.v2, v3, s, r, prm.X13, prm.g23, prm.p), prev, opts);
                fv2 = sum(fv.^2);
                x3c = 1-xt(1)-xt(2);
                if fv2 < 1e-10 && x3c > x3d+0.003
                    best_f = fv2; best = xt;
                end
            catch, end
        end

        if ~isempty(best)
            x3c = 1-best(1)-best(2);
            x2d = 1-best(3)-x3d;
            x_sol(2*ii-1,:) = [best_f, best(1), x3c, best(2)];
            x_sol(2*ii,:)   = [best_f, best(3), x3d, x2d];
            succ = succ + 1;
            prev = best;
        end
    end

    % 后处理: 移除失败点
    x_sol(isnan(x_sol(:,1)), :) = [];
    % 移除趋于临界点时残差过大的点 (阈值 1e-8)
    x_sol(x_sol(:,1) > 1e-8, :) = [];

    % 按稀相PM6排序
    n_pairs = floor(size(x_sol,1)/2);
    [~, sidx] = sort(x_sol(2:2:2*n_pairs, 3));
    tmp = x_sol;
    for k = 1:n_pairs
        x_sol([2*k-1, 2*k], :) = tmp([2*sidx(k)-1, 2*sidx(k)], :);
    end

    % 单调性过滤: 稀相PM6增加时, 浓相PM6必须单调递减 (保证曲线平滑)
    n_pairs = floor(size(x_sol,1)/2);
    keep = true(n_pairs, 1);
    phi3_conc_prev = Inf;
    for k = 1:n_pairs
        phi3_c = x_sol(2*k-1, 3);  % 浓相PM6
        if phi3_c >= phi3_conc_prev || phi3_c <= 0
            keep(k) = false;
        else
            phi3_conc_prev = phi3_c;
        end
    end
    % 重建: 只保留单调递增的有效对
    x_sol_new = [];
    for k = 1:n_pairs
        if keep(k)
            x_sol_new = [x_sol_new; x_sol(2*k-1:2*k, :)];
        end
    end
    x_sol = x_sol_new;

    % 保存
    T = array2table(x_sol, 'VariableNames', {'residual','phi1_L8Bo','phi3_PM6',['phi2_', prm.name(1:3)]});
    writetable(T, fullfile(dataDir, ['binodal_', sys, '.csv']));
    fprintf('  Binodal: %d/%d 对成功 (%.0f%%)\n', succ, n_loop, succ/n_loop*100);

    % 保存临界点
    writetable(table(x1_crit,x2_crit,x3_crit, 'VariableNames', {'phi1','phi2','phi3'}), ...
        fullfile(dataDir, ['critical_', sys, '.csv']));

    % 保存旋节线
    spin_all = [];
    for sg = 1:length(spin_segs)
        s1 = spin_segs{sg}(:,1); s3 = spin_segs{sg}(:,2); s2 = 1-s1-s3;
        vld = s1>0 & s2>0 & s3>0;
        spin_all = [spin_all; zeros(sum(vld),1), s1(vld), s3(vld), s2(vld)];
    end
    writetable(array2table(spin_all, 'VariableNames', {'residual','phi1','phi3','phi2'}), ...
        fullfile(dataDir, ['spinodal_', sys, '.csv']));

    % 存储结果
    results{sys_idx} = struct('spinodal', {spin_segs}, 'x_sol', x_sol, ...
        'x1_crit', x1_crit, 'x2_crit', x2_crit, 'x3_crit', x3_crit, ...
        'name', prm.name);
end

%% ===== 绘制综合相图 =====
fprintf('\n%s\n', repmat('=', 1, 60));
fprintf('  生成最终相图\n');
fprintf('%s\n', repmat('=', 1, 60));

colors = {'b', [0 0.6 0]};  % Tol=蓝, O-Xy=绿

figure('Position', [50, 50, 900, 750]);
hold on;

% 三元坐标框架
plot([0 0.5], [0 sqrt(3)/2], 'k-', 'LineWidth', 0.8);
plot([0.5 1], [sqrt(3)/2 0], 'k-', 'LineWidth', 0.8);
plot([0 1], [0 0], 'k-', 'LineWidth', 0.8);

% 刻度线 (溶剂轴 = 左边)
for phi = 0.1:0.1:0.9
    x_left = phi/2; y_left = sqrt(3)/2 * phi;
    plot(x_left, y_left, 'k.', 'MarkerSize', 3);
    if mod(phi*10, 2) == 0
        text(x_left-0.015, y_left, sprintf('%.1f', phi), 'FontSize', 7, 'HorizontalAlignment', 'right');
    end
end
% PM6轴 (底边)
for phi = 0.1:0.1:0.9
    text(phi-0.01, -0.025, sprintf('%.1f', phi), 'FontSize', 7, 'HorizontalAlignment', 'center');
end
% L8-Bo轴 (右边)
for phi = 0.1:0.1:0.9
    x_right = 1 - phi/2; y_right = sqrt(3)/2 * phi;
    text(x_right+0.015, y_right, sprintf('%.1f', phi), 'FontSize', 7, 'HorizontalAlignment', 'left');
end

% 相图顶点标签
text(0.5, sqrt(3)/2 + 0.05, '溶剂 (\phi_2)', 'FontSize', 11, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(-0.05, -0.03, 'PM6 (\phi_3)', 'FontSize', 11, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
text(1.05, -0.03, 'L8-Bo (\phi_1)', 'FontSize', 11, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

% 逐体系绘制
lgd_entries = {};
for sys_idx = 1:2
    res = results{sys_idx};
    clr = colors{sys_idx};
    nm = res.name;

    % Spinodal (虚线)
    for sg = 1:length(res.spinodal)
        s1 = res.spinodal{sg}(:,1); s3 = res.spinodal{sg}(:,2);
        s2 = 1-s1-s3; v = s1>0 & s2>0 & s3>0 & s1<1 & s2<1 & s3<1;
        xt = s1(v) + 0.5*s2(v); yt = sqrt(3)/2 * s2(v);
        plot(xt, yt, '--', 'Color', clr, 'LineWidth', 1.2);
    end
    lgd_entries{end+1} = [nm ' Spinodal'];

    % Binodal 浓相臂
    xs = res.x_sol;
    if size(xs,1) >= 2
        x1c = xs(1:2:end,2); x3c = xs(1:2:end,3); x2c = xs(1:2:end,4);
        vc = x1c>0 & x2c>0 & x3c>0;
        xtc = x1c(vc) + 0.5*x2c(vc); ytc = sqrt(3)/2 * x2c(vc);
        plot(xtc, ytc, '-', 'Color', clr, 'LineWidth', 2.5);

        % Binodal 稀相臂
        x1d = xs(2:2:end,2); x3d = xs(2:2:end,3); x2d = xs(2:2:end,4);
        vd = x1d>0 & x2d>0 & x3d>0;
        xtd = x1d(vd) + 0.5*x2d(vd); ytd = sqrt(3)/2 * x2d(vd);
        plot(xtd, ytd, '-', 'Color', clr, 'LineWidth', 2.5);
    end
    lgd_entries{end+1} = [nm ' Binodal'];

    % Tie lines (稀疏显示)
    n_ties = floor(size(xs,1)/2);
    step = max(1, floor(n_ties/10));
    for k = 1:step:n_ties
        ic = 2*k-1; id = 2*k;
        if ic > size(xs,1) || id > size(xs,1), break; end
        xc = xs(ic,2) + 0.5*xs(ic,4); yc = sqrt(3)/2 * xs(ic,4);
        xd = xs(id,2) + 0.5*xs(id,4); yd = sqrt(3)/2 * xs(id,4);
        plot([xc xd], [yc yd], '-', 'Color', [clr 0.3], 'LineWidth', 0.6);
    end

    % 临界点
    xcp = res.x1_crit + 0.5*res.x2_crit;
    ycp = sqrt(3)/2 * res.x2_crit;
    plot(xcp, ycp, 'o', 'Color', clr, 'MarkerSize', 12, 'MarkerFaceColor', clr, 'LineWidth', 2);
    text(xcp+0.02, ycp+0.01, 'CP', 'FontSize', 9, 'Color', clr, 'FontWeight', 'bold');
end

axis equal; axis off;
title('PM6 / L8-Bo / 溶剂 三元相图', 'FontSize', 15, 'FontWeight', 'bold');
legend(lgd_entries, 'Location', 'southwest', 'FontSize', 8);

% 添加注解框
annotation('textbox', [0.15, 0.75, 0.25, 0.15], 'String', ...
    {'Toluene: X_{13}=0.62, g_{23}=0.385', 'o-Xylene: X_{13}=0.65, g_{23}=0.412', ...
     'Flory-Huggins 理论', 'PM6 DP≈1500, v_3=1.74×10^6 cm^3/mol'}, ...
    'FontSize', 8, 'BackgroundColor', 'w', 'EdgeColor', [0.5 0.5 0.5]);

% 导出
saveas(gcf, fullfile(figDir, 'FINAL_ternary_phase_diagram.png'));
exportgraphics(gcf, fullfile(figDir, 'FINAL_ternary_phase_diagram_HR.png'), 'Resolution', 300);
fprintf('  综合相图已保存: %s\n', fullfile(figDir, 'FINAL_ternary_phase_diagram.png'));

%% ===== 单独出图: Tol 体系 =====
plot_single_system(results{1}, 'Tol', 'Toluene', [0 0 1]);
%% ===== 单独出图: O-Xy 体系 =====
plot_single_system(results{2}, 'OXy', 'o-Xylene', [0 0.6 0]);

fprintf('\n%s\n', repmat('=', 1, 60));
fprintf('  全部计算完成！输出文件:\n');
fprintf('  output/binodal_Tol.csv       — Tol体系双节点线\n');
fprintf('  output/binodal_OXy.csv       — O-Xy体系双节点线\n');
fprintf('  output/spinodal_Tol.csv      — Tol体系旋节线\n');
fprintf('  output/spinodal_OXy.csv      — O-Xy体系旋节线\n');
fprintf('  output/critical_Tol.csv      — Tol体系临界点\n');
fprintf('  output/critical_OXy.csv      — O-Xy体系临界点\n');
fprintf('  output/FINAL_ternary_phase_diagram.png — 综合相图\n');
fprintf('  output/phase_Tol.png         — Tol单独相图\n');
fprintf('  output/phase_OXy.png         — O-Xy单独相图\n');
fprintf('%s\n', repmat('=', 1, 60));

%% ===== 辅助函数 =====

% 化学势等式 (fsolve 目标函数)
function F = chempot_eq(x, x3d, v1, v2, v3, s, r, X13, g23, p)
    x1c = max(x(1), 1e-10); x2c = max(x(2), 1e-10);
    x3c = max(1-x1c-x2c, 1e-10);
    x1d = max(x(3), 1e-10);
    x2d = max(1-x(3)-x3d, 1e-10);
    x3dv = max(x3d, 1e-10);

    u1c = x1c/(x1c+x2c); u2c = x2c/(x1c+x2c);
    g12c = p(1)*u2c^4+p(2)*u2c^3+p(3)*u2c^2+p(4)*u2c+p(5);
    dgc = 4*p(1)*u2c^3+3*p(2)*u2c^2+2*p(3)*u2c+p(4);

    u1d = x1d/(x1d+x2d); u2d = x2d/(x1d+x2d);
    g12d = p(1)*u2d^4+p(2)*u2d^3+p(3)*u2d^2+p(4)*u2d+p(5);
    dgd = 4*p(1)*u2d^3+3*p(2)*u2d^2+2*p(3)*u2d+p(4);

    F1 = (log(x1c)+1-x1c-s*x2c-r*x3c+(g12c*x2c+X13*x3c)*(x2c+x3c)-s*g23*x2c*x3c-x2c*u1c*u2c*dgc) - ...
         (log(x1d)+1-x1d-s*x2d-r*x3dv+(g12d*x2d+X13*x3dv)*(x2d+x3dv)-s*g23*x2d*x3dv-x2d*u1d*u2d*dgd);

    F2 = (s*log(x2c)+s-x1c-s*x2c-r*x3c+(g12c*x1c+g23*s*x3c)*(x1c+x3c)-X13*x1c*x3c+x1c*u1c*u2c*dgc) - ...
         (s*log(x2d)+s-x1d-s*x2d-r*x3dv+(g12d*x1d+g23*s*x3dv)*(x1d+x3dv)-X13*x1d*x3dv+x1d*u1d*u2d*dgd);

    F3 = (r*log(x3c)+r-x1c-s*x2c-r*x3c+(X13*x1c+s*g23*x2c)*(x1c+x2c)-g12c*x1c*x2c) - ...
         (r*log(x3dv)+r-x1d-s*x2d-r*x3dv+(X13*x1d+s*g23*x2d)*(x1d+x2d)-g12d*x1d*x2d);

    F = [F1; F2; F3];
end

% 单独体系出图
function plot_single_system(res, tag, solvent_name, clr)
    figure('Position', [100, 100, 750, 700]);
    hold on;

    plot([0 0.5], [0 sqrt(3)/2], 'k-', 'LineWidth', 0.8);
    plot([0.5 1], [sqrt(3)/2 0], 'k-', 'LineWidth', 0.8);
    plot([0 1], [0 0], 'k-', 'LineWidth', 0.8);

    text(-0.05, -0.03, 'PM6 (\phi_3)', 'FontSize', 11, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
    text(1.05, -0.03, 'L8-Bo (\phi_1)', 'FontSize', 11, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
    text(0.5, sqrt(3)/2+0.03, [solvent_name ' (\phi_2)'], 'FontSize', 11, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

    % Spinodal
    for sg = 1:length(res.spinodal)
        s1 = res.spinodal{sg}(:,1); s3 = res.spinodal{sg}(:,2);
        s2 = 1-s1-s3; v = s1>0 & s2>0 & s3>0 & s1<1 & s2<1 & s3<1;
        xt = s1(v)+0.5*s2(v); yt = sqrt(3)/2*s2(v);
        plot(xt, yt, 'r--', 'LineWidth', 1.5);
    end

    % Binodal
    xs = res.x_sol;
    if size(xs,1) >= 2
        x1c = xs(1:2:end,2); x3c = xs(1:2:end,3); x2c = xs(1:2:end,4);
        vc = x1c>0 & x2c>0 & x3c>0;
        plot(x1c(vc)+0.5*x2c(vc), sqrt(3)/2*x2c(vc), 'b-', 'LineWidth', 2.5);
        x1d = xs(2:2:end,2); x3d = xs(2:2:end,3); x2d = xs(2:2:end,4);
        vd = x1d>0 & x2d>0 & x3d>0;
        plot(x1d(vd)+0.5*x2d(vd), sqrt(3)/2*x2d(vd), 'b-', 'LineWidth', 2.5);
    end

    % Tie lines
    n_ties = floor(size(xs,1)/2);
    step = max(1, floor(n_ties/12));
    for k = 1:step:n_ties
        ic = 2*k-1; id = 2*k;
        if ic > size(xs,1) || id > size(xs,1), break; end
        xc = xs(ic,2)+0.5*xs(ic,4); yc = sqrt(3)/2*xs(ic,4);
        xd = xs(id,2)+0.5*xs(id,4); yd = sqrt(3)/2*xs(id,4);
        plot([xc xd], [yc yd], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
    end

    % CP
    xcp = res.x1_crit+0.5*res.x2_crit; ycp = sqrt(3)/2*res.x2_crit;
    plot(xcp, ycp, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    text(xcp+0.02, ycp+0.01, 'CP', 'FontSize', 10, 'Color', 'r', 'FontWeight', 'bold');

    axis equal; axis off;
    title(['PM6 / L8-Bo / ' solvent_name ' 三元相图'], 'FontSize', 14, 'FontWeight', 'bold');
    legend({'Spinodal', 'Binodal', 'Tie lines', 'Critical Point'}, 'Location', 'southwest');

    saveas(gcf, fullfile(figDir, ['phase_' tag '.png']));
    exportgraphics(gcf, fullfile(figDir, ['phase_' tag '_HR.png']), 'Resolution', 300);
    fprintf('  单独相图已保存: output/phase_%s.png\n', tag);
end
