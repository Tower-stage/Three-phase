%% 从保存的扫描结果生成最终图表
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

load(fullfile(tempDir, 'param_sweep_v2.mat'));

v3_values = data2.v3_values;
X13_values = data2.X13_values;
g23_values = data2.g23_values;
succ_rate = data2.succ_rate;
shape_score = data2.shape_score;
x3_crit_v = data2.x3_crit_v;
x3c_range = data2.x3c_range;
nv = length(v3_values); nx = length(X13_values); ng = length(g23_values);

fprintf('已加载 %d 组合的参数扫描结果\n', nv*nx*ng);

%% ===== 找出最佳组合 (修正NaN排序) =====
fprintf('\n===== TOP 20 参数组合 (形状评分) =====\n');
fprintf('Rank  v3       X13   g23   Succ%%  Score   CP_phi3  Delta\n');
all_scores = shape_score(:);
all_scores(isnan(all_scores)) = -Inf;
[~, sidx] = sort(all_scores, 'descend');
rank = 0;
for k = 1:length(sidx)
    if all_scores(sidx(k)) <= 0, continue; end
    rank = rank + 1;
    if rank > 20, break; end
    [iv, ix, ig] = ind2sub([nv, nx, ng], sidx(k));
    fprintf('%2d    %-7.0f  %.2f  %.3f  %5.1f  %6.4f  %6.4f  %6.4f\n', ...
        rank, v3_values(iv), X13_values(ix), g23_values(ig), ...
        succ_rate(iv,ix,ig), shape_score(iv,ix,ig), ...
        x3_crit_v(iv,ix,ig), x3c_range(iv,ix,ig));
end

%% ===== 热力图: g23 效应 =====
figure('Position', [50, 50, 1400, 900]);
for ig = 1:ng
    subplot(2, 2, ig);
    score_slice = shape_score(:, :, ig);
    score_slice(isnan(score_slice)) = 0;
    imagesc(X13_values, log10(v3_values), score_slice);
    set(gca, 'YDir', 'normal');
    colormap(jet); colorbar;
    xlabel('X_{13} (L8-Bo/PM6)'); ylabel('log_{10}(v_3)');
    title(sprintf('Shape Score (g_{23}=%.3f)', g23_values(ig)));

    % 标注成功率
    for iv = 1:nv
        for ix = 1:nx
            sr = succ_rate(iv, ix, ig);
            if sr > 0
                text(X13_values(ix), log10(v3_values(iv)), ...
                    sprintf('%.0f', sr), 'HorizontalAlignment', 'center', ...
                    'FontSize', 7, 'Color', 'w');
            end
        end
    end
end
saveas(gcf, fullfile(figDir, 'heatmap_g23_effect.png'));
exportgraphics(gcf, fullfile(figDir, 'heatmap_g23_effect_HR.png'), 'Resolution', 300);
fprintf('g23效应热力图已保存\n');

%% ===== 挑选代表性组合生成相图 =====
% 3个代表性案例:
% 1) 全局最佳: v3=100000 X13=0.95 g23=0.385
% 2) 平衡型: v3=60000 X13=0.85 g23=0.55
% 3) 对比: 原始参数 v3=1743900 X13=0.62 g23=0.385

cases = {};

% Case 1: 全局最佳
[~, idx_best] = max(all_scores);
[iv1, ix1, ig1] = ind2sub([nv, nx, ng], idx_best);
cases{end+1} = [iv1, ix1, ig1];

% Case 2: 中v3最佳 (v3=20000-60000)
best_mid = -Inf; bi_mid = [];
for iv = 4:6
    for ix = 1:nx
        for ig = 1:ng
            sc = shape_score(iv,ix,ig);
            if ~isnan(sc) && sc > best_mid
                best_mid = sc; bi_mid = [iv, ix, ig];
            end
        end
    end
end
if ~isempty(bi_mid), cases{end+1} = bi_mid; end

% Case 3: 低v3还能成功的 (v3=20000)
best_low = -Inf; bi_low = [];
for iv = 4:4
    for ix = 1:nx
        for ig = 1:ng
            sc = shape_score(iv,ix,ig);
            if ~isnan(sc) && sc > best_low
                best_low = sc; bi_low = [iv, ix, ig];
            end
        end
    end
end
if ~isempty(bi_low), cases{end+1} = bi_low; end

fprintf('\n===== 代表性案例相图 =====\n');
for ci = 1:length(cases)
    iv = cases{ci}(1); ix = cases{ci}(2); ig = cases{ci}(3);
    res = all_results{iv, ix, ig};
    v3o = v3_values(iv); X13o = X13_values(ix); g23o = g23_values(ig);
    sc = shape_score(iv,ix,ig); sr = succ_rate(iv,ix,ig);

    fprintf('Case %d: v3=%.0f X13=%.2f g23=%.3f succ=%.0f%% score=%.4f\n', ...
        ci, v3o, X13o, g23o, sr, sc);

    xs = res.x_sol;
    if size(xs,1) < 4, continue; end

    c1 = xs(1:2:end,2); c3 = xs(1:2:end,3); c2 = xs(1:2:end,4);
    d1 = xs(2:2:end,2); d3 = xs(2:2:end,3); d2 = xs(2:2:end,4);
    rc = c1./(c1+c3); rd = d1./(d1+d3);

    nt = floor(size(xs,1)/2);
    step_ = max(1, floor(nt/10));

    % --- 子图1: 三元相图 ---
    figure('Position', [50+ci*300, 50, 600, 550]);
    hold on;
    plot([0 0.5],[0 sqrt(3)/2],'k-',[0.5 1],[sqrt(3)/2 0],'k-',[0 1],[0 0],'k-','LineWidth',0.8);

    for sg = 1:length(res.spin_segs)
        s1 = res.spin_segs{sg}(:,1); s3 = res.spin_segs{sg}(:,2);
        s2 = 1-s1-s3; v = s1>0 & s2>0 & s3>0 & s1<1 & s2<1 & s3<1;
        plot(s1(v)+0.5*s2(v), sqrt(3)/2*s2(v), 'r--', 'LineWidth', 1.2);
    end

    vc = c1>0 & c2>0 & c3>0;
    plot(c1(vc)+0.5*c2(vc), sqrt(3)/2*c2(vc), 'b-', 'LineWidth', 2.5);
    vd = d1>0 & d2>0 & d3>0;
    plot(d1(vd)+0.5*d2(vd), sqrt(3)/2*d2(vd), 'b-', 'LineWidth', 2.5);

    for k = 1:step_:nt
        ic = 2*k-1; id = 2*k;
        if ic>size(xs,1)||id>size(xs,1), break; end
        plot([xs(ic,2)+0.5*xs(ic,4) xs(id,2)+0.5*xs(id,4)], ...
             [sqrt(3)/2*xs(ic,4) sqrt(3)/2*xs(id,4)], ...
             'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    end

    if ~isnan(res.x1_crit)
        plot(res.x1_crit+0.5*(1-res.x1_crit-res.x3_crit), ...
             sqrt(3)/2*(1-res.x1_crit-res.x3_crit), ...
             'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    end

    text(-0.05,-0.03,'PM6','FontSize',10,'FontWeight','bold');
    text(1.05,-0.03,'L8-Bo','FontSize',10,'FontWeight','bold');
    text(0.5,sqrt(3)/2+0.03,'Tol','FontSize',10,'FontWeight','bold');
    axis equal off;
    title(sprintf('Ternary: v3=%.0f  X13=%.2f  g23=%.3f', v3o, X13o, g23o));
    legend({'Spinodal','Binodal','Tie lines','CP'}, 'Location', 'southwest');
    saveas(gcf, fullfile(figDir, sprintf('best_case%d_ternary.png', ci)));

    % --- 子图2: Ratio-Solvent 二维图 ---
    figure('Position', [350+ci*300, 50, 600, 550]);
    hold on;
    x_fill = [rc; flip(rd)];
    y_fill = [c2; flip(d2)];
    fill(x_fill, y_fill, [0.65 0.85 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    plot(rc, c2, 'b-o', 'MarkerSize', 5, 'LineWidth', 2);
    plot(rd, d2, 'r-s', 'MarkerSize', 5, 'LineWidth', 2);

    for k = 1:step_:nt
        ic = 2*k-1; id = 2*k;
        if ic>size(xs,1)||id>size(xs,1), break; end
        plot([rc(k) rd(k)], [c2(k) d2(k)], 'Color', [0.5 0.5 0.5], 'LineWidth', 0.4);
    end

    if ~isnan(res.x1_crit)
        cp_r = res.x1_crit/(res.x1_crit+res.x3_crit);
        plot(cp_r, 1-res.x1_crit-res.x3_crit, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    end

    xlabel('L8-Bo / (L8-Bo + PM6)  给/受体比', 'FontSize', 11);
    ylabel('Solvent volume fraction  \phi_2', 'FontSize', 11);
    title(sprintf('Ratio vs Solvent: v3=%.0f  X13=%.2f  g23=%.3f', v3o, X13o, g23o));
    xlim([0 1]); ylim([0 1]);
    legend({'Two-phase region', 'Conc. arm (PM6-rich)', 'Dilute arm (L8Bo-rich)', 'Tie lines', 'CP'}, ...
        'Location', 'best');
    grid on;
    saveas(gcf, fullfile(figDir, sprintf('best_case%d_ratio.png', ci)));
end

fprintf('各案例相图已保存\n');

%% ===== 综合对比: 最优案例的高质量图 =====
fprintf('\n===== 最终综合图 =====\n');

% 取全局最佳做高清图
[ivb, ixb, igb] = ind2sub([nv, nx, ng], idx_best);
res = all_results{ivb, ixb, igb};
v3b = v3_values(ivb); X13b = X13_values(ixb); g23b = g23_values(igb);
xs = res.x_sol;
c1 = xs(1:2:end,2); c3 = xs(1:2:end,3); c2 = xs(1:2:end,4);
d1 = xs(2:2:end,2); d3 = xs(2:2:end,3); d2 = xs(2:2:end,4);
rc = c1./(c1+c3); rd = d1./(d1+d3);
nt = floor(size(xs,1)/2);
step_ = max(1, floor(nt/12));

% 最终Ratio-Solvent图
figure('Position', [100, 100, 750, 650]);
hold on;

% 两相区填充
x_fill = [rc; flip(rd)];
y_fill = [c2; flip(d2)];
fill(x_fill, y_fill, [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.35);

% tie lines (背景)
for k = 1:step_:nt
    ic = 2*k-1; id = 2*k;
    if ic>size(xs,1)||id>size(xs,1), break; end
    plot([rc(k) rd(k)], [c2(k) d2(k)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
end

% Binodal 两臂
plot(rc, c2, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Binodal (conc. arm)');
plot(rd, d2, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Binodal (dilute arm)');

% 标记几根代表性 tie lines
for k = [1, floor(nt/3), floor(2*nt/3), nt]
    ic = 2*k-1; id = 2*k;
    if ic>size(xs,1)||id>size(xs,1), break; end
    plot([rc(k) rd(k)], [c2(k) d2(k)], 'r-', 'LineWidth', 1.2);
end

% Critical point
if ~isnan(res.x1_crit)
    cp_r = res.x1_crit/(res.x1_crit+res.x3_crit);
    cp_y = 1-res.x1_crit-res.x3_crit;
    plot(cp_r, cp_y, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k', 'DisplayName', 'Critical Point');
    text(cp_r+0.02, cp_y, sprintf('CP\n(%.3f, %.3f)', cp_r, cp_y), 'FontSize', 9);
end

xlabel('L8-Bo / (L8-Bo + PM6)   mass ratio in non-solvent fraction', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Solvent volume fraction  \phi_{solvent}', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('PM6/L8-Bo/Toluene Phase Diagram\nv3=%.0f  X13=%.2f  g23=%.3f', v3b, X13b, g23b), ...
    'FontSize', 14, 'FontWeight', 'bold');
xlim([0 1]); ylim([0 1]);
legend('Location', 'northeast');
grid on;
set(gca, 'FontSize', 11);

% 添加说明文本框
dim = [0.15, 0.15, 0.3, 0.12];
str = sprintf(['PM6 DP ≈ %.0f\n', ...
    'X_{13} = %.2f\n', ...
    'g_{23} = %.3f\n', ...
    'Two-phase region: blue shading\n', ...
    'Tie lines: gray/red'], ...
    v3b/1162.6, X13b, g23b);
annotation('textbox', dim, 'String', str, 'FontSize', 9, ...
    'BackgroundColor', 'w', 'EdgeColor', [0.4 0.4 0.4]);

saveas(gcf, fullfile(figDir, 'FINAL_ratio_solvent_diagram.png'));
exportgraphics(gcf, fullfile(figDir, 'FINAL_ratio_solvent_diagram_HR.png'), 'Resolution', 300);
fprintf('最终 Ratio-Solvent 图已保存\n');

%% ===== 参数调节指南总结 =====
fprintf('\n');
fprintf('========================================================\n');
fprintf(' 参数调节指南: 如何获得"规整美观"的三元相图\n');
fprintf('========================================================\n');
fprintf('\n');
fprintf('核心规律:\n');
fprintf('  1. v3 (PM6摩尔体积) 决定相图的"宽度"\n');
fprintf('     - v3 ↓ (PM6变小) → 相图更居中, 两臂更远离L8-Bo角\n');
fprintf('     - v3 ↑ (PM6变大) → 相图挤压到L8-Bo边, 狭窄\n');
fprintf('     - 推荐范围: v3 = 20,000 ~ 100,000 (DP ≈ 17~86)\n');
fprintf('\n');
fprintf('  2. X13 (L8-Bo/PM6 相互作用参数) 决定相分离强度\n');
fprintf('     - X13 ↑ → 两相区扩大, 临界点PM6增大\n');
fprintf('     - X13 ↓ → 两相区缩小甚至消失\n');
fprintf('     - 推荐范围: X13 = 0.65 ~ 0.95\n');
fprintf('\n');
fprintf('  3. g23 (溶剂/PM6 相互作用参数) 调节相容性\n');
fprintf('     - g23 ↑ → 溶剂与PM6更不相容, 拓宽两相区\n');
fprintf('     - g23 ↓ → 体系更相容, 需要更大v3才能分相\n');
fprintf('     - 推荐范围: g23 = 0.39 ~ 0.75\n');
fprintf('\n');
fprintf('推荐配方案例:\n');
fprintf('  A) 宽大规整型: v3≈100000, X13≈0.90, g23≈0.385\n');
fprintf('     → 两相区面积大, tie line 长, 适合展示全局形貌\n');
fprintf('  B) 平衡致密型: v3≈60000,  X13≈0.85, g23≈0.55\n');
fprintf('     → 两臂分布均匀, 略不对称, 形状更传统\n');
fprintf('  C) 原始参数型: v3≈1743900, X13≈0.62, g23≈0.385\n');
fprintf('     → 窄分布, 紧贴L8-Bo边, 适合实际体系模拟\n');
fprintf('\n');
fprintf('调节口诀:\n');
fprintf('  降v3 = 相图向中心拉宽\n');
fprintf('  升X13 = 相分离区扩大\n');
fprintf('  升g23 = 补偿v3降低时的分相驱动力\n');
fprintf('========================================================\n');
