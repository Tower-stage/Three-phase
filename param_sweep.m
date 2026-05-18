%% ============================================================
%% PM6 / L8-Bo / Toluene 参数网格扫描
%% 扫描变量: v3 (PM6摩尔体积) x X13 (给受体相互作用参数)
%% 目标: 找出产生"规整"三元相图的参数组合
%% ============================================================
clear; clc;

%% ===== 固定参数 (Tol体系) =====
v1 = 1132.1;      % L8-Bo 摩尔体积
v2 = 106.3;       % Toluene 摩尔体积
g23 = 0.3852;     % 溶剂/PM6
p  = [0, 0, 0, -0.2000, 0.6000];  % g12 线性依赖

%% ===== 扫描网格 =====
v3_values  = [2000, 3000, 5000, 8000, 12000, 20000, 35000, 60000, 100000, 200000, 500000, 1000000, 1743900];
X13_values = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.90, 1.00];

nv = length(v3_values);
nx = length(X13_values);

fprintf('参数扫描: %d x %d = %d 组合\n', nv, nx, nv*nx);
fprintf('v3 范围: %.0f ~ %.0f\n', min(v3_values), max(v3_values));
fprintf('X13 范围: %.2f ~ %.2f\n', min(X13_values), max(X13_values));

%% ===== 预分配结果存储 =====
results = cell(nv, nx);

% 质量指标矩阵
succ_rate  = NaN(nv, nx);  % 成功率
x3c_range  = NaN(nv, nx);  % 浓相PM6范围 (相图宽度)
x3_crit_v  = NaN(nv, nx);  % 临界点PM6
phi2_range = NaN(nv, nx);  % 溶剂浓度跨度
tie_len    = NaN(nv, nx);  % 平均tie line长度
shape_score = NaN(nv, nx); % 形状评分

%% ===== 并行选项 =====
opts = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 3000, ...
    'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10);

guesses = {[0.10, 0.80, 0.15]; [0.15, 0.70, 0.20]; ...
           [0.20, 0.60, 0.25]; [0.25, 0.50, 0.30]; ...
           [0.12, 0.75, 0.18]; [0.08, 0.85, 0.12]};

total_start = tic;
combo_count = 0;

for iv = 1:nv
for ix = 1:nx
    v3  = v3_values(iv);
    X13 = X13_values(ix);
    s = v1/v2;
    r = v1/v3;

    combo_count = combo_count + 1;
    t_combo = tic;

    % --- 1. Spinodal (快速) ---
    n_g = 400;  % 降分辨率加速
    x1v = linspace(1e-8, 0.999, n_g);
    x3v = linspace(1e-8, 0.999, n_g);
    [X1, X3] = meshgrid(x1v, x3v);
    X2 = 1 - X1 - X3;
    bad = (X1<=1e-8) | (X2<=1e-8) | (X3<=1e-8);

    u1 = X1./(X1+X2+eps); u2 = X2./(X1+X2+eps);
    g12  = p(1)*u2.^4+p(2)*u2.^3+p(3)*u2.^2+p(4)*u2+p(5);
    dg   = 4*p(1)*u2.^3+3*p(2)*u2.^2+2*p(3)*u2+p(4);
    d2g  = 12*p(1)*u2.^2+6*p(2)*u2+2*p(3);

    G22 = 1./X1+s./X2-2*g12+2*(1-2*u2).*dg+u1.*u2.*d2g;
    G23 = 1./X1-(g12+X13)+s*g23+u2.*(1-3*u2).*dg+u1.*u2.^2.*d2g;
    G33 = 1./X1+r./X3-2*X13-2*u2.^3.*dg+u2.^3.*u1.*d2g;
    detG = G22.*G33-G23.^2;
    detG(bad) = NaN;

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

    % --- 2. 临界点 ---
    x1_crit = NaN; x3_crit = NaN;
    for sg = 1:length(spin_segs)
        xs1 = spin_segs{sg}(:,1); xs3 = spin_segs{sg}(:,2);
        xs2 = 1-xs1-xs3;
        np = length(xs1);
        cv = NaN(np,1);
        for i = 1:np
            x1i=xs1(i); x3i=xs3(i); x2i=xs2(i);
            if x1i<=1e-8||x2i<=1e-8||x3i<=1e-8, continue; end
            u1i=x1i/(x1i+x2i); u2i=x2i/(x1i+x2i);
            gi=p(1)*u2i^4+p(2)*u2i^3+p(3)*u2i^2+p(4)*u2i+p(5);
            dgi=4*p(1)*u2i^3+3*p(2)*u2i^2+2*p(3)*u2i+p(4);
            d2gi=12*p(1)*u2i^2+6*p(2)*u2i+2*p(3);
            G22_i=1/x1i+s/x2i-2*gi+2*(1-2*u2i)*dgi+u1i*u2i*d2gi;
            G23_i=1/x1i-(gi+X13)+s*g23+u2i*(1-3*u2i)*dgi+u1i*u2i^2*d2gi;
            G33_i=1/x1i+r/x3i-2*X13-2*u2i^3*dgi+u2i^3*u1i*d2gi;
            G222=1/x1i^2-s/x2i^2-6*u2i/x2i*dgi+(3-6*u2i)*u2i/x2i*d2gi;
            G223=1/x1i^2-6*u2i^2/x2i*dgi+3*(1-2*u2i)*u2i^2/x2i*d2gi;
            G233=1/x1i^2-6*u2i^3/x2i*dgi+(3*u2i-6*u2i^2)*u2i^2/x2i*d2gi;
            G333=1/x1i^2-r/x3i^2+6*u2i^4/x2i*dgi+(3*u2i^2-2*u2i^3)*u2i^2/x2i*d2gi;
            cv(i)=G222*G33_i^2-3*G223*G23_i*G33_i+3*G233*G23_i^2-G22_i*G23_i*G333;
        end
        ok=~isnan(cv);
        if sum(ok)<2, continue; end
        sc=find(diff(sign(cv(ok)))~=0,1);
        if ~isempty(sc)
            x1v_ok=xs1(ok); x3v_ok=xs3(ok); cv_ok=cv(ok);
            t_=abs(cv_ok(sc))/(abs(cv_ok(sc))+abs(cv_ok(sc+1)));
            x1_crit=x1v_ok(sc)+t_*(x1v_ok(sc+1)-x1v_ok(sc));
            x3_crit=x3v_ok(sc)+t_*(x3v_ok(sc+1)-x3v_ok(sc));
            break;
        end
    end
    if isnan(x3_crit) && ~isempty(spin_segs)
        fallback_x3 = spin_segs{1}(:,2); fallback_x1 = spin_segs{1}(:,1);
        [x3_crit,idx]=max(fallback_x3); x1_crit=fallback_x1(idx);
    elseif isnan(x3_crit)
        x3_crit = 0.01; x1_crit = 0.2;  % 兜底值
    end
    x2_crit = 1-x1_crit-x3_crit;

    % --- 3. Binodal (精简版: 50点) ---
    x3d_max = max(0.005, x3_crit*0.88);
    n_loop = 50;
    x3d_vals = logspace(-5, log10(x3d_max), n_loop);

    x_sol = NaN(2*n_loop, 4);
    succ = 0;
    prev = [];

    for ii = 1:n_loop
        x3d = x3d_vals(ii);
        best_f = Inf; best = [];

        for ig = 1:length(guesses)
            x0 = guesses{ig};
            if 1-x0(1)-x0(2)<=0, continue; end
            try
                [xt, fv] = fsolve(@(x) chempot_sw(x, x3d, v1, v2, v3, s, r, X13, g23, p), x0, opts);
                fv2 = sum(fv.^2);
                x3c = 1-xt(1)-xt(2);
                if fv2 < 1e-8 && x3c > x3d+0.003 && abs(xt(1)-xt(3))>0.003
                    if fv2 < best_f, best_f = fv2; best = xt; end
                end
            catch, end
        end

        if isempty(best) && ~isempty(prev)
            try
                [xt, fv] = fsolve(@(x) chempot_sw(x, x3d, v1, v2, v3, s, r, X13, g23, p), prev, opts);
                fv2 = sum(fv.^2);
                x3c = 1-xt(1)-xt(2);
                if fv2 < 1e-8 && x3c > x3d+0.003
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

    % 后处理
    x_sol(isnan(x_sol(:,1)),:) = [];
    x_sol(x_sol(:,1)>1e-6,:) = [];

    % 排序 + 单调性过滤
    np_ = floor(size(x_sol,1)/2);
    if np_ >= 3
        [~, sidx] = sort(x_sol(2:2:2*np_,3));
        tmp = x_sol;
        for k = 1:np_
            x_sol([2*k-1,2*k],:) = tmp([2*sidx(k)-1,2*sidx(k)],:);
        end
        % 单调性
        keep = true(np_,1);
        pc = Inf;
        for k = 1:np_
            c3 = x_sol(2*k-1,3);
            if c3 >= pc || c3 <= 0, keep(k)=false;
            else pc = c3; end
        end
        x_new = [];
        for k = 1:np_
            if keep(k), x_new=[x_new; x_sol(2*k-1:2*k,:)]; end
        end
        x_sol = x_new;
    end

    np_final = floor(size(x_sol,1)/2);
    succ_rate(iv,ix) = np_final / n_loop * 100;
    x3_crit_v(iv,ix) = x3_crit;

    % --- 4. 质量指标 ---
    if np_final >= 5
        c3 = x_sol(1:2:end,3);  % 浓相PM6
        d3 = x_sol(2:2:end,3);  % 稀相PM6
        c2 = x_sol(1:2:end,4);  % 浓相溶剂
        d2 = x_sol(2:2:end,4);  % 稀相溶剂
        c1 = x_sol(1:2:end,2);  % 浓相L8Bo
        d1 = x_sol(2:2:end,2);  % 稀相L8Bo

        x3c_range(iv,ix) = c3(1) - c3(end);  % 浓相PM6跨度
        phi2_range(iv,ix) = max(abs(c2-d2)); % 溶剂浓度跨度

        % Tie line 平均长度
        tl = sqrt((c1-d1).^2 + (c3-d3).^2);
        tie_len(iv,ix) = mean(tl);

        % 形状评分: 平衡性越好分数越高
        % 好形状标准: 浓相臂远离L8Bo边 (>0.15), 稀相臂远离PM6边
        dist_from_L8Bo = min(c1);              % 浓相离L8Bo最近距离
        dist_from_PM6  = min(d3(d3>1e-4));     % 稀相离PM6最近距离
        balance = 1 - abs(dist_from_L8Bo - dist_from_PM6);  % 两臂平衡度

        % 形状分 = 跨度*平衡度*成功率因子
        shape_score(iv,ix) = x3c_range(iv,ix) * max(0,balance) * (succ_rate(iv,ix)/100);
    end

    % 存储精简结果 (只保留关键数据)
    results{iv,ix} = struct('v3',v3,'X13',X13,'x_sol',x_sol,...
        'x1_crit',x1_crit,'x3_crit',x3_crit,'spin_segs',{spin_segs});

    elapsed = toc(t_combo);
    if mod(combo_count, 20) == 0
        fprintf('  [%3d/%3d] v3=%.0f X13=%.2f | succ=%d/%d | score=%.3f | %.1fs\n', ...
            combo_count, nv*nx, v3, X13, np_final, n_loop, shape_score(iv,ix), elapsed);
    end
end
end

fprintf('\n总耗时: %.1f min\n', toc(total_start)/60);

%% ===== 保存扫描结果 =====
if ~exist('output','dir'), mkdir('output'); end
data = struct();
data.v3_values = v3_values;
data.X13_values = X13_values;
data.succ_rate = succ_rate;
data.shape_score = shape_score;
data.x3c_range = x3c_range;
data.phi2_range = phi2_range;
data.tie_len = tie_len;
data.x3_crit_v = x3_crit_v;
save('output/param_sweep_results.mat', 'data', 'results', '-v7.3');
fprintf('扫描结果已保存至 output/param_sweep_results.mat\n');

%% ===== 热力图: 形状评分 =====
figure('Position', [100, 100, 1000, 400]);

subplot(1,3,1);
imagesc(X13_values, log10(v3_values), shape_score);
set(gca, 'YDir', 'normal');
colormap(jet); colorbar;
xlabel('X_{13} (L8-Bo/PM6 interaction)');
ylabel('log_{10}(v_3) (PM6 molar volume)');
title('Shape Quality Score');
for iv = 1:nv
    for ix = 1:nx
        if succ_rate(iv,ix) > 0
            text(X13_values(ix), log10(v3_values(iv)), ...
                sprintf('%.0f', succ_rate(iv,ix)), ...
                'HorizontalAlignment', 'center', 'FontSize', 6);
        end
    end
end

subplot(1,3,2);
imagesc(X13_values, log10(v3_values), x3c_range);
set(gca, 'YDir', 'normal');
colormap(jet); colorbar;
xlabel('X_{13}');
ylabel('log_{10}(v_3)');
title('Conc. PM6 Span (\Delta\phi_3)');

subplot(1,3,3);
imagesc(X13_values, log10(v3_values), x3_crit_v);
set(gca, 'YDir', 'normal');
colormap(jet); colorbar;
xlabel('X_{13}');
ylabel('log_{10}(v_3)');
title('Critical Point \phi_3 (PM6)');

saveas(gcf, 'output/param_sweep_heatmaps.png');
exportgraphics(gcf, 'output/param_sweep_heatmaps_HR.png', 'Resolution', 300);
fprintf('热力图已保存\n');

%% ===== 找出最佳参数组合 =====
fprintf('\n===== TOP 10 参数组合 =====\n');
fprintf('Rank  v3        X13   Succ%%  Score   CP_phi3  Delta_phi3\n');
[sorted_score, idx] = sort(shape_score(:), 'descend');
for rank = 1:min(10, sum(~isnan(sorted_score)))
    [iv, ix] = ind2sub([nv, nx], idx(rank));
    fprintf('%2d    %-8.0f  %.2f  %5.1f  %6.4f  %6.4f  %6.4f\n', ...
        rank, v3_values(iv), X13_values(ix), succ_rate(iv,ix), ...
        shape_score(iv,ix), x3_crit_v(iv,ix), x3c_range(iv,ix));
end

%% ===== 绘制最佳几个参数组合的相图 =====
n_best = min(6, sum(~isnan(shape_score(:))));
best_idx = idx(1:n_best);

figure('Position', [50, 50, 1200, 800]);
for bi = 1:n_best
    [iv, ix] = ind2sub([nv, nx], best_idx(bi));
    res = results{iv,ix};
    if isempty(res), continue; end

    subplot(2, 3, bi);
    hold on;

    % Spinodal
    for sg = 1:length(res.spin_segs)
        s1 = res.spin_segs{sg}(:,1); s3 = res.spin_segs{sg}(:,2);
        s2 = 1-s1-s3; v = s1>0 & s2>0 & s3>0 & s1<1 & s2<1 & s3<1;
        xt = s1(v)+0.5*s2(v); yt = sqrt(3)/2*s2(v);
        plot(xt, yt, 'r--', 'LineWidth', 0.8);
    end

    % Binodal
    xs = res.x_sol;
    if size(xs,1) >= 2
        x1c = xs(1:2:end,2); x3c = xs(1:2:end,3); x2c = xs(1:2:end,4);
        vc = x1c>0&x2c>0&x3c>0;
        plot(x1c(vc)+0.5*x2c(vc), sqrt(3)/2*x2c(vc), 'b-', 'LineWidth', 1.5);
        x1d = xs(2:2:end,2); x3d = xs(2:2:end,3); x2d = xs(2:2:end,4);
        vd = x1d>0&x2d>0&x3d>0;
        plot(x1d(vd)+0.5*x2d(vd), sqrt(3)/2*x2d(vd), 'b-', 'LineWidth', 1.5);
    end

    % Tie lines
    nt = floor(size(xs,1)/2);
    step_ = max(1, floor(nt/8));
    for k = 1:step_:nt
        ic=2*k-1; id=2*k;
        if ic>size(xs,1)||id>size(xs,1), break; end
        plot([xs(ic,2)+0.5*xs(ic,4) xs(id,2)+0.5*xs(id,4)], ...
             [sqrt(3)/2*xs(ic,4) sqrt(3)/2*xs(id,4)], ...
             'Color', [0.6 0.6 0.6], 'LineWidth', 0.4);
    end

    % CP
    if ~isnan(res.x1_crit)
        plot(res.x1_crit+0.5*(1-res.x1_crit-res.x3_crit), sqrt(3)/2*(1-res.x1_crit-res.x3_crit), ...
            'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    end

    plot([0 0.5],[0 sqrt(3)/2],'k-',[0.5 1],[sqrt(3)/2 0],'k-',[0 1],[0 0],'k-');
    axis equal off;
    title(sprintf('v3=%.0f, X13=%.2f (score=%.3f)', res.v3, res.X13, shape_score(iv,ix)), 'FontSize', 9);
end
saveas(gcf, 'output/param_sweep_best_ternary.png');
exportgraphics(gcf, 'output/param_sweep_best_ternary_HR.png', 'Resolution', 300);
fprintf('最佳相图面板已保存\n');

%% ===== 新可视化: 给受体比 vs 溶剂百分比 =====
fprintf('\n===== 生成 ratio vs solvent 可视化 =====\n');

figure('Position', [50, 50, 1200, 800]);
for bi = 1:n_best
    [iv, ix] = ind2sub([nv, nx], best_idx(bi));
    res = results{iv,ix};
    if isempty(res), continue; end

    subplot(2, 3, bi);
    xs = res.x_sol;
    if size(xs,1) < 6, continue; end

    hold on;
    % 浓相臂: L8-Bo/(L8-Bo+PM6) vs 溶剂%
    c1 = xs(1:2:end,2); c3 = xs(1:2:end,3); c2 = xs(1:2:end,4);
    ratio_c = c1 ./ (c1 + c3);  % 浓相中L8-Bo占比
    plot(ratio_c, c2, 'b-o', 'MarkerSize', 4, 'LineWidth', 1.5);

    % 稀相臂
    d1 = xs(2:2:end,2); d3 = xs(2:2:end,3); d2 = xs(2:2:end,4);
    ratio_d = d1 ./ (d1 + d3);
    plot(ratio_d, d2, 'r-s', 'MarkerSize', 4, 'LineWidth', 1.5);

    % Tie lines
    nt = floor(size(xs,1)/2);
    step_ = max(1, floor(nt/8));
    for k = 1:step_:nt
        ic=2*k-1; id=2*k;
        if ic>size(xs,1)||id>size(xs,1), break; end
        plot([ratio_c(k) ratio_d(k)], [c2(k) d2(k)], ...
            'Color', [0.5 0.5 0.5], 'LineWidth', 0.4);
    end

    % Critical point
    if ~isnan(res.x1_crit)
        cp_r = res.x1_crit / (res.x1_crit + res.x3_crit);
        plot(cp_r, 1-res.x1_crit-res.x3_crit, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    end

    xlabel('L8-Bo / (L8-Bo+PM6)  ratio');
    ylabel('Solvent \phi_2');
    title(sprintf('v3=%.0f, X13=%.2f', res.v3, res.X13), 'FontSize', 9);
    xlim([0 1]); ylim([0 1]);
    legend({'Conc. arm', 'Dilute arm', 'Tie lines', 'CP'}, 'Location', 'best', 'FontSize', 6);
    grid on;
end
saveas(gcf, 'output/param_sweep_ratio_vs_solvent.png');
exportgraphics(gcf, 'output/param_sweep_ratio_vs_solvent_HR.png', 'Resolution', 300);
fprintf('Ratio vs Solvent 图已保存\n');

%% ===== 针对性分析: 最佳参数组合的详细相图 =====
fprintf('\n===== 推荐参数组合分析 =====\n');
[~, idx_best] = max(shape_score(:));
[iv_best, ix_best] = ind2sub([nv, nx], idx_best);
fprintf('最佳参数: v3 = %.0f, X13 = %.2f\n', v3_values(iv_best), X13_values(ix_best));
fprintf('形状评分: %.4f\n', shape_score(iv_best, ix_best));
fprintf('成功率: %.1f%%\n', succ_rate(iv_best, ix_best));
fprintf('临界点 phi3: %.4f\n', x3_crit_v(iv_best, ix_best));

% 找出几个不同风格的"好"参数
fprintf('\n推荐参数组合 (可根据需要选择):\n');
fprintf('  类型           v3      X13   特征\n');
for rank = 1:min(8, sum(~isnan(sorted_score)))
    [iv, ix] = ind2sub([nv, nx], idx(rank));
    style = '';
    if x3c_range(iv,ix) > 0.03, style = [style '大相区 ']; end
    if x3_crit_v(iv,ix) > 0.02, style = [style '高CP_PM6 ']; end
    if phi2_range(iv,ix) > 0.1, style = [style '宽溶剂跨 ']; end
    fprintf('  %2d. %-15s %-6.0f %.2f  %s\n', rank, style, v3_values(iv), X13_values(ix), style);
end

fprintf('\n===== 参数扫描完成 =====\n');

%% ===== 化学势等式 (精简版) =====
function F = chempot_sw(x, x3d, v1, v2, v3, s, r, X13, g23, p)
    x1c=max(x(1),1e-10); x2c=max(x(2),1e-10); x3c=max(1-x1c-x2c,1e-10);
    x1d=max(x(3),1e-10); x2d=max(1-x(3)-x3d,1e-10); x3dv=max(x3d,1e-10);

    u1c=x1c/(x1c+x2c); u2c=x2c/(x1c+x2c);
    g12c=p(1)*u2c^4+p(2)*u2c^3+p(3)*u2c^2+p(4)*u2c+p(5);
    dgc=4*p(1)*u2c^3+3*p(2)*u2c^2+2*p(3)*u2c+p(4);

    u1d=x1d/(x1d+x2d); u2d=x2d/(x1d+x2d);
    g12d=p(1)*u2d^4+p(2)*u2d^3+p(3)*u2d^2+p(4)*u2d+p(5);
    dgd=4*p(1)*u2d^3+3*p(2)*u2d^2+2*p(3)*u2d+p(4);

    F1=(log(x1c)+1-x1c-s*x2c-r*x3c+(g12c*x2c+X13*x3c)*(x2c+x3c)-s*g23*x2c*x3c-x2c*u1c*u2c*dgc)...
      -(log(x1d)+1-x1d-s*x2d-r*x3dv+(g12d*x2d+X13*x3dv)*(x2d+x3dv)-s*g23*x2d*x3dv-x2d*u1d*u2d*dgd);
    F2=(s*log(x2c)+s-x1c-s*x2c-r*x3c+(g12c*x1c+g23*s*x3c)*(x1c+x3c)-X13*x1c*x3c+x1c*u1c*u2c*dgc)...
      -(s*log(x2d)+s-x1d-s*x2d-r*x3dv+(g12d*x1d+g23*s*x3dv)*(x1d+x3dv)-X13*x1d*x3dv+x1d*u1d*u2d*dgd);
    F3=(r*log(x3c)+r-x1c-s*x2c-r*x3c+(X13*x1c+s*g23*x2c)*(x1c+x2c)-g12c*x1c*x2c)...
      -(r*log(x3dv)+r-x1d-s*x2d-r*x3dv+(X13*x1d+s*g23*x2d)*(x1d+x2d)-g12d*x1d*x2d);
    F=[F1;F2;F3];
end
