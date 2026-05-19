%% ============================================================
%% 四种三元体系最终相图 — [PM6/D18]:L8-Bo:[Tol/O-Xy]
%% ============================================================
%% 基于 Flory-Huggins 理论，计算并绘制四种体系的完整相图
%% PM6-Tol  |  PM6-OXy  |  D18-Tol  |  D18-OXy
%% 所有参数均通过三维网格扫描优化
%% ============================================================
clear; clc;

%% 自动定位项目根目录（脚本位于 scripts/main/）
scriptDir = fileparts(mfilename('fullpath'));
projectDir = fileparts(scriptDir);
dataDir = fullfile(projectDir, 'output', 'data');
figDir = fullfile(projectDir, 'output', 'figures');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
if ~exist(figDir, 'dir'), mkdir(figDir); end

%% ===== 四种体系配置 =====
% 固定参数
v1 = 1132.1;  % L8-Bo 摩尔体积 [cm^3/mol]

systems = {
    struct('tag','PM6_Tol','name','PM6 / L8-Bo / Toluene',    'donor','PM6','solvent','Toluene',    'v2',106.3,'X13',0.95,'g23',0.385,'v3',100000,'p',[0,0,0,-0.20,0.60]);
    struct('tag','PM6_OXy','name','PM6 / L8-Bo / o-Xylene',   'donor','PM6','solvent','o-Xylene',   'v2',120.6,'X13',0.80,'g23',0.412,'v3',120000,'p',[0,0,0,-0.20,0.75]);
    struct('tag','D18_Tol','name','D18 / L8-Bo / Toluene',    'donor','D18','solvent','Toluene',    'v2',106.3,'X13',0.95,'g23',0.450,'v3', 60000,'p',[0,0,0,-0.20,0.60]);
    struct('tag','D18_OXy','name','D18 / L8-Bo / o-Xylene',   'donor','D18','solvent','o-Xylene',   'v2',120.6,'X13',0.85,'g23',0.300,'v3',100000,'p',[0,0,0,-0.20,0.75]);
};

%% ===== 公共优化器设置 =====
opts = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 5000, ...
    'FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-12);

guesses = {[0.08,0.85,0.12]; [0.10,0.80,0.15]; [0.12,0.75,0.18]; ...
           [0.15,0.70,0.20]; [0.18,0.65,0.22]; [0.20,0.60,0.25]; ...
           [0.25,0.50,0.30]; [0.15,0.72,0.19]};

%% ===== 逐体系计算 =====
all_res = cell(1,4);

for sid = 1:4
    sys = systems{sid};
    s = v1/sys.v2; r = v1/sys.v3;
    fprintf('\n%s\n', repmat('=',1,70));
    fprintf('  %s\n', sys.name);
    fprintf('  v3=%.0f  X13=%.2f  g23=%.3f  s=%.4f  r=%.4e\n', sys.v3, sys.X13, sys.g23, s, r);
    fprintf('%s\n', repmat('=',1,70));

    % ---- 1. Spinodal (高分辨率) ----
    n_g = 1200;
    x1v = linspace(1e-8,0.999,n_g); x3v = linspace(1e-8,0.999,n_g);
    [X1,X3] = meshgrid(x1v,x3v); X2 = 1-X1-X3;
    bad = (X1<=1e-8)|(X2<=1e-8)|(X3<=1e-8);

    u1 = X1./(X1+X2+eps); u2 = X2./(X1+X2+eps);
    p_ = sys.p;
    g12  = p_(1)*u2.^4+p_(2)*u2.^3+p_(3)*u2.^2+p_(4)*u2+p_(5);
    dg   = 4*p_(1)*u2.^3+3*p_(2)*u2.^2+2*p_(3)*u2+p_(4);
    d2g  = 12*p_(1)*u2.^2+6*p_(2)*u2+2*p_(3);

    G22 = 1./X1+s./X2-2*g12+2*(1-2*u2).*dg+u1.*u2.*d2g;
    G23 = 1./X1-(g12+sys.X13)+s*sys.g23+u2.*(1-3*u2).*dg+u1.*u2.^2.*d2g;
    G33 = 1./X1+r./X3-2*sys.X13-2*u2.^3.*dg+u2.^3.*u1.*d2g;
    detG = G22.*G33-G23.^2; detG(bad) = NaN;

    fh = figure('Visible','off');
    [C_sp,~] = contour(X1,X3,detG,[0 0]); close(fh);

    spin_segs = {}; pos = 1;
    while pos <= size(C_sp,2)
        n = C_sp(2,pos);
        spin_segs{end+1} = [C_sp(1,pos+1:pos+n)', C_sp(2,pos+1:pos+n)'];
        pos = pos+n+1;
    end
    fprintf('  Spinodal: %d segments\n', length(spin_segs));

    % ---- 2. Critical Point ----
    x1_crit = NaN; x3_crit = NaN;
    for sg = 1:length(spin_segs)
        xs1 = spin_segs{sg}(:,1); xs3 = spin_segs{sg}(:,2); xs2 = 1-xs1-xs3;
        np = length(xs1); cv = NaN(np,1);
        for i = 1:np
            x1i=xs1(i);x3i=xs3(i);x2i=xs2(i);
            if x1i<=1e-8||x2i<=1e-8||x3i<=1e-8,continue;end
            u1i=x1i/(x1i+x2i);u2i=x2i/(x1i+x2i);
            gi=p_(1)*u2i^4+p_(2)*u2i^3+p_(3)*u2i^2+p_(4)*u2i+p_(5);
            dgi=4*p_(1)*u2i^3+3*p_(2)*u2i^2+2*p_(3)*u2i+p_(4);
            d2gi=12*p_(1)*u2i^2+6*p_(2)*u2i+2*p_(3);
            G22i=1/x1i+s/x2i-2*gi+2*(1-2*u2i)*dgi+u1i*u2i*d2gi;
            G23i=1/x1i-(gi+sys.X13)+s*sys.g23+u2i*(1-3*u2i)*dgi+u1i*u2i^2*d2gi;
            G33i=1/x1i+r/x3i-2*sys.X13-2*u2i^3*dgi+u2i^3*u1i*d2gi;
            G222=1/x1i^2-s/x2i^2-6*u2i/x2i*dgi+(3-6*u2i)*u2i/x2i*d2gi;
            G223=1/x1i^2-6*u2i^2/x2i*dgi+3*(1-2*u2i)*u2i^2/x2i*d2gi;
            G233=1/x1i^2-6*u2i^3/x2i*dgi+(3*u2i-6*u2i^2)*u2i^2/x2i*d2gi;
            G333=1/x1i^2-r/x3i^2+6*u2i^4/x2i*dgi+(3*u2i^2-2*u2i^3)*u2i^2/x2i*d2gi;
            cv(i)=G222*G33i^2-3*G223*G23i*G33i+3*G233*G23i^2-G22i*G23i*G333;
        end
        ok=~isnan(cv); if sum(ok)<2,continue;end
        sc=find(diff(sign(cv(ok)))~=0,1);
        if ~isempty(sc)
            x1o=xs1(ok);x3o=xs3(ok);co=cv(ok);
            t_=abs(co(sc))/(abs(co(sc))+abs(co(sc+1)));
            x1_crit=x1o(sc)+t_*(x1o(sc+1)-x1o(sc));
            x3_crit=x3o(sc)+t_*(x3o(sc+1)-x3o(sc));
            break;
        end
    end
    if isnan(x3_crit)&&~isempty(spin_segs)
        fb=spin_segs{1}(:,2);[x3_crit,idx]=max(fb);x1_crit=spin_segs{1}(idx,1);
    elseif isnan(x3_crit)
        x3_crit=0.01;x1_crit=0.2;
    end
    x2_crit = 1-x1_crit-x3_crit;
    fprintf('  Critical Point: phi1=%.4f phi2=%.4f phi3=%.4f\n', x1_crit, x2_crit, x3_crit);

    % ---- 3. Binodal (120点高分辨率) ----
    x3d_max = max(0.005, x3_crit*0.90);
    n_loop = 120;
    x3d_vals = logspace(-5, log10(x3d_max), n_loop);

    x_sol = NaN(2*n_loop, 4);
    succ = 0; prev = [];

    for ii = 1:n_loop
        x3d = x3d_vals(ii);
        best_f = Inf; best = [];

        for ig = 1:length(guesses)
            x0 = guesses{ig};
            if 1-x0(1)-x0(2) <= 0, continue; end
            try
                [xt,fv] = fsolve(@(x) chempot_eq(x, x3d, v1, sys.v2, sys.v3, s, r, sys.X13, sys.g23, sys.p), x0, opts);
                fv2 = sum(fv.^2); x3c = 1-xt(1)-xt(2);
                if fv2 < 1e-10 && x3c > x3d+0.003 && abs(xt(1)-xt(3)) > 0.005
                    if fv2 < best_f, best_f = fv2; best = xt; end
                end
            catch, end
        end

        if isempty(best) && ~isempty(prev)
            try
                [xt,fv] = fsolve(@(x) chempot_eq(x, x3d, v1, sys.v2, sys.v3, s, r, sys.X13, sys.g23, sys.p), prev, opts);
                fv2 = sum(fv.^2); x3c = 1-xt(1)-xt(2);
                if fv2 < 1e-10 && x3c > x3d+0.003
                    best_f = fv2; best = xt;
                end
            catch, end
        end

        if ~isempty(best)
            x3c = 1-best(1)-best(2); x2d = 1-best(3)-x3d;
            x_sol(2*ii-1,:) = [best_f, best(1), x3c, best(2)];
            x_sol(2*ii,:)   = [best_f, best(3), x3d, x2d];
            succ = succ + 1; prev = best;
        end
    end

    % === 后处理: 严格质量过滤 ===
    x_sol(isnan(x_sol(:,1)), :) = [];
    x_sol(x_sol(:,1) > 1e-8, :) = [];
    np_ = floor(size(x_sol,1)/2);
    if np_ >= 3
        [~, sidx] = sort(x_sol(2:2:2*np_,3));
        tmp = x_sol;
        for k = 1:np_, x_sol([2*k-1,2*k],:) = tmp([2*sidx(k)-1,2*sidx(k)],:); end
        keep = true(np_,1); pc = Inf; pd = -Inf;
        for k = 1:np_
            c3 = x_sol(2*k-1,3); d3 = x_sol(2*k,3);
            if c3 >= pc || c3 <= 0 || d3 <= pd, keep(k) = false;
            else pc = c3; pd = d3; end
        end
        filtered = [];
        for k = 1:np_
            if keep(k), filtered = [filtered; x_sol(2*k-1:2*k,:)]; end
        end
        x_sol = filtered;
    end

    np_final = floor(size(x_sol,1)/2);
    fprintf('  Binodal: %d/%d pairs (%.0f%%)\n', np_final, n_loop, np_final/n_loop*100);

    if np_final >= 3
        c3_f = x_sol(1:2:end,3); d3_f = x_sol(2:2:end,3);
        c1_f = x_sol(1:2:end,2); d1_f = x_sol(2:2:end,2);
        c2_f = x_sol(1:2:end,4); d2_f = x_sol(2:2:end,4);
        conc_mono = all(diff(c3_f) < 0); dil_mono = all(diff(d3_f) > 0);
        fprintf('  Quality: Conc monotonic=%d  Dil monotonic=%d\n', conc_mono, dil_mono);
        fprintf('  Conc phi3: %.4f -> %.4f  Dil phi3: %.6f -> %.4f\n', c3_f(1), c3_f(end), d3_f(1), d3_f(end));
        fprintf('  Residual: %.1e ~ %.1e\n', min(x_sol(:,1)), max(x_sol(:,1)));
    end

    % === 保存CSV ===
    cols = {'residual','phi1_L8Bo','phi3_Donor',['phi2_',sys.solvent]};
    T = array2table(x_sol, 'VariableNames', cols);
    writetable(T, fullfile(dataDir, ['binodal_',sys.tag,'.csv']));
    fprintf('  Saved: output/data/binodal_%s.csv\n', sys.tag);

    writetable(table(x1_crit,x2_crit,x3_crit,'VariableNames',{'phi1','phi2','phi3'}), ...
        fullfile(dataDir, ['critical_',sys.tag,'.csv']));
    fprintf('  Saved: output/data/critical_%s.csv\n', sys.tag);

    spin_all = [];
    for sg = 1:length(spin_segs)
        s1 = spin_segs{sg}(:,1); s3 = spin_segs{sg}(:,2); s2 = 1-s1-s3;
        v = s1>0 & s2>0 & s3>0;
        spin_all = [spin_all; zeros(sum(v),1), s1(v), s3(v), s2(v)];
    end
    writetable(array2table(spin_all,'VariableNames',{'residual','phi1','phi3','phi2'}), ...
        fullfile(dataDir, ['spinodal_',sys.tag,'.csv']));
    fprintf('  Saved: output/data/spinodal_%s.csv (%d pts)\n', sys.tag, size(spin_all,1));

    all_res{sid} = struct('spin_segs',{spin_segs}, 'x_sol',x_sol, ...
        'x1_crit',x1_crit, 'x3_crit',x3_crit, ...
        'tag',sys.tag, 'name',sys.name, 'donor',sys.donor, 'solvent',sys.solvent, ...
        'v3',sys.v3, 'X13',sys.X13, 'g23',sys.g23, 'v2',sys.v2);
end

%% ================================================================
%% 出图
%% ================================================================
fprintf('\n%s\n', repmat('=',1,70));
fprintf('  Generating figures\n');
fprintf('%s\n', repmat('=',1,70));

colors = {[0 0 0.9], [0 0.6 0], [0.9 0 0], [0.6 0 0.6]};

%% ---- 图1: 四体系综合三元相图 ----
figure('Position',[50,50,1200,900]); hold on;
plot([0 0.5],[0 sqrt(3)/2],'k-',[0.5 1],[sqrt(3)/2 0],'k-',[0 1],[0 0],'k-','LineWidth',0.8);

text(0.5, sqrt(3)/2+0.04, 'Solvent (\phi_2)', 'FontSize',12, 'HorizontalAlignment','center', 'FontWeight','bold');
text(-0.06,-0.03, 'Donor (\phi_3)', 'FontSize',12, 'HorizontalAlignment','right', 'FontWeight','bold');
text(1.06,-0.03, 'L8-Bo (\phi_1)', 'FontSize',12, 'HorizontalAlignment','left', 'FontWeight','bold');

lgd = {};
for sid = 1:4
    res = all_res{sid}; clr = colors{sid};

    for sg = 1:length(res.spin_segs)
        s1 = res.spin_segs{sg}(:,1); s3 = res.spin_segs{sg}(:,2);
        s2 = 1-s1-s3; v = s1>0 & s2>0 & s3>0 & s1<1 & s2<1 & s3<1;
        plot(s1(v)+0.5*s2(v), sqrt(3)/2*s2(v), '--', 'Color', clr, 'LineWidth', 1.0);
    end
    lgd{end+1} = [res.name ' Spinodal'];

    xs = res.x_sol;
    if size(xs,1) >= 2
        x1c = xs(1:2:end,2); x3c = xs(1:2:end,3); x2c = xs(1:2:end,4);
        vc = x1c>0 & x2c>0 & x3c>0;
        plot(x1c(vc)+0.5*x2c(vc), sqrt(3)/2*x2c(vc), '-', 'Color', clr, 'LineWidth', 2.5);
        x1d = xs(2:2:end,2); x3d = xs(2:2:end,3); x2d = xs(2:2:end,4);
        vd = x1d>0 & x2d>0 & x3d>0;
        plot(x1d(vd)+0.5*x2d(vd), sqrt(3)/2*x2d(vd), '-', 'Color', clr, 'LineWidth', 2.5);
    end
    lgd{end+1} = [res.name ' Binodal'];

    nt = floor(size(xs,1)/2); step_ = max(1, floor(nt/12));
    for k = 1:step_:nt
        ic=2*k-1; id=2*k;
        if ic>size(xs,1)||id>size(xs,1), break; end
        plot([xs(ic,2)+0.5*xs(ic,4), xs(id,2)+0.5*xs(id,4)], ...
             [sqrt(3)/2*xs(ic,4), sqrt(3)/2*xs(id,4)], ...
             '-', 'Color', [clr 0.20], 'LineWidth', 0.4);
    end

    cp_x = res.x1_crit+0.5*(1-res.x1_crit-res.x3_crit);
    cp_y = sqrt(3)/2*(1-res.x1_crit-res.x3_crit);
    plot(cp_x, cp_y, 'o', 'Color', clr, 'MarkerSize', 10, 'MarkerFaceColor', clr, 'LineWidth', 2);
end

axis equal off;
title('Four-System Ternary Phase Diagram Comparison', 'FontSize',15, 'FontWeight','bold');
legend(lgd, 'Location', 'southwest', 'FontSize', 7);

annotation('textbox',[0.15,0.70,0.32,0.22],'String',{...
    sprintf('PM6-Tol:  v3=%.0f X13=%.2f g23=%.3f',all_res{1}.v3,all_res{1}.X13,all_res{1}.g23),...
    sprintf('PM6-OXy: v3=%.0f X13=%.2f g23=%.3f',all_res{2}.v3,all_res{2}.X13,all_res{2}.g23),...
    sprintf('D18-Tol:  v3=%.0f X13=%.2f g23=%.3f',all_res{3}.v3,all_res{3}.X13,all_res{3}.g23),...
    sprintf('D18-OXy: v3=%.0f X13=%.2f g23=%.3f',all_res{4}.v3,all_res{4}.X13,all_res{4}.g23),...
    '','Flory-Huggins theory','Optimized via 3D grid sweep'},...
    'FontSize',8,'BackgroundColor','w','EdgeColor',[0.5 0.5 0.5]);

saveas(gcf, fullfile(figDir, 'ternary_four_systems.png'));
exportgraphics(gcf, fullfile(figDir, 'ternary_four_systems_HR.png'), 'Resolution', 300);
fprintf('  Saved: ternary_four_systems.png\n');

%% ---- 图2: 分面对比 (2x2 子图) ----
figure('Position',[50,50,1400,1200]);
for sid = 1:4
    res = all_res{sid}; xs = res.x_sol;
    subplot(2,2,sid); hold on;
    plot([0 0.5],[0 sqrt(3)/2],'k-',[0.5 1],[sqrt(3)/2 0],'k-',[0 1],[0 0],'k-','LineWidth',0.6);

    text(-0.04,-0.02,'Donor','FontSize',9,'FontWeight','bold');
    text(1.04,-0.02,'L8-Bo','FontSize',9,'FontWeight','bold');
    text(0.5,sqrt(3)/2+0.025,res.solvent,'FontSize',9,'FontWeight','bold');

    for sg = 1:length(res.spin_segs)
        s1=res.spin_segs{sg}(:,1); s3=res.spin_segs{sg}(:,2);
        s2=1-s1-s3; v=s1>0&s2>0&s3>0&s1<1&s2<1&s3<1;
        plot(s1(v)+0.5*s2(v), sqrt(3)/2*s2(v), 'r--', 'LineWidth', 1.0);
    end

    if size(xs,1)>=2
        x1c=xs(1:2:end,2);x3c=xs(1:2:end,3);x2c=xs(1:2:end,4);
        vc=x1c>0&x2c>0&x3c>0;
        plot(x1c(vc)+0.5*x2c(vc), sqrt(3)/2*x2c(vc), 'b-', 'LineWidth', 2.2);
        x1d=xs(2:2:end,2);x3d=xs(2:2:end,3);x2d=xs(2:2:end,4);
        vd=x1d>0&x2d>0&x3d>0;
        plot(x1d(vd)+0.5*x2d(vd), sqrt(3)/2*x2d(vd), 'b-', 'LineWidth', 2.2);
    end

    nt = floor(size(xs,1)/2); step_ = max(1, floor(nt/10));
    for k = 1:step_:nt
        ic=2*k-1;id=2*k;
        plot([xs(ic,2)+0.5*xs(ic,4),xs(id,2)+0.5*xs(id,4)],...
            [sqrt(3)/2*xs(ic,4),sqrt(3)/2*xs(id,4)],...
            'Color',[0.6 0.6 0.6],'LineWidth',0.4);
    end

    cp_x=res.x1_crit+0.5*(1-res.x1_crit-res.x3_crit);
    cp_y=sqrt(3)/2*(1-res.x1_crit-res.x3_crit);
    plot(cp_x,cp_y,'ro','MarkerSize',9,'MarkerFaceColor','r');

    axis equal off;
    title(sprintf('%s: %s\nv3=%.0f  X_{13}=%.2f  g_{23}=%.3f', ...
        res.donor, res.solvent, res.v3, res.X13, res.g23), ...
        'FontSize',11, 'FontWeight','bold');
    legend({'Spinodal','Binodal','Tie lines','CP'},'Location','southwest','FontSize',7);
end
sgtitle('Four-System Ternary Phase Diagrams (Grid View)', 'FontSize',14, 'FontWeight','bold');
saveas(gcf, fullfile(figDir, 'ternary_grid_2x2.png');
exportgraphics(gcf, fullfile(figDir, 'ternary_grid_2x2_HR.png'), 'Resolution', 300);
fprintf('  Saved: ternary_grid_2x2.png\n');

%% ---- 图3: 四种体系的 Ratio-Solvent 对比图 ----
figure('Position',[50,50,1400,1200]);
for sid = 1:4
    res = all_res{sid}; xs = res.x_sol;
    if size(xs,1) < 6, continue; end
    subplot(2,2,sid); hold on;

    c1=xs(1:2:end,2); c3=xs(1:2:end,3); c2=xs(1:2:end,4);
    d1=xs(2:2:end,2); d3=xs(2:2:end,3); d2=xs(2:2:end,4);
    rc=c1./(c1+c3); rd=d1./(d1+d3);
    nt = floor(size(xs,1)/2); step_ = max(1, floor(nt/10));

    fill([rc;flip(rd)], [c2;flip(d2)], [0.6 0.8 1], 'EdgeColor','none', 'FaceAlpha',0.35);

    for k = 1:step_:nt
        ic=2*k-1; id=2*k;
        if ic>size(xs,1)||id>size(xs,1), break; end
        plot([rc(k) rd(k)], [c2(k) d2(k)], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.6);
    end

    highlight = [1, floor(nt/3), floor(2*nt/3), nt];
    for h = 1:length(highlight)
        k = highlight(h);
        if 2*k-1>size(xs,1)||2*k>size(xs,1), continue; end
        plot([rc(k) rd(k)], [c2(k) d2(k)], 'r-', 'LineWidth', 1.2);
    end

    plot(rc, c2, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.8);
    plot(rd, d2, 'b-s', 'MarkerSize', 3, 'LineWidth', 1.8);

    cp_r = res.x1_crit/(res.x1_crit+res.x3_crit);
    cp_y = 1-res.x1_crit-res.x3_crit;
    plot(cp_r, cp_y, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');

    xlabel('L8-Bo / (L8-Bo + Donor)', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel('Solvent volume fraction  \phi_{solvent}', 'FontSize', 10, 'FontWeight', 'bold');
    title(sprintf('%s / %s\nv3=%.0f  X_{13}=%.2f  g_{23}=%.3f', ...
        res.donor, res.solvent, res.v3, res.X13, res.g23), 'FontSize', 11, 'FontWeight', 'bold');
    xlim([0 1]); ylim([0 1]);
    grid on;
end
sgtitle('Four-System Ratio-Solvent Phase Diagrams (Grid View)', 'FontSize',14, 'FontWeight','bold');
saveas(gcf, fullfile(figDir, 'ratio_grid_2x2.png');
exportgraphics(gcf, fullfile(figDir, 'ratio_grid_2x2_HR.png'), 'Resolution', 300);
fprintf('  Saved: ratio_grid_2x2.png\n');

%% ---- 图4: PM6 vs D18 同溶剂对比 (横向对比) ----
figure('Position',[50,50,1200,550]);

% Toluene 对比
subplot(1,2,1); hold on;
for sid = [1 3]
    res = all_res{sid}; xs = res.x_sol; clr = colors{sid};
    if size(xs,1) < 6, continue; end
    c1=xs(1:2:end,2); c3=xs(1:2:end,3); c2=xs(1:2:end,4);
    d1=xs(2:2:end,2); d3=xs(2:2:end,3); d2=xs(2:2:end,4);
    rc=c1./(c1+c3); rd=d1./(d1+d3);

    fill([rc;flip(rd)], [c2;flip(d2)], clr, 'EdgeColor','none', 'FaceAlpha',0.15);
    plot(rc, c2, '-', 'Color', clr, 'LineWidth', 2.5);
    plot(rd, d2, '-', 'Color', clr, 'LineWidth', 2.5);

    cp_r = res.x1_crit/(res.x1_crit+res.x3_crit);
    cp_y = 1-res.x1_crit-res.x3_crit;
    plot(cp_r, cp_y, 'o', 'Color', clr, 'MarkerSize', 10, 'MarkerFaceColor', clr);
end
xlabel('L8-Bo / (L8-Bo + Donor)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Solvent volume fraction  \phi_{solvent}', 'FontSize', 11, 'FontWeight', 'bold');
title('Toluene: PM6 vs D18', 'FontSize', 13, 'FontWeight', 'bold');
legend({'PM6-Tol','','D18-Tol',''}, 'Location', 'best', 'FontSize', 9);
xlim([0 1]); ylim([0 1]); grid on;

% o-Xylene 对比
subplot(1,2,2); hold on;
for sid = [2 4]
    res = all_res{sid}; xs = res.x_sol; clr = colors{sid};
    if size(xs,1) < 6, continue; end
    c1=xs(1:2:end,2); c3=xs(1:2:end,3); c2=xs(1:2:end,4);
    d1=xs(2:2:end,2); d3=xs(2:2:end,3); d2=xs(2:2:end,4);
    rc=c1./(c1+c3); rd=d1./(d1+d3);

    fill([rc;flip(rd)], [c2;flip(d2)], clr, 'EdgeColor','none', 'FaceAlpha',0.15);
    plot(rc, c2, '-', 'Color', clr, 'LineWidth', 2.5);
    plot(rd, d2, '-', 'Color', clr, 'LineWidth', 2.5);

    cp_r = res.x1_crit/(res.x1_crit+res.x3_crit);
    cp_y = 1-res.x1_crit-res.x3_crit;
    plot(cp_r, cp_y, 'o', 'Color', clr, 'MarkerSize', 10, 'MarkerFaceColor', clr);
end
xlabel('L8-Bo / (L8-Bo + Donor)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Solvent volume fraction  \phi_{solvent}', 'FontSize', 11, 'FontWeight', 'bold');
title('o-Xylene: PM6 vs D18', 'FontSize', 13, 'FontWeight', 'bold');
legend({'PM6-OXy','','D18-OXy',''}, 'Location', 'best', 'FontSize', 9);
xlim([0 1]); ylim([0 1]); grid on;

sgtitle('Donor Comparison under Same Solvent', 'FontSize',14, 'FontWeight','bold');
saveas(gcf, fullfile(figDir, 'ratio_donor_compare.png');
exportgraphics(gcf, fullfile(figDir, 'ratio_donor_compare_HR.png'), 'Resolution', 300);
fprintf('  Saved: ratio_donor_compare.png\n');

%% ---- 图5: 溶剂效应对比 (PM6/D18 分别对比 Tol vs OXy) ----
figure('Position',[50,50,1200,550]);

% PM6 溶剂对比
subplot(1,2,1); hold on;
for sid = [1 2]
    res = all_res{sid}; xs = res.x_sol; clr = colors{sid};
    if size(xs,1) < 6, continue; end
    c1=xs(1:2:end,2); c3=xs(1:2:end,3); c2=xs(1:2:end,4);
    d1=xs(2:2:end,2); d3=xs(2:2:end,3); d2=xs(2:2:end,4);
    rc=c1./(c1+c3); rd=d1./(d1+d3);

    fill([rc;flip(rd)], [c2;flip(d2)], clr, 'EdgeColor','none', 'FaceAlpha',0.15);
    plot(rc, c2, '-', 'Color', clr, 'LineWidth', 2.5);
    plot(rd, d2, '-', 'Color', clr, 'LineWidth', 2.5);

    cp_r = res.x1_crit/(res.x1_crit+res.x3_crit);
    cp_y = 1-res.x1_crit-res.x3_crit;
    plot(cp_r, cp_y, 'o', 'Color', clr, 'MarkerSize', 10, 'MarkerFaceColor', clr);
end
xlabel('L8-Bo / (L8-Bo + Donor)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Solvent volume fraction  \phi_{solvent}', 'FontSize', 11, 'FontWeight', 'bold');
title('PM6: Toluene vs o-Xylene', 'FontSize', 13, 'FontWeight', 'bold');
legend({'PM6-Tol','','PM6-OXy',''}, 'Location', 'best', 'FontSize', 9);
xlim([0 1]); ylim([0 1]); grid on;

% D18 溶剂对比
subplot(1,2,2); hold on;
for sid = [3 4]
    res = all_res{sid}; xs = res.x_sol; clr = colors{sid};
    if size(xs,1) < 6, continue; end
    c1=xs(1:2:end,2); c3=xs(1:2:end,3); c2=xs(1:2:end,4);
    d1=xs(2:2:end,2); d3=xs(2:2:end,3); d2=xs(2:2:end,4);
    rc=c1./(c1+c3); rd=d1./(d1+d3);

    fill([rc;flip(rd)], [c2;flip(d2)], clr, 'EdgeColor','none', 'FaceAlpha',0.15);
    plot(rc, c2, '-', 'Color', clr, 'LineWidth', 2.5);
    plot(rd, d2, '-', 'Color', clr, 'LineWidth', 2.5);

    cp_r = res.x1_crit/(res.x1_crit+res.x3_crit);
    cp_y = 1-res.x1_crit-res.x3_crit;
    plot(cp_r, cp_y, 'o', 'Color', clr, 'MarkerSize', 10, 'MarkerFaceColor', clr);
end
xlabel('L8-Bo / (L8-Bo + Donor)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Solvent volume fraction  \phi_{solvent}', 'FontSize', 11, 'FontWeight', 'bold');
title('D18: Toluene vs o-Xylene', 'FontSize', 13, 'FontWeight', 'bold');
legend({'D18-Tol','','D18-OXy',''}, 'Location', 'best', 'FontSize', 9);
xlim([0 1]); ylim([0 1]); grid on;

sgtitle('Solvent Effect Comparison for Each Donor', 'FontSize',14, 'FontWeight','bold');
saveas(gcf, fullfile(figDir, 'ratio_solvent_compare.png');
exportgraphics(gcf, fullfile(figDir, 'ratio_solvent_compare_HR.png'), 'Resolution', 300);
fprintf('  Saved: ratio_solvent_compare.png\n');

%% ===== 最终总结 =====
fprintf('\n%s\n', repmat('=',1,70));
fprintf('  FOUR-SYSTEM DELIVERY COMPLETE\n');
fprintf('%s\n', repmat('=',1,70));
fprintf('\nOutput files in output/data/ and output/figures/:\n');
for sid = 1:4
    res = all_res{sid};
    fprintf('  data/binodal_%s.csv / data/critical_%s.csv / data/spinodal_%s.csv\n', res.tag, res.tag, res.tag);
end
fprintf('  figures/ternary_four_systems.png     — 四体系综合三元相图\n');
fprintf('  figures/ternary_grid_2x2.png         — 2x2 分面三元相图\n');
fprintf('  figures/ratio_grid_2x2.png           — 2x2 分面 Ratio-Solvent 图\n');
fprintf('  figures/ratio_donor_compare.png      — 给体对比 (同溶剂)\n');
fprintf('  figures/ratio_solvent_compare.png    — 溶剂对比 (同给体)\n');
fprintf('  figures/*_HR.png                     — 300 DPI 高清版本\n');

fprintf('\nOptimized Parameters:\n');
for sid = 1:4
    res = all_res{sid};
    fprintf('  %-10s  v3=%6.0f  X13=%.2f  g23=%.3f  v2=%6.1f\n', ...
        res.tag, res.v3, res.X13, res.g23, res.v2);
end

%% ===== 化学势等式 =====
function F = chempot_eq(x, x3d, v1, v2, v3, s, r, X13, g23, p)
    x1c = max(x(1),1e-10); x2c = max(x(2),1e-10);
    x3c = max(1-x1c-x2c,1e-10);
    x1d = max(x(3),1e-10);
    x2d = max(1-x(3)-x3d,1e-10);
    x3dv = max(x3d,1e-10);

    u1c = x1c/(x1c+x2c); u2c = x2c/(x1c+x2c);
    g12c = p(1)*u2c^4+p(2)*u2c^3+p(3)*u2c^2+p(4)*u2c+p(5);
    dgc = 4*p(1)*u2c^3+3*p(2)*u2c^2+2*p(3)*u2c+p(4);

    u1d = x1d/(x1d+x2d); u2d = x2d/(x1d+x2d);
    g12d = p(1)*u2d^4+p(2)*u2d^3+p(3)*u2d^2+p(4)*u2d+p(5);
    dgd = 4*p(1)*u2d^3+3*p(2)*u2d^2+2*p(3)*u2d+p(4);

    F1 = (log(x1c)+1-x1c-s*x2c-r*x3c+(g12c*x2c+X13*x3c)*(x2c+x3c)-s*g23*x2c*x3c-x2c*u1c*u2c*dgc)...
       - (log(x1d)+1-x1d-s*x2d-r*x3dv+(g12d*x2d+X13*x3dv)*(x2d+x3dv)-s*g23*x2d*x3dv-x2d*u1d*u2d*dgd);
    F2 = (s*log(x2c)+s-x1c-s*x2c-r*x3c+(g12c*x1c+g23*s*x3c)*(x1c+x3c)-X13*x1c*x3c+x1c*u1c*u2c*dgc)...
       - (s*log(x2d)+s-x1d-s*x2d-r*x3dv+(g12d*x1d+g23*s*x3dv)*(x1d+x3dv)-X13*x1d*x3dv+x1d*u1d*u2d*dgd);
    F3 = (r*log(x3c)+r-x1c-s*x2c-r*x3c+(X13*x1c+s*g23*x2c)*(x1c+x2c)-g12c*x1c*x2c)...
       - (r*log(x3dv)+r-x1d-s*x2d-r*x3dv+(X13*x1d+s*g23*x2d)*(x1d+x2d)-g12d*x1d*x2d);
    F = [F1; F2; F3];
end
