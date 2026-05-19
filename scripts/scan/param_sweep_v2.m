%% ============================================================
%% 参数扫描 V2: 三维扫描 v3 x X13 x g23
%% 目标: 找出产生"规整美观"三元相图的参数组合规律
%% 新可视化: 给受体比 vs 溶剂百分比
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

%% ===== 固定参数 =====
v1 = 1132.1;
v2 = 106.3;       % Toluene
p  = [0, 0, 0, -0.2000, 0.6000];

%% ===== 扫描网格 (聚焦有效区) =====
v3_values  = [5000, 8000, 12000, 20000, 35000, 60000, 100000, 200000];
X13_values = [0.45, 0.55, 0.65, 0.75, 0.85, 0.95];
g23_values = [0.385, 0.55, 0.75, 1.00];

nv = length(v3_values);
nx = length(X13_values);
ng = length(g23_values);

fprintf('参数扫描 V2: %d x %d x %d = %d 组合\n', nv, nx, ng, nv*nx*ng);
fprintf('v3: %.0f ~ %.0f\n', min(v3_values), max(v3_values));
fprintf('X13: %.2f ~ %.2f\n', min(X13_values), max(X13_values));
fprintf('g23: %.3f ~ %.3f\n', min(g23_values), max(g23_values));

%% ===== 存储 =====
succ_rate  = NaN(nv, nx, ng);
shape_score = NaN(nv, nx, ng);
x3_crit_v  = NaN(nv, nx, ng);
x3c_range  = NaN(nv, nx, ng);
tie_len_v  = NaN(nv, nx, ng);
phi2_range = NaN(nv, nx, ng);
all_results = cell(nv, nx, ng);

%% ===== 优化器 =====
opts = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 3000, ...
    'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10);

guesses = {[0.10, 0.80, 0.15]; [0.15, 0.70, 0.20]; ...
           [0.20, 0.60, 0.25]; [0.25, 0.50, 0.30]; ...
           [0.12, 0.75, 0.18]; [0.08, 0.85, 0.12]};

total_t = tic;
cnt = 0;

for iv = 1:nv
for ix = 1:nx
for ig = 1:ng
    v3  = v3_values(iv);
    X13 = X13_values(ix);
    g23 = g23_values(ig);
    s = v1/v2;
    r = v1/v3;

    cnt = cnt + 1;

    % --- Spinodal ---
    ng_ = 400;
    x1v = linspace(1e-8, 0.999, ng_);
    x3v = linspace(1e-8, 0.999, ng_);
    [X1, X3] = meshgrid(x1v, x3v);
    X2 = 1-X1-X3;
    bad = (X1<=1e-8)|(X2<=1e-8)|(X3<=1e-8);

    u1=X1./(X1+X2+eps); u2=X2./(X1+X2+eps);
    g12=p(1)*u2.^4+p(2)*u2.^3+p(3)*u2.^2+p(4)*u2+p(5);
    dg=4*p(1)*u2.^3+3*p(2)*u2.^2+2*p(3)*u2+p(4);
    d2g=12*p(1)*u2.^2+6*p(2)*u2+2*p(3);

    G22=1./X1+s./X2-2*g12+2*(1-2*u2).*dg+u1.*u2.*d2g;
    G23=1./X1-(g12+X13)+s*g23+u2.*(1-3*u2).*dg+u1.*u2.^2.*d2g;
    G33=1./X1+r./X3-2*X13-2*u2.^3.*dg+u2.^3.*u1.*d2g;
    detG=G22.*G33-G23.^2;
    detG(bad)=NaN;

    fh=figure('Visible','off');
    [C_sp,~]=contour(X1,X3,detG,[0 0],'b-','LineWidth',1);
    close(fh);

    spin_segs={};
    pos=1;
    while pos<=size(C_sp,2)
        n=C_sp(2,pos);
        spin_segs{end+1}=[C_sp(1,pos+1:pos+n)',C_sp(2,pos+1:pos+n)'];
        pos=pos+n+1;
    end

    % --- Critical Point ---
    x1_crit=NaN; x3_crit=NaN;
    for sg=1:length(spin_segs)
        xs1=spin_segs{sg}(:,1); xs3=spin_segs{sg}(:,2); xs2=1-xs1-xs3;
        np=length(xs1); cv=NaN(np,1);
        for i=1:np
            x1i=xs1(i);x3i=xs3(i);x2i=xs2(i);
            if x1i<=1e-8||x2i<=1e-8||x3i<=1e-8,continue;end
            u1i=x1i/(x1i+x2i);u2i=x2i/(x1i+x2i);
            gi=p(1)*u2i^4+p(2)*u2i^3+p(3)*u2i^2+p(4)*u2i+p(5);
            dgi=4*p(1)*u2i^3+3*p(2)*u2i^2+2*p(3)*u2i+p(4);
            d2gi=12*p(1)*u2i^2+6*p(2)*u2i+2*p(3);
            G22i=1/x1i+s/x2i-2*gi+2*(1-2*u2i)*dgi+u1i*u2i*d2gi;
            G23i=1/x1i-(gi+X13)+s*g23+u2i*(1-3*u2i)*dgi+u1i*u2i^2*d2gi;
            G33i=1/x1i+r/x3i-2*X13-2*u2i^3*dgi+u2i^3*u1i*d2gi;
            G222=1/x1i^2-s/x2i^2-6*u2i/x2i*dgi+(3-6*u2i)*u2i/x2i*d2gi;
            G223=1/x1i^2-6*u2i^2/x2i*dgi+3*(1-2*u2i)*u2i^2/x2i*d2gi;
            G233=1/x1i^2-6*u2i^3/x2i*dgi+(3*u2i-6*u2i^2)*u2i^2/x2i*d2gi;
            G333=1/x1i^2-r/x3i^2+6*u2i^4/x2i*dgi+(3*u2i^2-2*u2i^3)*u2i^2/x2i*d2gi;
            cv(i)=G222*G33i^2-3*G223*G23i*G33i+3*G233*G23i^2-G22i*G23i*G333;
        end
        ok=~isnan(cv);
        if sum(ok)<2,continue;end
        sc=find(diff(sign(cv(ok)))~=0,1);
        if ~isempty(sc)
            x1o=xs1(ok);x3o=xs3(ok);co=cv(ok);
            t_=abs(co(sc))/(abs(co(sc))+abs(co(sc+1)));
            x1_crit=x1o(sc)+t_*(x1o(sc+1)-x1o(sc));
            x3_crit=x3o(sc)+t_*(x3o(sc+1)-x3o(sc));
            break;
        end
    end
    if isnan(x3_crit) && ~isempty(spin_segs)
        fb3=spin_segs{1}(:,2); fb1=spin_segs{1}(:,1);
        [x3_crit,idx]=max(fb3); x1_crit=fb1(idx);
    elseif isnan(x3_crit)
        x3_crit=0.01; x1_crit=0.2;
    end
    x2_crit=1-x1_crit-x3_crit;

    % --- Binodal ---
    x3d_max=max(0.005,x3_crit*0.88);
    n_loop=50;
    x3d_vals=logspace(-5,log10(x3d_max),n_loop);
    x_sol=NaN(2*n_loop,4);
    succ=0; prev=[];

    for ii=1:n_loop
        x3d=x3d_vals(ii);
        best_f=Inf; best=[];
        for igg=1:length(guesses)
            x0=guesses{igg};
            if 1-x0(1)-x0(2)<=0,continue;end
            try
                [xt,fv]=fsolve(@(x) chempot_v2(x,x3d,v1,v2,v3,s,r,X13,g23,p),x0,opts);
                fv2=sum(fv.^2); x3c=1-xt(1)-xt(2);
                if fv2<1e-8 && x3c>x3d+0.003 && abs(xt(1)-xt(3))>0.003
                    if fv2<best_f,best_f=fv2;best=xt;end
                end
            catch,end
        end
        if isempty(best) && ~isempty(prev)
            try
                [xt,fv]=fsolve(@(x) chempot_v2(x,x3d,v1,v2,v3,s,r,X13,g23,p),prev,opts);
                fv2=sum(fv.^2); x3c=1-xt(1)-xt(2);
                if fv2<1e-8 && x3c>x3d+0.003
                    best_f=fv2;best=xt;
                end
            catch,end
        end
        if ~isempty(best)
            x3c=1-best(1)-best(2); x2d=1-best(3)-x3d;
            x_sol(2*ii-1,:)=[best_f,best(1),x3c,best(2)];
            x_sol(2*ii,:)=[best_f,best(3),x3d,x2d];
            succ=succ+1; prev=best;
        end
    end

    % 后处理
    x_sol(isnan(x_sol(:,1)),:)=[];
    x_sol(x_sol(:,1)>1e-6,:)=[];
    np_=floor(size(x_sol,1)/2);
    if np_>=3
        [~,sidx]=sort(x_sol(2:2:2*np_,3));
        tmp=x_sol;
        for k=1:np_,x_sol([2*k-1,2*k],:)=tmp([2*sidx(k)-1,2*sidx(k)],:);end
        keep=true(np_,1); pc=Inf;
        for k=1:np_
            c3=x_sol(2*k-1,3);
            if c3>=pc||c3<=0,keep(k)=false;else pc=c3;end
        end
        xn=[];for k=1:np_,if keep(k),xn=[xn;x_sol(2*k-1:2*k,:)];end;end
        x_sol=xn;
    end

    np_final=floor(size(x_sol,1)/2);
    succ_rate(iv,ix,ig)=np_final/n_loop*100;
    x3_crit_v(iv,ix,ig)=x3_crit;

    if np_final>=5
        c3=x_sol(1:2:end,3); d3=x_sol(2:2:end,3);
        c2=x_sol(1:2:end,4); d2=x_sol(2:2:end,4);
        c1=x_sol(1:2:end,2); d1=x_sol(2:2:end,2);

        x3c_range(iv,ix,ig)=c3(1)-c3(end);
        phi2_range(iv,ix,ig)=max(abs(c2-d2));
        tl=sqrt((c1-d1).^2+(c3-d3).^2);
        tie_len_v(iv,ix,ig)=mean(tl);

        dist_L8Bo=min(c1); dist_PM6=min(d3(d3>1e-4));
        balance=1-abs(dist_L8Bo-dist_PM6);
        shape_score(iv,ix,ig)=x3c_range(iv,ix,ig)*max(0,balance)*(succ_rate(iv,ix,ig)/100);
    end

    all_results{iv,ix,ig}=struct('v3',v3,'X13',X13,'g23',g23,'x_sol',x_sol,...
        'x1_crit',x1_crit,'x3_crit',x3_crit,'spin_segs',{spin_segs});

    if mod(cnt,32)==0
        fprintf(' [%3d/%3d] v3=%.0f X13=%.2f g23=%.3f | succ=%d | %.1fmin\n',...
            cnt,nv*nx*ng,v3,X13,g23,np_final,toc(total_t)/60);
    end
end
end
end

fprintf('\n总耗时: %.1f min\n', toc(total_t)/60);

%% ===== 保存 =====
if ~exist('output','dir'),mkdir('output');end
data2=struct();
data2.v3_values=v3_values; data2.X13_values=X13_values;
data2.g23_values=g23_values; data2.succ_rate=succ_rate;
data2.shape_score=shape_score; data2.x3_crit_v=x3_crit_v;
data2.x3c_range=x3c_range; data2.phi2_range=phi2_range;
data2.tie_len_v=tie_len_v;
save(fullfile(tempDir, 'param_sweep_v2.mat'),'data2','all_results','-v7.3');
fprintf('结果已保存至 %s\n', fullfile(tempDir, 'param_sweep_v2.mat'));

%% ===== 找出最佳组合 =====
fprintf('\n===== TOP 20 参数组合 (按形状评分) =====\n');
fprintf('Rank  v3       X13   g23   Succ%%  Score   CP_phi3  Delta\n');
all_scores = shape_score(:);
all_scores(isnan(all_scores)) = -Inf; [~, sidx] = sort(all_scores, 'descend');
rank=0;
for k=1:length(sidx)
    if isnan(all_scores(sidx(k))), continue; end
    rank=rank+1;
    if rank>20, break; end
    [iv,ix,ig]=ind2sub([nv,nx,ng],sidx(k));
    fprintf('%2d    %-7.0f  %.2f  %.3f  %5.1f  %6.4f  %6.4f  %6.4f\n',...
        rank,v3_values(iv),X13_values(ix),g23_values(ig),...
        succ_rate(iv,ix,ig),shape_score(iv,ix,ig),...
        x3_crit_v(iv,ix,ig),x3c_range(iv,ix,ig));
end

%% ===== 选出每种 v3 下最佳的相图进行比较 =====
fprintf('\n===== 各 v3 最佳参数推荐 =====\n');
fprintf('v3        Best X13   Best g23  Succ%%  Score   CP_phi3\n');
for iv=1:nv
    best_score=-Inf; best_ix=1; best_ig=1;
    for ix=1:nx
        for ig=1:ng
            if ~isnan(shape_score(iv,ix,ig)) && shape_score(iv,ix,ig)>best_score
                best_score=shape_score(iv,ix,ig);
                best_ix=ix; best_ig=ig;
            end
        end
    end
    if best_score > -Inf
        fprintf('%-8.0f  %.2f       %.3f      %5.1f  %6.4f  %6.4f\n',...
            v3_values(iv),X13_values(best_ix),g23_values(best_ig),...
            succ_rate(iv,best_ix,best_ig),shape_score(iv,best_ix,best_ig),...
            x3_crit_v(iv,best_ix,best_ig));
    end
end

%% ===== 综合面板: 针对几个"最佳"参数做高质量相图 =====
% 选3个代表性组合: 低v3(大相区), 中v3(平衡), 高v3(接近原始)
candidates = {};
% 1. 找低v3中最好的
best_low = -Inf; bi_low = [];
for iv=1:3  % v3 <= 12000
    for ix=1:nx, for ig=1:ng
        if ~isnan(shape_score(iv,ix,ig)) && shape_score(iv,ix,ig)>best_low
            best_low=shape_score(iv,ix,ig);
            bi_low=[iv,ix,ig];
        end
    end; end
end
if ~isempty(bi_low), candidates{end+1}=bi_low; end

% 2. 中v3 (20000-60000)
best_mid=-Inf; bi_mid=[];
for iv=4:6
    for ix=1:nx, for ig=1:ng
        if ~isnan(shape_score(iv,ix,ig)) && shape_score(iv,ix,ig)>best_mid
            best_mid=shape_score(iv,ix,ig);
            bi_mid=[iv,ix,ig];
        end
    end; end
end
if ~isempty(bi_mid), candidates{end+1}=bi_mid; end

% 3. 全局最佳
tmp_scores = shape_score(:); tmp_scores(isnan(tmp_scores)) = -Inf; [~,gb]=max(tmp_scores);
[iv_gb,ix_gb,ig_gb]=ind2sub([nv,nx,ng],gb);
candidates{end+1}=[iv_gb,ix_gb,ig_gb];

fprintf('\n===== 代表性组合高质量相图 =====\n');

for ci=1:length(candidates)
    iv=candidates{ci}(1); ix=candidates{ci}(2); ig=candidates{ci}(3);
    res=all_results{iv,ix,ig};
    v3o=v3_values(iv); X13o=X13_values(ix); g23o=g23_values(ig);
    sc=shape_score(iv,ix,ig);
    sr=succ_rate(iv,ix,ig);

    fprintf('v3=%.0f  X13=%.2f  g23=%.3f  succ=%.0f%%  score=%.4f\n',v3o,X13o,g23o,sr,sc);

    % ===== 图A: 三元相图 (三角形) =====
    figure('Position',[ci*250,50,650,600]);
    hold on;
    plot([0 0.5],[0 sqrt(3)/2],'k-','LineWidth',0.8);
    plot([0.5 1],[sqrt(3)/2 0],'k-','LineWidth',0.8);
    plot([0 1],[0 0],'k-','LineWidth',0.8);

    % Spinodal
    for sg=1:length(res.spin_segs)
        s1=res.spin_segs{sg}(:,1); s3=res.spin_segs{sg}(:,2);
        s2=1-s1-s3; v=s1>0&s2>0&s3>0&s1<1&s2<1&s3<1;
        xt=s1(v)+0.5*s2(v); yt=sqrt(3)/2*s2(v);
        plot(xt,yt,'r--','LineWidth',1.2);
    end

    xs=res.x_sol;
    if size(xs,1)>=2
        % Binodal 浓相
        x1c=xs(1:2:end,2);x3c=xs(1:2:end,3);x2c=xs(1:2:end,4);
        vc=x1c>0&x2c>0&x3c>0;
        plot(x1c(vc)+0.5*x2c(vc),sqrt(3)/2*x2c(vc),'b-','LineWidth',2.5);
        % Binodal 稀相
        x1d=xs(2:2:end,2);x3d=xs(2:2:end,3);x2d=xs(2:2:end,4);
        vd=x1d>0&x2d>0&x3d>0;
        plot(x1d(vd)+0.5*x2d(vd),sqrt(3)/2*x2d(vd),'b-','LineWidth',2.5);
    end

    % Tie lines
    nt=floor(size(xs,1)/2); step_=max(1,floor(nt/10));
    for k=1:step_:nt
        ic=2*k-1;id=2*k;
        if ic>size(xs,1)||id>size(xs,1),break;end
        plot([xs(ic,2)+0.5*xs(ic,4),xs(id,2)+0.5*xs(id,4)],...
            [sqrt(3)/2*xs(ic,4),sqrt(3)/2*xs(id,4)],...
            'Color',[0.6 0.6 0.6],'LineWidth',0.6);
    end

    if ~isnan(res.x1_crit)
        plot(res.x1_crit+0.5*(1-res.x1_crit-res.x3_crit),...
            sqrt(3)/2*(1-res.x1_crit-res.x3_crit),...
            'ro','MarkerSize',10,'MarkerFaceColor','r');
    end

    text(-0.05,-0.03,'PM6','FontSize',11,'FontWeight','bold');
    text(1.05,-0.03,'L8-Bo','FontSize',11,'FontWeight','bold');
    text(0.5,sqrt(3)/2+0.03,'Toluene','FontSize',11,'FontWeight','bold');
    axis equal off;
    title(sprintf('Ternary: v3=%.0f X13=%.2f g23=%.3f',v3o,X13o,g23o),'FontSize',12);
    legend({'Spinodal','Binodal','Tie lines','CP'},'Location','southwest');
    saveas(gcf, fullfile(figDir, sprintf('ternary_v3%.0f_X13%.2f_g23%.3f.png',v3o,X13o,g23o)));

    % ===== 图B: 给受体比 vs 溶剂% (新型可视化) =====
    figure('Position',[ci*250+325,50,650,600]);
    hold on;

    c1=xs(1:2:end,2);c3=xs(1:2:end,3);c2=xs(1:2:end,4);
    d1=xs(2:2:end,2);d3=xs(2:2:end,3);d2=xs(2:2:end,4);

    rc=c1./(c1+c3);  % 浓相中 L8-Bo/(L8-Bo+PM6)
    rd=d1./(d1+d3);  % 稀相中 L8-Bo/(L8-Bo+PM6)

    % 填充两相区
    x_fill=[rc; flip(rd)];
    y_fill=[c2; flip(d2)];
    fill(x_fill,y_fill,[0.7 0.85 1],'EdgeColor','none','FaceAlpha',0.4);

    % 绘制两臂
    plot(rc,c2,'b-o','MarkerSize',5,'LineWidth',2,'DisplayName','Conc. phase (PM6-rich)');
    plot(rd,d2,'r-s','MarkerSize',5,'LineWidth',2,'DisplayName','Dilute phase (L8Bo-rich)');

    % Tie lines
    for k=1:step_:nt
        ic=2*k-1;id=2*k;
        if ic>size(xs,1)||id>size(xs,1),break;end
        plot([rc(k) rd(k)],[c2(k) d2(k)],'Color',[0.5 0.5 0.5],'LineWidth',0.4);
    end

    % Critical point
    if ~isnan(res.x1_crit)
        cp_r=res.x1_crit/(res.x1_crit+res.x3_crit);
        plot(cp_r,1-res.x1_crit-res.x3_crit,'ko','MarkerSize',10,'MarkerFaceColor','k');
    end

    xlabel('L8-Bo / (L8-Bo + PM6)  ratio  (给/受体比)','FontSize',12);
    ylabel('Solvent volume fraction  \phi_2','FontSize',12);
    title(sprintf('Ratio vs Solvent: v3=%.0f X13=%.2f g23=%.3f',v3o,X13o,g23o),'FontSize',12);
    xlim([0 1]); ylim([0 1]);
    legend('Location','best');
    grid on;
    saveas(gcf, fullfile(figDir, sprintf('ratio_v3%.0f_X13%.2f_g23%.3f.png',v3o,X13o,g23o)));
end

%% ===== 参数规律总结图: g23效应 =====
fprintf('\n===== 生成 g23 效应对比图 =====\n');

% 固定 v3=20000, X13=0.75, 变化 g23
iv_fixed=4; ix_fixed=4;  % v3=20000, X13=0.75
figure('Position',[50,50,1000,800]);
for ig=1:ng
    subplot(2,2,ig);
    res=all_results{iv_fixed,ix_fixed,ig};
    if isempty(res), title(sprintf('g23=%.3f: NO DATA',g23_values(ig))); continue; end
    xs=res.x_sol;
    if size(xs,1)<4, title(sprintf('g23=%.3f: INSUFFICIENT',g23_values(ig))); continue; end

    hold on;
    c1=xs(1:2:end,2);c3=xs(1:2:end,3);c2=xs(1:2:end,4);
    d1=xs(2:2:end,2);d3=xs(2:2:end,3);d2=xs(2:2:end,4);
    rc=c1./(c1+c3); rd=d1./(d1+d3);

    fill([rc;flip(rd)],[c2;flip(d2)],[0.7 0.85 1],'EdgeColor','none','FaceAlpha',0.4);
    plot(rc,c2,'b-o','MarkerSize',4,'LineWidth',1.5);
    plot(rd,d2,'r-s','MarkerSize',4,'LineWidth',1.5);

    nt=floor(size(xs,1)/2); step_=max(1,floor(nt/8));
    for k=1:step_:nt
        ic=2*k-1;id=2*k;
        if ic>size(xs,1)||id>size(xs,1),break;end
        plot([rc(k) rd(k)],[c2(k) d2(k)],'Color',[0.5 0.5 0.5],'LineWidth',0.3);
    end
    if ~isnan(res.x1_crit)
        cp_r=res.x1_crit/(res.x1_crit+res.x3_crit);
        plot(cp_r,1-res.x1_crit-res.x3_crit,'ko','MarkerSize',8,'MarkerFaceColor','k');
    end

    xlabel('L8-Bo/(L8-Bo+PM6)'); ylabel('Solvent \phi_2');
    xlim([0 1]); ylim([0 1]);
    title(sprintf('g_{23}=%.3f (v3=20000, X13=0.75)',g23_values(ig)));
    grid on;
end
saveas(gcf, fullfile(figDir, 'g23_effect_comparison.png'));
exportgraphics(gcf, fullfile(figDir, 'g23_effect_comparison_HR.png'),'Resolution',300);
fprintf('g23 效应对比图已保存\n');

fprintf('\n===== 全部参数扫描完成 =====\n');
fprintf('输出文件位于 %s 目录\n', fullfile(projectDir, 'output'));

%% ===== 化学势函数 =====
function F=chempot_v2(x,x3d,v1,v2,v3,s,r,X13,g23,p)
    x1c=max(x(1),1e-10);x2c=max(x(2),1e-10);x3c=max(1-x1c-x2c,1e-10);
    x1d=max(x(3),1e-10);x2d=max(1-x(3)-x3d,1e-10);x3dv=max(x3d,1e-10);
    u1c=x1c/(x1c+x2c);u2c=x2c/(x1c+x2c);
    g12c=p(1)*u2c^4+p(2)*u2c^3+p(3)*u2c^2+p(4)*u2c+p(5);
    dgc=4*p(1)*u2c^3+3*p(2)*u2c^2+2*p(3)*u2c+p(4);
    u1d=x1d/(x1d+x2d);u2d=x2d/(x1d+x2d);
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
