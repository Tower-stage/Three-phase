%% ============================================================
%% D18-o-Xylene 快速参数扫描
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

v1 = 1132.1;
v2 = 120.6;       % o-Xylene
p  = [0, 0, 0, -0.2000, 0.7500];

%% D18: 分子量更高，与L8-Bo相容性略好
%% 扫描范围基于PM6-OXy最佳参数(v3=120000, X13=0.80, g23=0.412)向外扩展
v3_values  = [60000, 80000, 100000, 150000, 200000, 300000];
X13_values = [0.45, 0.55, 0.65, 0.75, 0.85, 0.95];
g23_values = [0.30, 0.412, 0.50, 0.60];

nv = length(v3_values); nx = length(X13_values); ng = length(g23_values);
fprintf('D18-OXy 扫描: %d 组合\n', nv*nx*ng);

opts = optimoptions('fsolve', 'Display', 'off', 'MaxIterations', 3000, ...
    'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-10);

guesses = {[0.10, 0.80, 0.15]; [0.15, 0.70, 0.20]; ...
           [0.20, 0.60, 0.25]; [0.25, 0.50, 0.30]; ...
           [0.12, 0.75, 0.18]; [0.08, 0.85, 0.12]};

succ_rate = NaN(nv, nx, ng);
shape_score = NaN(nv, nx, ng);
all_res = cell(nv, nx, ng);

tic;
cnt = 0;
for iv = 1:nv
for ix = 1:nx
for ig = 1:ng
    v3 = v3_values(iv); X13 = X13_values(ix); g23 = g23_values(ig);
    s = v1/v2; r = v1/v3;
    cnt = cnt + 1;

    % Spinodal (低分辨率加速)
    ng_ = 300;
    x1v = linspace(1e-8, 0.999, ng_); x3v = linspace(1e-8, 0.999, ng_);
    [X1, X3] = meshgrid(x1v, x3v); X2 = 1-X1-X3;
    bad = (X1<=1e-8)|(X2<=1e-8)|(X3<=1e-8);
    u1=X1./(X1+X2+eps); u2=X2./(X1+X2+eps);
    g12=p(1)*u2.^4+p(2)*u2.^3+p(3)*u2.^2+p(4)*u2+p(5);
    dg=4*p(1)*u2.^3+3*p(2)*u2.^2+2*p(3)*u2+p(4);
    d2g=12*p(1)*u2.^2+6*p(2)*u2+2*p(3);
    G22=1./X1+s./X2-2*g12+2*(1-2*u2).*dg+u1.*u2.*d2g;
    G23=1./X1-(g12+X13)+s*g23+u2.*(1-3*u2).*dg+u1.*u2.^2.*d2g;
    G33=1./X1+r./X3-2*X13-2*u2.^3.*dg+u2.^3.*u1.*d2g;
    detG=G22.*G33-G23.^2; detG(bad)=NaN;
    fh=figure('Visible','off'); [C_sp,~]=contour(X1,X3,detG,[0 0]); close(fh);
    spin_segs={}; pos=1;
    while pos<=size(C_sp,2), n=C_sp(2,pos); spin_segs{end+1}=[C_sp(1,pos+1:pos+n)',C_sp(2,pos+1:pos+n)']; pos=pos+n+1; end

    % Critical Point
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
    if isnan(x3_crit)&&~isempty(spin_segs), fb3=spin_segs{1}(:,2);[x3_crit,idx]=max(fb3);x1_crit=spin_segs{1}(idx,1); end
    if isnan(x3_crit), x3_crit=0.01; x1_crit=0.2; end

    % Binodal
    x3d_max=max(0.005,x3_crit*0.88);
    n_loop=40;
    x3d_vals=logspace(-5,log10(x3d_max),n_loop);
    x_sol=NaN(2*n_loop,4); succ=0; prev=[];
    for ii=1:n_loop
        x3d=x3d_vals(ii); best_f=Inf; best=[];
        for igg=1:length(guesses)
            x0=guesses{igg};
            if 1-x0(1)-x0(2)<=0,continue;end
            try
                [xt,fv]=fsolve(@(x) chempot_D18_OXy(x,x3d,v1,v2,v3,s,r,X13,g23,p),x0,opts);
                fv2=sum(fv.^2); x3c=1-xt(1)-xt(2);
                if fv2<1e-8 && x3c>x3d+0.003 && abs(xt(1)-xt(3))>0.003
                    if fv2<best_f, best_f=fv2; best=xt; end
                end
            catch, end
        end
        if isempty(best)&&~isempty(prev)
            try
                [xt,fv]=fsolve(@(x) chempot_D18_OXy(x,x3d,v1,v2,v3,s,r,X13,g23,p),prev,opts);
                fv2=sum(fv.^2); x3c=1-xt(1)-xt(2);
                if fv2<1e-8 && x3c>x3d+0.003, best_f=fv2; best=xt; end
            catch, end
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

    if np_final>=5
        c3=x_sol(1:2:end,3); c1=x_sol(1:2:end,2); c2=x_sol(1:2:end,4);
        d1=x_sol(2:2:end,2); d3=x_sol(2:2:end,3); d2=x_sol(2:2:end,4);
        x3c_range=c3(1)-c3(end);
        dist_L8Bo=min(c1); dist_PM6=min(d3(d3>1e-4));
        balance=1-abs(dist_L8Bo-dist_PM6);
        shape_score(iv,ix,ig)=x3c_range*max(0,balance)*(succ_rate(iv,ix,ig)/100);
    end
    all_res{iv,ix,ig}=struct('x_sol',x_sol,'x1_crit',x1_crit,'x3_crit',x3_crit);

    if mod(cnt,16)==0
        fprintf(' [%3d/%3d] v3=%.0f X13=%.2f g23=%.3f | succ=%d/40 | %.1fmin\n',...
            cnt,nv*nx*ng,v3,X13,g23,np_final,toc/60);
    end
end
end
end

fprintf('\nD18-OXy 总耗时: %.1f min\n', toc/60);

% 输出最佳参数
fprintf('\n===== D18-OXy TOP 参数 =====\n');
all_s = shape_score(:); all_s(isnan(all_s))=-Inf; [~,sidx]=sort(all_s,'descend');
for k=1:min(10,length(sidx))
    if isnan(all_s(sidx(k))), continue; end
    [iv,ix,ig]=ind2sub([nv,nx,ng],sidx(k));
    fprintf('Rank %d: v3=%.0f  X13=%.2f  g23=%.3f  succ=%.1f%%  score=%.4f\n',...
        k,v3_values(iv),X13_values(ix),g23_values(ig),succ_rate(iv,ix,ig),shape_score(iv,ix,ig));
end

% 保存
save(fullfile(tempDir, 'scan_D18_OXy.mat'), 'v3_values','X13_values','g23_values','succ_rate','shape_score','all_res');

function F=chempot_D18_OXy(x,x3d,v1,v2,v3,s,r,X13,g23,p)
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
