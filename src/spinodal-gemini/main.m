clear
clc

% 1. 生成 Toluene(x1) 和 PIM-1(x3) 的高密度网格点
x1_val = linspace(1e-8, 0.999, 800); % Toluene从0到1
x3_val = linspace(1e-8, 0.999, 800); % PIM1从0到1
[X1, X3] = meshgrid(x1_val, x3_val);

% 2. 一次性计算所有点的 F 值
F = fun(X1, X3);

% 3. 提取 F = 0 的曲线坐标 (这就是你要的触底的那条线)
% 使用 contour 函数寻找 F=0 的等高线
figure(1);
[C, h] = contour(X1, X3, F, [0 0], 'b-', 'LineWidth', 2);

xlabel('Volume fraction of Toluene (x1)');
ylabel('Volume fraction of PIM-1 (x3)');
title('Spinodal Curve in Cartesian coordinates');
grid on;

% ---- 提取数据点以便后续画三元相图 ----
% C矩阵包含了提取出来的曲线(x1, x3)坐标数据
if size(C, 2) > 1
    % 解析 contour 返回的矩阵 C 
    % (C的第一列是等高线数值和点的数量，后面是具体的x和y坐标)
    num_points = C(2,1);
    Toluene_curve = C(1, 2:1+num_points)';
    PIM1_curve = C(2, 2:1+num_points)';
    THF_curve = 1 - Toluene_curve - PIM1_curve;

    % 合并到 xsol 变量中 (与你原来的数据格式类似: [fval, Tol, PIM, THF])
    % 这里 fval 必定为 0，因为是等高线提取的
    xsol =[zeros(num_points,1), Toluene_curve, PIM1_curve, THF_curve];

    disp('数据计算完成，曲线已经触底！你可以用 xsol 去画你的三元相图了。');
else
    disp('未找到相图边界，请检查相互作用参数！');
end