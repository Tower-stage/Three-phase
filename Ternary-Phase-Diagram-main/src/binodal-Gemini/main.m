clear
clc

% 优化器设置：关闭平时刷屏的输出(Display=off)，提高运行速度
options = optimoptions('fmincon','MaxIterations',10000,'OptimalityTolerance',1e-10, 'Display', 'off');

% 初始猜测值
x0_default =[0.2; 0.5; 0.5];
x0 = x0_default;

% 【重要修复】上下界不能为绝对的 0 或 1，否则会导致 log(0) 报错！
lb =[1e-6; 1e-6; 1e-6];
ub =[1-1e-6; 1-1e-6; 1-1e-6];

% 【重要修复】线性约束 A*x <= b 
% 对应: x1c + x2c <= 1 - 1e-6 (保证 x3c 大于 0)
A = [1 1 0]; 
b = 1 - 1e-6;

% 输入参数
iX13 = 1.7;
nloop = 60;
ix3d_values = logspace(-5, log10(0.3), nloop);

% 预分配一个全是 NaN 的矩阵，避免没有解的时候存入零或者废数据
xsol = NaN(2*nloop, 4); 

for i = 1:nloop
    ix3d = ix3d_values(i);   

    % 动态更新第三个变量(x1d)的上限，确保 x2d = 1 - x1d - ix3d > 0
    ub(3) = 1 - ix3d - 1e-6;

    % 使用匿名函数传参，彻底抛弃 global
    fun_anon = @(x) fun(x, ix3d, iX13);

    % 执行优化
    [x, f] = fmincon(fun_anon, x0, A, b, [], [], lb, ub,[], options);

    % 检测平凡解：两相组成差异是否足够大  
    x2d_calc = 1 - x(3) - ix3d;  
    phase_diff = abs(x(1) - x(3)) + abs(x(2) - x2d_calc);  

    % 如果误差足够小，且两相差距足够大 (视为找到了真实的Binodal点)
    if f < 1e-7 && phase_diff > 1e-3  
        x0 = x;  % 延续到下一个点，帮助收敛

        % 【重要修复】只有解算成功，才把数据录入结果矩阵
        j = i*2;
        % 浓相 (Concentrated phase)
        xsol(j-1,1) = f;
        xsol(j-1,2) = x(1);            % phi1 (nonsolvent)
        xsol(j-1,3) = 1 - x(1) - x(2); % phi3 (polymer)
        xsol(j-1,4) = x(2);            % phi2 (solvent)

        % 稀相 (Dilute phase)
        xsol(j,1) = f;
        xsol(j,2) = x(3);              % phi1 (nonsolvent)
        xsol(j,3) = ix3d;              % phi3 (polymer)
        xsol(j,4) = x2d_calc;          % phi2 (solvent)

        fprintf('迭代 %d 成功: PIM-1稀相浓度 = %.4e\n', i, ix3d);
    else  
        % 寻找失败或找到平凡解，重置初始值，跳过数据记录（保留NaN）
        x0 = x0_default;  
        fprintf('迭代 %d 失败/平凡解: 尝试重置初值...\n', i);
    end 
end

% 提取结果：将所有未成功计算的 NaN 行删除，得到纯净的相图数据
xsol(isnan(xsol(:,1)), :) =[];

% 【重要修复】导出文件放到循环外面！
writetable(array2table(xsol, 'VariableNames', {'residual','phi1','phi3','phi2'}), 'binodal_results.csv');
disp('Binodal 解算完成，结果已保存！');