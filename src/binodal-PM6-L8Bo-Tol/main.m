clear
clc

% 优化器设置
options = optimoptions('fmincon','MaxIterations',10000,'OptimalityTolerance',1e-10, 'Display', 'off');

% 初始猜测值 [浓相phi1(L8Bo); 浓相phi2(Tol); 稀相phi1(L8Bo)]
x0_default = [0.2; 0.5; 0.5];
x0 = x0_default;

% 边界保护
lb = [1e-6; 1e-6; 1e-6];
ub = [1-1e-6; 1-1e-6; 1-1e-6];

% 线性约束：x1c + x2c <= 1 - 1e-6 (保证 x3c > 0)
A = [1 1 0];
b = 1 - 1e-6;

% L8-Bo / PM6 相互作用参数
% 注意：readme.md 推荐值为 0.5500；spinodal-gemini 调试中使用 0.447
% 如果相图异常（如两相区过大/过小），可尝试调整此值
iX13 = 0.55;

nloop = 60;
ix3d_values = logspace(-5, log10(0.3), nloop);

xsol = NaN(2*nloop, 4);

for i = 1:nloop
    ix3d = ix3d_values(i);

    % 动态更新 x1d 上限，确保 x2d = 1 - x1d - ix3d > 0
    ub(3) = 1 - ix3d - 1e-6;

    fun_anon = @(x) fun(x, ix3d, iX13);
    [x, f] = fmincon(fun_anon, x0, A, b, [], [], lb, ub, [], options);

    x2d_calc = 1 - x(3) - ix3d;
    phase_diff = abs(x(1) - x(3)) + abs(x(2) - x2d_calc);

    if f < 1e-7 && phase_diff > 1e-3
        x0 = x;
        j = i*2;
        xsol(j-1,1) = f;
        xsol(j-1,2) = x(1);            % phi1 (L8-Bo)
        xsol(j-1,3) = 1 - x(1) - x(2); % phi3 (PM6)
        xsol(j-1,4) = x(2);            % phi2 (Toluene)

        xsol(j,1) = f;
        xsol(j,2) = x(3);              % phi1 (L8-Bo)
        xsol(j,3) = ix3d;              % phi3 (PM6)
        xsol(j,4) = x2d_calc;          % phi2 (Toluene)

        fprintf('迭代 %d 成功: PM6稀相浓度 = %.4e\n', i, ix3d);
    else
        x0 = x0_default;
        fprintf('迭代 %d 失败/平凡解: 尝试重置初值...\n', i);
    end
end

xsol(isnan(xsol(:,1)), :) = [];
writetable(array2table(xsol, 'VariableNames', {'residual','phi1_L8Bo','phi3_PM6','phi2_Tol'}), 'binodal_PM6_L8Bo_Tol.csv');
disp('Binodal 解算完成，结果已保存至 binodal_PM6_L8Bo_Tol.csv');
