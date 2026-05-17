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
% 调试提示：若相图异常，可尝试 [0.40, 0.45, 0.447, 0.50, 0.55]
iX13 = 0.55;

nloop = 60;
ix3d_values = logspace(-5, log10(0.3), nloop);

xsol = NaN(2*nloop, 4);
fail_log = {};  % 记录失败详情以便分析

success_count = 0;
consecutive_fail = 0;

for i = 1:nloop
    ix3d = ix3d_values(i);

    % 动态更新 x1d 上限，确保 x2d = 1 - x1d - ix3d > 0
    ub(3) = 1 - ix3d - 1e-6;

    fun_anon = @(x) fun(x, ix3d, iX13);
    [x, f] = fmincon(fun_anon, x0, A, b, [], [], lb, ub, [], options);

    x2d_calc = 1 - x(3) - ix3d;
    x3c = 1 - x(1) - x(2);
    phase_diff = abs(x(1) - x(3)) + abs(x(2) - x2d_calc);
    comp_diff = abs(x3c - ix3d);  % 浓相与稀相的 PM6 浓度差

    % 判定条件：
    % 1. 残差足够小（化学势真正相等）
    % 2. 两相组成差异足够大（防止假相分离）
    % 3. 浓相 PM6 浓度显著高于稀相（核心物理判据）
    is_valid = (f < 1e-7) && (phase_diff > 1e-3) && (comp_diff > 0.01);

    if is_valid
        x0 = x;
        j = i*2;
        xsol(j-1,1) = f;
        xsol(j-1,2) = x(1);            % phi1 (L8-Bo)
        xsol(j-1,3) = x3c;             % phi3 (PM6)
        xsol(j-1,4) = x(2);            % phi2 (Toluene)

        xsol(j,1) = f;
        xsol(j,2) = x(3);              % phi1 (L8-Bo)
        xsol(j,3) = ix3d;              % phi3 (PM6)
        xsol(j,4) = x2d_calc;          % phi2 (Toluene)

        success_count = success_count + 1;
        consecutive_fail = 0;
        fprintf('迭代 %3d 成功 | ix3d=%.4e | f=%.2e | pd=%.3f | cd=%.3f | 浓相=[%.3f,%.3f,%.3f] 稀相=[%.3f,%.3f,%.3f]\n', ...
            i, ix3d, f, phase_diff, comp_diff, x(1), x(2), x3c, x(3), x2d_calc, ix3d);
    else
        consecutive_fail = consecutive_fail + 1;
        fail_log{end+1} = sprintf(...
            '迭代 %3d 失败 | ix3d=%.4e | f=%.2e | pd=%.3f | cd=%.3f | 原因: %s | x=[%.3f,%.3f,%.3f]', ...
            i, ix3d, f, phase_diff, comp_diff, ...
            get_fail_reason(f, phase_diff, comp_diff), x(1), x(2), x(3));
        fprintf('迭代 %3d 失败 | ix3d=%.4e | f=%.2e | pd=%.3f | cd=%.3f | 原因: %s\n', ...
            i, ix3d, f, phase_diff, comp_diff, get_fail_reason(f, phase_diff, comp_diff));

        % 若连续失败超过 5 次，尝试回退到默认初值
        if consecutive_fail >= 5
            x0 = x0_default;
            consecutive_fail = 0;
            fprintf('  -> 连续失败超限，重置初值为默认值\n');
        end
    end
end

xsol(isnan(xsol(:,1)), :) = [];
writetable(array2table(xsol, 'VariableNames', {'residual','phi1_L8Bo','phi3_PM6','phi2_Tol'}), 'binodal_PM6_L8Bo_Tol.csv');

fprintf('\n========================================\n');
fprintf('Binodal 解算完成\n');
fprintf('成功率: %d/%d (%.1f%%)\n', success_count, nloop, success_count/nloop*100);
fprintf('有效数据点: %d 对\n', size(xsol,1)/2);
fprintf('结果已保存至 binodal_PM6_L8Bo_Tol.csv\n');
fprintf('========================================\n');

% 若成功率过低，给出提示
if success_count < nloop * 0.3
    fprintf('\n警告: 成功率低于 30%%，可能原因:\n');
    fprintf('  1. iX13 过大或过小，可尝试 [0.40, 0.45, 0.447, 0.50, 0.55] 对比\n');
    fprintf('  2. x0_default 不合适，可尝试 [0.5; 0.1; 0.8] 或 [0.3; 0.1; 0.9]\n');
    fprintf('  3. g12 多项式系数导致相容性过高，相分离区本身很小\n');
    fprintf('  4. v3 值不正确（当前 v3=%.1f）\n', 1743900);
end

% 辅助函数：生成失败原因字符串
function reason = get_fail_reason(f, pd, cd)
    reasons = {};
    if f >= 1e-7,      reasons{end+1} = '残差大'; end
    if pd <= 1e-3,     reasons{end+1} = '组成差小'; end
    if cd <= 0.01,     reasons{end+1} = 'PM6差小(假相分离)'; end
    reason = strjoin(reasons, '+');
    if isempty(reason), reason = '未知'; end
end
