clear
clc
global ix3d iX13

options = optimoptions('fmincon','MaxIterations',100000,'OptimalityTolerance',1e-10);
% Might be good choice for initial guess
% Initial guess for DMF-THF-PIM1 system
x0_default = [0.2; 0.5; 0.5];
x0 = [0.2;0.5;0.5];
lb = [0;0;0];
ub = [1;1;1];

A = [0 0 0;0 0 0;1 1 0];
b = [0;0;1];

% interaction parameter nonsolvent(1)-polymer(3)

iX13 =1.7;
nloop = 60;

% polymer volume fraction in dilute phase (no need to change except want to 
% cover larger range)
% ix3d = 1e-70;
ix3d_values = logspace(-5, log10(0.3), nloop);
% Input parameters

xsol = zeros(2*nloop,4);

for i = 1:nloop
    ix3d = ix3d_values(i);   % 直接赋值，不再用 base 累加

    [x,f] = fmincon(@fun,x0,A,b,[],[],lb,ub,[],options);

    % 检测平凡解：两相组成差异是否足够大  
    x2d_calc = 1 - x(3) - ix3d;  
    phase_diff = abs(x(1) - x(3)) + abs(x(2) - x2d_calc);  

    if f < 1e-8 && phase_diff > 1e-4  
        % 真正的两相平衡解，延续  
        x0 = x;  
    else  
        % 平凡解或收敛失败，重置初始猜测  
        x0 = x0_default;  
    end 

    j = i*2;
    xsol(j-1,1) = fun(x);
    xsol(j-1,2) = x(1);
    xsol(j-1,4) = x(2);
    xsol(j-1,3) = 1-x(1)-x(2);writetable(array2table(xsol, 'VariableNames', {'residual','phi1','phi3','phi2'}), 'binodal_results.csv')
    xsol(j,1) = fun(x);
    xsol(j,2) = x(3);
    xsol(j,4) = 1- x(3)-ix3d;
    xsol(j,3) = ix3d;
end
