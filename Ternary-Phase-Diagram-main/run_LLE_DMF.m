function xsol = run_LLE_DMF(script_dir)  
cd(script_dir);  
run('main.m');  
% clear 只清了本函数的局部变量，xsol 在 run 后仍存在于本函数工作区  
end