original_dir = pwd;  
output_dir   = fullfile(original_dir, 'output');  
if ~exist(output_dir, 'dir'), mkdir(output_dir); end  
  
%% 1. Binodal (DMF)  
xsol = run_LLE_DMF(fullfile(original_dir, 'src'));  
cd(original_dir);  % run 内部 cd 了，需要切回来  
writematrix(xsol, fullfile(output_dir, 'LLE_DMF.csv'));  
fprintf('LLE_DMF.csv saved. Size: %dx%d\n', size(xsol,1), size(xsol,2));  
  
%% 2. Spinodal  
xsol = run_Spinodal(fullfile(original_dir, 'src', 'Spinodal'));  
cd(original_dir);  
writematrix(xsol, fullfile(output_dir, 'Spinodal.csv'));  
fprintf('Spinodal.csv saved. Size: %dx%d\n', size(xsol,1), size(xsol,2));  
  
%% 3. Critical Point  
xsol = run_CriticalPoint(fullfile(original_dir, 'src', 'Critical_point'));  
cd(original_dir);  
writematrix(xsol, fullfile(output_dir, 'CriticalPoint.csv'));  
fprintf('CriticalPoint.csv saved. Size: %dx%d\n', size(xsol,1), size(xsol,2));  
  
disp('All done.');