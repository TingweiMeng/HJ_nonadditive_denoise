%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the test file for parameter selection s.t. residual norms are close
% case_num is the number of the example: 
%           1 for 4box, 2 for wfc3_uvis_full_field, 3 for abell_2744
% alp_additive is alp for our model, which is selected by hand
% the section "try upper bound" need to be run for several times, until an
%           upper bound for alp_literature is selected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
%% load data
case_num = 1; % choose which example
% generate/load data
if ~exist(sprintf('case%d/data.mat',case_num), 'file')
    t = 20; % parameter for the noise
    generate_data(case_num, t);
end    
load(sprintf('case%d/data.mat',case_num));
%% choose alp for additive model
alp_additive = 4.0;
filename = sprintf('case%d/alp_%.1f/result.mat', case_num,alp_additive);
if exist(filename, 'file')
    load(filename);
else
    v_additive = ADMM_dual(x_noisy*t, t, alp_additive); 
end
resnorm_additive = norm(x_noisy(:)-v_additive(:));
v0_tmp = v_additive; % initialization for ADMM_literature
%% try upper bound
alp_tmp = 8.0; % candidate for alp_upper
v0_tmp = ADMM_literature(x_noisy*t, t, alp_tmp, v0_tmp);
resnorm_tmp = norm(x_noisy(:)-v0_tmp(:));
if resnorm_tmp > resnorm_additive
    alp_upper = alp_tmp;
    fprintf('alp_upper selected.\n');
end
% in general, alp_additive is smaller than alp_litertuare
alp_lower = alp_additive; 
%% choose alp_literature to make resnorm of two models close
epsl = resnorm_additive * 0.001; % stopping criteria
N = 10; % number of attempts
for j = 1:N
    fprintf('\niteration %d\n\n', j);
    alp_tmp  = (alp_lower + alp_upper) * 0.5; % binary search
    v0_tmp = ADMM_literature(x_noisy*t, t, alp_tmp, log(v0_tmp));
    resnorm_tmp = norm(x_noisy(:)-v0_tmp(:));
    if abs(resnorm_tmp - resnorm_additive) < epsl
        break;
    end
    if resnorm_tmp > resnorm_additive
        alp_upper = alp_tmp;
    else
        alp_lower = alp_tmp;
    end
end
v_literature = v0_tmp;
alp_literature = alp_tmp;

folder_name = sprintf('case%d/same_resnorm_alpadditive_%.1f', case_num, alp_additive);
if ~exist(folder_name, 'dir')
   mkdir(folder_name)
end
save(sprintf('%s/result.mat', folder_name), 'v_additive','v_literature',...
    'alp_additive', 'alp_literature', 'x_noisy', 'x_ori', 't');
