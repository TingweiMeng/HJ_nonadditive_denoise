%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the test file for solving the two models with a range of alp
% case_num is the number of the example: 
%           1 for 4box, 2 for wfc3_uvis_full_field, 3 for abell_2744
% gap is the gap between different alp's 
% num_test is the total number of alp's to try
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
%% choose which example and which alp
case_num = 1;
gap = 1;
num_test = 10; % solve the two models with alp= gap,..., num_test*gap

% generate/load data
str = sprintf('./case%d', case_num);
if ~exist(sprintf('%s/data.mat',str), 'file')
    t = 20.0; % parameter for the noise
    generate_data(case_num, t);
end
load(sprintf('%s/data.mat',str));

err = zeros(num_test, 2);
v_literature = 0*x_noisy + 1.0; % init for literature model

for i = 1:num_test
    fprintf('\nexample %d\n\n', i);
    alp = i;
    % additive model
    v_additive = ADMM_dual(x_noisy*t, t, alp);
    err(i, 1) = sum((v_additive(:)-x_ori(:)).^2);
    % literature model
    v_literature = ADMM_literature(x_noisy*t, t, alp, v_literature);
    err(i, 2) = sum((v_literature(:)-x_ori(:)).^2);
    
    folder_name = sprintf('%s/alp_%.1f',str,i); % note: change this if gap is smaller
    if ~exist(folder_name, 'dir')
       mkdir(folder_name)
    end
    save(sprintf('%s/result.mat', folder_name), 'v_additive','v_literature','alp');
end

save(sprintf('%s/errors.mat',str), 'err');