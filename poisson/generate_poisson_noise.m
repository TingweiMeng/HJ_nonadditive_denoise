function x = generate_poisson_noise(x0, t)
% this function add Poisson noise on the original image x0
% % INPUTS:
%       x0:     n1*n2 matrix, original image
%       t:      positive scalar, the paramter in Poisson distribution (see our paper)
% % OUTPUT:
%       x:      n1*n2 matrix, noisy image
%
rng(1); % fix random seed 1
x = poissrnd(t*x0)/t;
end