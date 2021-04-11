function x = generate_multiplicative_noise(x0, L)
% this function add multiplicative noise on the original image x0
% % INPUTS:
%       x0:     n1*n2 matrix, original image
%       L:      positive scalar, the paramter in gamma distribution (t in our paper)
% % OUTPUT:
%       x:      n1*n2 matrix, noisy image
%
rng(1); % fix random seed 1
noise = gamrnd(L,1/L,size(x0));
x = x0 .* noise;
end