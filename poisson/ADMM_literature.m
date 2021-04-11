function v = ADMM_literature(x, t, alp, v0)
% this function solves min_v sum_i (tv_i - x_i log v_i) + TV(v)*alp
% % INPUTS: 
%       x:      n1*n2 matrix, noisy image
%       t:      positive scalar, the parameter in noise (see our paper)
%       alp:    positive scalar, the parameter in front of TV
%       v0:     n1*n2 matrix, initialization
% % OUTPUT:
%       v:      n1*n2 matrix, the restored image
%
N = 100000; epsl = 1e-4; % stopping criteria
lam = 1.0;
y0 = 0.*x; w0 = v0;
for i=1:N
	% v1 solves min_v tv-xlogv + lam/2(v-w+y)^2
    s = (w0-y0-t/lam)/2;
    v1 = s + sqrt(s.^2 + x/lam);
    w1 = TVL2(v1+y0, alp/lam / 0.78539816339744830962, 4, 0);
    y1 = y0 + v1 - w1;
    err = max([sum((y1(:)-y0(:)).^2); sum((v1(:)-v0(:)).^2); sum((w1(:)-w0(:)).^2); sum((v1(:)-w1(:)).^2)]);
    if err < epsl
        break;
    end
    fprintf('iter %d: error: %.4f\n', i, err);
    y0 = y1; v0=v1; w0 = w1;
end
v = v1;
end