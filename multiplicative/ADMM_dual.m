function v = ADMM_dual(x, t, alp)
% this function solves min_v -sum_i t log v_i + J((x-tv)/alp); J*=TV
% % INPUTS: 
%       x:      n1*n2 matrix, noisy image
%       t:      positive scalar, the parameter in noise (see our paper)
%       alp:    positive scalar, the parameter in front of TV
% % OUTPUT:
%       v:      n1*n2 matrix, the restored image
%
N = 100000; epsl = 1e-4; % stopping criteria
lam = 1.0;
v0 = x + .2; w0 = x-t*v0; y0 = 0.*x;
for i=1:N
	% v1 solves min_v -t log v + lam/2 |u+tv-x+y|^2
    vec = w0 - x + y0;
	v1 = -vec / 2/t + sqrt(lam^2 *vec.^2 +4*lam*t)/(2*lam*t);
    z0 = x - t*v1 - y0;
    w1 = z0 - TVL2(z0, alp/ 0.78539816339744830962, 4, 0);
    y1 = y0 + w1 + t*v1 - x;
    err = max([sum((y1(:)-y0(:)).^2); sum((v1(:)-v0(:)).^2); sum((w1(:)-w0(:)).^2); sum((t*v1(:)+w1(:)-x(:)).^2)]);
    if err < epsl
        break;
    end
    fprintf('iter %d: error: %.4f\n', i, err);
    y0 = y1; w0=w1; v0 = v1;
end
v = v1;
end




   