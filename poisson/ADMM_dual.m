function v = ADMM_dual(x, t, alp)
% this function solves min_v sum_i(tv_i log v_i - tv_i) + J((x-tv)/alp); J*=TV
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
	v1 = Newton_step(v0, t, lam*t^2, lam*t*(x-w0-y0), epsl/10);
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

function v = Newton_step(v0, t, lam, c, epsl)
% this function implements Newton's method for tlog v + lam*v = c (pointwisely)
% % INPUTS: 
%       v0:     n1*n2 matrix, initialization
%       t:      scalar
%       lam:    scalar
%       c:      n1*n2 matrix
%       epsl:   stopping tolerate
% % OUTPUT:
%       v:      n1*n2 matrix, solution
%
N = 100000;
for i=1:N
    v1 = v0 - (t*log(v0) + lam*v0 - c) ./ (t./v0 + lam);
    v1 = max(v1, epsl);
    err = sum((v1(:)-v0(:)).^2);
    if err < epsl
        break;
    end
    v0 = v1;
end
if i == N
    fprintf('Newton does not converge.\n');
end
v = v1;
end