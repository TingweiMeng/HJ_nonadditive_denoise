function v = ADMM_literature(x, t, alp, v0)
% this function solves min_v sum_i (t log v_i + xi/vi) + alp*TV(log v) 
% after change of variable(logv=w): min_w sum_i t wi + xi exp(-wi) + alp*TV(w)
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
	% v1 solves min_v t vi + xi exp(-vi) + lam/2(v-w+y)^2
	% by Newton's method solving t - lam(wi-yi) = -lam vi + xi exp(-vi)
	v1 = -Newton_step(v0, x, lam, t - lam*(w0-y0), epsl * 0.1);
    w1 = TVL2(v1+y0, alp/lam / 0.78539816339744830962, 4, 0);
    y1 = y0 + v1 - w1;
    err = max([sum((y1(:)-y0(:)).^2); sum((v1(:)-v0(:)).^2); sum((w1(:)-w0(:)).^2); sum((v1(:)-w1(:)).^2)]);
    if err < epsl
        break;
    end
    fprintf('iter %d: error: %.4f\n', i, err);
    y0 = y1; v0=v1; w0 = w1;
end
v = exp(v1);
end

function p = Newton_step(p0, t, lam, c, epsl)
% this function implements Newton's method for te^p + lam*p = c (pointwisely)
% % INPUTS: 
%       p0:     n1*n2 matrix, initialization
%       t:      n1*n2 matrix
%       lam:    scalar
%       c:      n1*n2 matrix
%       epsl:   stopping tolerate
% % OUTPUT:
%       p:      n1*n2 matrix, solution
%
N = 100000;
for i=1:N
    p1 = p0 - (t.*exp(p0) + lam.*p0 - c) ./ (t.*exp(p0) + lam);
    err = sum((p1(:)-p0(:)).^2);
    if err < epsl
        break;
    end
    p0 = p1;
end
if i == N
    fprintf('Newton does not converge.\n');
end
p = p1;
end