function [x] = conjGrad(A, b)
% Conjugate gradient optimization -- direct method
% Solves Ax = b for symmetric semidefinite A
n = size(A,1);
alpha_e = zeros(n,n);
% choose p to be unit vectors e_i
% alpha coefficients

for i = 1:n
    e = zeros(n,1);
    e(i) = 1;
    alpha_coef = dot(e,b)/dot(e,A*e);
    alpha_e(:,i) = alpha_coef*e;
end
x = sum(alpha_e,2);
end

