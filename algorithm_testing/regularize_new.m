function [x, L] = regularize_new(A,b,lambda)
% New regularization matrix method for Tikhonov regularization
% A NEW TIKHONOV REGULARIZATION METHOD -MARTIN FUHRY AND LOTHAR REICHEL

n = size(A,1);
[U S V] = svd(A);

D2 = zeros(n,n);
for i = 1:n
    D2(i,i) = max((lambda^2 - S(i,i)^2), 0);
end
b_tilde = U'*b;

% solution:
x = V*(inv(S'*S + D2)*S'*b_tilde);
