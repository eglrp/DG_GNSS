function [x,L, past] = regularize(A,b)
% provide Tikhonov regularization solution

%[U S V] = svd(A);
n = size(A,1);
lim_low = 0;
lim_hi = max(max(A));
res_p = norm(b);
iter = 0;
past = [];

lambda1 = (lim_hi - lim_low)/4;
lambda2 = 3*(lim_hi - lim_low)/4;
L1 = lambda1*eye(n);
x1 = inv(A'*A + L1'*L1)*A'*b;
L2 = lambda2*eye(n);
x2 = inv(A'*A + L2'*L2)*A'*b;
res1 = norm(A*x1 - b);
res2 = norm(A*x2 - b);
if res1 < res2
    lambda = lambda1;
    res_p = res1;
    lim_hi = lambda;
else
    lambda = lambda2;
    res_p = res2;
    lim_low = lambda;
end

while iter < 50    
    lambda = (lim_hi - lim_low)/2 + lim_low;
    L = lambda*eye(n);
    x = inv(A'*A + L'*L)*A'*b;
    res = norm(A*x - b);
    if res < res_p
        lim_low = lambda;
    else
        r = lim_hi - lim_low;
        lim_hi = lim_low;
        lim_low = lim_low - r;
    end
    res_p = res;
    past = [past res];
    iter = iter + 1;
end



