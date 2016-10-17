function Ainv = genInv(A)
% Ben Noble Generalized Inverse

n = size(A,1); % number of satellites
A11 = A(1:4, 1:4); % non-singular matrix
A12 = A(1:4, 5:n);
A21 = A(5:n, 1:4);
A22 = A(5:n, 5:n);

% A12 = A11*Q so:
Q = A11\A12;

% A12 = A21' and A21 = P*A11 --> P = Q'*A11'*All^-1
P = Q'*A11*pinv(A11);

Ainv = [eye(4); Q']*inv((eye(4) + P'*P)*A11*(eye(4) + Q*Q'))*[eye(4) P'];

end