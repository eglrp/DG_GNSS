function [U_SV, V_SV, W_SV] = velocitySum(XS1, XS2)
% This finds the sum velocity vector for each satellite in the two previous
% frames (assumes satellite index remains constant)

for i = 1:size(XS1,2)
    U_SV(i) = XS2(1,i) - XS1(1,i);
    V_SV(i) = XS2(2,i) - XS1(2,i);
    W_SV(i) = XS2(3,i) - XS1(3,i);
end
