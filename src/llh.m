function [phi, lambda, h] = llh(X, Y, Z)
% converts WGS84 ECEF coordinates to the geometric reference
% frame.

% Geometry of earth ellipsoid based on WGS84:
% semimajor axis radius (m)
a = 6378137.0;

% flattening (s2/m3)
f = (1/298.257223563);

% semi-minor axis radius (m)
b = a*(1-f);

% first eccentricity (squared):
e2 = (a^2 - b^2)/(a^2);

% second eccentricity:
%e_p = sqrt((a/b)^2 -1);

epsilon = 1;
p = sqrt(X^2 + Y^2);
phi_0 = atan((Z/p)*((1 - e2)^(-1)));
phi_1 = phi_0; 
k = 0;
while epsilon > 0.000000000000000000000000000000000000000000001
    k = k+1;
    phi_0 = phi_1;
    N_0 = (a^2)/sqrt((a*cos(phi_0))^2 + (b*sin(phi_0))^2);
    h = p/(cos(phi_0)) - N_0;
    phi_1 = atan((Z/p)*((1 - (e2*(N_0)/(N_0+h)))^(-1)));
    epsilon = abs(phi_1 - phi_0);
end
phi = radtodeg(phi_0);
if sign(X) == sign(Y)
    if sign(X) == 1
        lambda = radtodeg(atan(Y/X));
    else
        lambda = radtodeg(atan(Y/X)) - 180;
    end
elseif sign(X) == -1
    lambda = radtodeg(atan(Y/X)) + 180;
else
    lambda = radtodeg(atan(Y/X));
end
end