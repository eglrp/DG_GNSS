function [X, Y, Z] = ecef(phi, lambda, h)
% converts the coordinates of a point expressed in the geographic coordinate system to the coordinates of the WGS84 ECEF.
% phi is latitude and lambda is longitude.

% Geometry of earth ellipsoid based on WGS84:
% semimajor axis radius (m)
a = 6378137.0;

% flattening (s2/m3)
f = (1/298.257223563);

% semi-minor axis radius (m)
b = a*(1-f);

% first eccentricity:
%e = sqrt(2*f - f^2);

% second eccentricity:
%e_p = sqrt((a/b)^2 -1);

% radius of curvature in prime vertical
N = a^2/sqrt((a*cos(degtorad(phi)))^2 + (b*sin(degtorad(phi)))^2);

% ECEF coordinates in WGS84 frame
X = (N + h)*cos(degtorad(phi))*cos(degtorad(lambda));
Y = (N + h)*cos(degtorad(phi))*sin(degtorad(lambda));
Z = (((b/a)^2)*N + h)*sin(degtorad(phi));

end