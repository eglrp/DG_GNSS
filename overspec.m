% Over-specified clock bias and receiver coordinates calculation

clear;
format long g;
n = 6;
hr = [1; 1; 1; 1; 1; 1];
h = [hr; 0];
e7_unit = [0; 0; 0; 0; 0; 0; 1];
k = 100;

% Note: scale all distances by 1e-7:
disp ' For proper scaling all distances are divided by 1e7.'

% pseudoranges vector:
r_tilda = [21119278.32; 22527064.18; 23674159.88; 20951647.38; ...
		   20155401.42; 24222110.91; 0] * 1e-7;

% Reduced r_tilda:
r_tilda_r = r_tilda(1:6);

% Squared pseudo-ranges vector:
u_tilda = diag(r_tilda) * r_tilda +e7_unit;

u_tilda_r = u_tilda(1:6);

% Coordinates of satellites S1, S2, .., S6:

s1 = [14177553.47, -18814768.09, 12243866.38] * 1e-7;
s2 = [15097199.81,  -4636088.67, 21326706.55] * 1e-7;
s3 = [23460342.33,  -9433518.58,  8174941.25] * 1e-7;
s4 = [-8206488.95, -18217989.14, 17605231.99] * 1e-7;
s5 = [ 1399988.07, -17563734.90, 19705591.18] * 1e-7;
s6 = [ 6995655.48, -23537808.26, -9927906.48] * 1e-7;

% Satellite Coordinate matrix Sr:
Sr = [s1', s2', s3', s4', s5', s6'] * 1e7;

s12 = (s1-s2) * (s1-s2)';
s13 = (s1-s3) * (s1-s3)';
s14 = (s1-s4) * (s1-s4)';
s15 = (s1-s5) * (s1-s5)';
s16 = (s1-s6) * (s1-s6)';

s23 = (s2-s3) * (s2-s3)';
s24 = (s2-s4) * (s2-s4)';
s25 = (s2-s5) * (s2-s5)';
s26 = (s2-s6) * (s2-s6)';

s34 = (s3-s4) * (s3-s4)';
s35 = (s3-s5) * (s3-s5)';
s36 = (s3-s6) * (s3-s6)';

s45 = (s4-s5) * (s4-s5)';
s46 = (s4-s6) * (s4-s6)';

s56 = (s5-s6) * (s5-s6)';

% Reduced A matrix:

Ar = [0 s12 s13 s14 s15 s16;
	  s12 0 s23 s24 s25 s26;
	  s13 s23 0 s34 s35 s36;
	  s14 s24 s34 0 s45 s46;
	  s15 s25 s35 s45 0 s56;
	  s16 s26 s36 s46 s56 0];

% A = Satellite reference matrix:

A = [Ar, hr; hr', 0];

P = eye(n) - hr*hr'/n; % P matrix

% Calculate minimal norm reduced xu_tilda_r:
xu_tilda_r = pinv(P * Ar + k * hr * hr') * (P * u_tilda_r + k * hr);

% Find its last element:
xu_tilda_last = hr' * (r_tilda_r - Ar * xu_tilda_r) / n;

% Put the two together:
xu_tilda = [xu_tilda_r; xu_tilda_last];

% Calculate the minimal norm xr_tilda_r (no need to find the last element):
xr_tilda_r = pinv(P * Ar + k * hr * hr') * P * r_tilda_r;

% Find its last element:
xr_tilda_last = hr' * (r_tilda_r - Ar * xr_tilda_r)/n;
xr_tilda = [xr_tilda_r; xr_tilda_last];

% Intermediate results to calculate clock bias:
rAr = xr_tilda' * r_tilda;
rAu = xr_tilda' * u_tilda; %xu_tilda' * r_tilda;
uAu = xu_tilda' * u_tilda;

% Clock Bias with + sqrt(discriminant)
bias_plus = (rAu + sqrt(rAu^2 - uAu * (.5 + rAr)))/(.5 + rAr)/2;

% Clock Bias with - sqrt(discriminant)
% bias_minus is admissible here since xu_tilda' * r_tilda > 0.

bias_minus = (rAu - sqrt(rAu^2 - uAu * (.5 + rAr)))/(.5 + rAr)/2;
bias_small = uAu / rAu /4;

% Calculate xr:
xr = xu_tilda_r - 2 * bias_minus * xr_tilda(1:6);

disp 'Receiver Clock Bias (meters): '
bias = bias_minus * 1e7

% Receiver Position Coordinates:
disp ' Receiver Position Coordinates (meters): '
P_xyz = Sr * xr

