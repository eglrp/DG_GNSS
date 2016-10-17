% 2016  3 23 21  0  0.00000000
% Approx. XYZ : -2693674.7660 -4273831.4480  3880383.2100
pos = [-2693674.7660; -4273831.4480; 3880383.2100];

% Satellite PRN's
sat = [3; 5; 7; 8; 9; 16; 23; 27; 28; 30];
 
% PR  =  distance +
%        c * (receiver clock offset - satellite clock offset +
%                 other biases)



% L1 pseudorange (m)
PR = [24796342.320; 25190885.640; 21051150.540; 23421564.560; ...
    20185617.780; 23130311.660; 20937631.020; 23236950.980; ...
    24745795.000; 22689126.620]; 

% Satellite positions in ECEF (m)
XS = [-9853659.976 -20962339.661 -12986790.632;
      -9266243.080  12439812.263  21487001.848;
      -17396304.120  -8149173.727  18639706.887;
      5862868.051 -25873417.197    980370.260;
      -9701382.340 -16809642.699  18105841.341;
      10061719.918 -10449386.830  22149509.102;
      -2767399.415 -24088669.699  10488753.345;
      11825676.781 -21058403.483  10842200.774;
      -20982809.734 -10574143.312 -11796671.097;
      -24139911.085  -1901619.513  11002104.420];
 
 % Satellite clock bias (microsec) 
 dtS = [-39.492412; -135.077317; 463.797220; -28.615227; 95.719609; ...
      -29.944114; -166.059303; 74.938506; 511.165127; 109.725965];
  
 % phase(corr) = phase(r) -  dT(r)*freq 
  
 % L1 Carrier phase
 Ph1 = [130305654.375; 132378983.006; 110624620.259; 123081202.207; ...
     106076176.016; 121550668.474; 110028009.908; 122111060.630; ...
     130040039.748; 119232196.073];
 
 % L2 Carrier phase
 Ph2 = [101536910.34844; 103152461.29945; 86201008.51149; 95907479.42146; ...
     82656791.17049; 94714828.27745; 85736089.91848; 95151519.20447; ...
     101329927.53244; 92908238.72547];
 
 % P2
 P2 = [24796350.780; 25190888.160; 21051150.260; 23421569.400; ...
     20185619.020; 23130312.240; 20937628.860; 23236954.460; ...
     24745801.000; 22689129.940];

 %
 % LX:    >= 25dBHz -> 1; 26-27dBHz -> 2; 28-31dBHz -> 3       COMMENT
 %      32-35dBHz -> 4; 36-38dBHz -> 5; 39-41dBHz -> 6       COMMENT
 %     42-44dBHz -> 7; 45-48dBHz -> 8; >= 49dBHz -> 9 
 
 SNR1 = [7; 6; 9; 8; 9; 8; 9; 8; 5; 8];
 SNR2 = [44; 45; 49; 46; 49; 45; 48; 47; 44; 47];
  
 % Calculating elevations based on pos and XS:
%  for i = 1:numel(PR)
%      [Az(i,1), El(i,1)] = topocent(pos, XS(i,:)' - pos);
%  end
%  
 %% DISTANCE GEOMETRY SOLUTION -- ALL SATELLITES
% This solution has not removed the clock bias.

r = [PR(1:10);0] ;
u_bar_r = diag(PR(1:10)*PR(1:10)');
u_bar = [u_bar_r; 1];
hr = ones(numel(PR), 1);
h = [hr; 0];
e_n_unit = zeros(numel(PR)+1,1);
k = 100;

P = eye(numel(PR)) - 1/(numel(PR))*hr*hr';

Ar = zeros(numel(PR));
 
 for i = 1:numel(PR)
     for j = 1:numel(PR)
         Ar(i,j) = (XS(i,:) - XS(j,:))*(XS(i,:) - XS(j,:))';
     end
 end

 A = [Ar, hr; hr', 0];
 % X matrix
 x_r_bar = pinv(P*Ar)*P*u_bar_r;
 x_f_bar = (1/numel(PR))*hr'*u_bar_r - (1/numel(PR))*hr'*Ar*x_r_bar;
 
 x_bar = [x_r_bar; x_f_bar];
 
 %XR_m = XS'*x_r_bar;
 
 % Calculate minimal norm reduced xu_tilda_r:
% xu_bar_r = pinv(P * Ar + k * hr * hr') * (P * u_bar_r + k * hr);

% Find its last element:
%xu_tilda_last = hr' * (r(1:10) - Ar * xu_tilda_r) / n;

% Put the two together:
%xu_tilda = [xu_tilda_r; xu_tilda_last];

% Calculate the minimal norm xr_tilda_r (no need to find the last element):
%xr_tilda_r = pinv(P * Ar + k * hr * hr') * P * r(1:10);

% Find its last element:
%xr_tilda_last = hr' * (r_tilda_r - Ar * xr_tilda_r)/n;
%xr_tilda = [xr_tilda_r; xr_tilda_last];

% Intermediate results to calculate clock bias:
rAr = x_r_bar' * r(1:10);
rAu = x_r_bar' * u_bar_r; %xu_tilda' * r_tilda;
uAu = x_bar' * u_bar;

if x_bar'*r < 0
    bias = (rAu + sqrt(rAu^2 - uAu * (.5 + rAr)))/(.5 + rAr)/2;
else
    bias = (rAu - sqrt(rAu^2 - uAu * (.5 + rAr)))/(.5 + rAr)/2;
end

%bias_small = uAu / rAu /4;

% Calculate xr:
xr = x_r_bar - 2 * bias * x_r_bar;


% Receiver Position Coordinates:
disp ' Receiver Position Coordinates (meters): '
P_xyz = XS' * xr
 
 
 
%----------------------------------------------------------

% guess b:
n = numel(PR);
b = 0;
b_pre = 0;
iter = 0
r = PR(1:10);
u = u_bar_r;
tol = 25;
%k_const = (1/n)*norm(Ar, 1);
k_const = 100;

while abs(b - b_pre) > 0.0001
    iter = iter + 1;
    b_pre = b;
    x = pinv(P*Ar + k_const*hr*hr', tol)*(P*(u-2*b*r) + k_const*hr);
    xu = pinv(P*Ar + k_const*hr*hr', tol)*(P*u + k_const*hr);
    xr = pinv(P*Ar + k_const*hr*hr', tol)*(P*r + k_const*hr);
    
    if xr'*r < 0
        b = ((xu'*r + xr'*u) + sqrt(((xu'*r + xr'*u)^2) -2*(1 + 2*xr'*r)*(xu'*u)))/(2 + 4*xr'*r);
    else
        b = ((xu'*r + xr'*u) - sqrt(((xu'*r + xr'*u)^2) -2*(1 + 2*xr'*r)*(xu'*u)))/(2 + 4*xr'*r);
    end
end
 
P_new = XS' * x;
%%

x_ur = pinv(P*Ar + k_const*hr*hr', tol)*(P*u + k_const*hr);
x_ur_last = (1/n)*hr'*(u - Ar*x_ur);

x_u = [x_ur; x_ur_last];

x_rr = pinv(P*Ar + k_const*hr*hr', tol)*(P*r+ k_const*hr);
x_rr_last = (1/n)*hr'*(r - Ar*x_rr);

x_r = [x_rr; x_rr_last];

cdt_plus = (x_r'*[u;1] + sqrt((x_r'*[u;1])^2 - x_u'*[u;1]*(0.5 + x_r'*[r;0])))/(0.5*(0.5 + x_r'*[r;0]));

cdt_minus = (x_r'*[u;1] - sqrt((x_r'*[u;1])^2 - x_u'*[u;1]*(0.5 + x_r'*[r;0])))/(0.5*(0.5 + x_r'*[r;0]));

if x_r'*[r;0] < 0
    cdt = cdt_plus;
else
    cdt = cdt_minus;
end

x_calc_r = x_ur - 2*cdt*x_rr;

P_new2 = XS' * x_calc_r;
%%

x_ur = pinv(P*Ar + k_const*hr*hr', tol)*(P*u + k_const*hr);
x_ur_last = (1/n)*hr'*(u - Ar*x_ur);

x_u = [x_ur; x_ur_last];

x_rr = pinv(P*Ar + k_const*hr*hr', tol)*(P*r+ k_const*hr);
x_rr_last = (1/n)*hr'*(r - Ar*x_rr);

x_r = [x_rr; x_rr_last];


if x_r'*[r;0] < 0
    dbt = ((x_u'*[r;0] + x_r'*[u;1]) + sqrt(((x_u'*[r;0] + x_r'*[u;1])^2) -2*(1 + 2*x_r'*[r;0])*(x_u'*[u;1])))/(2 + 4*x_r'*[r;0]);
else
    dbt = ((x_u'*[r;0] + x_r'*[u;1]) - sqrt(((x_u'*[r;0] + x_r'*[u;1])^2) -2*(1 + 2*x_r'*[r;0])*(x_u'*[u;1])))/(2 + 4*x_r'*[r;0]);
end
    
x_c_r = x_ur + 2*dbt*x_rr;

P_new3 = XS' * x_c_r;