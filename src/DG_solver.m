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

PR1 = PR * 1e-7;

% Satellite positions in ECEF (m)
S = [-9853659.976 -20962339.661 -12986790.632;
      -9266243.080  12439812.263  21487001.848;
      -17396304.120  -8149173.727  18639706.887;
      5862868.051 -25873417.197    980370.260;
      -9701382.340 -16809642.699  18105841.341;
      10061719.918 -10449386.830  22149509.102;
      -2767399.415 -24088669.699  10488753.345;
      11825676.781 -21058403.483  10842200.774;
      -20982809.734 -10574143.312 -11796671.097;
      -24139911.085  -1901619.513  11002104.420];
  
XS = S .* 1e-7;
 
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

n = numel(PR1);

r_b_r = PR1(1:10);
r_b = [r_b_r; 0];
u_b_r = diag(PR1(1:10)*PR1(1:10)');
u_b = [u_b_r; 1];
h_r = ones(n, 1);
h = [h_r; 0];
e_n_unit = zeros(n+1,1);


P = eye(n) - 1/n * h_r*h_r';

A_r = zeros(n);
 
for i = 1:n
    for j = 1:n
        A_r(i,j) = (XS(i,:) - XS(j,:))*(XS(i,:) - XS(j,:))';
    end
end
A = [A_r, h_r; h_r', 0];

A_r_norm_1 = norm(A_r,1);

%k = 100; 
k = 1/n * norm(A_r,1);
 
x_ub_r = pinv(A_r)*u_b_r;
%x_ub_r = pinv(P * A_r + k * h_r * h_r') * (P * u_b_r + k * h_r);
%x_ub_last = h_r' * (r_b_r - A_r * x_ub_r) / n;
%x_ub = [x_ub_r; x_ub_last];
 
x_rb_r = pinv(P*A_r + A_r_norm_1 * h_r*h_r') * P*r_b_r;
%x_rb_r = pinv(P * A_r + k * h_r * h_r') * P * r_b_r;
%x_rb_last = h_r' * (r_b_r - A_r * x_rb_r)/n;
%x_rb = [x_rb_r; x_rb_last];
 
 
% if x_ub' * r_b < 0
%    dtR_dist = ((x_ub' * r_b) + sqrt((x_ub' * r_b)^2 - (1 + x_rb' * r_b) * (x_ub' * u_b)))/(1 + 2 * x_rb' * r_b);
% else
%    dtR_dist = ((x_ub' * r_b) - sqrt((x_ub' * r_b)^2 - (1 + x_rb' * r_b) * (x_ub' * u_b)))/(1 + 2 * x_rb' * r_b);
% end

if x_ub_r' * r_b_r < 0
   dtR_dist = ((x_ub_r' * r_b_r) + sqrt((x_ub_r' * r_b_r)^2 - (1 + x_rb_r' * r_b_r) * (x_ub_r' * u_b_r)))/(1 + 2 * x_rb_r' * r_b_r);
else
   dtR_dist = ((x_ub_r' * r_b_r) - sqrt((x_ub_r' * r_b_r)^2 - (1 + x_rb_r' * r_b_r) * (x_ub_r' * u_b_r)))/(1 + 2 * x_rb_r' * r_b_r);
end
 
 
%x_r = x_ub_r - 2 * dtR_dist * x_rb_r;
    
 
x_r = pinv(P*A_r + 1/n * A_r_norm_1 * h_r*h_r') * (P*u_b_r + 1/n * A_r_norm_1*h_r);
 
%dtR_dist = -.5*(x_r - x_ub_r)./x_rb_r;
 
%dtR = dtR_dist*1e7/goGNSS.V_LIGHT;
 
XR = S'*x_r;
 
 
 
r = [PR(1:10);0] ;
u_bar_r = diag(PR(1:10)*PR(1:10)');
u_bar = [u_bar_r; 1];
h_r = ones(numel(PR), 1);

P = eye(numel(PR)) - 1/(numel(PR))*h_r*h_r';

A = zeros(numel(PR)+1);
 
 for i = 1:numel(PR)
     for j = 1:numel(PR)
         A(i,j) = (S(i,:) - S(j,:))*(S(i,:) - S(j,:))';
     end
     A(i,numel(PR)+1) = 1;
     A(numel(PR)+1, i) = 1;
 end
 
 A_r1 = A(1:numel(PR), 1:numel(PR));
 
 % X matrix
 x_r_bar = pinv(P*A_r1)*P*u_bar_r;
 x_f_bar = (1/numel(PR))*h_r'*u_bar_r - (1/numel(PR))*h_r'*A_r1*x_r_bar;
 
 x_bar = [x_r_bar; x_f_bar];
 
 XR_m = S'*x_r_bar;