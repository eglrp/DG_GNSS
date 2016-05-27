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
 for i = 1:numel(PR)
     [Az(i,1), El(i,1)] = topocent(pos, XS(i,:)' - pos);
 end
 
  
  %% DISTANCE GEOMETRY METHOD -- FUNCTIONS
  [dtR_t, A_t, Ainv_t] = DG_SA_code_clock(XS(1:4,:), dtS(1:4), [], [], PR, [1;2;3;4]);
  [XR_t] = DG_SA_code(XS(1:4,:), PR, dtR_t, dtS(1:4), [0;0;0;0], [0;0;0;0], A_t, Ainv_t, [1;2;3;4]);
  
  
  
  %% DISTANCE GEOMETRY METHOD -- FOUR SATELLITES ONLY
  
 r = [PR(1:4);0] ;
 u_bar = [diag(PR(1:4)*PR(1:4)'); 1];
 
 QR2 = (XS(1,:) - XS(2,:))*(XS(1,:) - XS(2,:))';
 QS2 = (XS(1,:) - XS(3,:))*(XS(1,:) - XS(3,:))'; 
 QT2 = (XS(1,:) - XS(4,:))*(XS(1,:) - XS(4,:))';
 RS2 = (XS(2,:) - XS(3,:))*(XS(2,:) - XS(3,:))';
 RT2 = (XS(2,:) - XS(4,:))*(XS(2,:) - XS(4,:))';
 ST2 = (XS(3,:) - XS(4,:))*(XS(3,:) - XS(4,:))';

 A = [0, QR2, QS2, QT2, 1;
      QR2, 0, RS2, RT2, 1;
      QS2, RS2, 0, ST2, 1;
      QT2, RT2, ST2, 0, 1;
      1, 1, 1, 1, 0];
  
 Ainv = inv(A);
 
 if r'*Ainv*u_bar > 0;
    dtR = (r'*Ainv*u_bar - sqrt((r'*Ainv*u_bar)^2 - u_bar'*Ainv*u_bar*(0.5 + r'*Ainv*r)))/(2*(0.5 + r'*Ainv*r));
 else
    dtR = (r'*Ainv*u_bar + sqrt((r'*Ainv*u_bar)^2 - u_bar'*Ainv*u_bar*(0.5 + r'*Ainv*r)))/(2*(0.5 + r'*Ainv*r));
 end
 
 PR_corr = PR(1:4) - dtR*[1;1;1;1] - (10^-6)*dtS(1:4);
 
 u = [diag(PR_corr(1:4)*PR_corr(1:4)'); 1];
 
 S = [XS(1:4, :)' , [0;0;0]];
 
 XR = S*Ainv*u;
 
 
 % New choice of satellites:
 
 r = [PR(3:6);0] ;
 u_bar = [diag(PR(3:6)*PR(3:6)'); 1];
 
 QR2 = (XS(3,:) - XS(4,:))*(XS(3,:) - XS(4,:))';
 QS2 = (XS(3,:) - XS(5,:))*(XS(3,:) - XS(5,:))'; 
 QT2 = (XS(3,:) - XS(6,:))*(XS(3,:) - XS(6,:))';
 RS2 = (XS(4,:) - XS(5,:))*(XS(4,:) - XS(5,:))';
 RT2 = (XS(4,:) - XS(6,:))*(XS(4,:) - XS(6,:))';
 ST2 = (XS(5,:) - XS(6,:))*(XS(5,:) - XS(6,:))';

 A = [0, QR2, QS2, QT2, 1;
      QR2, 0, RS2, RT2, 1;
      QS2, RS2, 0, ST2, 1;
      QT2, RT2, ST2, 0, 1;
      1, 1, 1, 1, 0];
  
 Ainv = inv(A);
 
 if r'*Ainv*u_bar > 0;
    dtR1 = (r'*Ainv*u_bar - sqrt((r'*Ainv*u_bar)^2 - u_bar'*Ainv*u_bar*(0.5 + r'*Ainv*r)))/(2*(0.5 + r'*Ainv*r));
 else
    dtR1 = (r'*Ainv*u_bar + sqrt((r'*Ainv*u_bar)^2 - u_bar'*Ainv*u_bar*(0.5 + r'*Ainv*r)))/(2*(0.5 + r'*Ainv*r));
 end
 
 PR_corr = PR(3:6) - dtR1*[1;1;1;1] - (10^-6)*dtS(3:6);
 
 u = [diag(PR_corr(1:4)*PR_corr(1:4)'); 1];
 
 S = [XS(3:6, :)' , [0;0;0]];
 
 XR1 = S*Ainv*u;
 
 % A third choice
 r = [PR(7:10);0] ;
 u_bar = [diag(PR(7:10)*PR(7:10)'); 1];
 
 QR2 = (XS(7,:) - XS(8,:))*(XS(7,:) - XS(8,:))';
 QS2 = (XS(7,:) - XS(9,:))*(XS(7,:) - XS(9,:))'; 
 QT2 = (XS(7,:) - XS(10,:))*(XS(7,:) - XS(10,:))';
 RS2 = (XS(8,:) - XS(9,:))*(XS(8,:) - XS(9,:))';
 RT2 = (XS(8,:) - XS(10,:))*(XS(8,:) - XS(10,:))';
 ST2 = (XS(9,:) - XS(10,:))*(XS(9,:) - XS(10,:))';

 A = [0, QR2, QS2, QT2, 1;
      QR2, 0, RS2, RT2, 1;
      QS2, RS2, 0, ST2, 1;
      QT2, RT2, ST2, 0, 1;
      1, 1, 1, 1, 0];
  
 Ainv = inv(A);
 
 if r'*Ainv*u_bar > 0;
    dtR2 = (r'*Ainv*u_bar - sqrt((r'*Ainv*u_bar)^2 - u_bar'*Ainv*u_bar*(0.5 + r'*Ainv*r)))/(2*(0.5 + r'*Ainv*r));
 else
    dtR2 = (r'*Ainv*u_bar + sqrt((r'*Ainv*u_bar)^2 - u_bar'*Ainv*u_bar*(0.5 + r'*Ainv*r)))/(2*(0.5 + r'*Ainv*r));
 end
 
 PR_corr = PR(7:10) - dtR2*[1;1;1;1] - (10^-6)*dtS(7:10);
 
 u = [diag(PR_corr(1:4)*PR_corr(1:4)'); 1];
 
 S = [XS(7:10, :)' , [0;0;0]];
 
 XR2 = S*Ainv*u;
 
 % A fourth choice
 r = [PR([2,4,6,8]);0] ;
 u_bar = [diag(PR([2,4,6,8])*PR([2,4,6,8])'); 1];
 
 QR2 = (XS(2,:) - XS(4,:))*(XS(2,:) - XS(4,:))';
 QS2 = (XS(2,:) - XS(6,:))*(XS(2,:) - XS(6,:))'; 
 QT2 = (XS(2,:) - XS(8,:))*(XS(2,:) - XS(8,:))';
 RS2 = (XS(4,:) - XS(6,:))*(XS(4,:) - XS(6,:))';
 RT2 = (XS(4,:) - XS(8,:))*(XS(4,:) - XS(8,:))';
 ST2 = (XS(6,:) - XS(8,:))*(XS(6,:) - XS(8,:))';

 A = [0, QR2, QS2, QT2, 1;
      QR2, 0, RS2, RT2, 1;
      QS2, RS2, 0, ST2, 1;
      QT2, RT2, ST2, 0, 1;
      1, 1, 1, 1, 0];
  
 Ainv = inv(A);
 
 if r'*Ainv*u_bar > 0;
    dtR3 = (r'*Ainv*u_bar - sqrt((r'*Ainv*u_bar)^2 - u_bar'*Ainv*u_bar*(0.5 + r'*Ainv*r)))/(2*(0.5 + r'*Ainv*r));
 else
    dtR3 = (r'*Ainv*u_bar + sqrt((r'*Ainv*u_bar)^2 - u_bar'*Ainv*u_bar*(0.5 + r'*Ainv*r)))/(2*(0.5 + r'*Ainv*r));
 end
 
 PR_corr = PR([2,4,6,8]) - dtR3*[1;1;1;1] - (10^-6)*dtS([2,4,6,8]);
 
 u = [diag(PR_corr(1:4)*PR_corr(1:4)'); 1];
 
 S = [XS([2,4,6,8], :)' , [0;0;0]];
 
 XR3 = S*Ainv*u;

% Averaging the four estimates
XR_avg = XR + XR1 + XR2 + XR3;
XR_avg = XR_avg./4;
 
%% LEAST-SQUARES SOLUTION -- FOUR SATELLITES

[XR_ls, dtR_ls, A_ls, b_ls] = leastSquare([1;2;3;4], XS(1:4,:), PR(1:4));

%% LEAST-SQUARES SOLUTION -- ALL SATELLITES

[XR_ls2, dtR_ls2, A_ls2, b_ls2] = leastSquare([1;2;3;4;5;6;7;8;9;10], XS, PR);


%% DISTANCE GEOMETRY SOLUTION -- ALL SATELLITES
% This solution has not removed the clock bias.

r = [PR(1:10);0] ;
u_bar_r = diag(PR(1:10)*PR(1:10)');
u_bar = [u_bar_r; 1];
h_r = ones(numel(PR), 1);

P = eye(numel(PR)) - 1/(numel(PR))*h_r*h_r';

A = zeros(numel(PR)+1);
 
 for i = 1:numel(PR)
     for j = 1:numel(PR)
         A(i,j) = (XS(i,:) - XS(j,:))*(XS(i,:) - XS(j,:))';
     end
     A(i,numel(PR)+1) = 1;
     A(numel(PR)+1, i) = 1;
 end
 
 A_r = A(1:numel(PR), 1:numel(PR));
 
 % X matrix
 x_r_bar = pinv(P*A_r)*P*u_bar_r;
 x_f_bar = (1/numel(PR))*h_r'*u_bar_r - (1/numel(PR))*h_r'*A_r*x_r_bar;
 
 x_bar = [x_r_bar; x_f_bar];
 
 XR_m = XS'*x_r_bar;
 
 %% Multi-sat select based on elevation and SNR
snr_ind = find(SNR1 >= 7);
el_ind = find(El >= 20);

for i = 1:numel(el_ind)
    ind(i) = find(el_ind(i) == snr_ind);
end

r = [PR(ind);0] ;
u_bar_r = diag(PR(ind)*PR(ind)');
u_bar = [u_bar_r; 1];
h_r = ones(numel(ind), 1);

P = eye(numel(ind)) - 1/(numel(ind))*h_r*h_r';

A = zeros(numel(ind)+1);
 
 for i = 1:numel(ind)
     for j = 1:numel(ind)
         A(i,j) = (XS(ind(i),:) - XS(ind(j),:))*(XS(ind(i),:) - XS(ind(j),:))';
     end
     A(i,numel(ind)+1) = 1;
     A(numel(ind)+1, i) = 1;
 end
 
 A_r = A(1:numel(ind), 1:numel(ind));
 
 % X matrix
 x_r_bar2 = pinv(P*A_r)*P*u_bar_r;
% x_f_bar2 = (1/numel(ind))*h_r'*u_bar_r - (1/numel(ind))*h_r'*A_r*x_r_bar;
 
 %x_bar2 = [x_r_bar; x_f_bar];
 
 XR_m2 = XS(ind,:)'*x_r_bar2;
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
   
  
  
 
  
  