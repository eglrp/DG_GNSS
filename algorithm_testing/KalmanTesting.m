% This script tests the Kalman filter function.
% Copyright 2016, Hadi Tabatabaee, all rights reserved.

load('sat_final.mat');
load('signal.mat');

time_steps = numel(time);

global state state_cov

pr = PR(k,:);
nSats = size(pr,1);
ph_l1 = PH_L1(k,:);
ph_l2 = PH_L2(k,:);
snr = SNR(k,:);
Ar = zeros(nSats);
hr = ones(nSats,1);
h = [hr; 0];

for i = 1:time_steps

    [XR(:,i), XR_cov(:,:,i), dtR(i)] = kalmanFilter_P(XS(:,:,i), pr(:,i));

end