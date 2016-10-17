load('sat_final.mat');
load('signal.mat');

time_steps = numel(time);
pr = PR(k,:);
nSats = size(pr,1);
ph_l1 = PH_L1(k,:);
ph_l2 = PH_L2(k,:);
snr = SNR(k,:);
Ar = zeros(nSats);
hr = ones(nSats,1);
h = [hr; 0];

for i = 1:time_steps
    %[XR(:,i), XR_cov(:,:,i), dtR(i)] = KF_simple(XS(:,:,i), pr(:,i));
    [XR2(:,i), XR_cov2(:,:,i), dtR2(i)] = DG_sol(XS(:,:,i), pr(:,i));
    [XRg(:,i), dtRg(i)] = DG_sol2(XS(:,:,i), pr(:,i));
    [XR3(:,i), XR_cov3(:,:,i), dtR3(i)] = leastSquare2(XS(:,:,i), pr(:,i));
end

figure
hold on
plot(1:400, XRg(1,:))
plot(1:400, XR2(1,:))
plot(1:400, XR3(1,:))
