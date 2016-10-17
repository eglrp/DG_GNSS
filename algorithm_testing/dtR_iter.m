load('sat_final.mat');
load('signal.mat');

time_steps = numel(time);

pr = PR(k,:);
pr1 = PR(k,:);
n = numel(pr(:,1));
c = 299792458;

for t = 1:time_steps
    for i = 1:50
        [XR(:,t), XR_cov(:,:,t), dtR(t)] = DG_sol(XS(:,:,t), pr1(:,t));
        pr1(:,t) = pr1(:,t) - c*dtR(t)*ones(n,1);
    end
    [XR_LS(:,t), XR_cov_LS(:,:,t), dtR_LS(t)] = LS_dtR(XS(:,:,t), pr(:,t), dtR(t), XR(:,t));
end

    