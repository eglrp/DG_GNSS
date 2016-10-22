function [indx, param] = pickSat(S)
% Picks the 'best' set of four satellite based on the largest tetrahedron volume.
% Copyright 2016, Hadi Tabatabaee. All rights reserved.
% -------------------------------------------------------------------------------
n = size(S,2);
ind = pick(1:n, 4,'');
param = [];
for k = 1:size(ind,1)
    XS = S(ind(k,:));
    for i = 1:4
        for j = 1:4
            tmp = XS(:,i) - XS(:,j);
            Ar(i,j) = tmp'*tmp;
        end
    end
    param = [param, abs(det(Ar))];
end

[~, l] = max(param);
indx = ind(l,:);