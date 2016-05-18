function [k] = max_snr(pr1, snr1)
% SYNTAX:
%   k = max_snr(pr1, snr1)
% INPUT:
%   pr1  = L1 code pseudorange for the current epoch (all satellites in the
%          constellation)
%   snr1 = SNR for L1 signals for all signals for the current epoch
%
% OUTPUT:
%   k    = index number for the selected satellites
%
% DESCRIPTION:
% This function selects the four satellites with the highest SNR to be used
% for positioning using the distance geometry method.
%
% -------------------------------------------------------------------------
%
%   Copyright 2016, Hadi Tabatabaee, all rights reserved.
%
% -------------------------------------------------------------------------


i = find(~(pr1(:)==0)); % find pseudorange indexes that are non-zero
SNR = snr1(i);

% sort based on SNR (low to high):
[~, I] = sort(SNR);

k = flipud(I); % flip order to get k
end
