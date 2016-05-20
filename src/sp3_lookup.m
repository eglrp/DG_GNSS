function [XS, dtS] = sp3_lookup(sp3file, sat, time)
% SYNTAX: [XS, dtS] = sp3_lookup(sp3file, sat)
%
% INPUT:
%   sp3file = Precise ephemeride SP3 file
%   sat     = PRN for satellites of interest
%   time    = time of epoch of interest
%
% OUTPUT:
%   XS      = matrix containing satellite coordinates (ECEF) for the PRN's
%             marked in sat
%   dtS     = satellite clock error
% 
% DESCRIPTION:
%   Reads satellite positions and clock errors at a single epoch using the
%   precise ephemerides.
%
% -------------------------------------------------------------------------
%
%   Copyright 2016, Hadi Tabatabaee, all rights reserved.
%
% -------------------------------------------------------------------------
