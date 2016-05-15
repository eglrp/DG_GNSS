function [XS_p] = satellite_projection(XS, VS_tx, XS_tx, dtS, time_tx, time_interval, traveltime)

% INPUT:
%   XS            = satellite position at transmission time in ECEF(time_rx) (X,Y,Z)
%   dtS           = satellite clock error (vector)
%   XS_tx         = satellite position at transmission time in ECEF(time_tx) (X,Y,Z)
%   VS_tx         = satellite velocity at transmission time in ECEF(time_tx) (X,Y,Z)
%   time_tx       = transmission time (vector)
%   time_interval = time difference between current and previous epochs
%   traveltime    = signal travel time (from previous epoch)
%
% OUTPUT:
%   XS_p    = Projected satellite position in ECEF(time_rx) (X,Y,Z)
%
% Description:
%   This function estimates the satellite position at the next epoch. These
%   positions can be used for receiver clock bias estimation. The goal is
%   to use these to recalculate the satellite positions once the clock bias
%   has been estimated.
%
%
% * Conversion from ECEF(time_tx) to ECEF(time_rx) is done using code
% adapted from the goGPS_MATLAB open-source project 
% (https://gogps-project.org).
%
% -------------------------------------------------------------------------
%
%   Copyright 2016, Hadi Tabatabaee, All rights reserved.
%
% -------------------------------------------------------------------------

% Projected satellite position in ECEF(time_tx)
XS_p_tx = XS_tx + time_interval*(VS_tx);

% Convert to ECEF(time_rx) -- adapted from goGPS_MATLAB (for GPS only)
Omegae_dot = goGNSS.OMEGAE_DOT_GPS;
XS_p(i,:) = earth_rotation_correction(traveltime, XS_p_tx(i,:), Omegae_dot);

