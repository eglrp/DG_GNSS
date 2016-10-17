function [XS, dtS] = sp3_linear(sp3file, sat, date)
% SYNTAX: [XS, dtS] = sp3_linear(sp3file, sat)
%
% INPUT:
%   sp3file = Precise ephemeride SP3 file name
%   sat     = PRN for satellites of interest
%   date    = date of epoch (YYYY, MM, DD, HH, MM, SS) 
%
% OUTPUT:
%   XS      = matrix containing satellite coordinates (ECEF) for the PRN's
%             marked in sat (m)
%   dtS     = satellite clock error (sec)
% 
% DESCRIPTION:
%   Reads satellite positions and clock errors at a single epoch using the
%   precise ephemerides. Uses linear interpolation to estimate position
%   between epochs.
%
% -------------------------------------------------------------------------
%
%   Copyright 2016, Hadi Tabatabaee, all rights reserved.
%
% -------------------------------------------------------------------------



XS = zeros(numel(sat), 3);
dtS = zeros(numel(sat), 1);

XS_l = zeros(numel(sat), 3);
dtS_l = zeros(numel(sat), 1);

XS_h = zeros(numel(sat), 3);
dtS_h = zeros(numel(sat), 1);


fid=fopen(sp3file,'r');

% Find first line
line = fgetl(fid);
while line(1) ~= '#'
    line = fgetl(fid);
end
line_num = 1;
sp3_version = line(2);

% Check if velocity and drift rate are available:
V_toggle = (line(3) == 'V');

% Start date and time
st_year = str2num(line(4:7));
st_month = str2num(line(9:10));
st_day = str2num(line(12:13));
st_hr = str2num(line(15:16));
st_min = str2num(line(18:19));
st_sec = str2num(line(21:31));

% Total number of epochs recorded
nEpoch = str2num(line(33:39));

% Finding our epoch of interest
if ~(st_year == date(1,1) & st_month == date(1,2));
   warning('Month or year of observations and ephemerides do not match!');
end

day_diff = date(1,3) - st_day;
hr_diff = date(1,4) - st_hr;
min_diff = date(1,5) - st_min;
sec_diff = date(1,6) - st_sec;

if day_diff & hr_diff & min_diff < 0
    warning('Epoch not available in ephemeride.');
end

if sec_diff <= 0
    epoch = 1; % if the epoch is seconds before the ephemeride starts 
               %just use first recorded satellite positions.
end

% Rounding to closest lower epoch
epoch = floor(day_diff*96 + hr_diff*4 + (min_diff + (sec_diff/60))/15) + 1;

if epoch > nEpoch
    warning('Epoch of interest takes place after last recorded epoch.');
end

while line(1) ~= '+'
    line = fgetl(fid);
end
nSat = str2num(line(5:6));
lines_per_epoch = nSat + 1;

epoch_counter = 0;

while 1
    
    while line(1) ~= '*'
        line = fgetl(fid);
    end
    epoch_counter = epoch_counter + 1;
    if epoch_counter == epoch
        break
    end
    line = fgetl(fid);
end

[t] = sscanf(line, '%c %d %d %d %d %d %lf'); % time of epoch
time_diff = 60*(date(1, 5) - t(6, 1)) + (date(1, 6) - t(7, 1));
interval = 60*15;

line = fgetl(fid);

% Searching through the epoch for satellites of interest
while line(1) ~= '*'
    if (line(1)=='P')
        [p] = sscanf(line, '%c %c %d %lf %lf %lf %f');
        sat_ind = find(p(3)==sat);
        if ~isempty(sat_ind)
            XS_l(sat_ind,:) = 1000*p(4:6);
            dtS_l(sat_ind,1) = (10^-6)*p(7);
        end
    end
    line = fgetl(fid);
end

% going to next epoch
while line(1) ~= '*'
    line = fgetl(fid);
end
epoch_counter = epoch_counter + 1;
line = fgetl(fid);
% Searching through the epoch for satellites of interest
while line(1) ~= '*'
    if (line(1)=='P')
        [p_h] = sscanf(line, '%c %c %d %lf %lf %lf %f');
        sat_ind = find(p_h(3)==sat);
        if ~isempty(sat_ind)
            XS_h(sat_ind,:) = 1000*p_h(4:6);
            dtS_h(sat_ind,1) = (10^-6)*p_h(7);
        end
    end
    line = fgetl(fid);
end

XS = (XS_h - XS_l)*time_diff/interval + XS_l;
dtS = (dtS_h - dtS_l)*time_diff/interval + dtS_l;

fclose(fid);
end
    




    