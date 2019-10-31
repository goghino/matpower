function [profile, ratios] = createLoadProfile(N, mpc)
%CREATELOADPROFILE Creates model of the load scaling for the planning
%horizon. Load at each bus is scaled by profile(i) for each time period.
%If ratios is set, it splits the load at each bus to
%residendtial/commertial and each can be than scaled using different
%profile
close all;
%% 
% f=@(x)(3+sin(2*x-pi/2)+2*sin(4*x-pi/2))/5 ;
% f=@(x)(3+sin(2.*x-pi/2)+2*sin(4.*x-pi/2).*x/10)/5;
% fplot(f, [0 2*pi]); hold on;
% hours = linspace(0, 2*pi, 24);
% profile = [f(hours)]'; %needs to be a column vector

%% 
%hours = linspace(0, 2*pi, 48);
%profile = FitTicinoData('Ticino303_390.dat',303,390,hours);

% profile24 = [0.0335    0.0532    0.0806    0.1176    0.1972    0.4933    0.8275    0.8183    0.8635 ...
%            0.7868    0.5832    0.4466    0.3895    0.3765    0.4756    0.7566    0.9793    0.9768 ...
%            0.8524    0.3665    0.1255    0.1549    0.0886    0.0184 ]';
% profile = profile24;

% profile48 = [    0.0335    0.0423    0.0527    0.0650    0.0793    0.0956    0.1147    0.1402    0.1853    0.2796 ...
%     0.4484    0.6575    0.8078    0.8430    0.8218    0.8253    0.8558    0.8636    0.8181 ...
%     0.7286    0.6252    0.5348    0.4693    0.4269    0.3998    0.3821    0.3740    0.3838    0.4267    0.5171 ...
%     0.6543    0.8098    0.9336    0.9882    0.9868    0.9733    0.9402    0.8142    0.5773 ...
%     0.3278    0.1710    0.1543    0.1396    0.1300    0.10    0.0864    0.0427    0.0400]';   
%
% hours = linspace(0, 2*pi, length(profile));

%%
addpath('/Users/Juraj/Documents/Optimization/matpower/lib/mpopf/ticinoData');
profile = csvread('TI240hrs.dat');
hours = linspace(0, 2*pi, length(profile));
% normalize data to 0-2pi and 0-1
[~, profile] = ScaleData(hours, profile);
ratios = ones(size(mpc.bus, 1),1);

%% DSO profile for case15cigre
if size(mpc.bus,1) == 15
    % Residential profile over time
    rp  = [0.25 0.20 0.19 0.195 0.19 0.25 0.4 0.60 0.65 0.65 0.65 0.6 0.75 0.65 0.55 0.5 0.45 0.60 0.75 0.90 0.8 0.7 0.6 0.40 0.3]';
    % Commercial and Industrial profile over time
    cip = [0.35 0.35 0.30 0.375 0.40 0.50 0.6 0.85 1.00 0.95 1.00 0.8 0.85 0.90 0.90 0.9 0.80 0.55 0.50 0.45 0.4 0.4 0.3 0.35 0.3]';

    % 1st row is percentage of resedential loads and 2nd row is percentage of
    % commercial/industrial loads at each bus
    rcip = [0 0.75	0.00	0.51	1.00	1.00	1.00	0.00	1.00	0.00	0.86	1.00	0.74	0.00	0.35;
            0 0.25	0.00	0.49	0.00	0.00	0.00	1.00	0.00	1.00	0.14	0.00	0.26	1.00	0.65]';

    profile = [rp, cip];
    ratios = rcip;
end
%% Create load of the size N that was requested (repeat if data above are smaller)
Nprofile = size(profile,1);
if (N <= Nprofile)
    profile = profile(1:N,:);
else
    profile = repmat(profile, ceil(N/Nprofile), 1);
    profile = profile(1:N,:);
end

%storage_load = abs(p_storage.E_storage_max_MWh * p_storage.rPminEmax_MW_per_MWh);
%storage_injection = p_storage.E_storage_max_MWh * p_storage.rPmaxEmax_MW_per_MWh;
%load_scaling_profile = scaleLoadProfile(load_scaling_profile, mpc, storage_injection, storage_load);
profile = scaleLoadProfile(profile, mpc, 0, 0);
end
