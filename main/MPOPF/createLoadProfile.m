function profile_scaled = createLoadProfile(mpc, storage_injection, storage_load)
%CREATELOADPROFILE Creates model of the load scaling for the planning horizon
close all;
%% 
% f=@(x)(3+sin(2*x-pi/2)+2*sin(4*x-pi/2))/5 ;
% f=@(x)(3+sin(2.*x-pi/2)+2*sin(4.*x-pi/2).*x/10)/5;
% fplot(f, [0 2*pi]); hold on;
% hours = linspace(0, 2*pi, 24);
% profile = [f(hours)]'; %needs to be a column vector

%% 
addpath('/home/i1042002/matpower-fork/main/QuadraticCosts/ticinoData');
%hours = linspace(0, 2*pi, 48);
%profile = FitTicinoData('Ticino303_390.dat',303,390,hours);

% profile24 = [0.0335    0.0532    0.0806    0.1176    0.1972    0.4933    0.8275    0.8183    0.8635 ...
%            0.7868    0.5832    0.4466    0.3895    0.3765    0.4756    0.7566    0.9793    0.9768 ...
%            0.8524    0.3665    0.1255    0.1549    0.0886    0.0184 ]';

% profile = [    0.0335    0.0423    0.0527    0.0650    0.0793    0.0956    0.1147    0.1402    0.1853    0.2796 ...
%     0.4484    0.6575    0.8078    0.8430    0.8218    0.8253    0.8558    0.8636    0.8181 ...
%     0.7286    0.6252    0.5348    0.4693    0.4269    0.3998    0.3821    0.3740    0.3838    0.4267    0.5171 ...
%     0.6543    0.8098    0.9336    0.9882    0.9868    0.9733    0.9402    0.8142    0.5773 ...
%     0.3278    0.1710    0.1543    0.1396    0.1300    0.10    0.0864    0.0427    0.0400]';   
% hours = linspace(0, 2*pi, length(profile));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile = csvread('ticinoData/TI240hrs.dat');
hours = linspace(0, 2*pi, length(profile));
[hours, profile] = ScaleData(hours, profile);
%profile = profile(1:24);
%profile = repmat(profile, 10, 1);
%fprintf('Running MPOPF with N=%d periods\n', size(profile,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute min and max generation capabilities
%  so that we can scale the load profile

%generator limits idx
GEN_STATUS = 8;
PMAX = 9;
PMIN = 10;
QMIN = 5;
QMAX = 4;

%nominal load idx
PD = 3;
QD = 4;

%identify ON generators
generatorsON = find(mpc.gen(:,GEN_STATUS) > 0);

PGmin_sum = sum(mpc.gen(generatorsON,PMIN));
PGmax_sum = sum(mpc.gen(generatorsON,PMAX)) + sum(storage_injection);
QGmin_sum = sum(mpc.gen(generatorsON,QMIN));
QGmax_sum = sum(mpc.gen(generatorsON,QMAX));

PD_sum = sum(mpc.bus(:,PD)) + sum(storage_load);
QD_sum = sum(mpc.bus(:,QD));

% PGmin_sum < prof * PD_sum < PGmax_sum 
% QGmin_sum < prof * QD_sum < QGmax_sum
% which means that bounds on prof are:
% PGmin_sum / PD_sum < prof < PGmax_sum / PD_sum 
% QGmin_sum / QD_sum < prof < QGmax_sum / QD_sum

fprintf('P: %.2f <= alpha * %.2f <= %.2f\n', PGmin_sum, PD_sum, PGmax_sum);
fprintf('Q: %.2f <= alpha * %.2f <= %.2f\n', QGmin_sum, QD_sum, QGmax_sum);

% scaling of the nominal load is withing bounds 30-100% or depending on
% PGmin and PGmax, not considering transmission losses
prof_min = PGmin_sum / PD_sum;
prof_max = PGmax_sum / PD_sum;
fprintf('%.2f <= alpha <= %.2f\n\n', prof_min, prof_max);

ALPHA_MIN = 0.8; 
ALPHA_MAX = 1.0;
fprintf('Used fixed bounds (%.2f, %.2f) on alfa:\n', ALPHA_MIN, ALPHA_MAX);

prof_min = max( prof_min, ALPHA_MIN); 
prof_max = min( prof_max, ALPHA_MAX);
fprintf('%.2f <= alpha <= %.2f\n', prof_min, prof_max);

%% 
%scale the profile prof so that it is in bounds prof_min, prof_max 
%which makes sure that load does not exceed generation
profile_scaled = (prof_max - prof_min) * profile + prof_min;

%%

%plot(hours, profile); hold on;
%plot(hours, profile_scaled);
%legend('Approximated load', 'Approximated and scaled load')
end
