function profile = createLoadProfile(mpc)
%CREATELOADPROFILE Creates model of the load scaling for the 24hrs horizon
close all;
%% 
% f=@(x)(3+sin(2*x-pi/2)+2*sin(4*x-pi/2))/5 ;
% f=@(x)(3+sin(2.*x-pi/2)+2*sin(4.*x-pi/2).*x/10)/5;
% fplot(f, [0 2*pi]); hold on;
% hours = linspace(0, 2*pi, 24);
% profile = [f(hours)]'; %needs to be a column vector

%% 
addpath('ticinoData');
hours = linspace(0, 2*pi, 24);
profile = FitTicinoData('Ticino303_390.dat',303,390,hours);

%% compute min and max generation capabilities
%  so that we can scale the load profile

%generator limits idx
PMAX = 9;
PMIN = 10;
QMIN = 5;
QMAX = 4;

%nominal load idx
PD = 3;
QD = 4;

% take only 80% of max generation (account for e.g. transmillion losses)
PGmin_sum = sum(mpc.gen(:,PMIN));
PGmax_sum = sum(mpc.gen(:,PMAX)) * 0.8;
QGmin_sum = sum(mpc.gen(:,QMIN));
QGmax_sum = sum(mpc.gen(:,QMAX));
PD_sum = sum(mpc.bus(:,PD));
QD_sum = sum(mpc.bus(:,QD));

% PGmin_sum < prof * PD_sum < PGmax_sum 
% QGmin_sum < prof * QD_sum < QGmax_sum
% which means that bounds on prof are:
% PGmin_sum / PD_sum < prof < PGmax_sum / PD_sum 
% QGmin_sum / QD_sum < prof < QGmax_sum / QD_sum

fprintf('P: %.2f <= alpha * %.2f <= %.2f\n', PGmin_sum, PD_sum, PGmax_sum);
fprintf('Q: %.2f <= alpha * %.2f <= %.2f\n', QGmin_sum, QD_sum, QGmax_sum);

% scaling of the nominal load is withing bounds 30-100% or depending on
% PGmin and PGmax
prof_min = max( PGmin_sum / PD_sum, 0.3); 
prof_max = min( PGmax_sum / PD_sum, 1.0);

fprintf('%.2f <= alpha <= %.2f\n', prof_min, prof_max);

%% curtail

%% 
plot(hours, profile); hold on;

%scale the profile prof so that load does not exceed generation
profile = (prof_max - prof_min) * profile + prof_min;

plot(hours, profile);
legend('Approximated', 'Approximated and scaled')
end