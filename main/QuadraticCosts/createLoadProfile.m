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
%profile = FitTicinoData('Ticino303_390.dat',303,390,hours);
profile = [0.0335    0.0532    0.0806    0.1176    0.1972    0.4933    0.8275    0.8183    0.8635 ...
           0.7868    0.5832    0.4466    0.3895    0.3765    0.4756    0.7566    0.9793    0.9768 ...
           0.8524    0.3665    0.1255    0.1549    0.0886    0.0184 ]';

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
PGmax_sum = sum(mpc.gen(generatorsON,PMAX)); % ??? take only 80% of max generation (account for e.g. transmillion losses)
QGmin_sum = sum(mpc.gen(generatorsON,QMIN));
QGmax_sum = sum(mpc.gen(generatorsON,QMAX));

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

%% 
plot(hours, profile); hold on;

%scale the profile prof so that it is in bounds prof_min, prof_max 
%which makes sure that load does not exceed generation
profile = (prof_max - prof_min) * profile + prof_min;

plot(hours, profile);
legend('Approximated load', 'Approximated and scaled load')
end