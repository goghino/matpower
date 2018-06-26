function profile = createLoadProfile()
%CREATELOADPROFILE Creates model of the load scaling for the planning horizon
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

% profile48 = [    0.0335    0.0423    0.0527    0.0650    0.0793    0.0956    0.1147    0.1402    0.1853    0.2796 ...
%     0.4484    0.6575    0.8078    0.8430    0.8218    0.8253    0.8558    0.8636    0.8181 ...
%     0.7286    0.6252    0.5348    0.4693    0.4269    0.3998    0.3821    0.3740    0.3838    0.4267    0.5171 ...
%     0.6543    0.8098    0.9336    0.9882    0.9868    0.9733    0.9402    0.8142    0.5773 ...
%     0.3278    0.1710    0.1543    0.1396    0.1300    0.10    0.0864    0.0427    0.0400]';   
%
% hours = linspace(0, 2*pi, length(profile));

%%
profile = csvread('ticinoData/TI240hrs.dat');
hours = linspace(0, 2*pi, length(profile));
% normalize data to 0-2pi and 0-1
[~, profile] = ScaleData(hours, profile);

profile = profile(1:10);

end