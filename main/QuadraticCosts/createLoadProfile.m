function profile = createLoadProfile(mpc)
%CREATELOADPROFILE Creates model of the load scaling for the 24hrs horizon

%f=@(x)(3+sin(2*x-pi/2)+2*sin(4*x-pi/2))/5 ;
f=@(x)(3+sin(2.*x-pi/2)+2*sin(4.*x-pi/2).*x/10)/5;

hours = linspace(0, 2*pi, 24);

profile = [f(hours)]'; %needs to be a column vector
%plot(hours, profile); hold on; fplot(f, [0 2*pi])

%% scale the profile so that load does not exceed generation
% PMAX = 9;
% PMIN = 10;
% PD = 3; %PD and QD is scaled
% QD = 4;
% 
% PGmin_sum = sum(mpc.gen(:,PMIN));
% PGmax_sum = sum(mpc.gen(:,PMAX));
% PD_sum = sum(mpc.bus(:,PD));
% QD_sum = sum(mpc.bus(:,QD));
end