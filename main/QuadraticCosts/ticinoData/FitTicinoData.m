function profile = FitTicinoData(filename, t0, t1, tnew)
% filename - name of the load data
% t0,t1 - time of the load data
% tnew - new time of the rescaled data
n=7;
data        = 1e-5*load(filename);
t           = t0:t1; t = t';
[t, data]   = ScaleData(t, data);
params.t    = t;
params.data = data;

options     = optimoptions('fmincon','Algorithm','interior-point','Display','iter','TolFun', 1e-13, 'TolX', 1e-14, ...
	                       'ScaleProblem','obj-and-constr','MaxFunEvals',10^5);

options     = optimoptions('fmincon','Algorithm','sqp','Display','iter','TolFun', 1e-14, 'TolX', 1e-15, ...
	                       'ScaleProblem','obj-and-constr','MaxFunEvals',10^5);
                       
%ScaleProblem = jacobian                       

LB = [repmat(.0, 2*n, 1);  .7; 1.4;  2.2;  2.8; 4.0;  4.8; 5.6 ];
UB = [repmat(10, 2*n, 1);  1.2; 1.6;  2.5;  3.0; 4.5;  5.1; 6.2 ];
x0 = 0.5*(LB+UB);

%x = fminunc(@(x) misfit(x,params),x0,options);
x = fmincon(@(x) misfit(x,params),x0,[],[],[],[],LB, UB,[],options);
%fprintf('%25.16e\n', x);

%% sample approximated data on hourly basis
hourly.t = tnew';
profile = proxy(x, hourly);

%%
% f = proxy(x, params);
% close all;
% plot(t, data, 'k--'); hold on;
% plot(t, f); hold on;
% plot(hourly.t, profile, 'x');
% legend('Original data', 'Approximation', 'Profile');
% 
% figure;
% plot(t, data, 'k--');
% hold on;
% for i=1:n
% wavelets(:,i) = Wavelet(i, x, params);
% end
% plot(t, wavelets);
% grid on;
% grid minor;
% legend('Original data', 'n=1', 'n=2', 'n=3', 'n=4', 'n=5', 'n=6', 'n=7')

%f=@(x)(3+sin(2*x-pi/2)+2*sin(4*x-pi/2)*x/10)/5 ; fplot(f, [0 2*pi])
end


function y = misfit(x, params)
   t    = params.t;
   data = params.data;

   f = proxy(x, params);
   %y = sqrt((f-data).' * (f - data));
   y = (f-data).' * (f - data);
end


function f = proxy(x, params)
   t = params.t;

  % f = (x(1)+sin(x(2)*t-x(3))  +  x(4)*sin(x(5)*t-x(6)).*t)/5;
   n = length(x)/3;
   f = zeros(length(t),1);
   for i=1:n
   	 %f = f + x(i)*exp(-x(n+i)*(t-x(2*n+i)).^2);
     f = f + Wavelet(i, x, params);
   end

 end


function f = Wavelet(i, x, params)
   t = params.t;
  % f = (x(1)+sin(x(2)*t-x(3))  +  x(4)*sin(x(5)*t-x(6)).*t)/5;
   n = length(x)/3;
   
   f = x(i)*exp(-x(n+i)*(t-x(2*n+i)).^2);
 end


function [t, y] = ScaleData(t0, data)

	data_min = min(data);
	data_max = max(data);

	y        = (data - data_min) / (data_max - data_min);
    
    t0_min   = min(t0);
    t0_max   = max(t0);

    t        = 2*pi * (t0 - t0_min)/(t0_max - t0_min); 

end

% f=@(x)(3+sin(2*x-pi/2)+2*sin(4*x-pi/2)*x/10)/5 ; fplot(f, [0 2*pi])