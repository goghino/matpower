function profile_scaled = scaleLoadProfile(profile, mpc, storage_injection, storage_load)
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