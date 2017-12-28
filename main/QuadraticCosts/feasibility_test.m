addpath( ...
    '/home/juraj/matpower-3.12.8', ...
    '/home/juraj/matpower-3.12.8/lib', ...
    '/home/juraj/matpower-3.12.8/lib/t', ...
    '/home/juraj/matpower-3.12.8/data', ...
    '/home/juraj/matpower-3.12.8/mips/lib', ...
    '/home/juraj/matpower-3.12.8/mips/lib/t', ...
    '/home/juraj/matpower-3.12.8/most/lib', ...
    '/home/juraj/matpower-3.12.8/most/lib/t', ...
    '/home/juraj/matpower-3.12.8/mptest/lib', ...
    '/home/juraj/matpower-3.12.8/mptest/lib/t', ...
    '-end' );

setenv('OMP_NUM_THREADS', '1')

mpopt0 = mpoption('verbose', 2, 'out.all', 0, 'opf.start', 3);
mpopt = mpoption(mpopt0, 'opf.ac.solver', 'IPOPT');

%% create scaling profile
mpc0 = case1888rte;
profile = createLoadProfile(mpc0);

%% test each step of load profile and test if it is feasible
for i = 1:length(profile)
   
    alpha = profile(i);
    
    %scale Pd, Qd by load profile
    mpc = mpc0;
    mpc.bus(:,3:4) = mpc.bus(:,3:4) .* alpha;
    
    [RESULTS, SUCCESS] = runopf(mpc, mpopt);
    
    if(~SUCCESS)
       fprintf('Not feasible for alpha %d, hrs = %d.\n', alpha, i); 
    end
end




