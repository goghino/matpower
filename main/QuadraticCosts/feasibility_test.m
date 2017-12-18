%% create scaling profile
load_scaling_profile0 = [0.4544 0.3570 0.2860 0.2783 0.3795 0.5822 0.8086 0.9633 1.0086 0.9883 0.9761 1.0000 1.0193 0.9773 0.8772 0.7991 0.8359 1.0023 1.2063 1.3123 1.2438 1.0343 0.7873 0.5885]';
pv_scaling_profile =    [0 0 0  0 0 0.0046 0.0548 0.1686 0.3457 0.5100 0.6687 0.7496 0.8175 0.8305 0.8026 0.7212 0.5988 0.4453 0.2718 0.1203 0.0350 0.0019 0 0 ]';

load_scaling_profile = load_scaling_profile0- 1*pv_scaling_profile;
load_scaling_profile = max(load_scaling_profile, 0.0779);
load_scaling_profile = min(load_scaling_profile, 1.0324);


mpopt = mpoption('verbose', 2, 'out.all', 0);
mpopt = mpoption(mpopt, 'opf.ac.solver', 'IPOPT', 'verbose', 2);

%% test each step of load profile and test if it is feasible
for i = 1:length(load_scaling_profile)
   
    alpha = load_scaling_profile(i);
    
    mpc = case300;
    
    %scale Pd, Qd by load profile
    mpc.bus(:,3:4) = mpc.bus(:,3:4) .* alpha;
    
    [RESULTS, SUCCESS] = runopf(mpc, mpopt);
    
    if(~SUCCESS)
       error('Not feasible for alpha %d, hrs = %d.', alpha, i); 
    end
end




