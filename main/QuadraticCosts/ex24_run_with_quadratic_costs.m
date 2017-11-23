function opt_solution = ex24_run_with_quadratic_costs()
clc
close all
define_constants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set options

opt = mpoption('verbose',1,'out.all',0, 'opf.ac.solver','IPOPT');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set time step data
%nstorage_ref = 10;  %% 1 ... 100
factor_timesteps = 2;  %% 1...365 (number of time periods)
Kpv = 1;  %% 0 ...1  (size of PV penetration)


load_scaling_profile0 = [0.4544 0.3570 0.2860 0.2783 0.3795 0.5822 0.8086 0.9633 1.0086 0.9883 0.9761 1.0000 1.0193 0.9773 0.8772 0.7991 0.8359 1.0023 1.2063 1.3123 1.2438 1.0343 0.7873 0.5885]';
pv_scaling_profile =  [  0 0 0  0 0 0.0046 0.0548 0.1686 0.3457 0.5100 0.6687 0.7496 0.8175 0.8305 0.8026 0.7212 0.5988 0.4453 0.2718 0.1203 0.0350 0.0019 0 0 ]';

%load_scaling_profile0=load_scaling_profile0(1:5);
%pv_scaling_profile=pv_scaling_profile(1:5);

save -ascii -double 'load.dat' load_scaling_profile0;
save -ascii -double 'pv.dat' pv_scaling_profile;

load_scaling_profile = load_scaling_profile0- 1*pv_scaling_profile;
save -ascii -double 'loadpv.dat' load_scaling_profile;

figure;
plot(load_scaling_profile0,'b--');
hold on;
plot(load_scaling_profile,'b');
legend('load','net load with PV injection')


load_scaling_profile       = kron(ones(factor_timesteps,1), load_scaling_profile  );
%load_scaling_profile       = load_scaling_profile(1:5); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create base case file
mpc= case118;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fixes:
mpc = ext2int(mpc);
mpc.branch(mpc.branch(:,RATE_A)==0,RATE_A) = 9900;
mpc.gen(:,PMIN) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% prepare storage data
% id_load_log = abs(mpc.bus(:,PD))>0;
% id_load = find(id_load_log);
% nload = length(id_load);
% nstorage_applied = min(nstorage_ref,nload);
% id_storage_location = id_load(1:nstorage_applied);

[sortedPD, id_load] = sort(mpc.bus(:,PD), 'descend'); %% sort buses w.r.t PD
nload = length(find(abs(mpc.bus(:,PD))>0));
nstorage_ref = round(size(mpc.bus,1)*0.02); %% apply storage to top 2% of buses with highest load
nstorage_applied = min(nstorage_ref,nload);
id_storage_location = id_load(1:nstorage_applied); %% apply storages to top 2% loaded buses

% remove storage if applied to REF bus, add it to different bus
if find(mpc.bus(id_storage_location, BUS_TYPE) == 3)
    nstorage_applied = min(nstorage_ref+1,nload);
    id_storage_location = id_load(1:nstorage_applied);
    id_storage_location(mpc.bus(id_storage_location,BUS_TYPE) ==3) = [];
end

p_storage.id_storage_location = id_storage_location;
p_storage.E_storage_max_MWh  = 2*abs(mpc.bus(id_storage_location,PD));
p_storage.E_storage_init_MWh = p_storage.E_storage_max_MWh*.7;
p_storage.rPmaxEmax_MW_per_MWh = 1/3;
p_storage.rPminEmax_MW_per_MWh = -1/2;
p_storage.c_discharge        = .97;
p_storage.c_charge           = .95;

mpcN_opf_storage = create_storage_case_file3(mpc,load_scaling_profile, p_storage)
opt_solution     = runopf(mpcN_opf_storage, opt)


end
