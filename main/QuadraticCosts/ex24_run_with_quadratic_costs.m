function opt_solution = ex24_run_with_quadratic_costs()
clc
close all

addpath( ...
    '/Users/Juraj/Documents/Optimization/matpower/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower/data', ...
    '/Users/Juraj/Documents/Optimization/matpower/mips/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower/mips/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower/most/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower/most/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower/mptest/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower/mptest/lib/t', ...
    '-end' );

setenv('OMP_NUM_THREADS', '1')

define_constants;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set options

opt = mpoption('verbose',2,'out.all',0, 'opf.ac.solver','IPOPT');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set time step data
factor_timesteps = 1;  %% 1...365 (number of time periods)
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
mpc = case118;
%mpc = case30; % works with top 2%, not with static = 10 or top 1%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fixes:
mpc = ext2int(mpc);
mpc.branch(mpc.branch(:,RATE_A)==0,RATE_A) = 9900;
mpc.gen(:,PMIN) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% prepare storage data

%do we place static number of storages to all grids, or we pick top 2% of
%the load buses
static_placement = 0;

if (static_placement)
    id_load_log = abs(mpc.bus(:,PD))>0;
    id_load = find(id_load_log);
    nload = length(id_load);
    nstorage_ref = 3;  %% 1 ... 100
    nstorage_applied = min(nstorage_ref,nload);
    id_storage_location = id_load(1:nstorage_applied);
else
    [sortedPD, id_load] = sort(mpc.bus(:,PD), 'descend'); %% sort buses w.r.t PD
    nload = length(find(abs(mpc.bus(:,PD))>0));
    nstorage_ref = round(size(mpc.bus,1)*0.02); %% apply storage to top 2% of buses with highest load
    nstorage_applied = min(nstorage_ref,nload);
    id_storage_location = id_load(1:nstorage_applied); %% apply storages to top 2% loaded buses
end

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

%% run OPF
mpcN_opf_storage = create_storage_case_file3(mpc,load_scaling_profile, p_storage);
opt_solution     = runopf(mpcN_opf_storage, opt);

%% verify constraints
N = length(load_scaling_profile) * factor_timesteps;

%linear constraints
E = opt_solution.A * opt_solution.x;
Emin = -p_storage.E_storage_init_MWh;
Emax = p_storage.E_storage_max_MWh - p_storage.E_storage_init_MWh;

%PF constraints and branch PF
[h, g] = opf_consfcn(opt_solution.x, opt_solution.om);
gn = 2*size(mpc.bus,1);
hn = 2*size(mpc.branch,1);


for i = 1:N
    E_i = E((1:nstorage_applied) + (i-1)*nstorage_applied);
    idx = find(E_i > Emax );
    nviol = length(idx);
    if (nviol > 0)
       fprintf('Vioated MAX energy limit for %d storages in period %d with error %e\n', nviol, i, max(E_i(idx) - Emax(idx))); 
    end
    
    idx = find(E_i < Emin );
    nviol = length(idx);
    if (nviol > 0)
       fprintf('Vioated MIN energy limit for %d storages in period %d with error %e\n', nviol, i, max(Emin(idx) - E_i(idx))); 
    end
    
    g_local = g((1:gn) + (i-1)*gn);
    idx = find(g_local > 1e-5 );
    nviol = length(idx);
    if (nviol > 0)
       fprintf('Vioated %d PF equations in period %d with error %e\n', nviol, i, max(g_local(idx))); 
    end
    
    h_local = h((1:hn) + (i-1)*hn);
    idx = find(h_local > 0 );
    nviol = length(idx);
    if (nviol > 0)
       fprintf('Vioated %d branch flow limits in period %d with error %e\n', nviol, i, max(h_local(idx))); 
    end
end

end
