clc
close all

addpath( ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/data', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/mips/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/mips/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/most/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/most/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/mptest/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/mptest/lib/t', ...
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
%plot_storage_results(opt_solution)


%% verify constraints graphically
nstorage = opt_solution.nstorage;
N = opt_solution.horizon;
Pgen_discharge = reshape(opt_solution.gen(opt_solution.id_gen_storages_discharge,PG),[nstorage, N ]);
Pgen_charge    = reshape(opt_solution.gen(opt_solution.id_gen_storages_charge   ,PG),[nstorage, N ]);
Pgen_storage = Pgen_discharge+Pgen_charge;
Pgen_storage_max = repmat(opt_solution.P_storage_max_MW,[1,N]);
Pgen_storage_min = repmat(opt_solution.P_storage_min_MW,[1,N]);

E_storage = [opt_solution.E_storage_init_MWh, repmat(opt_solution.E_storage_init_MWh,[1,N]) - cumsum(Pgen_discharge./repmat(opt_solution.c_discharge,[1,N]) + Pgen_charge.*repmat(opt_solution.c_charge,[1,N]),2)];
E_storage_max = repmat(opt_solution.E_storage_max_MWh,[1,N+1]);

stairs(0:N, [Pgen_storage Pgen_storage(:,end)  ]', 'LineWidth',2);
hold on;
stairs(0:N, [Pgen_storage_max Pgen_storage_max(:,end)  ]', '--');
stairs(0:N, [Pgen_storage_min Pgen_storage_min(:,end)  ]', '--');
grid on
title('storage network injection power : trajectory (solid) and maximum (dashed)')
xlabel('t [hours]')
ylabel('storage power [MW]')

figure();
plot(0:N, E_storage', 'LineWidth',2);
hold on;
plot(0:N, E_storage_max','--');
grid on
title('storage energy level: trajectory (solid) and maximum (dashed)')
xlabel('t [hours]')
ylabel('storage energy [MWh]')

%% verify constraints

%linear constraints
% 0 < Einit + T*P1              < Emax 
% 0 < Einit + T*P1 + T*P2       < Emax 
% 0 < Einit + T*P1 + T*P2 + ... < Emax 
% E = opt_solution.A * opt_solution.x;
% Emin = -p_storage.E_storage_init_MWh;
% Emax = p_storage.E_storage_max_MWh - p_storage.E_storage_init_MWh;

E = E_storage;
Emax = opt_solution.E_storage_max_MWh;
Emin = 0*opt_solution.E_storage_max_MWh;

%PF constraints and branch PF
[h, g] = opf_consfcn(opt_solution.x, opt_solution.om);
gn = 2*size(mpc.bus,1);
hn = 2*size(mpc.branch,1);


for i = 1:N
    E_i = E((1:nstorage_applied), i);
    idx = find(E_i > Emax );
    nviol = length(idx);
    if (nviol > 0)
       fprintf('[N=%2d] Vioated MAX energy limit for %d storages with error %e\n', i, nviol, max(E_i(idx) - Emax(idx))); 
    end
    
    idx = find(E_i < Emin );
    nviol = length(idx);
    if (nviol > 0)
       fprintf('[N=%2d] Vioated MIN energy limit for %d storages with error %e\n', i, nviol, max(Emin(idx) - E_i(idx))); 
    end
    
    g_local = g((1:gn) + (i-1)*gn);
    idx = find(g_local > 1e-5 );
    nviol = length(idx);
    if (nviol > 0)
       fprintf('[N=%2d] Vioated %d PF equations with error %e\n', i, nviol, max(g_local(idx))); 
    end
    
    h_local = h((1:hn) + (i-1)*hn);
    idx = find(h_local > 0 );
    nviol = length(idx);
    if (nviol > 0)
       fprintf('[N=%2d] Vioated %d branch flow limits with error %e\n', i, nviol, max(h_local(idx))); 
    end
end

%figure; plot(E(1:nstorage_applied:end)); title('Energy level of 1st storage over time'); xlabel('T [h]'); ylabel('E [MWh]')
%figure; plot(E(2:nstorage_applied:end)); title('Energy level of 2st storage over time'); xlabel('T [h]'); ylabel('E [MWh]')

%% parse x vector
% nb = size(mpc.bus,1);
% ng = size(mpc.gen,1);
% 
% offset = nb*N + nb*N; %skip Va Vm
% 
% Pg_discharge = opt_solution.x(offset + ng*N + (1:nstorage_applied:nstorage_applied*N)) * opt_solution.baseMVA;
% offset = offset + nstorage_applied*N;
% Pg_charge = opt_solution.x(offset + (1:nstorage_applied:nstorage_applied*N)) * opt_solution.baseMVA;
% offset = offset + nstorage_applied*N + ng*N;
% Qg_discharge = opt_solution.x(offset + (1:nstorage_applied*N)) * opt_solution.baseMVA;
% offset = offset + nstorage_applied*N;
% Qg_charge = opt_solution.x(offset + (1:nstorage_applied*N)) * opt_solution.baseMVA;
% 
% figure; plot(Pg_discharge(1:nstorage_applied:end)); title('Pg discharge over time'); xlabel('T [h]'); ylabel('P [MW]')
% figure; plot(Pg_charge(1:nstorage_applied:end));title('Pg charge over time'); xlabel('T [h]'); ylabel('P [MW]')
