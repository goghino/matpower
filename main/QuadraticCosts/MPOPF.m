function opt_solution = MPOPF(mpc, mpopt, N, Emax, Rcount, Rfirst)
%MPOPF Run a multiperiod OPF problem
%   mpc  - Matpower case
%   N    - Number of time periods
%   Emax - max capacity of the storage relative to a PD at given bus
%   Rcount - fraction of the PD buses that will have storage
%   Rfirst - relative position of the storage placement when sorted by PD from largest first

%% set options
define_constants;

%% set time step data
factor_timesteps = N;  %% 1...365 (number of time periods)
Kpv = 1;  %% 0 ...1  (size of PV penetration)

load_scaling_profile0 = [0.4544 0.3570 0.2860 0.2783 0.3795 0.5822 0.8086 0.9633 1.0086 0.9883 0.9761 1.0000 1.0193 0.9773 0.8772 0.7991 0.8359 1.0023 1.2063 1.3123 1.2438 1.0343 0.7873 0.5885]';
pv_scaling_profile =  [  0 0 0  0 0 0.0046 0.0548 0.1686 0.3457 0.5100 0.6687 0.7496 0.8175 0.8305 0.8026 0.7212 0.5988 0.4453 0.2718 0.1203 0.0350 0.0019 0 0 ]';

load_scaling_profile = load_scaling_profile0- 1*pv_scaling_profile;

% curtail load scaling profile
MIN = 0.0779;
MAX = 1.0324;

figure;
plot(load_scaling_profile0,'b--');
hold on;
plot(load_scaling_profile,'b');
hold on; plot(MIN*ones(length(load_scaling_profile),1),'r--');
hold on; plot(MAX*ones(length(load_scaling_profile),1),'r--');
legend('load','net load with PV injection', 'Load curtailment')


load_scaling_profile = max(load_scaling_profile, MIN);
load_scaling_profile = min(load_scaling_profile, MAX);

load_scaling_profile       = kron(ones(factor_timesteps,1), load_scaling_profile  );

%% mpc fixes
mpc = ext2int(mpc);
mpc.branch(mpc.branch(:,RATE_A)==0,RATE_A) = 9900;
mpc.gen(:,PMIN) = 0;

%% prepare storage data

%do we place static number of storages to all grids, or we pick top 2% of the load buses?
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
    nstorage_ref = round(size(mpc.bus,1)* Rcount); %% apply storage N% of buses
    first = round(size(mpc.bus,1)* Rfirst); %% place storages starting at 'first' bus
    nstorage_applied = min(nstorage_ref,nload);
    id_storage_location = id_load(first + (1:nstorage_applied)); %% apply storages to top 2% loaded buses
end

if(nstorage_applied < 1)
   error('Number of storateges has to be > 0'); 
end

% remove storage if applied to REF bus, add it to different bus
if find(mpc.bus(id_storage_location, BUS_TYPE) == 3)
    nstorage_applied = min(nstorage_ref+1,nload);
    id_storage_location = id_load(1:nstorage_applied);
    id_storage_location(mpc.bus(id_storage_location,BUS_TYPE) ==3) = [];
end

p_storage.id_storage_location = id_storage_location;
p_storage.E_storage_max_MWh  = Emax*abs(mpc.bus(id_storage_location,PD));
p_storage.E_storage_init_MWh = p_storage.E_storage_max_MWh*.7;
p_storage.rPmaxEmax_MW_per_MWh = 1/3;
p_storage.rPminEmax_MW_per_MWh = -1/2;
p_storage.c_discharge        = .97;
p_storage.c_charge           = .95;

%% run OPF
mpcN_opf_storage = create_storage_case_file3(mpc,load_scaling_profile, p_storage);
opt_solution     = runopf(mpcN_opf_storage, mpopt);
%plot_storage_results(opt_solution)

end