function [RESULTS, SUCCESS] = runmpopf_dummy(mpc, mpopt, N, Emax, Rcount, Rfirst)
%MPOPF Run a multiperiod OPF problem
%   mpc  - Matpower case
%   N    - Number of time periods
%   Emax - max capacity of the storage relative to a PD at given bus
%   Rcount - fraction of the PD buses that will have storage
%   Rfirst - relative position of the storage placement when sorted by PD from largest first

%% set options
define_constants;

%% mpc fixes
mpc = ext2int(mpc);
mpc.branch(mpc.branch(:,RATE_A)==0,RATE_A) = 9900;
mpc.gen(:,PMIN) = 0;

%% prepare storage data

p_storage = createStorage(mpc, Rfirst, Rcount, Emax);

%% set load scaling profile and scale it to respect PG/QG generation limits
load_scaling_profile = createLoadProfile();

%storage_load = abs(p_storage.E_storage_max_MWh * p_storage.rPminEmax_MW_per_MWh);
%storage_injection = p_storage.E_storage_max_MWh * p_storage.rPmaxEmax_MW_per_MWh;
%load_scaling_profile = scaleLoadProfile(load_scaling_profile, mpc, storage_injection, storage_load);

%repeat data for required no. of days
load_scaling_profile       = kron(ones(N,1), load_scaling_profile);

%% run OPF
mpcN_opf_storage = create_storage_case_file3(mpc,load_scaling_profile, p_storage);
[RESULTS, SUCCESS] = runopf(mpcN_opf_storage, mpopt);

% plot_storage_results(opt_solution)

end