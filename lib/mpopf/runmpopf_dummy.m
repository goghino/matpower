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
Nprofile = size(load_scaling_profile,1);
if (N <= Nprofile)
    load_scaling_profile = load_scaling_profile(1:N);
else
    load_scaling_profile = repmat(load_scaling_profile, ceil(N/Nprofile), 1);
    load_scaling_profile = load_scaling_profile(1:N);
end

storage_load = abs(p_storage.E_storage_max_MWh * p_storage.rPminEmax_MW_per_MWh);
storage_injection = p_storage.E_storage_max_MWh * p_storage.rPmaxEmax_MW_per_MWh;
%load_scaling_profile = scaleLoadProfile(load_scaling_profile, mpc, storage_injection, storage_load);
load_scaling_profile = scaleLoadProfile(load_scaling_profile, mpc, 0, 0);

%% run OPF
mpcN_opf_storage = create_storage_case_file3(mpc,load_scaling_profile, p_storage);
[RESULTS, SUCCESS] = runopf(mpcN_opf_storage, mpopt);

% plot_storage_results(opt_solution)

end