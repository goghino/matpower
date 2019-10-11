function [RESULTS, SUCCESS] = runmpopf_dummy(mpc, mpopt, N)
%MPOPF Run a multiperiod OPF problem
%   mpc   - Matpower case
%   mpopt - Matpower options struct
%   N     - Number of time periods

%%
%   Emax - max capacity of the storage relative to a PD at given bus
Emax = 2;

%   Rcount - fraction of the PD buses that will have storage, negative is absolute count
Rcount = -2;

%   Rfirst - relative position of the storage placement when sorted by PD from largest first
Rfirst = 0.00;

%% set options
define_constants;

%% mpc fixes
mpc = ext2int(mpc);
mpc.branch(mpc.branch(:,RATE_A)==0,RATE_A) = 9900;
mpc.gen(:,PMIN) = 0;

%% prepare storage data
p_storage = createStorage(mpc, Rfirst, Rcount, Emax);

%% set load scaling profile and scale it to respect PG/QG generation limits
[load_scaling_profile, residential_ratios] = createLoadProfile(N, mpc);

mpc.load_profile = load_scaling_profile;
mpc.load_ratios = residential_ratios;

%% run OPF
mpcN_opf_storage = create_storage_case_file3(mpc, p_storage);
[RESULTS, SUCCESS] = runopf(mpcN_opf_storage, mpopt);

plot_storage_results(RESULTS)

end