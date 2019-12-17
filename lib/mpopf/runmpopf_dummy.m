function [RESULTS, SUCCESS] = runmpopf_dummy(mpc, mpopt, N)
%MPOPF Run a multiperiod OPF problem
%   mpc   - Matpower case
%   mpopt - Matpower options struct
%   N     - Number of time periods

%% set options
define_constants;

%% mpc fixes
mpc = ext2int(mpc);
mpc.branch(mpc.branch(:,RATE_A)==0,RATE_A) = 9900;
mpc.gen(:,PMIN) = 0;

%% prepare storage data
%   Emax - max capacity of the storage relative to a PD at given bus
Emax = 2;
%   Rcount - fraction of the PD buses that will have storage, negative is absolute count
Rcount = -2;
%   Rfirst - relative position of the storage placement when sorted by PD from largest first
Rfirst = 0.00;
%   enable_flexibility - flag if the storage flexibility should be provided
enable_flexibility = 1;

%prepares storage data: location, capacity, efficiency, ...
p_storage = createStorage(mpc, Rfirst, Rcount, Emax);

%% set load scaling profile and scale it to respect PG/QG generation limits
[load_scaling_profile, residential_ratios] = createLoadProfile(N, mpc);

mpc.load_profile = load_scaling_profile;
mpc.load_ratios = residential_ratios;

%% 
if(enable_flexibility == 1)
    mpc.storageFlexibility = 1;
    
    %Specify storage flexibility parameters and requirements
    mpc.FlexibilityReq.up = [0 3 2 0 0 4 0 0 0 0 zeros(1, N-10)]'; 
    mpc.FlexibilityReq.down = [-1 0 0 -2 -1 0 -15 -15 -15 -15 zeros(1, N-10)]';
    %in case N is very small, N < size(mpc.FlexibilityReq.XXX)
    mpc.FlexibilityReq.up = mpc.FlexibilityReq.up(1:N);
    mpc.FlexibilityReq.down = mpc.FlexibilityReq.down(1:N);
    assert(length(mpc.FlexibilityReq.up) == N);
    assert(length(mpc.FlexibilityReq.down) == N);
    
    mpc.enableDemandShift = 1;

    %Specify demand shift and flexibility parameters and requirements
    mpc.demandShift.busesID = [2,13]; %which buses offer demand shift
    mpc.demandShift.responsePowerMW = [0.15, 0.15]; %demand shift offer MW
    mpc.demandShift.responseTimeH = [3, 3]; %demand shift duration
    mpc.demandShift.reboundPowerMW = [0.15, 0.15]; %demand shift offer MW
    mpc.demandShift.reboundTimeH = [3, 3]; %demand shift duration
    mpc.demandShift.recoveryTimeH = [2, 2]; %demand shift duration
    assert(length(mpc.demandShift.busesID) == length(mpc.demandShift.responsePowerMW));
    assert(length(mpc.demandShift.busesID) == length(mpc.demandShift.responseTimeH));
    assert(length(mpc.demandShift.busesID) == length(mpc.demandShift.reboundPowerMW));
    assert(length(mpc.demandShift.busesID) == length(mpc.demandShift.reboundTimeH));
    assert(length(mpc.demandShift.busesID) == length(mpc.demandShift.recoveryTimeH));
end

%% run OPF
mpcN_opf_storage = create_storage_case_file3(mpc, p_storage);
[RESULTS, SUCCESS] = runopf(mpcN_opf_storage, mpopt);

if mpcN_opf_storage.storageFlexibility == 0
    plot_storage_results(RESULTS)
else
    plot_flexibility_results(RESULTS)
end

end