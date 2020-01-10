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
%   Einit - initial energy of the storage device at t0 (of the max capacity)
Emax = 2;
Einit = 0.7;
%   Rcount - fraction of the PD buses that will have storage, negative is absolute count
Rcount = -2;
%   Rfirst - relative position of the storage placement when sorted by PD from largest first
Rfirst = 0.00;
%   enable_flexibility - flag if the storage flexibility should be provided
enable_flexibility = 1;

%prepares storage data: location, capacity, efficiency, ...
p_storage = createStorage(mpc, Rfirst, Rcount, Emax, Einit);

%% set load scaling profile and scale it to respect PG/QG generation limits
[load_scaling_profile, residential_ratios] = createLoadProfile(N, mpc);

mpc.load_profile = load_scaling_profile;
mpc.load_ratios = residential_ratios;

%% 
if(enable_flexibility == 1)
    mpc.storageFlexibility = 1;
    
    %Specify demand shift and flexibility parameters and requirements
    mpc.enableDemandShift = 0;
    mpc.demandShift.busesID = [2,13]'; %which buses offer demand shift
    mpc.demandShift.responsePowerMW = 100*[0.1, 0.1]'; %demand shift offer MW
    mpc.demandShift.responseTimeH = [2, 2]'; %demand shift duration
    mpc.demandShift.reboundPowerMW = 100*[0.1, 0.1]'; %demand shift offer MW
    mpc.demandShift.reboundTimeH = [1, 1]'; %demand shift duration
    mpc.demandShift.recoveryTimeH = [2, 2]'; %demand shift duration
    assert(length(mpc.demandShift.busesID) == length(mpc.demandShift.responsePowerMW));
    assert(length(mpc.demandShift.busesID) == length(mpc.demandShift.responseTimeH));
    assert(length(mpc.demandShift.busesID) == length(mpc.demandShift.reboundPowerMW));
    assert(length(mpc.demandShift.busesID) == length(mpc.demandShift.reboundTimeH));
    assert(length(mpc.demandShift.busesID) == length(mpc.demandShift.recoveryTimeH));
    %binary variable is set in add_storage2bNoPeriodic: mpc_storage.z0=...
    
    
    %Specify requirements
    %mpc.FlexibilityReq.up =   100*[0.1 0.1  0       0   0  0.15 0.05 0    0 0 zeros(1, N-10)]'; 
    %mpc.FlexibilityReq.down = 100*[0   0    -0.01   0   0  0    0    -0.1 0 0 zeros(1, N-10)]';
    mpc.FlexibilityReq.up =   [0  3 2  0  0 4  0  0  0  0]';
    mpc.FlexibilityReq.down = [-1 0 0 -2 -1 0 -5 -5 -5 -5]';
    
    %in case N is very small, N < size(mpc.FlexibilityReq.XXX)
    assert(size(mpc.FlexibilityReq.up,2)==1);%must be col. vector
    assert(size(mpc.FlexibilityReq.down,2)==1);%must be col. vector
    mpc.FlexibilityReq.up = mpc.FlexibilityReq.up(1:N);
    mpc.FlexibilityReq.down = mpc.FlexibilityReq.down(1:N);
    assert(length(mpc.FlexibilityReq.up) == N);
    assert(length(mpc.FlexibilityReq.down) == N);
    assert(sum(mpc.FlexibilityReq.up .* mpc.FlexibilityReq.down) == 0 );
end

%% run OPF
mpcN_opf_storage = create_storage_case_file3(mpc, p_storage);
[RESULTS, SUCCESS] = runopf(mpcN_opf_storage, mpopt);

if ~isfield(mpcN_opf_storage, 'storageFlexibility') || mpcN_opf_storage.storageFlexibility == 0
    plot_storage_results(RESULTS)
end

plot_flexibility_results(RESULTS)
plot_demand_shift_results(RESULTS)

end