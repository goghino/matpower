function [RESULTS, SUCCESS] = MPOPF(mpc, mpopt, N, Emax, Rcount, Rfirst)
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
%do we place static number of storages to all grids, or we pick top 2% of the load buses?
%Rcount = -50; %if negative, use this number of storages abs(Rcoung), otherwise
%use Rcount % top of the load buses 
if (Rcount < 0)
    load_sorted = find(abs(mpc.bus(:,PD)) > 0);
    nload = length(load_sorted);
    nstorage = abs(Rcount);  %% 1 ... 100
    nstorage_applied = min(nstorage,nload);
    first = 0;
    id_storage_location = load_sorted(first + (1:nstorage_applied));
else
    [sortedPD, load_sorted] = sort(mpc.bus(:,PD), 'descend'); %% sort buses w.r.t PD
    nload = length(find(abs(mpc.bus(:,PD)) > 0));
    nstorage = round(nload * Rcount); %% apply storage N% of load buses
    first = round(nload * Rfirst); %% place storages starting at 'first' load bus
    nstorage_applied = max(min(nstorage,nload), 1); %% use at least 1 storage, use maximum nload storages
    id_storage_location = load_sorted(first + (1:nstorage_applied)); %% apply storages to top 2% loaded buses
end

% Do not place storage to REF bus
ref_idx = find(mpc.bus(id_storage_location, BUS_TYPE) == 3);
if(ref_idx)
    id_storage_location(ref_idx) = [];
    id_storage_location = [id_storage_location; load_sorted(first + nstorage_applied+1)];
end

if(nstorage_applied < 1)
   error('Number of storateges has to be > 0'); 
end

p_storage.id_storage_location = id_storage_location;
p_storage.E_storage_max_MWh  = Emax*abs(mpc.bus(id_storage_location,PD));
p_storage.E_storage_init_MWh = p_storage.E_storage_max_MWh*.7;
p_storage.rPmaxEmax_MW_per_MWh = 1/3;
p_storage.rPminEmax_MW_per_MWh = -1/2;
p_storage.c_discharge        = .97;
p_storage.c_charge           = .95;

%% set time step data
SelectProfile = 1;

if (SelectProfile == 0)
    Kpv = 1;  %% 0 ...1  (size of PV penetration)

    load_scaling_profile0 = [0.4544 0.3570 0.2860 0.2783 0.3795 0.5822 0.8086 0.9633 1.0086 0.9883 0.9761 1.0000 1.0193 0.9773 0.8772 0.7991 0.8359 1.0023 1.2063 1.3123 1.2438 1.0343 0.7873 0.5885]';
    pv_scaling_profile    = [  0 0 0  0 0 0.0046 0.0548 0.1686 0.3457 0.5100 0.6687 0.7496 0.8175 0.8305 0.8026 0.7212 0.5988 0.4453 0.2718 0.1203 0.0350 0.0019 0 0 ]';
    load_scaling_profile  = load_scaling_profile0 - Kpv*pv_scaling_profile;
elseif (SelectProfile == 1)
    storage_load = abs(p_storage.E_storage_max_MWh * p_storage.rPminEmax_MW_per_MWh);
    storage_injection = p_storage.E_storage_max_MWh * p_storage.rPmaxEmax_MW_per_MWh;
    load_scaling_profile       = createLoadProfile(mpc, storage_injection, storage_load);
    Nprofile = size(load_scaling_profile,1);
    if (N <= Nprofile)
        load_scaling_profile = load_scaling_profile(1:N);
    else    
        load_scaling_profile = repmat(load_scaling_profile, ceil(N/Nprofile), 1);
        load_scaling_profile = load_scaling_profile(1:N);
    end
    fprintf('Running MPOPF with N=%d periods\n', size(load_scaling_profile,1));
end

%repeat data for required no. of days
%load_scaling_profile       = kron(ones(N,1), load_scaling_profile);

%% run OPF
mpcN_opf_storage = create_storage_case_file3(mpc,load_scaling_profile, p_storage);
[RESULTS, SUCCESS] = runopf(mpcN_opf_storage, mpopt);

% plot_storage_results(opt_solution)

end
