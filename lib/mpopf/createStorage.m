function storage = createStorage(mpc, Rfirst, Rcount, Emax)
define_constants;

if (Rcount < 0)
    [load_sorted, load_sorted_Idx] = sort(mpc.bus(:,PD), 'descend');
    nload = length(load_sorted);
    nstorage = abs(Rcount);  %% 1 ... 100
    nstorage_applied = min(nstorage,nload);
    first = 0;
    storage_locations = load_sorted_Idx(first + (1:nstorage_applied));
else
    [~, load_sorted] = sort(mpc.bus(:,PD), 'descend'); %% sort buses w.r.t PD
    nload = length(find(abs(mpc.bus(:,PD)) > 0));
    nstorage = round(nload * Rcount); %% apply storage N% of load buses
    first = round(nload * Rfirst); %% place storages starting at 'first' load bus
    nstorage_applied = max(min(nstorage,nload), 1); %% use at least 1 storage, use maximum nload storages
    storage_locations = load_sorted(first + (1:nstorage_applied)); %% apply storages to top 2% loaded buses
end

% Do not place storage to REF bus
ref_idx = find(mpc.bus(storage_locations, BUS_TYPE) == 3);
if(ref_idx)
    storage_locations(ref_idx) = [];
    storage_locations = [storage_locations; load_sorted(first + nstorage_applied+1)];
end

if(nstorage_applied < 1)
   error('Number of storateges has to be > 0'); 
end

storage.id_storage_location = storage_locations; %storage placement
storage.E_storage_max_MWh  = Emax*abs(mpc.bus(storage_locations,PD)); %max storage level
storage.E_storage_init_MWh = storage.E_storage_max_MWh*.7; %initial charge level
storage.rPmaxEmax_MW_per_MWh = 1/3;
storage.rPminEmax_MW_per_MWh = -1/2;
storage.P_storage_max_MW   =  storage.rPmaxEmax_MW_per_MWh*storage.E_storage_max_MWh;
storage.P_storage_min_MW   =  storage.rPminEmax_MW_per_MWh*storage.E_storage_max_MWh;
storage.c_discharge        = .97;
storage.c_charge           = .95;

end