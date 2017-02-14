clear all;
clc
close all
% addpath('/home/alex/checkout/gitFEN/matlab/matpower51')
% addpath('/home/alex/checkout/gitFEN/matlab/matpower51/t')
define_constants

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TBD 
%% write down problem structure on paper
%% create case
%%   where long horizon is needed hydro storage and aggregated flexibility?
%%   case file with 118 and EU
%%   scalable storages#

%% set options

%opt = mpoption('verbose',2,'out.all',0, 'opf.ac.solver','IPOPT');
opt = mpoption('verbose',2,'out.all',0);

%% set time step data
nstorage_ref = 3;  %% 1 ... 100
factor_timesteps = 1;  %% 1...365
Kpv = 1;  %% 0 ...1  (size of PV penetration)

%load scaling profile for each period
load_scaling_profile0 = [0.4544 0.3570 0.2860 0.2783 0.3795 0.5822 0.8086 0.9633 1.0086 0.9883 0.9761 1.0000 1.0193 0.9773 0.8772 0.7991 0.8359 1.0023 1.2063 1.3123 1.2438 1.0343 0.7873 0.5885]';
pv_scaling_profile =  [  0 0 0  0 0 0.0046 0.0548 0.1686 0.3457 0.5100 0.6687 0.7496 0.8175 0.8305 0.8026 0.7212 0.5988 0.4453 0.2718 0.1203 0.0350 0.0019 0 0 ]'  ;
save -ascii -double 'load.dat' load_scaling_profile0;
save -ascii -double 'pv.dat' pv_scaling_profile;

load_scaling_profile = load_scaling_profile0 - 1*pv_scaling_profile;
save -ascii -double 'loadpv.dat' load_scaling_profile;

figure;
plot(load_scaling_profile0,'b--');
hold on;
plot(load_scaling_profile,'b');
legend('load','net load with PV injection')


load_scaling_profile       = kron(ones(factor_timesteps,1), load_scaling_profile  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create base case file

mpc = case118;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fixes:
mpc = ext2int(mpc); %remove any isolated buses or branches and reorder remaining buses in inc order
mpc.branch(mpc.branch(:,RATE_A)==0,RATE_A) = 9900; %redefine zero MVA rating A (long term rating) branches 
mpc.gen(:,PMIN) = 0;
%mpc.branch(:,RATE_A) = 9900;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% prepare storage data
id_load_log = abs(mpc.bus(:,PD))>0; %find real power demand (MW) buses
id_load = find(id_load_log);
nload = length(id_load); %how many load buses do we have
nstorage_applied = min(nstorage_ref,nload); %how many storages we use (max nstorage_ref or less if there are less loads)
id_storage_location = id_load(1:nstorage_applied); %place storages to load buses
if find(mpc.bus(id_storage_location,2) == 3) %if BUS_TYPE == 3 (aka ref)
    nstorage_applied = min(nstorage_ref+1,nload); %add single additional storage (if aplicable)
    id_storage_location = id_load(1:nstorage_applied);
    id_storage_location(mpc.bus(id_storage_location,2) == 3) = [];
end

%set some network propertis related to storages
p_storage.id_storage_location = id_storage_location;
p_storage.E_storage_max_MWh  = 2*abs(mpc.bus(id_storage_location,PD));
p_storage.E_storage_init_MWh = p_storage.E_storage_max_MWh*.7;
p_storage.rPmaxEmax_MW_per_MWh = 1/3;
p_storage.rPminEmax_MW_per_MWh = -1/2;
p_storage.c_discharge        = .97;
p_storage.c_charge           = .95;

mpcN_opf_storage = create_storage_case_file3(mpc,load_scaling_profile, p_storage)
resN_opf_storage = runopf(mpcN_opf_storage, opt)
resN_opf_storage.success

%% Plot actual ramping behavior of PG - real power output (MW)
%in gen we have firs data for generators for all time periods, only then storages
ng = size(mpc.gen,1);
NG = ng * length(load_scaling_profile);
ri = 76:ng:1317; %TODO: hardcoded number for the specific testcase
ci = 5665:ng:6960;
figure; hold on;
plot(1:size(ri,2), repmat(150, [1, size(ri,2)]), 'r-');
plot(1:size(ri,2), repmat(-150, [1, size(ri,2)]), 'r-');
for i = 1:ng
    ramps = mpcN_opf_storage.A(ri,ci)*resN_opf_storage.gen(i:ng:NG, PG);
    plot(ramps);
end
title('Generator ramping over time horizon')




