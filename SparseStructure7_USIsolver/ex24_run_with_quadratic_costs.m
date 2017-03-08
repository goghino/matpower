%close all
% addpath('/home/alex/checkout/gitFEN/matlab/matpower51')
% addpath('/home/alex/checkout/gitFEN/matlab/matpower51/t')
define_constants

%defines whether to enable explicit generator ramping constraint and specifies
%its limits (actual limit is value given multiplied by BaseMVA=100). If ramp
%is not specified ramping limits are enforced in between IPM iterations
%by checking iterate solutions x_k
RAMP = 1;

ramp_max = 120;
ramp_min = -80;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TBD 
%% write down problem structure on paper
%% create case
%%   where long horizon is needed hydro storage and aggregated flexibility?
%%   case file with 118 and EU
%%   scalable storages#

%% set options

opt = mpoption('verbose',2,'out.all',0, 'opf.ac.solver','MIPS'); %'opf.ac.solver','IPOPT'
opt = mpoption(opt, 'ipopt.opts', struct('tol', 1e-04));
setenv('OMP_NUM_THREADS', '1');
%opt = mpoption('verbose',2,'out.all',0);

%% set time step data
nstorage_ref = 3;  %% number of storages in network
factor_timesteps = 1;  %% 1...365
Kpv = 1;  %% 0 ...1  (size of PV penetration)

%load scaling profile for each period
load_scaling_profile0 = [0.0544 0.0544 0.0544 0.1544 0.2544 0.4544 0.3570 0.2860 0.2783 0.3795 0.5822 0.8086 0.9633 1.0086 0.9883 0.9761 1.0000 1.0193 0.9773 0.8772 0.7991 0.8359 1.0023 1.2063]';
pv_scaling_profile =  [0 0 0 0 0 0 0 0  0 0 0.0046 0.0548 0.1686 0.3457 0.5100 0.6687 0.7496 0.8175 0.8305 0.8026 0.7212 0.5988 0.4453 0.2718 ]'  ;
save -ascii -double 'load.dat' load_scaling_profile0;
save -ascii -double 'pv.dat' pv_scaling_profile;

load_scaling_profile = load_scaling_profile0 - 1*pv_scaling_profile;
save -ascii -double 'loadpv.dat' load_scaling_profile;

% figure;
% plot(load_scaling_profile0,'b--');
% hold on;
% plot(load_scaling_profile,'b');
% legend('load','net load with PV injection')

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
%network properties related to generator ramps
p_storage.RAMP               = RAMP;
p_storage.ramp_max           = ramp_max;
p_storage.ramp_min           = ramp_min;

mpcN_opf_storage = create_storage_case_file3(mpc,load_scaling_profile, p_storage);
resN_opf_storage = runopf(mpcN_opf_storage, opt);

disp('[success, f]');
disp([resN_opf_storage.success, resN_opf_storage.f]);

%% Plot actual ramping behavior of PG - real power output (MW)
%in gen we have firs data for generators for all time periods, only then storages
ng = size(mpc.gen,1);
NG = ng * length(load_scaling_profile);

%number of preceeding constraints in A for storages (N+1)*nstorages
r_offset = (length(load_scaling_profile)+1) * length(id_storage_location);
%rows of A matrix for ramping constraints of g1, constraints are always
%same g1=g2=g3=... ([-I I]) so we extract just for first generator and reuse it
%for all the rest 1..ng
ri = r_offset+(1:ng:(NG-ng));

%skip comuns in A corresponding to bus variables (angle, magnitude)
c_offset = 2*size(mpcN_opf_storage.bus,1);
%columns of A matrix, same for all generators
ci = c_offset + (1:ng:NG);

figure; subplot(2,1,1); hold on; title('Generator ramping over time horizon');
xlabel('N'); ylabel('Real power ramp [MW]');
plot(1:size(ri,2), repmat(p_storage.ramp_max, [1, size(ri,2)]), 'r-');
plot(1:size(ri,2), repmat(p_storage.ramp_min, [1, size(ri,2)]), 'r-');
for i = 1:ng
    %ramp rates of i-th generator, PG_t+1 - PG_t
    if RAMP
        ramps = mpcN_opf_storage.A(ri,ci)*resN_opf_storage.gen(i:ng:NG, PG);
    else
        ramps = resN_opf_storage.gen((i+ng):ng:NG, PG) - resN_opf_storage.gen(i:ng:(NG-ng), PG);
    end
    plot(ramps);
end

%plot generator profile
subplot(2,1,2); hold on; title('Generator power profile');
xlabel('N'); ylabel('Real power output [MW]');
for i = 1:ng
    plot(resN_opf_storage.gen(i:ng:NG, PG))
end
