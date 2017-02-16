function mpc_storage = add_storage(mpc,nnodes,storage_nodes,P_storage_max_MW,P_storage_min_MW, E_storage_max_MWh, E_storage_init_MWh,c_discharge, c_charge,T_timestep_hours,ramp_max,ramp_min)
    
    %defines whether to enable generator ramping constraint and specifies
    %its limits (actual limit is value given multiplied by 100)
    if (ramp_min == -Inf && ramp_max == Inf)
        RAMPS = 0; %do not even create ramping constraints
    else
        RAMPS = 1;
    end
    
    %%
    mpc_storage      = mpc;
    nnodes_total     = size(mpc.bus,1); %buses for all time periods
    if mod(nnodes_total,nnodes)
        error('wrong number of nodes')
    end
    N                = nnodes_total/nnodes; %number of time periods
    nstorage         = length(storage_nodes); %length of placement of storages (load buses)
    storage_nodesN   = repmat(storage_nodes,[N,1])+kron( (0:(N-1))' , ones(nstorage,1)*nnodes ); %replicate storage nodes and add offsets for each timestep
    
    if find(mpc_storage.bus(storage_nodesN,2)==3)
        error('can not add storage to slack bus')
    end
    mpc_storage.bus(storage_nodesN,2) = 2; %Bus type 2 = PV


    n_gen_cols = size(mpc_storage.gen,2); %number of columns in generator table
%% generator data
% Place storages into generator table, and replicate N times for each time
% period. Do this two times, for charge and discharge.
% Put proper bus ID, ones into VG, GEN_STATUS columns and set proper
% PMAX/PMIN - otherwise set zeros.
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
    mpc_storage.gen  = [mpc_storage.gen;
                       [storage_nodesN, repmat( [ones(nstorage,1)*[0 0 0 0 1 0 1],  P_storage_max_MW,  0*P_storage_max_MW, zeros(nstorage,n_gen_cols-10)], [N,1] )   ];   % discharger
                       [storage_nodesN, repmat( [ones(nstorage,1)*[0 0 0 0 1 0 1],0*P_storage_max_MW,   P_storage_min_MW, zeros(nstorage,n_gen_cols-10)], [N,1] )   ] ]; %  charger
                   
% %% generator cost data
% %   1   startup shutdown    n   x1  y1  ... xn  yn
% %   2   startup shutdown    n   c(n-1)  ... c0
%  mpc.gencost = ones(ngens,1)*[2 0 0 3 0 200 0];
mpc_storage.gencost = [mpc_storage.gencost;
                       ones(nstorage*N*2,1)*[2 0 0 3 0 0 0]];


%% add userconstraints
%% x                 = [theta_bus, Vm_bus, P_gen, Q_gen]
%% P_gen             = [ Pgen_1  ...  Pgen_N, Pdischarge_1 ... Pdischarge_N , Pcharge_1 ... Pcharge_N]
%% Pdischarge_i      = [ Pdischarge_1_i  ...   Pdischarge_nstorage_i]
%% Pcharge_i         = [ Pcharge_1_i  ...   Pcharge_nstorage_i]
ngen_nostorage_total = size(mpc.gen,1);
ngen_total           = size(mpc_storage.gen,1);
% size of A is [NR,NC]
% where NR is:
%   (N+1)*number of storages (start=end state and each period limits) for storage based constraints
%   ng*(N-1) for ramping limits
NR = nstorage+N*nstorage;
if RAMPS
    ng = ngen_nostorage_total/N;  NR = NR + ng*(N-1);
end
% NC is size of x - see above; nnodes = [theta,Vm] and ngen = [P Q]
NC = 2*(nnodes_total+ngen_total);
A                    = sparse(NR,NC);
l                    = zeros(NR,1); 
u                    = l;
offset               = 0; %count of already constructed linear constraints


%diagonal matrices nstorage x nstorage with charging/discharging coeffs.
M_diag_discharge = sparse(1:nstorage,1:nstorage,1./c_discharge);
M_diag_charge = sparse(1:nstorage,1:nstorage,c_charge);


%% sum(Pdischarge+Pcharge)_1^N = 0 
%ordering of the variables:
% theta_bus, Vm_bus, generators -- storages (discharge/charge)-- Q gen
A(offset+(1:nstorage),2*nnodes_total+ngen_nostorage_total+(1:(2*N*nstorage)) ) = ...
    [repmat(M_diag_discharge,[1,N]),...
     repmat(M_diag_charge   ,[1,N]) ];
l(offset+(1:nstorage))   = -Inf; 
u(offset+(1:nstorage))   = Inf;
offset = offset + nstorage;

%% 0 < Einit + T*P1 < Emax
%% 0 < Einit + T*(P1+P2) < Emax ...
%couples charge/discarge of the single storage in current period and all
%previous periods -> Pcharge_i_t - Pdischarge_i_t
A(offset+(1:(N*nstorage)),2*nnodes_total+ngen_nostorage_total+(1:(2*N*nstorage))  ) = ...
    [ -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_discharge), ...
      -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_charge)    ]; %TODO should be +kron() here???
                                                                              
l(offset+(1:N*nstorage)) = repmat(-E_storage_init_MWh,[N,1]);
u(offset+(1:N*nstorage)) = repmat(E_storage_max_MWh-E_storage_init_MWh,[N,1]);
offset = offset + N*nstorage;

%% Ramping limits on real power
%% dp_min < p_g_t+1 - p_g_t < dp_max for all gen g and periods N
if RAMPS
    A(offset+(1:ng*(N-1)), 2*nnodes_total+(1:ngen_nostorage_total)) = ... 
        [kron(diag(ones(N-1,1)), -eye(ng)), zeros((N-1)*ng,ng)] + [zeros((N-1)*ng,ng), kron(diag(ones(N-1,1)), eye(ng))];
    l(offset+(1:ng*(N-1))) = repmat( ramp_min * ones(ng,1), [N-1,1]);
    u(offset+(1:ng*(N-1))) = repmat( ramp_max * ones(ng,1), [N-1,1]);
    offset = offset + ng*(N-1);
end

assert(offset == NR); %check if we have all constraints

if 0
    figure; spy(A(:,2*nnodes_total+(1:ngen_total))); title('Storage constraints');
    figure; spy(A(:,2*nnodes_total+(1:ngen_total))); title('Storage+ramping constraints');
end

%TODO ramp limits for reactive power


%Matpower uses additional A,l,u fields as user specified constraints to standard OPF
mpc_storage.A                         = sparse(A);
clear A;
mpc_storage.l                         = l;
mpc_storage.u                         = u;

mpc_storage.c_charge                  = c_charge;
mpc_storage.c_discharge               = c_discharge;
mpc_storage.nstorage                  = nstorage;
mpc_storage.storage_nodes             = storage_nodes;

mpc_storage.P_storage_max_MW          = P_storage_max_MW;
mpc_storage.P_storage_min_MW          = P_storage_min_MW;
mpc_storage.E_storage_init_MWh        = E_storage_init_MWh;
mpc_storage.E_storage_max_MWh         = E_storage_max_MWh;

mpc_storage.id_gen_storages_discharge = ngen_nostorage_total+(1:(N*nstorage));
mpc_storage.id_gen_storages_charge    = ngen_nostorage_total+N*nstorage+(1:(N*nstorage));
mpc_storage.id_x_gen_storages_discharge = 2*nnodes_total+mpc_storage.id_gen_storages_discharge;
mpc_storage.id_x_gen_storages_charge = 2*nnodes_total+mpc_storage.id_gen_storages_charge;


end
