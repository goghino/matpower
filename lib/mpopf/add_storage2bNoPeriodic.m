function mpc_storage = add_storage(mpc,nnodes,N,storage_nodes_idx,P_discharge_max_MW,P_charge_min_MW, E_storage_max_MWh, E_storage_init_MWh,c_discharge, c_charge,T_timestep_hours, storageFlexibility, storageFlexibilityReq)
    
    mpc_storage      = mpc;
    nstorage         = length(storage_nodes_idx);
    nstorageN        = N*nstorage;
    if mod(size(mpc.bus,1),nnodes)
        error('wrong number of nodes, MPOPF bus count not a multiple of base OPF')
    end
    
    %distinguish between replicated mpc and single period mpc
    % if provided with mpc replicated N times, than replicate also
    % storages, otherwise keep only single period.
    % Lin. constraints are always expanded to full horizon and replicated.
    if(N == size(mpc.bus,1)/nnodes)
        Nkron = N;
        ngen_nostorageN = size(mpc.gen, 1);
        ngenerators      = size(mpc.gen,1)/N;
    else
        Nkron = 1;
        ngen_nostorageN = size(mpc.gen,1) * N;
        ngenerators      = size(mpc.gen,1);
    end
    nnodesN = nnodes * N;
    storage_nodes_idxN = repmat(storage_nodes_idx,[Nkron,1])+kron( (0:(Nkron-1))' , ones(nstorage,1)*nnodes ); %proper bus number for different periods N
    
    ngen_storageN =  nstorageN*2; %charge/discharge
    ngen_flexibilityN = 0;
    if storageFlexibility
        ngen_flexibilityN = nstorageN*4; %up/down charg, up/down discharg.
    end
    ngenN_total = ngen_nostorageN + ngen_storageN + ngen_flexibilityN;
    
    if find(mpc_storage.bus(storage_nodes_idxN,2)==3)
        error('can not add storage to slack bus')
    end
    mpc_storage.bus(storage_nodes_idxN,2) = 2; %BUS_TYPE = PV


    
%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
    Qmin = 0;
    Qmax = 0;%1e-10;
    n_gen_cols = size(mpc_storage.gen,2);
    mpc_storage.gen  = [mpc_storage.gen; %	Pg	Qg	Qmax	Qmin	Vg	mBase	status	    Pmax	          Pmin	            Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf   
                       [storage_nodes_idxN, repmat( [ones(nstorage,1)*[0 0 Qmax Qmin 1 0 1],  P_discharge_max_MW, 0*P_discharge_max_MW, zeros(nstorage,n_gen_cols-10)], [Nkron,1] )   ];   % discharger
                       [storage_nodes_idxN, repmat( [ones(nstorage,1)*[0 0 Qmax Qmin 1 0 1],0*P_charge_min_MW,   P_charge_min_MW, zeros(nstorage,n_gen_cols-10)], [Nkron,1] )   ]; %  charger
                       ];
    %discharging 0 < P < Pmax 
    %charging Pmin < P < 0   (Pmin is negative)
    
    if storageFlexibility                   
    mpc_storage.gen  = [mpc_storage.gen; %	Pg	Qg	Qmax	Qmin	Vg	mBase	status	  Pmax	               Pmin	              Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf   
                       [storage_nodes_idxN, repmat( [ones(nstorage,1)*[0 0 Qmax Qmin 1 0 1],  P_discharge_max_MW,   0*P_discharge_max_MW, zeros(nstorage,n_gen_cols-10)], [Nkron,1] )   ];   % discharger up-flexibility
                       [storage_nodes_idxN, repmat( [ones(nstorage,1)*[0 0 Qmax Qmin 1 0 1],  0*P_discharge_max_MW,  -P_discharge_max_MW, zeros(nstorage,n_gen_cols-10)], [Nkron,1] )   ];   % discharger down-flexibility
                       [storage_nodes_idxN, repmat( [ones(nstorage,1)*[0 0 Qmax Qmin 1 0 1],  -P_charge_min_MW,   0*P_charge_min_MW, zeros(nstorage,n_gen_cols-10)], [Nkron,1] )   ]; %  charger up-flexibility
                       [storage_nodes_idxN, repmat( [ones(nstorage,1)*[0 0 Qmax Qmin 1 0 1],  0*P_charge_min_MW,  P_charge_min_MW, zeros(nstorage,n_gen_cols-10)], [Nkron,1] )   ]; %  charger down-flexibility
                       ];    
    end
    %up-discharge    0 < P < Pmax   
    %down-discharge  -Pmax < P < 0
    %up-charge       0 < P < -Pmin    (Pmin is negative)      
    %down-charge     Pmin < P < 0     (Pmin is negative)
    
                   
% %% generator cost data
% %   polynomial cost function
% %   1   startup shutdown    n   x1  y1  ... xn  yn
% %   2   startup shutdown    n   c(n-1)  ... c0
%  mpc.gencost = ones(ngens,1)*[2 0 0 3 0 200 0];
% % piecewise linear cost function
% %   1   startup shutdown    n   p1  f1  ... pn  fn
MODEL = 1;
PWLINEAR = 1; POLY = 2;
if(mpc_storage.gencost(1,MODEL) == POLY) %POLYNOMIAL MODEL
    NCOST = 4;
    NC = mpc_storage.gencost(1, NCOST);
    gencost = [POLY 0 0 NC zeros(1, NC)]; %polynomial cost, NCOST=3, c0=c1=c2=0 (no cost)

    mpc_storage.gencost = [mpc_storage.gencost;
                           ones(nstorage*Nkron*2,1)*gencost]; 
                       
   if storageFlexibility
    mpc_storage.gencost = [mpc_storage.gencost;
                           ones(nstorage*Nkron*4,1)*gencost];     
   end
else
    NCOST = 4;
    NC = mpc_storage.gencost(1, NCOST);
    coeffs = mpc_storage.gencost(1,NCOST+1:end);
    coeffs(2:2:end) = 0; %set 0 to even entries, parameter f_i
    gencost = [PWLINEAR 0 0 NC coeffs]; %piecewise linear cost, f0=f1=f2=fn=0 (no cost)

    mpc_storage.gencost = [mpc_storage.gencost;
                           ones(nstorage*Nkron*2,1)*gencost]; 
                       
   if storageFlexibility
    mpc_storage.gencost = [mpc_storage.gencost;
                           ones(nstorage*Nkron*4,1)*gencost];     
   end                       
end

%% add user constraints
%% x                 = [theta_bus, Vm_bus, P_gen, P_flex, Q_gen, Q_flex]
%% P_gen             = [ Pgen_1  ...  Pgen_N, Pdischarge_1 ... Pdischarge_N , Pcharge_1 ... Pcharge_N]
%% Pdischarge_i      = [ Pdischarge_1_i  ...   Pdischarge_nstorage_i]
%% Pcharge_i         = [ Pcharge_1_i  ...   Pcharge_nstorage_i]
%% P_flex            = [ Pdis_up_1...N Pdis_down_1...N  Pch_up_1...N Pch_down_1...N]
%nA = (N+1)*nstorage;

nA = nstorageN; %number of lin. constr. (energy of each storage at each time interval)
nX = 2*(nnodesN+ngenN_total); %number of variables, |x|

if storageFlexibility
    nA_flexibility_limits = nstorageN*2; %storage (discharge + up/down) and (charge + up/down) flexibility is 0 < ... < Pmax and Pmin < ... < 0
    nA_flexibility_requirements = 2*N; %sum of up flexibility > requirement or sum of down flex < requirement
    nA = nA + nA_flexibility_limits + nA_flexibility_requirements; 
end

A                    = sparse(nA,nX);
l                    = zeros(nA,1); 
u                    = l;


%% 0 < Einit + T*P1 < Emax
%% 0 < Einit + T*(P1+P2) < Emax ...
% T_timestep_hours = 1
M_diag_discharge = sparse(1:nstorage,1:nstorage,1./c_discharge);
M_diag_charge = sparse(1:nstorage,1:nstorage,c_charge);

A(1:(nstorageN),2*nnodesN+ngen_nostorageN+(1:(ngen_storageN))  ) = ...
    [ -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_discharge), ...
      -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_charge)    ];

%overall energy needs to account for the flexibility as well  
if storageFlexibility
    A(1:(nstorageN),2*nnodesN+ngen_nostorageN+ngen_storageN+(1:(ngen_flexibilityN))  ) = ...
        [-kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_discharge), ...
        -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_discharge), ...
        -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_charge), ...
        -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_charge)];
end

l(1:nstorageN) = repmat(-E_storage_init_MWh,[N,1]);
u(1:nstorageN) = repmat(E_storage_max_MWh-E_storage_init_MWh,[N,1]);



%% 0 < p_sd + u_d < p_max
%% 0 < p_sd + d_d < p_max
%% pmin < p_sc + u_c < 0 (pmin is negative)
%% pmin < p_sc + d_c < 0 (pmin is negative)
if storageFlexibility
    I = mpc.baseMVA*speye(nstorageN);
    z = sparse(nstorageN, nstorageN);

    A(nstorageN+(1:nA_flexibility_limits), 2*nnodesN+ngen_nostorageN+(1:(ngen_storageN+ngen_flexibilityN))  ) = ...
        [ I  z  I  I  z  z; ...
          z  I  z  z  I  I; ]; %assumes complementarity of ud*dd and uc*dc
       % Pd Pc  ud dd uc dc
       
    l(nstorageN+(1:nA_flexibility_limits)) = [zeros(nstorageN, 1); repmat(P_charge_min_MW,[N,1]) ];    
    u(nstorageN+(1:nA_flexibility_limits)) = [repmat(P_discharge_max_MW,[N,1]); zeros(nstorageN, 1) ];  
      
    %no assumptions on complementarity    
    %[ I  z  I  z  z  z; ...
    %  I  z  z  I  z  z; ...
    %  z  I  z  z  I  z; ...
    %  z  I  z  z  z  I; ];
    %l((nstorageN+1):(nstorageN+ngen_flexibilityN)) = [zeros(2*nstorageN, 1); repmat(P_charge_min_MW,[2*N,1]) ];    
    %u((nstorageN+1):(nstorageN+ngen_flexibilityN)) = [repmat(P_discharge_max_MW,[2*N,1]); zeros(2*nstorageN, 1) ];      
end


%% u     < uc_1 + ud_1 < u+10%
%% u     < uc_2 + ud_2 < u+10%
%% d-10% < dc_1 + dd_1 < d
%% d-10% < dc_2 + dd_2 < d
if storageFlexibility 
    e = mpc.baseMVA*ones(1,nstorage);   
    A(nstorageN+nA_flexibility_limits+(1:nA_flexibility_requirements), 2*nnodesN+ngen_nostorageN+ngen_storageN+(1:ngen_flexibilityN)  ) = ...
        [kron(eye(N),e),        sparse(N,nstorageN), kron(eye(N),e),      sparse(N,nstorageN); ...
         sparse(N,nstorageN),   kron(eye(N),e),      sparse(N,nstorageN), kron(eye(N),e)];
     %   ud1 ud2...udN                                 uc1  uc2...ucN        
     %                           dd1 dd2...ddN                             dd1   dd2...ddN
       
    assert(length(storageFlexibilityReq.up) == N);
    assert(length(storageFlexibilityReq.down) == N);
    assert(all(storageFlexibilityReq.up >= 0)); % up flex must be non-negative
    assert(all(storageFlexibilityReq.down <= 0)); % down flex must be non-positive
     
    l(nstorageN+nA_flexibility_limits+(1:nA_flexibility_requirements)) = [storageFlexibilityReq.up; 1.1.*storageFlexibilityReq.down];
    u(nstorageN+nA_flexibility_limits+(1:nA_flexibility_requirements)) = [1.1.*storageFlexibilityReq.up; storageFlexibilityReq.down];
   
end
%%
mpc_storage.A                           = sparse(A);
clear A;
mpc_storage.l                           = l;
mpc_storage.u                           = u;

mpc_storage.c_charge                    = c_charge;
mpc_storage.c_discharge                 = c_discharge;
mpc_storage.nstorage                    = nstorage;
mpc_storage.ngenerators                 = ngenerators;
mpc_storage.storage_nodes               = storage_nodes_idx;
mpc_storage.storageFlexibility          = storageFlexibility;
mpc_storage.storageFlexibilityReq       = storageFlexibilityReq;

mpc_storage.P_storage_max_MW            = P_discharge_max_MW;
mpc_storage.P_storage_min_MW            = P_charge_min_MW;
mpc_storage.E_storage_init_MWh          = E_storage_init_MWh;
mpc_storage.E_storage_max_MWh           = E_storage_max_MWh;

mpc_storage.id_gen_storages_discharge   = ngen_nostorageN+(1:(Nkron*nstorage));
mpc_storage.id_gen_storages_charge      = ngen_nostorageN+Nkron*nstorage+(1:(Nkron*nstorage));
mpc_storage.id_x_gen_storages_discharge = 2*nnodesN+mpc_storage.id_gen_storages_discharge;
mpc_storage.id_x_gen_storages_charge    = 2*nnodesN+mpc_storage.id_gen_storages_charge;


end
