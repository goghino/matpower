function mpc_storage = add_storage(mpc,nbus,N,storage_nodes_idx,P_discharge_max_MW,P_charge_min_MW, E_storage_max_MWh, E_storage_init_MWh,c_discharge, c_charge,T_timestep_hours)
    
    mpc_storage      = mpc;
    nstorage         = length(storage_nodes_idx);
    nstorageN        = N*nstorage;
    if mod(size(mpc.bus,1),nbus)
        error('wrong number of nodes, MPOPF bus count is not a multiple of base OPF')
    end
    
    enableDemandShift = 0; storageFlexibility=0;
    if isfield(mpc, 'storageFlexibility')
        storageFlexibility = mpc.storageFlexibility;
    end
    if isfield(mpc, 'enableDemandShift')
        enableDemandShift = mpc.enableDemandShift;
        demandShift = mpc.demandShift;
    end
    
    if enableDemandShift || storageFlexibility
        FlexibilityReq = mpc.FlexibilityReq;
    else
        FlexibilityReq = [];
    end
    
    
    %distinguish between replicated mpc and single period mpc
    % if provided with mpc replicated N times, than replicate also
    % storages, otherwise keep only single period.
    % Lin. constraints are always expanded to full horizon and replicated.
    if(N == size(mpc.bus,1)/nbus)
        Nkron = N;
        ngen_nostorageN = size(mpc.gen, 1);
        ngenerators      = size(mpc.gen,1)/N;
    else
        Nkron = 1;
        ngen_nostorageN = size(mpc.gen,1) * N;
        ngenerators      = size(mpc.gen,1);
    end
    nbusN = nbus * N;
    
    
    ngen_storageN =  nstorageN*2; %charge/discharge
    ngen_storageFlexibilityN = 0;
    if storageFlexibility
        ngen_storageFlexibilityN = nstorageN*4; %up/down charg, up/down discharg.
    end
    ngen_demandShiftN = 0;
    if enableDemandShift
        nbus_demandShift = length(mpc.demandShift.busesID); %buses which provide DS
        nbus_demandShiftN = nbus_demandShift*N; %buses which provide DS, overall count
        ngen_demandShiftN = nbus_demandShift*4*N; %response+upFlex and rebound+downFlex
    end
    
    ngenN_total = ngen_nostorageN + ngen_storageN + ngen_storageFlexibilityN + ngen_demandShiftN;
    
    
    storage_nodes_idxN = repmat(storage_nodes_idx, [Nkron,1]) + kron( (0:(Nkron-1))' , ones(nstorage,1)*nbus ); %proper bus number for different periods N
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
    
    if enableDemandShift
    demandShift_nodes_idxN = repmat(mpc.demandShift.busesID, [Nkron,1])+kron( (0:(Nkron-1))' , ones(nbus_demandShift,1)*nbus ); %proper bus number for different periods N
    Prsp = mpc.demandShift.responsePowerMW;
    Prb = mpc.demandShift.reboundPowerMW;
    %Qmax = max(Prsp)*10;
    %Qmin = -max(Prsp)*10;
    
    mpc_storage.gen  = [mpc_storage.gen; %	Pg	Qg	Qmax	Qmin	Vg	mBase	status	  Pmax	               Pmin	              Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf   
                       [demandShift_nodes_idxN, repmat( [ones(nbus_demandShift,1)*[0 0 Qmax Qmin 1 0 1],  0*Prsp,   -Prsp, zeros(nbus_demandShift,n_gen_cols-10)], [Nkron,1] )   ];   % response up-shift p^du
                       [demandShift_nodes_idxN, repmat( [ones(nbus_demandShift,1)*[0 0 Qmax Qmin 1 0 1],  Prsp,  0*Prsp, zeros(nbus_demandShift,n_gen_cols-10)], [Nkron,1] )   ];   % up flexibility u^d
                       [demandShift_nodes_idxN, repmat( [ones(nbus_demandShift,1)*[0 0 Qmax Qmin 1 0 1],  Prb,   0*Prb, zeros(nbus_demandShift,n_gen_cols-10)], [Nkron,1] )   ]; % response down-shift p^dd
                       [demandShift_nodes_idxN, repmat( [ones(nbus_demandShift,1)*[0 0 Qmax Qmin 1 0 1],  0*Prb,  -Prb, zeros(nbus_demandShift,n_gen_cols-10)], [Nkron,1] )   ];   % down flexibility d^d
                       ];    
    end
    %up-shift          -Prsp < P < 0
    %up flexibility    0     < P < Prsp
    %down-shift        0     < P < Prb
    %down flexibility  -Prb  < P < 0
    
                   
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
   
   if enableDemandShift
    mpc_storage.gencost = [mpc_storage.gencost;
                           ones(nbus_demandShift*Nkron*4,1)*gencost];     
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
   
   if enableDemandShift
    mpc_storage.gencost = [mpc_storage.gencost;
                           ones(nbus_demandShift*Nkron*4,1)*gencost];     
   end
end

%% add user constraints
%% x                 = [theta_bus, Vm_bus, P_gen, P_flex, Q_gen, Q_flex]
%% P_gen             = [ Pgen_1  ...  Pgen_N, Pdischarge_1 ... Pdischarge_N , Pcharge_1 ... Pcharge_N]
%% Pdischarge_i      = [ Pdischarge_1_i  ...   Pdischarge_nstorage_i]
%% Pcharge_i         = [ Pcharge_1_i  ...   Pcharge_nstorage_i]
%% P_flex            = [ Pdis_up_1...N Pdis_down_1...N  Pch_up_1...N Pch_down_1...N]
%nA = (N+1)*nstorage;

nA = nstorageN; %number of lin. constr. (energy of each storage at each time interval, no init=end constraint)
nX = 2*(nbusN+ngenN_total); %number of variables, |x|

if storageFlexibility || enableDemandShift
    nA_flexibility_requirements = 2*N; %sum of up flexibility > requirement or sum of down flex < requirement
    nA = nA + nA_flexibility_requirements; 
end

if storageFlexibility
    nA_flexibility = nstorageN; %(energy of each storage at each time interval with ES flexibility, no init=end constraint)
    nA_flexibility_limits = nstorageN*2; %storage (discharge + up/down) and (charge + up/down) flexibility is 0 < ... < Pmax and Pmin < ... < 0
    nA = nA + nA_flexibility + nA_flexibility_limits;
end

if enableDemandShift
    Trsp = mpc.demandShift.responseTimeH;
    Trb = mpc.demandShift.reboundTimeH;
    
    nA_demandShiftLimits = nbus_demandShiftN*2; %limits fo the (-pdu+ud) and (pdd-dd)
    nA_demandShift = 0;
    for i = 1:nbus_demandShift
        %t->t+Trsp for each bus for each N
        nA_demandShift = nA_demandShift + 2*N;
        %t+Trsp->t+Trsp+Trb for each bus for each N such that we don't insert empty rows
        nA_demandShift = nA_demandShift + 2*(N-Trsp(i));
    end
    
    nA = nA + nA_demandShiftLimits + nA_demandShift;
    nX = nX + nbus_demandShiftN; %additional binary variables
    
    %limits and initial value of the extra binary variables b=(0,1)
    mpc_storage.z0 = [1;1;1;1;1;1;0;0;0;0;1;1;1;1;1;1;0;0;0;0;0;0;0;0;1;1;1;1;1;1;zeros(2*9,1)]; %0.5*ones(nbus_demandShiftN, 1);
    mpc_storage.zl = zeros(nbus_demandShiftN, 1);
    mpc_storage.zu = ones(nbus_demandShiftN, 1);
end

%variables ordering
%first all Va, then Vm, then Pg from all gens/storages/flex, then all Qg
%[Va Vm] P-Q[gens] P-Q[discharge,charge] P-Q[discharge up/down charge-up/down] P-Q[pdu ud pdd dd] [binary] 
%Pg = [Pgen -- Ped Pec -- ued ded uec dec -- Pdu ud Pdd dd]
%[2nb 2ng 4ns 8ns 8nf nf] where nf is number of demand buses providing demand shift

%constraints ordering
%storage energy levels (with or without flexibility) - N*ns
%storage flexibility limits (assuming complementarity) - 2*N*ns
%....
%flexibility provision requirements 2*N

A                    = sparse(nA,nX);
l                    = zeros(nA,1); 
u                    = l;

%offsets of A in row (offset1) and columns (offset2)
offset1 = 0;
offset2storages    = 2*nbusN+ngen_nostorageN;
offset2flexibility = 2*nbusN+ngen_nostorageN + ngen_storageN;
offset2demandShift = 2*nbusN+ngen_nostorageN + ngen_storageN + ngen_storageFlexibilityN;
offset2binaryVars  = 2*nbusN+2*(ngen_nostorageN + ngen_storageN + ngen_storageFlexibilityN + ngen_demandShiftN);

%% SoC levels for the storages over the time horizont
%% 0 < Einit + T*(P1c+P1d) < Emax
%% 0 < Einit + T*(P1c+P1d+P2c+P2d) < Emax ...
% T_timestep_hours = 1
M_diag_discharge = sparse(1:nstorage,1:nstorage,1./c_discharge);
M_diag_charge = sparse(1:nstorage,1:nstorage,c_charge);

% 0 < Einit + T*(P1c+P1d) < Emax
% 0 < Einit + T*(P1c+P1d+P2c+P2d) < Emax ...
A(offset1+(1:nstorageN), offset2storages+(1:(ngen_storageN))  ) = ...
    [ -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_discharge), ...
      -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_charge)    ];
  
l(offset1+(1:nstorageN)) = repmat(-E_storage_init_MWh,[N,1]);
u(offset1+(1:nstorageN)) = repmat(E_storage_max_MWh-E_storage_init_MWh,[N,1]);

offset1 = offset1 + nstorageN;  

% 0 < Einit + T*(P1c+P1d+u1_d+d1_d+u1_c+d1_c) < Emax
% 0 < Einit + T*(P1c+P1d+u1_d+d1_d+u1_c+d1_c +P2c+P2d+u2_d+d2_d+u2_c+d2_c) < Emax ...
%overall energy needs to account for the flexibility as well  
if storageFlexibility
    A(offset1+(1:nA_flexibility), offset2storages+(1:(ngen_storageN))  ) = ...
    [ -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_discharge), ...
      -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_charge)    ];
    A(offset1+(1:nA_flexibility), offset2flexibility+(1:(ngen_storageFlexibilityN))  ) = ...
        [-kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_discharge), ...
        -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_discharge), ...
        -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_charge), ...
        -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_charge)];

    l(offset1+(1:nA_flexibility)) = repmat(-E_storage_init_MWh,[N,1]);
    u(offset1+(1:nA_flexibility)) = repmat(E_storage_max_MWh-E_storage_init_MWh,[N,1]);

    offset1 = offset1 + nA_flexibility;
end


%% Respect max. charge/discharge rate considering also the flexibility provision
%% 0 < p_sd + ue_d + de_d < p_max
%% pmin < p_sc + ue_c + de_c < 0 (pmin is negative)
if storageFlexibility
    I = mpc.baseMVA*speye(nstorageN);
    z = sparse(nstorageN, nstorageN);

    A(offset1+(1:nA_flexibility_limits), offset2storages+(1:(ngen_storageN+ngen_storageFlexibilityN))  ) = ...
        [ I  z  I  I  z  z; ...
          z  I  z  z  I  I; ]; %assumes complementarity of ued*ded and uec*dec
       % Ped Pec  ued ded uec dec
       
    l(offset1+(1:nA_flexibility_limits)) = [zeros(nstorageN, 1); repmat(P_charge_min_MW,[N,1]) ];    
    u(offset1+(1:nA_flexibility_limits)) = [repmat(P_discharge_max_MW,[N,1]); zeros(nstorageN, 1) ];  
      
    offset1 = offset1 + nA_flexibility_limits;
    
    %alternative set of constraints with no assumptions on complementarity    
    %[ I  z  I  z  z  z; ...
    %  I  z  z  I  z  z; ...
    %  z  I  z  z  I  z; ...
    %  z  I  z  z  z  I; ];
    %l((nstorageN+1):(nstorageN+ngen_flexibilityN)) = [zeros(2*nstorageN, 1); repmat(P_charge_min_MW,[2*N,1]) ];    
    %u((nstorageN+1):(nstorageN+ngen_flexibilityN)) = [repmat(P_discharge_max_MW,[2*N,1]); zeros(2*nstorageN, 1) ];      
    %offset1 = offset1 + ...
end

%% Overall flexibility provisions in each time period
%% u     < ud_1 + uec_1 + ued_1 + ... + ud_ns + uec_ns + ued_ns < u+0.1% @t1
%% u     < ud_1 + uec_1 + ued_1 + ... + ud_ns + uec_ns + ued_ns < u+0.1% @t2
%% d-0.1% < dd_1 + dec_1 + ded_1 + ... + dd_ns + dec_ns + ded_ns < d  @t1
%% d-0.1% < dd_1 + dec_1 + ded_1 + ... + dd_ns + dec_ns + ded_ns < d  @t2
if storageFlexibility 
    e = mpc.baseMVA*ones(1,nstorage);   
    A(offset1+(1:nA_flexibility_requirements), offset2flexibility+(1:ngen_storageFlexibilityN)  ) = ...
        [kron(eye(N),e),        sparse(N,nstorageN), kron(eye(N),e),      sparse(N,nstorageN); ...
         sparse(N,nstorageN),   kron(eye(N),e),      sparse(N,nstorageN), kron(eye(N),e)];
     %   ued1 ued2...uedN                            uec1  uec2...uecN        
     %                           ded1 ded2...dedN                         ded1   ded2...dedN
end

if enableDemandShift
    e = mpc.baseMVA*ones(1,nbus_demandShift);   
    A(offset1+(1:nA_flexibility_requirements), offset2demandShift+(1:(4*nbus_demandShiftN))  ) = ...
        [sparse(N,nbus_demandShiftN), kron(eye(N),e),    sparse(N,2*nbus_demandShiftN); ...
         sparse(N,3*nbus_demandShiftN),                                     kron(eye(N),e)];
     %   Pdu1...PduN (skip)         ud1 ud2...udN          
     %                                               Pdd1...PddN (skip)    dd1 dd2...ddN    
end

if storageFlexibility || enableDemandShift
    assert(length(FlexibilityReq.up) == N);
    assert(length(FlexibilityReq.down) == N);
    assert(all(FlexibilityReq.up >= 0)); % up flex must be non-negative
    assert(all(FlexibilityReq.down <= 0)); % down flex must be non-positive 
    
    l(offset1+(1:nA_flexibility_requirements)) = [FlexibilityReq.up; 1.01.*FlexibilityReq.down];
    u(offset1+(1:nA_flexibility_requirements)) = [1.01.*FlexibilityReq.up; FlexibilityReq.down];
    offset1 = offset1 + nA_flexibility_requirements;
end

%% Demand shift limits (-Pdu+ud = Prsp) and (Pdd-dd = Prb)
%% 0 < -Pdu1 + ud1 - Prsp*b <= 0
%% 0 < -Pdu_nf_N + ud_nf_N - Prsp*b <= 0
%%
%% 0 < Pdd1 - dd1 - Prb*b <= 0
%% 0 < Pdd_nf_N - dd_nf_N - Prb*b <= 0
if enableDemandShift  
    e = mpc.baseMVA*speye(nbus_demandShift);
    
    A(offset1+(1:nA_demandShiftLimits), offset2demandShift + (1:4*nbus_demandShiftN)) = ...
    [kron(eye(N),-e) kron(eye(N),e)   sparse(nbus_demandShiftN, 2*nbus_demandShiftN); ...
     sparse(nbus_demandShiftN, 2*nbus_demandShiftN),  kron(eye(N),e), kron(eye(N),-e)];
    % -Pdu1...-PduN  ud1  ud2...udN                                                  
    %              
    %                                               Pdu1 Pdu2...PduN  -dd1 -dd2..ddN
 
    A(offset1+(1:nA_demandShiftLimits), offset2binaryVars + (1:nbus_demandShiftN)) = ...
    [ kron(eye(N), diag(-Prsp)); ...
      kron(eye(N), diag(-Prb))];
    %b1 b2 bN              
    %b1 b2 bN      
    
    l(offset1+(1:nA_demandShiftLimits)) = [-Inf*ones(nA_demandShiftLimits, 1)];
    u(offset1+(1:nA_demandShiftLimits)) = [zeros(nA_demandShiftLimits, 1)];
 
    offset1 = offset1 + nA_demandShiftLimits;
end

%% Demand shift provision and complementarity
%% sum_{t->t+Trsp} -Pdu+ud >= Trsp*Prsp*(bt - bt-1), where t+Trsp<=N
%% sum_{t+Trsp->t+Trsp+Trb} Pdd-dd >= Trb*Prb*(bt - bt-1), where t+Trsp+Trb<=N
%% complementarity
%% sum_{t->t+Trsp} Pdd-dd <= Trb*Prb*[1 - (bt - bt-1)], where t+Trsp<=N
%% sum_{t+Trsp->t+Trsp+Trb} -Pdu+ud <= Trsp*Prsp*[1- (bt - bt-1)], where t+Trsp+Trb<=N
if enableDemandShift
    
    %prepare the blocks to be inserted into the 'big' A matrix
    %each block represents sum over Trsp/Trsp+Trb for a single bus (one sum for each bus)
    sumBlockTrsp=sparse(nbus_demandShift, nbus_demandShift*max(Trsp));%t->t+Trsp
    sumBlockTrb=sparse(nbus_demandShift, nbus_demandShift*max(Trsp+Trb));%t->t+Trsp+Trb
    binaryBlock = sparse(nbus_demandShift, 2*nbus_demandShift);%[t-1, t]
    for i = 1:nbus_demandShift
       % sum_{t->t+Trsp} XXX
       sumBlockTrsp(i,(i-1)+(1:nbus_demandShift:nbus_demandShift*Trsp(i))) = 1*mpc.baseMVA;
       
       % sum_{t+Trb->t+Trsp+Trb} XXX
       skipTrspVars = nbus_demandShift*Trsp(i);
       sumBlockTrb(i,(skipTrspVars+i-1)+(1:nbus_demandShift:nbus_demandShift*Trb(i))) = 1*mpc.baseMVA; 
       
       % -b(t-1) + b(t)
       binaryBlock(i,:) = [sparse(1,i-1) -1 sparse(1,nbus_demandShift-1) 1 sparse(1,nbus_demandShift-i)];
    end
    
    % insert the 'sum' blocks into A for each time period t
    % sum_{t->t+Trsp} -Pdu+ud >= Trsp*Prsp*(bt - bt-1), where t+Trsp<=N
    % sum_{t+Trsp->t+Trsp+Trb} Pdd-dd >= Trb*Prb*(bt - bt-1), where t+Trsp+Trb<=N
    rowsInsertedRsp = 0;
    rowsInsertedRb = 0;
    colsInserted = 0;
    for t = 1:N
        %make sure that t+Trsp<=N so that we do not overflow to the next set of variables
        sumBlockCurtailedTrsp = sumBlockTrsp(:,1:nbus_demandShift*min(max(Trsp), N-t+1));
        sz_Trsp_row = size(sumBlockCurtailedTrsp,1);
        sz_Trsp_col = size(sumBlockCurtailedTrsp,2);
        %make sure that t+Trsp+Trb<=N so that we do not overflow to the next set of variables
        sumBlockCurtailedTrb = sumBlockTrb(:,1:nbus_demandShift*min(max(Trsp+Trb), N-t+1));
        sz_Trb_row = size(sumBlockCurtailedTrb,1);
        sz_Trb_col = size(sumBlockCurtailedTrb,2);
        assert(sz_Trsp_row==sz_Trb_row);
        assert(sz_Trsp_row==nbus_demandShift);
        
        %and for t==1 we need to curtail the binary block, since there is no t0
        if (t==1)
            binaryBlockCurtailed = -binaryBlock(:,1:nbus_demandShift);
            extraSkip = 0;
        else
            binaryBlockCurtailed = binaryBlock; 
            extraSkip = -nbus_demandShift; %offsets for binary blocks are: 0, nbus_demandShift, 3nbus_demandShift, 5nbus_demandShift,7nbus_demandShift,...
        end
        b2 = size(binaryBlockCurtailed,2); %2*nbus_demandShift
        
        
        %sum_{t->t+Trsp} -Pdu+ud >= Trsp*Prsp*(bt - bt-1)
        assert(max(colsInserted+(1:sz_Trsp_col)) <= nbus_demandShiftN);%no overflow to another set of variables
        assert(max(colsInserted+extraSkip+(1:b2)) <= nbus_demandShiftN);%no overflow to another set of variables
        %nonzeros for -Pdu: additional row and col offset (nbus_demandShift for each N)
        A(offset1+rowsInsertedRsp+(1:sz_Trsp_row), offset2demandShift+colsInserted+(1:sz_Trsp_col)) = -sumBlockCurtailedTrsp;
        %nonzeros for ud: additional col offset to skip Pdu variables
        A(offset1+rowsInsertedRsp+(1:sz_Trsp_row), offset2demandShift+nbus_demandShiftN+colsInserted+(1:sz_Trsp_col)) = sumBlockCurtailedTrsp;
        %nonzeros for the binary vars: additional col offset to skip Pdu,ud,Pdd,dd variables
        A(offset1+rowsInsertedRsp+(1:sz_Trsp_row), offset2binaryVars+colsInserted+extraSkip+(1:b2)) = -(Prsp.*Trsp).*binaryBlockCurtailed;
        rowsInsertedRsp = rowsInsertedRsp + sz_Trsp_row;
        
        % sum_{t+Trsp->t+Trsp+Trb} Pdd-dd >= Trb*Prb*(bt - bt-1), where t+Trsp+Trb<=N
        if (t+min(Trsp) <= N)
        assert(max(colsInserted+(1:sz_Trb_col)) <= nbus_demandShiftN);%no overflow to another set of variables
            
        %Pdd
        A(offset1+nbus_demandShiftN+rowsInsertedRb+(1:sz_Trb_row), offset2demandShift+2*nbus_demandShiftN+colsInserted+(1:sz_Trb_col)) = sumBlockCurtailedTrb;
        %-dd
        A(offset1+nbus_demandShiftN+rowsInsertedRb+(1:sz_Trb_row), offset2demandShift+3*nbus_demandShiftN+colsInserted+(1:sz_Trb_col)) = -sumBlockCurtailedTrb;
        %binary var from RHS of the equation in the paper
        A(offset1+nbus_demandShiftN+rowsInsertedRb+(1:sz_Trb_row), offset2binaryVars+colsInserted+extraSkip+(1:b2)) = -(Prb.*Trb).*binaryBlockCurtailed;
        rowsInsertedRb = rowsInsertedRb + sz_Trb_row;
        assert(sz_Trb_row > 0);
        end
        
        colsInserted = colsInserted + nbus_demandShift;
    end
    %%%%if trsp(i)!=trsp(j) then rows inserted is not a multiple of nbus_demandShift
    %%%%since we insert different number of trailing Trb constraints for each bus
    rowsInserted = rowsInsertedRsp + rowsInsertedRb;
    
    l(offset1+(1:rowsInserted)) = zeros(rowsInserted,1);
    u(offset1+(1:rowsInserted)) = Inf*ones(rowsInserted,1);
    offset1 = offset1 + rowsInserted;
        
    % complementarity
    % sum_{t->t+Trsp} Pdd-dd <= Trb*Prb*[1 - (bt - bt-1)], where t+Trsp<=N
    % sum_{t+Trsp->t+Trsp+Trb} -Pdu+ud <= Trsp*Prsp*[1- (bt - bt-1)], where t+Trsp+Trb<=N
    %insert the 'sum' blocks into A for each time period t
    rowsInsertedRsp = 0;
    rowsInsertedRb = 0;
    colsInserted = 0;
    for t = 1:N
        %make sure that t+Trsp<=N so that we do not overflow to the next set of variables
        sumBlockCurtailedTrsp = sumBlockTrsp(:,1:nbus_demandShift*min(max(Trsp), N-t+1));
        sz_Trsp_row = size(sumBlockCurtailedTrsp,1);
        sz_Trsp_col = size(sumBlockCurtailedTrsp,2);
        %make sure that t+Trsp+Trb<=N so that we do not overflow to the next set of variables
        sumBlockCurtailedTrb = sumBlockTrb(:,1:nbus_demandShift*min(max(Trsp+Trb), N-t+1));
        sz_Trb_row = size(sumBlockCurtailedTrb,1);
        sz_Trb_col = size(sumBlockCurtailedTrb,2);
        assert(sz_Trsp_row==sz_Trb_row);
        assert(sz_Trsp_row==nbus_demandShift);
        
        %and for t==1 we need to curtail the binary block, since there is no t0
        if (t==1)
            binaryBlockCurtailed = -binaryBlock(:,1:nbus_demandShift);
            extraSkip = 0;
        else
            binaryBlockCurtailed = binaryBlock; 
            extraSkip = -nbus_demandShift; %offsets for binary blocks are 0, nbus_demandShift, 3nbus_demandShift, 5nbus_demandShift,7nbus_demandShift,...
        end
        b2 = size(binaryBlockCurtailed,2);
        
        
        %sum_{t->t+Trsp} Pdd-dd <= Trb*Prb*[1 - (bt - bt-1)], where t+Trsp<=N
        assert(max(colsInserted+(1:sz_Trsp_col)) <= nbus_demandShiftN);%no overflow to another set of variables
        assert(max(colsInserted+extraSkip+(1:b2)) <= nbus_demandShiftN);%no overflow to another set of variables
        %nonzeros for Pdd: additional row and col offset (nbus_demandShift for each N)
        A(offset1+rowsInsertedRsp+(1:sz_Trsp_row), offset2demandShift+2*nbus_demandShiftN+colsInserted+(1:sz_Trsp_col)) = sumBlockCurtailedTrsp;
        %nonzeros for -dd: additional col offset to skip Pdu variables
        A(offset1+rowsInsertedRsp+(1:sz_Trsp_row), offset2demandShift+3*nbus_demandShiftN+colsInserted+(1:sz_Trsp_col)) = -sumBlockCurtailedTrsp;
        %nonzeros for the binary vars: additional col offset to skip Pdu,ud,Pdd,dd variables
        A(offset1+rowsInsertedRsp+(1:sz_Trsp_row), offset2binaryVars+colsInserted+extraSkip+(1:b2)) = (Prb.*Trb).*binaryBlockCurtailed;
        u(offset1+rowsInsertedRsp+(1:sz_Trsp_row)) = Prb.*Trb;
        rowsInsertedRsp = rowsInsertedRsp + sz_Trsp_row;
        
        % sum_{t+Trsp->t+Trsp+Trb} -Pdu+ud <= Trsp*Prsp*[1- (bt - bt-1)], where t+Trsp+Trb<=N
        if (t+min(Trsp) <= N)
          assert(sz_Trb_row ==  size(sumBlockCurtailedTrb,1));
          assert(max(colsInserted+(1:sz_Trb_col)) <= nbus_demandShiftN);%no overflow to another set of variables

          %-Pdu
          A(offset1+nbus_demandShiftN+rowsInsertedRb+(1:sz_Trb_row), offset2demandShift+colsInserted+(1:sz_Trb_col)) = -sumBlockCurtailedTrb;
          %du
          A(offset1+nbus_demandShiftN+rowsInsertedRb+(1:sz_Trb_row), offset2demandShift+nbus_demandShiftN+colsInserted+(1:sz_Trb_col)) = sumBlockCurtailedTrb;
          %binary var from RHS of the equation in the paper
          A(offset1+nbus_demandShiftN+rowsInsertedRb+(1:sz_Trb_row), offset2binaryVars+colsInserted+extraSkip+(1:b2)) = (Prsp.*Trsp).*binaryBlockCurtailed;
          u(offset1+nbus_demandShiftN+rowsInsertedRb+(1:sz_Trb_row)) = Prsp.*Trsp; %Prsp(t+Trsp <= N).*Trsp(t+Trsp <= N)
          rowsInsertedRb = rowsInsertedRb + sz_Trb_row;
        end
        
        colsInserted = colsInserted + nbus_demandShift;
    end
    rowsInserted = rowsInsertedRsp + rowsInsertedRb;
    
    l(offset1+(1:rowsInserted)) = -Inf*ones(rowsInserted,1);
    %u() is set in the for loop above
    offset1 = offset1 + rowsInserted;
    assert(nA == offset1);
end

%% Demand Shift complementarity constraint for the binary var b*(b-1)

if isfield(mpc_storage, 'enableDemandShift') && mpc_storage.enableDemandShift
    mpc_storage.user_constraints.nli = {
     {'binary_complementarity_constr', nbus_demandShiftN, 'binary_complementarity', 'binary_complementarity_hess', {'z'}, {}}
    };
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
mpc_storage.FlexibilityReq              = FlexibilityReq;

mpc_storage.P_storage_max_MW            = P_discharge_max_MW;
mpc_storage.P_storage_min_MW            = P_charge_min_MW;
mpc_storage.E_storage_init_MWh          = E_storage_init_MWh;
mpc_storage.E_storage_max_MWh           = E_storage_max_MWh;

mpc_storage.id_gen_storages_discharge   = ngen_nostorageN+(1:(Nkron*nstorage));
mpc_storage.id_gen_storages_charge      = ngen_nostorageN+Nkron*nstorage+(1:(Nkron*nstorage));
mpc_storage.id_x_gen_storages_discharge = 2*nbusN+mpc_storage.id_gen_storages_discharge;
mpc_storage.id_x_gen_storages_charge    = 2*nbusN+mpc_storage.id_gen_storages_charge;


end
