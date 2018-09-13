function mpc_storage = add_storage(mpc,nnodes,N,storage_nodes,P_storage_max_MW,P_storage_min_MW, E_storage_max_MWh, E_storage_init_MWh,c_discharge, c_charge,T_timestep_hours)
    
    mpc_storage      = mpc;
    nnodes_total     = size(mpc.bus,1);
    if mod(nnodes_total,nnodes)
        error('wrong number of nodes')
    end
    nstorage         = length(storage_nodes);
    
    %distinguish between replicated mpc and single mpc
    % if provided with mpc replicated N times, than replicate also
    % storages, otherwise keep only single period.
    % Lin. constraints are always expanded to full horizon and replicated.
    if(N == nnodes_total/nnodes)
        Nkron = N;
        ngen_nostorage_total = size(mpc.gen, 1);
        ngen_total = ngen_nostorage_total + N*nstorage*2;
    else
        Nkron = 1;
        nnodes_total = nnodes_total * N;
        ngen_nostorage_total = size(mpc.gen,1) * N;
        ngen_total = ngen_nostorage_total + N*nstorage*2;
    end
    storage_nodesN = repmat(storage_nodes,[Nkron,1])+kron( (0:(Nkron-1))' , ones(nstorage,1)*nnodes ); %proper bus number for different periods N
    
    if find(mpc_storage.bus(storage_nodesN,2)==3)
        error('can not add storage to slack bus')
    end
    mpc_storage.bus(storage_nodesN,2) = 2;


    n_gen_cols = size(mpc_storage.gen,2);
%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
    Qmin = 0;
    Qmax = 0;%1e-10;
    mpc_storage.gen  = [mpc_storage.gen;
                       [storage_nodesN, repmat( [ones(nstorage,1)*[0 0 Qmax Qmin 1 0 1],  P_storage_max_MW, 0*P_storage_max_MW, zeros(nstorage,n_gen_cols-10)], [Nkron,1] )   ];   % discharger
                       [storage_nodesN, repmat( [ones(nstorage,1)*[0 0 Qmax Qmin 1 0 1],0*P_storage_max_MW,   P_storage_min_MW, zeros(nstorage,n_gen_cols-10)], [Nkron,1] )   ] ]; %  charger                   

% %% generator cost data
% %   1   startup shutdown    n   x1  y1  ... xn  yn
% %   2   startup shutdown    n   c(n-1)  ... c0
%  mpc.gencost = ones(ngens,1)*[2 0 0 3 0 200 0];
NCOST = 4;
NC = mpc_storage.gencost(1, NCOST);
gencost = [2 0 0 NC zeros(1, NC)]; %polynomial cost, NCOST=3, c0=c1=c2=0 (no cost)

mpc_storage.gencost = [mpc_storage.gencost;
                       ones(nstorage*Nkron*2,1)*gencost];                    


%% add userconstraints
%% x                 = [theta_bus, Vm_bus, P_gen, Q_gen]
%% P_gen             = [ Pgen_1  ...  Pgen_N, Pdischarge_1 ... Pdischarge_N , Pcharge_1 ... Pcharge_N]
%% Pdischarge_i      = [ Pdischarge_1_i  ...   Pdischarge_nstorage_i]
%% Pcharge_i         = [ Pcharge_1_i  ...   Pcharge_nstorage_i]
ngen_nostorage_total = size(mpc.gen,1);
ngen_total           = size(mpc_storage.gen,1);
A                    = sparse(nstorage+N*nstorage,2*(nnodes_total+ngen_total) );
l                    = zeros(nstorage+N*nstorage,1); 
u                    = l;



M_diag_discharge = sparse(1:nstorage,1:nstorage,1./c_discharge);
M_diag_charge = sparse(1:nstorage,1:nstorage,c_charge);


%% sum(Pdischarge+Pcharge)_1^N = 0 
A(1:nstorage,2*nnodes_total+ngen_nostorage_total+(1:(2*N*nstorage)) ) = ...
    [repmat(M_diag_discharge,[1,N]),...
     repmat(M_diag_charge   ,[1,N]) ];
%     l(1:nstorage)   = 0; 
%     u(1:nstorage)   = 0;
      l(1:nstorage)   = -Inf; 
      u(1:nstorage)   = Inf;

%% 0 < Einit + T*P1 < Emax
%% 0 < Einit + T*(P1+P2) < Emax ...
A(nstorage+(1:(N*nstorage)),2*nnodes_total+ngen_nostorage_total+(1:(2*N*nstorage))  ) = ...
    [ -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_discharge), ...
      -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_charge)    ];
                                                                              
l((nstorage+1):(nstorage+N*nstorage)) = repmat(-E_storage_init_MWh,[N,1]);
u((nstorage+1):(nstorage+N*nstorage)) = repmat(E_storage_max_MWh-E_storage_init_MWh,[N,1]);


%%
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

mpc_storage.id_gen_storages_discharge   = ngen_nostorage_total+(1:(Nkron*nstorage));
mpc_storage.id_gen_storages_charge      = ngen_nostorage_total+Nkron*nstorage+(1:(Nkron*nstorage));
mpc_storage.id_x_gen_storages_discharge = 2*nnodes_total+mpc_storage.id_gen_storages_discharge;
mpc_storage.id_x_gen_storages_charge    = 2*nnodes_total+mpc_storage.id_gen_storages_charge;


end
