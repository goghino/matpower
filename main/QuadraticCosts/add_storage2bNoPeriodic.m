function mpc_storage = add_storage(mpc,nnodes,storage_nodes,P_storage_max_MW,P_storage_min_MW, E_storage_max_MWh, E_storage_init_MWh,c_discharge, c_charge,T_timestep_hours)
    
    mpc_storage      = mpc;
    nnodes_total     = size(mpc.bus,1);
    if mod(nnodes_total,nnodes)
        error('wrong number of nodes')
    end
    N                = nnodes_total/nnodes;
    nstorage         = length(storage_nodes);
    storage_nodesN   = repmat(storage_nodes,[N,1])+kron( (0:(N-1))' , ones(nstorage,1)*nnodes ); %proper bus number for different periods N
    
    if find(mpc_storage.bus(storage_nodesN,2)==3)
        error('can not add storage to slack bus')
    end
    mpc_storage.bus(storage_nodesN,2) = 2; %BUS_TYPE = PV


    n_gen_cols = size(mpc_storage.gen,2);
%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
    Qmin = 0;
    Qmax = 0;%1e-10;
    mpc_storage.gen  = [mpc_storage.gen;
                       [storage_nodesN, repmat( [ones(nstorage,1)*[0 0 Qmax Qmin 1 0 1],  P_storage_max_MW, 0*P_storage_max_MW, zeros(nstorage,n_gen_cols-10)], [N,1] )   ];   % discharger
                       [storage_nodesN, repmat( [ones(nstorage,1)*[0 0 Qmax Qmin 1 0 1],0*P_storage_max_MW,   P_storage_min_MW, zeros(nstorage,n_gen_cols-10)], [N,1] )   ] ]; %  charger
    %discharging 0 < P < Pmax 
    %charging Pmin < P < 0
    
                   
% %% generator cost data
% %   1   startup shutdown    n   x1  y1  ... xn  yn
% %   2   startup shutdown    n   c(n-1)  ... c0
%  mpc.gencost = ones(ngens,1)*[2 0 0 3 0 200 0];
mpc_storage.gencost = [mpc_storage.gencost;
                       ones(nstorage*N*2,1)*[2 0 0 3 0 0 0]]; %polynomial cost, NCOST=3, c0=c1=c2=0 (no cost)


%% add user constraints
%% x                 = [theta_bus, Vm_bus, P_gen, Q_gen]
%% P_gen             = [ Pgen_1  ...  Pgen_N, Pdischarge_1 ... Pdischarge_N , Pcharge_1 ... Pcharge_N]
%% Pdischarge_i      = [ Pdischarge_1_i  ...   Pdischarge_nstorage_i]
%% Pcharge_i         = [ Pcharge_1_i  ...   Pcharge_nstorage_i]
ngen_nostorage_total = size(mpc.gen,1);
ngen_total           = size(mpc_storage.gen,1);
%A                    = sparse(nstorage+N*nstorage,2*(nnodes_total+ngen_total) );
%l                    = zeros(nstorage+N*nstorage,1); 
%u                    = l;
 
A                    = sparse(N*nstorage,2*(nnodes_total+ngen_total) );
l                    = zeros(N*nstorage,1); 
u                    = l;


M_diag_discharge = sparse(1:nstorage,1:nstorage,1./c_discharge);
M_diag_charge = sparse(1:nstorage,1:nstorage,c_charge);


%% 0 < Einit + T*P1 < Emax
%% 0 < Einit + T*(P1+P2) < Emax ...
% T_timestep_hours = 1
A(1:(N*nstorage),2*nnodes_total+ngen_nostorage_total+(1:(2*N*nstorage))  ) = ...
    [ -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_discharge), ...
      -kron(tril(ones(N)), mpc.baseMVA*T_timestep_hours*M_diag_charge)    ];
                                                                              

l(1:N*nstorage) = repmat(-E_storage_init_MWh,[N,1]);
u(1:N*nstorage) = repmat(E_storage_max_MWh-E_storage_init_MWh,[N,1]);

mpc_storage.A                           = sparse(A);
clear A;
mpc_storage.l                           = l;
mpc_storage.u                           = u;

mpc_storage.c_charge                    = c_charge;
mpc_storage.c_discharge                 = c_discharge;
mpc_storage.nstorage                    = nstorage;
mpc_storage.storage_nodes               = storage_nodes;

mpc_storage.P_storage_max_MW            = P_storage_max_MW;
mpc_storage.P_storage_min_MW            = P_storage_min_MW;
mpc_storage.E_storage_init_MWh          = E_storage_init_MWh;
mpc_storage.E_storage_max_MWh           = E_storage_max_MWh;

mpc_storage.id_gen_storages_discharge   = ngen_nostorage_total+(1:(N*nstorage));
mpc_storage.id_gen_storages_charge      = ngen_nostorage_total+N*nstorage+(1:(N*nstorage));
mpc_storage.id_x_gen_storages_discharge = 2*nnodes_total+mpc_storage.id_gen_storages_discharge;
mpc_storage.id_x_gen_storages_charge    = 2*nnodes_total+mpc_storage.id_gen_storages_charge;


end
