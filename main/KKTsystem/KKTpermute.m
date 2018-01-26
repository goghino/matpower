function [P Pinv npart] = KKTpermute(mpc, ns)

    nbus = size(mpc.bus,1);
    nbrch = size(mpc.branch,1);
    ngen = size(mpc.gen,1);
    
    [PVbus_idx, nPVbus_idx] = getXbuses(mpc, 2); %2==PV
    [REFgen_idx, nREFgen_idx] = getREFgens(mpc);
    
    %local and global primal variables
    nprimal_l = (nbus - 1) + length(nPVbus_idx) + ngen + 1;
    nprimal_g = length(PVbus_idx) + (ngen - 1);
    
    %primal variables, slack variables, equality and inequality constraints
    npart = (nprimal_l) + (2*nbrch) + (2*nbus) + (2*nbrch);
    N = ns*npart + nprimal_g;
    
    P = zeros(N,1);

    offset = 1;
    for i = 0:ns-1
        %[Va_i Vm_i(nPV) Qg_i Pg_i(REF)]
        P(offset:offset-1+nprimal_l) = (1:nprimal_l) + i*nprimal_l;
        offset = offset + nprimal_l;
        %[lam_i]
        P(offset:offset-1+2*nbus) = (1:2*nbus)   + (ns*nprimal_l) + (nprimal_g)  + (ns*2*nbrch) + (i*2*nbus);
        offset = offset + 2*nbus;
        %[mu_i]
        P(offset:offset-1+2*nbrch) = (1:2*nbrch) + (ns*nprimal_l) + (nprimal_g)  + (ns*2*nbrch) + (ns*2*nbus) + (i*2*nbrch);
        offset = offset + 2*nbrch;
        %[Z_i]
        P(offset:offset-1+2*nbrch) = (1:2*nbrch) + (ns*nprimal_l) + (nprimal_g) + (i*2*nbrch);
        offset = offset + 2*nbrch;
    end

    %[Vm_PV Pg_nREF] global
    P(offset:offset-1+nprimal_g) = (1:nprimal_g) + ns*nprimal_l;
    offset = offset + nprimal_g;
    
    Pinv = zeros(N,1);
    for i = 1:N
       Pinv(P(i)) = i;
    end
end

function [Xbus_idx, nXbus_idx] = getXbuses(mpc, type)
    %returns indices of buses with specified type and its complement to the full bus set
    BUS_TYPE = 2;
    Xbus_idx = find(mpc.bus(:,BUS_TYPE) == type);
    nXbus_idx = find(mpc.bus(:,BUS_TYPE) ~= type);
end

function [REFgen_idx, nREFgen_idx] = getREFgens(mpc)
    %returns indices of generators connected to reference bus and its
    %complement to the full generator set
    BUS_TYPE = 2;
    REF = 3;
    GEN_BUS = 1;
    REFbus_idx = find(mpc.bus(:,BUS_TYPE) == REF);
    REFgen_idx = find(mpc.gen(:,GEN_BUS) == REFbus_idx); %index of gen connected to ref_bus
    nREFgen_idx = find(mpc.gen(:,GEN_BUS) ~= REFbus_idx); %index of gens not connected to ref_bus
end