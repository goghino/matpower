function [P Pinv] = KKTpermute(mpc, nc)

    nbus = size(mpc.bus,1);
    nbrch = size(mpc.branch,1);
    ngen = size(mpc.gen,1);
    
    BUS_TYPE = 2;
    PV = 2;
    busPV = find(mpc.bus(:,BUS_TYPE) == PV); %find PV buses
    nPV = size(busPV,1); %number of PV buses
    
    npart = (2*nbus + ngen + 1) + (2*nbrch) + (2*nbus) + (2*nbrch); %#equations in hess/jac of bus power flow/line limits
    N = nc*npart + ngen-1 + 2*(nc-1)*nPV;
    
    P = zeros(N,1);

    offset = 1;
    for i = 0:nc-1
        %[Va_i Vm_i Qg_i Pg_ref]
        P(offset:offset-1+2*nbus+ngen+1) = (1:2*nbus+ngen+1) + i*(2*nbus+ngen+1);
        offset = offset + 2*nbus+ngen+1;
        %[lam_i]
        P(offset:offset-1+2*nbus) = (1:2*nbus)   + (nc*(2*nbus+ngen+1)) + (ngen-1) + ((nc-1)*nPV) + (nc*2*nbrch) + (i*2*nbus);
        offset = offset + 2*nbus;
        %[mu_i]
        P(offset:offset-1+2*nbrch) = (1:2*nbrch) + (nc*(2*nbus+ngen+1)) + (ngen-1) + ((nc-1)*nPV) + (nc*2*nbrch) + (nc*2*nbus) + (i*2*nbrch);
        offset = offset + 2*nbrch;
        %[Z_i]
        P(offset:offset-1+2*nbrch) = (1:2*nbrch) + (nc*(2*nbus+ngen+1)) + (ngen-1) + (i*2*nbrch);
        offset = offset + 2*nbrch;
    end

    %[Pg] global
    P(offset:offset-1+ngen-1) = (1:ngen-1) + nc*(2*nbus+ngen+1);
    offset = offset + ngen-1;
    
    %[Z_i of A]
    P(offset:offset-1+(nc-1)*nPV) = (1:(nc-1)*nPV) + nc*(2*nbus+ngen+1)+ngen-1+nc*2*nbrch;
    offset = offset + (nc-1)*nPV;
    
    %[A] treated as inequality constr
    P(offset:offset-1+(nc-1)*nPV) = (1:(nc-1)*nPV) + (nc*(2*nbus+ngen+1)) + (ngen-1) + ((nc-1)*nPV) + (nc*2*nbrch) + (nc*2*nbus) + (nc*2*nbrch);
    offset = offset + (nc-1)*nPV;
    
    Pinv = zeros(N,1);
    for i = 1:N
       Pinv(P(i)) = i;
    end
end