function [P Pinv] = KKTpermute(nbus, nbrch, ngen, nc)

    npart = 2*(2*nbus + 2*nbrch); %hess/jac of bus power flow/line limits
    N = nc*npart + 2*ngen;
    
    P = zeros(N,1);

    offset = 1;
    for i = 0:nc-1
        %[Va_i Vm_i]
        P(offset:offset-1+2*nbus) = (1:2*nbus)   + i*2*nbus;
        offset = offset + 2*nbus;
        %[lam_i]
        P(offset:offset-1+2*nbus) = (1:2*nbus)   + (nc*2*nbus)+2*ngen + nc*2*nbrch + i*2*nbus;
        offset = offset + 2*nbus;
        %[mu_i]
        P(offset:offset-1+2*nbrch) = (1:2*nbrch) + (nc*2*nbus)+2*ngen + nc*2*nbrch + nc*2*nbus + i*2*nbrch;
        offset = offset + 2*nbrch;
        %[Z_i]
        P(offset:offset-1+2*nbrch) = (1:2*nbrch) + (nc*2*nbus)+2*ngen + i*2*nbrch;
        offset = offset + 2*nbrch;
    end

    %[Pg Qg] global
    P(offset:offset-1+2*ngen) = (1:2*ngen) + nc*2*nbus;
    
    Pinv = zeros(N,1);
    for i = 1:N
       Pinv(P(i)) = i;
    end
end