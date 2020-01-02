function h = discharge_charge_hess(x, lam, mpc)

ns = mpc.nstorage;
ng = mpc.ngenerators;
N = mpc.horizon;

%Here we are provided with variables x using the internal ordering
%reorder the variables back so that we can compute the offsets
%for the flexibility variables ourselves.
%We assume ordering as defined in add_storage2bNoPeriodic.m
%Pg = [Pgen -- Ped Pec -- ued ded uec dec -- Pdu ud Pdd dd]

%ordering of discharge/charge complementarity is as follows
%t=1 [ped11*pec11, ped21*pec21, ped_ns1*pec_ns1]
%t=2 [ped12*pec12, ped22*pec22, ped_ns2*pec_ns2]
%t=N [ped11*pec1N, ped21*pec2N, ped_ns1*pec_nsN]
%f = [t1, t2, .. tN]

L=diag(lam);
z=sparse(N*ns,N*ns);

%% simple complementarity Ped*Pec
if mpc.storageFlexibility == 0
    h_=[z  L;
        L  z];

    %add zero blocks for (Pg), (--, --), (ued, ded, uec, dec)
    n11 = N*ng;
    n22 = size(h_, 1);
    n33 = 4*N*ns;
    h = [sparse(n11,n11) sparse(n11,n22) sparse(n11, n33);
        sparse(n22,n11)  h_ sparse(n22, n33);
        sparse(n33, n11) sparse(n33, n22) sparse(n33, n33)];
else
    %% overall complementarity (Pd+ud-dd)(-Pc+uc-dc)
    h_a = [L z L -L z z ];
    h_b = [z -L z z L -L]; 

    h_ = [h_b; -h_a; h_b; -h_b; h_a; -h_a];

    %add zero blocks for (Pg)
    n11 = N*ng;
    n22 = size(h_, 1);
    h = [sparse(n11,n11) sparse(n11,n22);
        sparse(n22,n11)  h_ ];
end

if mpc.enableDemandShift
    %                         [Pdu Pdd ud dd]
    n11 = size(h,1);
    n12 = N*4*length(mpc.demandShift.busesID);
    h = [h sparse(n11,n12);
        sparse(n12,n11) sparse(n12,n12)];
end

%% revert back to the internal ordering
h = h(mpc.order.gen.e2i, mpc.order.gen.e2i);

end