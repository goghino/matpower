function h = discharge_charge_hess(x, lam, mpc)

ns = mpc.nstorage;
ng = mpc.ngenerators;
N = mpc.horizon;

%Here we are provided with variables x using the internal ordering
%reorder the variables back so that we can compute the offsets
%for the flexibility variables ourselves.
%We assume ordering as defined in add_storage2bNoPeriodic.m
%Pg = [Pgen Pd Pc ud dd uc dc]

%ordering of discharge/charge complementarity is as follows
%t=1 [pd11*pc11, pd21*pc21, pd_ns1*pc_ns1]
%t=2 [pd12*pc12, pd22*pc22, pd_ns2*pc_ns2]
%t=N [pd11*pc1N, pd21*pc2N, pd_ns1*pc_nsN]
%f = [t1, t2, .. tN]

L=diag(lam);
z=sparse(N*ns,N*ns);

%% simple complementarity Pd*Pc
if mpc.storageFlexibility == 0
    h_=[z  L;
        L  z];

    %add zero blocks for (Pg), (--, --), (ud, dd, uc, dc)
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

%% revert back to the internal ordering
h = h(mpc.order.gen.e2i, mpc.order.gen.e2i);

end