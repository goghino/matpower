function h = flexibility_up_down_hess(x, lam, addDischargeChargeFlexibilityPairs, mpc)

ns = mpc.nstorage;
ng = mpc.ngenerators;
N = mpc.horizon;

%Here we are provided with variables x using the internal ordering
%reorder the variables back so that we can compute the offsets
%for the flexibility variables ourselves.
%We assume ordering as defined in add_storage2bNoPeriodic.m
%Pg = [Pgen -- Ped Pec -- ued ded uec dec -- Pdu ud Pdd dd]

%ordering of flexibility up/down complementarity is as follows
%first discharging complementarity
%t=1 [ued11*ded11, ued21*ded21, ued_ns1*ded_ns1]
%t=2 [ued12*ded12, ued22*ded22, ued_ns2*ded_ns2]
%t=N [ued1N*ded1N, ued2N*ded2N, ued_nsN*ded_nsN]
%f = [t1, t2, .. tN]
lam1 = lam(1:(N*ns)); %for the discharging complementarity
lam2 = lam((N*ns) + (1:(N*ns))); %for the discharging complementarity

if addDischargeChargeFlexibilityPairs
    lam3 = lam((2*N*ns) + (1:(N*ns))); %for the discharging/charging up/down complementarity
    lam4 = lam((3*N*ns) + (1:(N*ns))); %for the discharging/charging up/down complementarity
else
    lam3 = zeros(N*ns,1);
    lam4 = zeros(N*ns,1);
end

L1=diag(lam1);
L2=diag(lam2);
L3=diag(lam3);
L4=diag(lam4);
z=sparse(N*ns,N*ns);

h_=[z   L1 z   L4;
    L1  z  L3  z;
    z   L3 z   L2;
    L4  z  L2  z];

%add zero blocks for Pg, Ped, Pec
n11 = N*ng + 2*N*ns;
n22 = size(h_, 1);
h = [sparse(n11,n11) sparse(n11,n22);
    sparse(n22,n11)  h_];


if mpc.enableDemandShift
    %                         [Pdu Pdd ud dd]
    n11 = size(h,1);
    n12 = N*4*length(mpc.demandShift.busesID);
    h = [h sparse(n11,n12);
        sparse(n12,n11) sparse(n12,n12)];
end

%revert back to the internal ordering
h = h(mpc.order.gen.e2i, mpc.order.gen.e2i);

end