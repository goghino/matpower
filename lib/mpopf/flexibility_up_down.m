function [f, df] = flexibility_up_down(x, addDischargeChargeFlexibilityPairs, mpc)

ns = mpc.nstorage;
ng = mpc.ngenerators;
N = mpc.horizon;


%Here we are provided with variables x using the internal ordering
%reorder the variables back so that we can compute the offsets
%for the flexibility variables ourselves.
%We assume ordering as defined in add_storage2bNoPeriodic.m
%Pg = [Pgen -- Ped Pec -- ued ded uec dec -- Pdu ud Pdd dd]
Pg = x{1};
Pg = Pg(mpc.order.gen.i2e);

offset = N*ng + 2*N*ns; %skip Pg, Ped, Pec
ued = Pg(offset + (1:N*ns)); offset = offset + N*ns;
ded = Pg(offset + (1:N*ns)); offset = offset + N*ns;
uec = Pg(offset + (1:N*ns)); offset = offset + N*ns;
dec = Pg(offset + (1:N*ns)); offset = offset + N*ns;


%ordering of flexibility up/down complementarity is as follows
%first discharging up-down complementarity
%t=1 [ued11*ded11, ued21*ded21, ued_ns1*ded_ns1]
%t=2 [ued12*ded12, ued22*ded22, ued_ns2*ded_ns2]
%t=N [ued1N*ded1N, ued2N*ded2N, ued_nsN*ded_nsN]
%f = [t1, t2, .. tN]
%second charging up-down complementarity
f = [ued.*ded; uec.*dec];

if (addDischargeChargeFlexibilityPairs)
    %third charging up - discharging down complementarity
    %fourth discharging up - charging down complementarity
    f = [f; uec.*ded; ued.*dec];
end

if nargout > 1
    z1 = sparse(N*ns, N*ng + 2*N*ns); %skip Pg, Ped, Pec
    z2 = sparse(N*ns, N*ns); %skip one of the flexibility ued, ded, uec, dec
    %   *   ued        ded        uec        dec
    df=[z1, diag(ded), diag(ued), z2,       z2;
        z1, z2,       z2,       diag(dec), diag(uec)];
    
    if (addDischargeChargeFlexibilityPairs)
    df=[ df ;
        z1, z2,       diag(uec), diag(ded), z2; 
        z1, diag(dec), z2,       z2,       diag(ued) ];
    end
    
    if mpc.enableDemandShift
        %                         [Pdu Pdd ud dd]
        df = [df,   sparse(2*N*ns,N*4*length(mpc.demandShift.busesID)) ];
    end
    
    %revert back to the internal ordering
    df = df(:, mpc.order.gen.e2i);
end
end