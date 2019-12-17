function [f, df] = flexibility_up_down(x, addDischargeChargeFlexibilityPairs, mpc)

ns = mpc.nstorage;
ng = mpc.ngenerators;
N = mpc.horizon;


%Here we are provided with variables x using the internal ordering
%reorder the variables back so that we can compute the offsets
%for the flexibility variables ourselves.
%We assume ordering as defined in add_storage2bNoPeriodic.m
%Pg = [Pgen Pd Pc ud dd uc dc]
Pg = x{1};
Pg = Pg(mpc.order.gen.i2e);

offset = N*ng + 2*N*ns; %skip Pg, Pd, Pc
ud = Pg(offset + (1:N*ns)); offset = offset + N*ns;
dd = Pg(offset + (1:N*ns)); offset = offset + N*ns;
uc = Pg(offset + (1:N*ns)); offset = offset + N*ns;
dc = Pg(offset + (1:N*ns)); offset = offset + N*ns;


%ordering of flexibility up/down complementarity is as follows
%first discharging up-down complementarity
%t=1 [ud11*dd11, ud21*dd21, ud_ns1*dd_ns1]
%t=2 [ud12*dd12, ud22*dd22, ud_ns2*dd_ns2]
%t=N [ud1N*dd1N, ud2N*dd2N, ud_nsN*dd_nsN]
%f = [t1, t2, .. tN]
%second charging up-down complementarity
f = [ud.*dd; uc.*dc];

if (addDischargeChargeFlexibilityPairs)
    %third charging up - discharging down complementarity
    %fourth discharging up - charging down complementarity
    f = [f; uc.*dd; ud.*dc];
end

if nargout > 1
    z1 = sparse(N*ns, N*ng + 2*N*ns); %skip Pg, Pd, Pc
    z2 = sparse(N*ns, N*ns); %skip one of the flexibility ud, dd, uc, dc
    %   *   ud        dd        uc        dc
    df=[z1, diag(dd), diag(ud), z2,       z2;
        z1, z2,       z2,       diag(dc), diag(uc)];
    
    if (addDischargeChargeFlexibilityPairs)
    df=[ df ;
        z1, z2,       diag(uc), diag(dd), z2; 
        z1, diag(dc), z2,       z2,       diag(ud) ];
    end
    
    %revert back to the internal ordering
    df = df(:, mpc.order.gen.e2i);
end
end