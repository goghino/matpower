function [Ybus, Yf, Yt] = updateYbus(branch, Ybus, Yf, Yt, cont)
%MAKEYBUS   Builds the bus admittance matrix and branch admittance matrices.
%   [YBUS, YF, YT] = MAKEYBUS(MPC)
%   [YBUS, YF, YT] = MAKEYBUS(BRANCH, Ybus, Yf, Yt, C)
%   
%   Updates admittance matrices for given branch contingency C


%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;


%% build connection matrices
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
 
%% modifications of admittance matrices for contingency
if ( cont > 0)

%update Ybus    
Ybus(f(cont), f(cont)) = Ybus(f(cont), f(cont)) - Yf(cont, f(cont));
Ybus(f(cont), t(cont)) = Ybus(f(cont), t(cont)) - Yf(cont, t(cont));

Ybus(t(cont), f(cont)) = Ybus(t(cont), f(cont)) - Yt(cont, f(cont));
Ybus(t(cont), t(cont)) = Ybus(t(cont), t(cont)) - Yt(cont, t(cont));

%update Yf, Yt
Yt(cont, :) = 0;
Yf(cont, :) = 0;
end
