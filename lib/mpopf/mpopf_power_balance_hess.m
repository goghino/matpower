function d2G = mpopf_power_balance_hess(x, lambda, mpc, mpopf_aux, Ybus, mpopt)
%MPOPF_POWER_BALANCE_HESS  Evaluates Hessian of power balance constraints.
%The function iterates for each timeperiod and puts the resulting matrices
%in the global Jacobian/Hessian with proper offsets
%ACCEPT FULL X AND LAMBDA HERE. PROCESS EACH PERIOD INDIVIDUALLY
%AND ASSEMBLE HESSIAN WITH CORRECT ORDERING
%
%   Hessian evaluation function for AC active and reactive power balance
%   constraints.
%
%   Inputs:
%     X : optimization vector 
%     LAMBDA : column vector of Lagrange multipliers on active and reactive
%              power balance constraints
%     MPOPF_AUX:
%           .profile
%           .index.getLocalIndicesOPF
%           .index.getGlobalIndices
%     MPC : MATPOWER case struct
%     YBUS : bus admittance matrix
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     D2G : Hessian of power balance constraints.
%
%   Example:
%       d2G = opf_power_balance_hess(x, lambda, mpc, Ybus, mpopt);
%
%   See also OPF_POWER_BALANCE_FCN.

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
[VAopf, VMopf, PGopf, QGopf] = mpopf_aux.index.getLocalIndicesOPF(mpc);

profile = mpopf_aux.profile;
Nt = length(profile);

%nominal grid load
mpcBusLoadNominal = mpc.bus(:,3:4);

%% problem dimensions
nb = size(mpc.bus,1);            %% number of buses
ng = size(mpc.gen,1);            %% number of dispatchable injections
nx = length(x); 

d2G = sparse(nx,nx);

for i = 1:Nt
    %% update mpc by load scaling profile for this period by scaling PD, QD
    load_factor = profile(i);
    mpc.bus(:,3:4) = mpc.bus(:,3:4) * load_factor;
    
    %% manage variable and constraints indexing
    
    % extract x indexes for the current time period
    idx = mpopf_aux.index.getGlobalIndices(mpc, Nt, i-1);
    
    % compute indexes of re/im parts of the mismatch of current period
    idxRe = (i-1)*nb + (1:nb); 
    idxIm = Nt*nb + (i-1)*nb + (1:nb);

    %% unpack data
    Va = x(idx(VAopf));
    Vm = x(idx(VMopf));
    Pg = x(idx(PGopf));
    Qg = x(idx(QGopf));

    %% reconstruct V
    V = Vm .* exp(1j * Va);

    %%----- evaluate Hessian of power balance constraints -----
    lamP = lambda(idxRe);
    lamQ = lambda(idxIm);
    
    % real part
    [Gpaa, Gpav, Gpva, Gpvv] = d2Sbus_dV2(Ybus, V, lamP);

    % adjust for voltage dependent loads (constant impedance part of ZIP loads)
    diaglam = sparse(1:nb, 1:nb, lamP, nb, nb);
    Sd = makeSdzip(mpc.baseMVA, mpc.bus, mpopt);
    diagSdz = sparse(1:nb, 1:nb, Sd.z, nb, nb);
    Gpvv = Gpvv + 2 * diaglam * diagSdz;

    % imaginary part 
    [Gqaa, Gqav, Gqva, Gqvv] = d2Sbus_dV2(Ybus, V, lamQ);

    %% construct Hessian for period i with proper offsets
    d2G(idx([VAopf VMopf]), idx([VAopf VMopf])) = real([Gpaa Gpav; Gpva Gpvv]) + imag([Gqaa Gqav; Gqva Gqvv]);
end

