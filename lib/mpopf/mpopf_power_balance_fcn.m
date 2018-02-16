function [g, dg] = mpopf_power_balance_fcn(x, mpc, mpopf_aux, Ybus, mpopt)
%OPF_POWER_BALANCE_FCN  Evaluates AC power balance constraints and their gradients.
%The function iterates for each timeperiod and puts the resulting matrices
%in the global Jacobian/Hessian with proper offsets
%ACCEPT FULL X AND LAMBDA HERE. PROCESS EACH PERIOD INDIVIDUALLY
%AND ASSEMBLE HESSIAN WITH CORRECT ORDERING
%
%   Computes the active or reactive power balance equality constraints for
%   AC optimal power flow. Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
%     MPC : MATPOWER case struct
%     MPOPF_AUX:
%           .profile
%           .index.getLocalIndicesOPF
%           .index.getGlobalIndices
%     YBUS : bus admittance matrix
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     G  : vector of equality constraint values (active/reactive power balances)
%     DG : (optional) equality constraint gradients
%
%   Examples:
%       g = opf_power_balance_fcn(x, mpc, Ybus, mpopt);
%       [g, dg] = opf_power_balance_fcn(x, mpc, Ybus, mpopt);
%
%   See also OPF_POWER_BALANCE_HESS

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

[VAopf, VMopf, PGopf, QGopf] = mpopf_aux.index.getLocalIndicesOPF(mpc);

profile = mpopf_aux.profile;
Nt = length(profile);

%nominal grid load
mpcBusLoadNominal = mpc.bus(:,3:4);

%% problem dimensions
nb = size(mpc.bus,1);            %% number of buses
ng = size(mpc.gen,1);            %% number of dispatchable injections
nx = length(x);                  %% number of optimization variables

g = zeros(Nt*2*nb,1);
if nargout > 1
    dg = sparse(Nt*2*nb, nx);
end

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
    [baseMVA, bus, gen] = deal(mpc.baseMVA, mpc.bus, mpc.gen);
    Va = x(idx(VAopf));
    Vm = x(idx(VMopf));
    Pg = x(idx(PGopf));
    Qg = x(idx(QGopf));

    %% ----- evaluate constraints -----

    %% put Pg, Qg back in gen
    gen(:, PG) = Pg * baseMVA;  %% active generation in MW
    gen(:, QG) = Qg * baseMVA;  %% reactive generation in MVAr

    %% reconstruct V
    V = Vm .* exp(1j * Va);

    %% rebuild Sbus
    Sbus = makeSbus(baseMVA, bus, gen, mpopt, Vm);  %% net injected power in p.u.

    %% evaluate complex power balance mismatches
    mis = V .* conj(Ybus * V) - Sbus;

    %% assemble active and reactive power balance constraints    
    g(idxRe,1) = real(mis);    %% active power mismatch
    g(idxIm,1) = imag(mis);    %% reactive power mismatch

    %%----- evaluate constraint gradients -----
    if nargout > 1
        %% compute partials of injected bus powers
        [dSbus_dVm, dSbus_dVa] = dSbus_dV(Ybus, V);         %% w.r.t. V
        neg_Cg = sparse(gen(:, GEN_BUS), 1:ng, -1, nb, ng); %% Pbus w.r.t. Pg

        %% adjust for voltage dependent loads
        [dummy, neg_dSd_dVm] = makeSbus(baseMVA, bus, gen, mpopt, Vm);
        dSbus_dVm = dSbus_dVm - neg_dSd_dVm;

        %% create proper offsets in Jacobian
        dg(idxRe, idx([VAopf VMopf])) = real([dSbus_dVa dSbus_dVm]);  %% P mismatch w.r.t Va, Vm
        dg(idxRe, idx(PGopf)) = neg_Cg ;                              %% P mismatch w.r.t Pg
        
        dg(idxIm, idx([VAopf VMopf])) = imag([dSbus_dVa dSbus_dVm]);  %% Q mismatch w.r.t Va, Vm
        dg(idxIm, idx(QGopf)) = neg_Cg ;                              %% Q mismatch w.r.t Qg
    end
    
    %% return the load scaling for the next iteration to nominal state
    mpc.bus(:,3:4) = mpcBusLoadNominal;
end