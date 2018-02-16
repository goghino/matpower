function d2H = mpopf_branch_flow_hess(x, lambda, mpc, mpopf_aux, Yf, Yt, il, mpopt)
%OPF_BRANCH_FLOW_HESS  Evaluates Hessian of branch flow constraints.
%The function iterates for each timeperiod and puts the resulting matrices
%in the global Jacobian/Hessian with proper offsets
%ACCEPT FULL X AND LAMBDA HERE. PROCESS EACH PERIOD INDIVIDUALLY
%AND ASSEMBLE HESSIAN WITH CORRECT ORDERING
%
%   Hessian evaluation function for AC branch flow constraints.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA : column vector of Kuhn-Tucker multipliers on constrained
%              branch flows
%     MPC : MATPOWER case struct
%     MPOPF_AUX:
%           .profile
%           .index.getLocalIndicesOPF
%           .index.getGlobalIndices
%     YF : admittance matrix for "from" end of constrained branches
%     YT : admittance matrix for "to" end of constrained branches
%     IL : vector of branch indices corresponding to branches with
%          flow limits (all others are assumed to be unconstrained).
%          YF and YT contain only the rows corresponding to IL.
%     MPOPT : MATPOWER options struct
%
%   Outputs:
%     D2H : Hessian of AC branch flow constraints.
%
%   Example:
%       d2H = opf_branch_flow_hess(x, lambda, mpc, Yf, Yt, il, mpopt);
%
%   See also OPF_BRANCH_FLOW_FCN.

%   MATPOWER
%   Copyright (c) 1996-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialize -----
%% define named indices into data matrices
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

[VAopf, VMopf, PGopf, QGopf] = mpopf_aux.index.getLocalIndicesOPF(mpc);

profile = mpopf_aux.profile;
Nt = length(profile);

%% problem dimensions
nb = size(mpc.bus,1);            %% number of buses
nl2 = size(mpc.branch,1);        %% number of constrained lines
nx = length(x);                  %% number of optimization variables

d2H = sparse(nx,nx);

for i = 1:Nt
    %% manage variable and constraints indexing
    
    % extract x indexes for the current time period
    idx = mpopf_aux.index.getGlobalIndices(mpc, Nt, i-1);
    
    % compute indexes of from/to parts of the current period
    idxF = (i-1)*nl2 + (1:nl2); 
    idxT = Nt*nl2 + (i-1)*nl2 + (1:nl2);
    
    %% unpack data
    lim_type = upper(mpopt.opf.flow_lim(1));
    Va = x(idx(VAopf));
    Vm = x(idx(VMopf));

    %% reconstruct V
    V = Vm .* exp(1j * Va);

    %%----- evaluate Hessian of flow constraints -----
    %% keep dimensions of empty matrices/vectors compatible
    %% (required to avoid problems when using Knitro
    %%  on cases with all lines unconstrained)

    muF = lambda(idxF);
    muT = lambda(idxT);

    if lim_type == 'I'          %% square of current
        [dIf_dVa, dIf_dVm, dIt_dVa, dIt_dVm, If, It] = dIbr_dV(mpc.branch(il,:), Yf, Yt, V);
        [Hfaa, Hfav, Hfva, Hfvv] = d2AIbr_dV2(dIf_dVa, dIf_dVm, If, Yf, V, muF);
        [Htaa, Htav, Htva, Htvv] = d2AIbr_dV2(dIt_dVa, dIt_dVm, It, Yt, V, muT);
    else
        f = mpc.branch(il, F_BUS);    %% list of "from" buses
        t = mpc.branch(il, T_BUS);    %% list of "to" buses
        Cf = sparse(1:nl2, f, ones(nl2, 1), nl2, nb);   %% connection matrix for line & from buses
        Ct = sparse(1:nl2, t, ones(nl2, 1), nl2, nb);   %% connection matrix for line & to buses
        [dSf_dVa, dSf_dVm, dSt_dVa, dSt_dVm, Sf, St] = dSbr_dV(mpc.branch(il,:), Yf, Yt, V);
        if lim_type == '2'        %% square of real power
            [Hfaa, Hfav, Hfva, Hfvv] = d2ASbr_dV2(real(dSf_dVa), real(dSf_dVm), real(Sf), Cf, Yf, V, muF);
            [Htaa, Htav, Htva, Htvv] = d2ASbr_dV2(real(dSt_dVa), real(dSt_dVm), real(St), Ct, Yt, V, muT);
        elseif lim_type == 'P'    %% real power                                 
            [Hfaa, Hfav, Hfva, Hfvv] = d2Sbr_dV2(Cf, Yf, V, muF);
            [Htaa, Htav, Htva, Htvv] = d2Sbr_dV2(Ct, Yt, V, muT);
            [Hfaa, Hfav, Hfva, Hfvv] = deal(real(Hfaa), real(Hfav), real(Hfva), real(Hfvv));
            [Htaa, Htav, Htva, Htvv] = deal(real(Htaa), real(Htav), real(Htva), real(Htvv));
        else                      %% square of apparent power
            [Hfaa, Hfav, Hfva, Hfvv] = d2ASbr_dV2(dSf_dVa, dSf_dVm, Sf, Cf, Yf, V, muF);
            [Htaa, Htav, Htva, Htvv] = d2ASbr_dV2(dSt_dVa, dSt_dVm, St, Ct, Yt, V, muT);
        end
    end

    %% construct Hessian for period i with proper offsets
    d2H(idx([VAopf VMopf]), idx([VAopf VMopf])) = [Hfaa Hfav; Hfva Hfvv] + [Htaa Htav; Htva Htvv];
end