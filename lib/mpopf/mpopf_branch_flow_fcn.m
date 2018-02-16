function [h, dh] = mpopf_branch_flow_fcn(x, mpc, mpopf_aux, Yf, Yt, il, mpopt)
%MPOPF_BRANCH_FLOW_FCN  Evaluates AC branch flow constraints and Jacobian.
%The function iterates for each timeperiod and puts the resulting matrices
%in the global Jacobian/Hessian with proper offsets
%ACCEPT FULL X AND LAMBDA HERE. PROCESS EACH PERIOD INDIVIDUALLY
%AND ASSEMBLE HESSIAN WITH CORRECT ORDERING
%
%   Active power balance equality constraints for AC optimal power flow.
%   Computes constraint vectors and their gradients.
%
%   Inputs:
%     X : optimization vector
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
%     H  : vector of inequality constraint values (flow limits)
%          where the flow can be apparent power, real power, or
%          current, depending on the value of opf.flow_lim in MPOPT
%          (only for constrained lines), normally expressed as
%          (limit^2 - flow^2), except when opf.flow_lim == 'P',
%          in which case it is simply (limit - flow).
%     DH : (optional) inequality constraint gradients, column j is
%          gradient of H(j)
%
%   Examples:
%       h = opf_branch_flow_fcn(x, mpc, Yf, Yt, il, mpopt);
%       [h, dh] = opf_branch_flow_fcn(x, mpc, Yf, Yt, il, mpopt);
%
%   See also OPF_BRANCH_FLOW_HESS.

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

h = zeros(Nt*2*nl2,1);
if nargout > 1
    dh = sparse(Nt*2*nl2, nx);
end

for i = 1:Nt    
    %% manage variable and constraints indexing
    
    % extract x indexes for the current time period
    idx = mpopf_aux.index.getGlobalIndices(mpc, Nt, i-1);
    
    % compute indexes of from/to parts of the current period
    idxF = (i-1)*nl2 + (1:nl2); 
    idxT = Nt*nl2 + (i-1)*nl2 + (1:nl2);
    
    %% unpack data
    lim_type = upper(mpopt.opf.flow_lim(1));
    branch = mpc.branch;
    Va = x(idx(VAopf));
    Vm = x(idx(VMopf));

    %% ----- evaluate constraints -----
    if nl2 > 0
        %% reconstruct V
        V = Vm .* exp(1j * Va);

        flow_max = branch(il, RATE_A) / mpc.baseMVA;
        if lim_type ~= 'P'      %% typically use square of flow
            flow_max = flow_max.^2;
        end
        if lim_type == 'I'      %% current magnitude limit, |I|
            If = Yf * V;
            It = Yt * V;
            flow = [ If .* conj(If) - flow_max;    %% branch current limits (from bus)
                  It .* conj(It) - flow_max ];  %% branch current limits (to bus)
        else
            %% compute branch power flows
            Sf = V(branch(il, F_BUS)) .* conj(Yf * V);  %% complex power injected at "from" bus (p.u.)
            St = V(branch(il, T_BUS)) .* conj(Yt * V);  %% complex power injected at "to" bus (p.u.)
            if lim_type == '2'                      %% active power limit, P squared (Pan Wei)
                flow = [ real(Sf).^2 - flow_max;       %% branch real power limits (from bus)
                      real(St).^2 - flow_max ];     %% branch real power limits (to bus)
            elseif lim_type == 'P'                  %% active power limit, P
                flow = [ real(Sf) - flow_max;          %% branch real power limits (from bus)
                      real(St) - flow_max ];        %% branch real power limits (to bus
            else                                    %% apparent power limit, |S|
                flow = [ Sf .* conj(Sf) - flow_max;    %% branch apparent power limits (from bus)
                      St .* conj(St) - flow_max ];  %% branch apparent power limits (to bus)
            end
        end

        % put the constraints to proper positions
        h(idxF) = flow(1:nl2);
        h(idxT) = flow((nl2+1):end);
    else
        h = zeros(0,1);
    end

    %%----- evaluate partials of constraints -----
    if nargout > 1
        if nl2 > 0
            %% compute partials of Flows w.r.t. V
            if lim_type == 'I'                      %% current
                [dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft] = dIbr_dV(branch(il,:), Yf, Yt, V);
            else                                    %% power
                [dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft] = dSbr_dV(branch(il,:), Yf, Yt, V);
            end
            if lim_type == 'P' || lim_type == '2'   %% real part of flow (active power)
                dFf_dVa = real(dFf_dVa);
                dFf_dVm = real(dFf_dVm);
                dFt_dVa = real(dFt_dVa);
                dFt_dVm = real(dFt_dVm);
                Ff = real(Ff);
                Ft = real(Ft);
            end

            if lim_type == 'P'
                %% active power
                [df_dVa, df_dVm, dt_dVa, dt_dVm] = deal(dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm);
            else
                %% squared magnitude of flow (of complex power or current, or real power)
                [df_dVa, df_dVm, dt_dVa, dt_dVm] = ...
                      dAbr_dV(dFf_dVa, dFf_dVm, dFt_dVa, dFt_dVm, Ff, Ft);
            end
            %% construct Jacobian of "from" branch flow ineq constraints
            dh(idxF, idx([VAopf VMopf])) = [df_dVa df_dVm]; %% "from" flow limit
            dh(idxT, idx([VAopf VMopf])) = [dt_dVa dt_dVm]; %% "to" flow limit

        else
            dh = sparse(0, 2*nb);
        end
    end
end