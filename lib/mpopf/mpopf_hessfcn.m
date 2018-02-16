function Lxx = mpopf_hessfcn(x, lambda, cost_mult, mpopf_aux)
%MPOPF_HESSFCN  Evaluates Hessian of Lagrangian for AC MPOPF.
%   LXX = MPOPF_HESSFCN(X, LAMBDA, COST_MULT, SCOPF_AUX)
%
%   Hessian evaluation function for AC optimal power flow, suitable
%   for use with MIPS or FMINCON's interior-point algorithm.
%
%   Inputs:
%     X : optimization vector
%     LAMBDA (struct)
%       .eqnonlin : Lagrange multipliers on power balance equations
%       .ineqnonlin : Kuhn-Tucker multipliers on constrained branch flows
%     COST_MULT : (optional) Scale factor to be applied to the cost
%          (default = 1).
%     SCOPF_AUX : auxiliary SCOPF data
%           .om
%           .cont
%           .index
%
%   Outputs:
%     LXX : Hessian of the Lagrangian.
%
%   Examples:
%       Lxx = opf_hessfcn(x, lambda, cost_mult, d);
%
%   See also OPF_COSTFCN, OPF_CONSFCN.

%%----- initialize -----
%% define named indices into data matrices
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

%% default args
if isempty(cost_mult)
    cost_mult = 1;
end

%% unpack data
om = mpopf_aux.om;
mpc = om.get_mpc();
[baseMVA, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.gen, mpc.branch, mpc.gencost);
cp = om.get_cost_params();
[N, Cw, H, dd, rh, kk, mm] = deal(cp.N, cp.Cw, cp.H, cp.dd, ...
                                    cp.rh, cp.kk, cp.mm);
vv = om.get_idx();


%% MPOPF related parameters
profile = mpopf_aux.profile;
Nt = length(profile);

[VAopf, VMopf, PGopf, QGopf] = mpopf_aux.index.getLocalIndicesOPF(mpc);

%% unpack needed parameters
nl = size(branch, 1);       %% number of branches
ng = size(gen, 1);          %% number of dispatchable injections
nxyz = length(x);           %% total number of control vars of all types

%% set default constrained lines
if nargin < 8
    il = (1:nl);            %% all lines have limits by default
end

%% reconstruct V
pcost = gencost(1:ng, :);
if size(gencost, 1) > ng
    qcost = gencost(ng+1:2*ng, :);
else
    qcost = [];
end

%% ----- evaluate d2f -----
d2f = sparse(nxyz, nxyz);
for i = 1:Nt
    % grab Pg & Qg
    idx = mpopf_aux.index.getGlobalIndices(mpc, Nt, i-1);
    Pg = x(idx(PGopf));  %% active generation in p.u.
    Qg = x(idx(QGopf));  %% reactive generation in p.u.


    d2f_dPg2 = sparse(ng, 1);               %% w.r.t. p.u. Pg
    d2f_dQg2 = sparse(ng, 1);               %% w.r.t. p.u. Qg
    ipolp = find(pcost(:, MODEL) == POLYNOMIAL);
    d2f_dPg2(ipolp) = baseMVA^2 * polycost(pcost(ipolp, :), Pg(ipolp)*baseMVA, 2);
    if ~isempty(qcost)          %% Qg is not free
        ipolq = find(qcost(:, MODEL) == POLYNOMIAL);
        d2f_dQg2(ipolq) = baseMVA^2 * polycost(qcost(ipolq, :), Qg(ipolq)*baseMVA, 2);
    end

    %get index of the nominal case variables and use only Pg, Qg
    PGidx = idx(PGopf);
    QGidx = idx(QGopf);
    d2f(PGidx, PGidx) = sparse(1:ng, 1:ng, d2f_dPg2, ng, ng);
    d2f(QGidx, QGidx) = sparse(1:ng, 1:ng, d2f_dQg2, ng, ng);

    %% generalized cost
    if ~isempty(N)
        error('Generalized cost not supported in MPOPF');
        % ... deleted code
    end
end
d2f = d2f * cost_mult;

%% ----- evaluate Hessian of power balance constraints -----
%calls our callbacks hess_miss(x, lam) from mpopf_setup.m
%where all the MPOPF specific indexing is handled
d2G = om.eval_nln_constraint_hess(x, lambda.eqnonlin, 1);

%% ----- evaluate Hessian of flow constraints -----
%calls our callbacks hess_flow(x, lam) from mpopf_setup.m
%where all the MPOPF specific indexing is handled
d2H = om.eval_nln_constraint_hess(x, lambda.ineqnonlin, 0);

Lxx = tril(d2f + d2G + d2H);
