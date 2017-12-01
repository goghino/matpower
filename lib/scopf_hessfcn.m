function Lxx = scopf_hessfcn(x, lambda, cost_mult, scopf_aux)
%SCOPF_HESSFCN  Evaluates Hessian of Lagrangian for AC SCOPF.
%   LXX = SCOPF_HESSFCN(X, LAMBDA, COST_MULT, SCOPF_AUX)
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
om = scopf_aux.om;
mpc = om.get_mpc();
[baseMVA, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.gen, mpc.branch, mpc.gencost);
cp = om.get_cost_params();
[N, Cw, H, dd, rh, kk, mm] = deal(cp.N, cp.Cw, cp.H, cp.dd, ...
                                    cp.rh, cp.kk, cp.mm);
vv = om.get_idx();


%% SCOPF related parameters
cont = scopf_aux.cont;
NS = length(cont);

idx0 = scopf_aux.index.getGlobalIndices(mpc, NS, 0);

% get indices of REF gen and PV bus
[REFgen_idx, nREFgen_idx] = scopf_aux.index.getREFgens(mpc);
[PVbus_idx, nPVbus_idx] = scopf_aux.index.getXbuses(mpc,2);%2==PV

[VAscopf, VMscopf, PGscopf, QGscopf] = scopf_aux.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = scopf_aux.index.getLocalIndicesOPF(mpc);

%% unpack needed parameters
nl = size(branch, 1);       %% number of branches
ng = size(gen, 1);          %% number of dispatchable injections
nxyz = length(x);           %% total number of control vars of all types

%% set default constrained lines
if nargin < 8
    il = (1:nl);            %% all lines have limits by default
end

%% grab Pg & Qg
Pg = x(idx0(PGscopf));  %% active generation in p.u.
Qg = x(idx0(QGscopf));  %% reactive generation in p.u.

%% reconstruct V
pcost = gencost(1:ng, :);
if size(gencost, 1) > ng
    qcost = gencost(ng+1:2*ng, :);
else
    qcost = [];
end

%% ----- evaluate d2f -----
d2f_dPg2 = sparse(ng, 1);               %% w.r.t. p.u. Pg
d2f_dQg2 = sparse(ng, 1);               %% w.r.t. p.u. Qg
ipolp = find(pcost(:, MODEL) == POLYNOMIAL);
d2f_dPg2(ipolp) = baseMVA^2 * polycost(pcost(ipolp, :), Pg(ipolp)*baseMVA, 2);
if ~isempty(qcost)          %% Qg is not free
    ipolq = find(qcost(:, MODEL) == POLYNOMIAL);
    d2f_dQg2(ipolq) = baseMVA^2 * polycost(qcost(ipolq, :), Qg(ipolq)*baseMVA, 2);
end
i = [idx0(PGscopf) idx0(QGscopf)]';
d2f = sparse(i, i, [d2f_dPg2; d2f_dQg2], nxyz, nxyz);

%% generalized cost
if ~isempty(N)
    nw = size(N, 1);
    r = N * x - rh;                 %% Nx - rhat
    iLT = find(r < -kk);            %% below dead zone
    iEQ = find(r == 0 & kk == 0);   %% dead zone doesn't exist
    iGT = find(r > kk);             %% above dead zone
    iND = [iLT; iEQ; iGT];          %% rows that are Not in the Dead region
    iL = find(dd == 1);             %% rows using linear function
    iQ = find(dd == 2);             %% rows using quadratic function
    LL = sparse(iL, iL, 1, nw, nw);
    QQ = sparse(iQ, iQ, 1, nw, nw);
    kbar = sparse(iND, iND, [   ones(length(iLT), 1);
                                zeros(length(iEQ), 1);
                                -ones(length(iGT), 1)], nw, nw) * kk;
    rr = r + kbar;                  %% apply non-dead zone shift
    M = sparse(iND, iND, mm(iND), nw, nw);  %% dead zone or scale
    diagrr = sparse(1:nw, 1:nw, rr, nw, nw);
    
    %% linear rows multiplied by rr(i), quadratic rows by rr(i)^2
    w = M * (LL + QQ * diagrr) * rr;
    HwC = H * w + Cw;
    AA = N' * M * (LL + 2 * QQ * diagrr);
    d2f = d2f + AA * H * AA' + 2 * N' * M * QQ * sparse(1:nw, 1:nw, HwC, nw, nw) * N;
end
d2f = d2f * cost_mult;

%%----- evaluate Hessian of power balance constraints -----
d2G = om.eval_nln_constraint_hess(x, lambda.eqnonlin, 1);

%%----- evaluate Hessian of flow constraints -----
d2H = om.eval_nln_constraint_hess(x, lambda.ineqnonlin, 0);

Lxx = tril(d2f + d2G + d2H);
