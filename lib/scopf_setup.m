function om = scopf_setup(mpc, scopf_aux, mpopt)
%SCOPF  Constructs an OPF model object from a MATPOWER case struct.
%   OM = SCOPF_SETUP(MPC, MPOPT)
%
%   Assumes that MPC is a MATPOWER case struct with internal indexing,
%   all equipment in-service, etc.
%
%   See also OPF, EXT2INT, OPF_EXECUTE.

%% options
dc  = strcmp(upper(mpopt.model), 'DC');
if dc==1
   error('DC not supported in SCOPF');
end

use_vg = mpopt.opf.use_vg;
if use_vg ~= 0
   error('opf.use_vg not supported in SCOPF');
end

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

% get indices of REF gen and of REF/PV buses
[REFgen_idx, nREFgen_idx] = scopf_aux.index.getREFgens(mpc);
[REFbus_idx,nREFbus_idx] = scopf_aux.index.getXbuses(mpc,3);%3==REF
[PVbus_idx, nPVbus_idx] = scopf_aux.index.getXbuses(mpc,2);%2==PV

% indices of local OPF solution vector x = [VA VM PG QG]
[VAscopf, VMscopf, PGscopf, QGscopf] = scopf_aux.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = scopf_aux.index.getLocalIndicesOPF(mpc);

%% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections
ns   = length(scopf_aux.cont);

if isfield(mpc, 'A')
  error('User constraints not supported in SCOPF');
end

if isfield(mpc, 'N')
    error('general cost vars not supported in SCOPF');
end

%% create (read-only) copies of individual fields for convenience
[baseMVA, bus, gen, branch, gencost, Au, lbu, ubu, mpopt, ...
    N, fparm, H, Cw, z0, zl, zu, userfcn] = opf_args(mpc, mpopt);

%% warn if there is more than one reference bus
refs = find(bus(:, BUS_TYPE) == REF);
if length(refs) > 1 && mpopt.verbose > 0
  errstr = ['\nopf_setup: Warning: Multiple reference buses.\n', ...
              '           For a system with islands, a reference bus in each island\n', ...
              '           may help convergence, but in a fully connected system such\n', ...
              '           a situation is probably not reasonable.\n\n' ];
  fprintf(errstr);
end

%% set up initial variables and bounds
Va   = bus(:, VA) * (pi/180);
Vm   = bus(:, VM);
Pg   = gen(:, PG) / baseMVA;
Qg   = gen(:, QG) / baseMVA;
Pmin = gen(:, PMIN) / baseMVA;
Pmax = gen(:, PMAX) / baseMVA;
Qmin = gen(:, QMIN) / baseMVA;
Qmax = gen(:, QMAX) / baseMVA;


%% voltage angle reference constraints
Vau = Inf(nb, 1);
Val = -Vau;
Vau(refs) = Va(refs);
Val(refs) = Va(refs);

if any(gencost(:, MODEL) ~= POLYNOMIAL & gencost(:, MODEL) ~= PW_LINEAR)
    error('opf_setup: some generator cost rows have invalid MODEL value');
end

%% find branches with flow limits
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
nl2 = length(il);         %% number of constrained lines

%% construct OPF model object
om = opf_model(mpc);

for i = 1:ns
    om.add_var(strcat('Va',int2str(i)), length(Va), Va, Val, Vau); %adds OPF variables for voltage angles and set initial values
    om.add_var(strcat('Vm',int2str(i)), length(nPVbus_idx), Vm(nPVbus_idx), bus(nPVbus_idx, VMIN), bus(nPVbus_idx, VMAX)); %voltage magnitude
    om.add_var(strcat('Qg',int2str(i)), length(Qg), Qg, Qmin, Qmax); %generator imaginary power
    om.add_var(strcat('Pg',int2str(i)), length(REFgen_idx), Pg(REFgen_idx), Pmin(REFgen_idx), Pmax(REFgen_idx)); %generator real power
end

om.add_var('Vmg', length(PVbus_idx), Vm(PVbus_idx), bus(PVbus_idx, VMIN), bus(PVbus_idx, VMAX)); %voltage magnitude
om.add_var('Pgg', length(nREFgen_idx), Pg(nREFgen_idx), Pmin(nREFgen_idx), Pmax(nREFgen_idx)); %generator real power

for i = 1:ns
    %make Ybus reflect the i-th contingency
    cont = scopf_aux.cont(i);
    [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
    
    %extract local OPF vector for i-th contingency
    idx = scopf_aux.index.getGlobalIndices(mpc, ns, i-1);
    
    % Functions below extract local OPF vector from SCOPF x
    % and format is as a {Va, Vm, Pg, Qg} cell array for PF equations
    % and {Va, Vm} cell array for branch power limits.
    % There is also different Ybus, Yf, Yt for each contingency.
    % Matlab somehow remembers it when creating the handles (so we don't
    % have to create an array of handles with references to elements of
    % admittance matrix array)
    fcn_mis = @(x)opf_power_balance_fcn({x(idx(VAscopf)), x(idx(VMscopf)), x(idx(PGscopf)), x(idx(QGscopf))}, mpc, Ybus, mpopt);
    hess_mis = @(x, lam)scopf_power_balance_hess(x, lam, i-1, mpc, scopf_aux, Ybus, mpopt);
    
    fcn_flow = @(x)opf_branch_flow_fcn({x(idx(VAscopf)), x(idx(VMscopf))}, mpc, Yf(il, :), Yt(il, :), il, mpopt);
    hess_flow = @(x, lam)scopf_branch_flow_hess(x, lam, i-1, mpc, scopf_aux, Yf, Yt, il, mpopt);
    
    % mismatch equations (eq. constr.)
    om.add_nln_constraint({strcat('Pmis',int2str(i)), strcat('Qmis',int2str(i))}, [nb;nb], 1, fcn_mis, hess_mis);
    % branch power flow limits (ineq. constr.)
    om.add_nln_constraint({strcat('Sf',int2str(i)), strcat('St',int2str(i))}, [nl2;nl2], 0, fcn_flow, hess_flow);
end

%% find/prepare polynomial generator costs
cpg = [];
cqg = [];
[pcost qcost] = pqcost(mpc.gencost, ng);
ip0 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 1);   %% constant
ip1 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 2);   %% linear
ip2 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 3);   %% quadratic
ip3 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) > 3);    %% cubic or greater
if ~isempty(ip2) || ~isempty(ip1) || ~isempty(ip0)
    kpg = zeros(ng, 1);
    cpg = zeros(ng, 1);
    if ~isempty(ip2)
        Qpg = zeros(ng, 1);
        Qpg(ip2) = 2 * pcost(ip2, COST) * baseMVA^2;
        cpg(ip2) = cpg(ip2) + pcost(ip2, COST+1) * baseMVA;
        kpg(ip2) = kpg(ip2) + pcost(ip2, COST+2);
    else
        Qpg = [];   %% no quadratic terms
    end
    cpg(ip1) = cpg(ip1) + pcost(ip1, COST) * baseMVA;
    kpg(ip1) = kpg(ip1) + pcost(ip1, COST+1);
    kpg(ip0) = kpg(ip0) + pcost(ip0, COST);
end
if ~isempty(qcost)
    iq0 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) == 1);   %% constant
    iq1 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) == 2);   %% linear
    iq2 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) == 3);   %% quadratic
    iq3 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) > 3);    %% cubic or greater
    if ~isempty(iq2) || ~isempty(iq1) || ~isempty(iq0)
        kqg = zeros(ng, 1);
        cqg = zeros(ng, 1);
        if ~isempty(iq2)
            Qqg = zeros(ng, 1);
            Qqg(iq2) = 2 * qcost(iq2, COST) * baseMVA^2;
            cqg(iq2) = cqg(iq2) + qcost(iq2, COST+1) * baseMVA;
            kqg(iq2) = kqg(iq2) + qcost(iq2, COST+2);
        else
            Qqg = [];   %% no quadratic terms
        end
        cqg(iq1) = cqg(iq1) + qcost(iq1, COST) * baseMVA;
        kqg(iq1) = kqg(iq1) + qcost(iq1, COST+1);
        kqg(iq0) = kqg(iq0) + qcost(iq0, COST);
    end
end

%% quadratic/linear generator costs
if ~isempty(cpg)
  %consider REF generator first, then the rest - permute cost coefficients
  %accordingly, only nominal case Pg_ref is considered
  om.add_quad_cost('polPg', Qpg([REFgen_idx; nREFgen_idx]), cpg([REFgen_idx; nREFgen_idx]), kpg([REFgen_idx; nREFgen_idx]), {'Pg1', 'Pgg'});
end
if ~isempty(cqg)
  %only nominal scenario Qg is considered
  om.add_quad_cost('polQg', Qqg, cqg, kqg, {'Qg1'});
end

%% higher order polynomial generator costs
if ~isempty(ip3)
  om.add_nln_cost('polPg', 1, cost_Pg([REFgen_idx; nREFgen_idx]), {'Pg1', 'Pgg'});
end
if ~isempty(qcost) && ~isempty(iq3)
  om.add_nln_cost('polQg', 1, cost_Qg, {'Qg1'});
end

%calls opf_power_balance_hess(x, lam).
%Before passing x to evaluate hessian it is permuted to OPF ordering
%and each set of constraints for given contingency extract different part 
%of x containing local variables for the contingency. Returned hessian
%is in OPF ordering, therefore we perform appropriate permutation
%with indices corresponding to i-th contingency
function H_scopf = scopf_power_balance_hess(x, lam, cont_i, mpc, scopf_aux, Ybus, mpopt)
n = length(x);
ns = length(scopf_aux.cont);

% get indices of REF gen and of REF/PV buses
[REFgen_idx, nREFgen_idx] = scopf_aux.index.getREFgens(mpc);
[PVbus_idx, nPVbus_idx] = scopf_aux.index.getXbuses(mpc,2);%2==PV

% indices of local OPF solution vector x = [VA VM PG QG]
[VAscopf, VMscopf, PGscopf, QGscopf] = scopf_aux.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = scopf_aux.index.getLocalIndicesOPF(mpc);

idx = scopf_aux.index.getGlobalIndices(mpc, ns, cont_i);

H_opf = opf_power_balance_hess({x(idx(VAscopf)), x(idx(VMscopf)), x(idx(PGscopf)), x(idx(QGscopf))}, lam, mpc, Ybus, mpopt);

H_scopf = sparse(n,n);
H_scopf(idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx) VMscopf(PVbus_idx) PGscopf(nREFgen_idx)]),...
   idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx) VMscopf(PVbus_idx) PGscopf(nREFgen_idx)])) ...
   = ...
H_opf([VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx) VMopf(PVbus_idx) PGopf(nREFgen_idx)],...
  [VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx) VMopf(PVbus_idx) PGopf(nREFgen_idx)]);

%calls opf_branch_flow_hess(x, lam).
%Before passing x to evaluate hessian it is permuted to OPF ordering
%and each set of constraints for given contingency extract different part 
%of x containing local variables for the contingency. Returned hessian
%is in OPF ordering, therefore we perform appropriate permutation
%with indices corresponding to i-th contingency
function H_scopf = scopf_branch_flow_hess(x, lam, cont_i, mpc, scopf_aux, Yf, Yt, il, mpopt)
n = length(x);
ns = length(scopf_aux.cont);

% get indices of REF gen and of REF/PV buses
[PVbus_idx, nPVbus_idx] = scopf_aux.index.getXbuses(mpc,2);%2==PV

% indices of local OPF solution vector x = [VA VM PG QG]
[VAscopf, VMscopf, PGscopf, QGscopf] = scopf_aux.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = scopf_aux.index.getLocalIndicesOPF(mpc);

idx = scopf_aux.index.getGlobalIndices(mpc, ns, cont_i);

H_opf = opf_branch_flow_hess({x(idx(VAscopf)), x(idx(VMscopf))}, lam, mpc, Yf(il, :), Yt(il, :), il, mpopt);

H_scopf = sparse(n,n);
H_scopf(idx([VAscopf VMscopf(nPVbus_idx) VMscopf(PVbus_idx)]),...
        idx([VAscopf VMscopf(nPVbus_idx) VMscopf(PVbus_idx)])) ...
   = ...
H_opf([VAopf VMopf(nPVbus_idx) VMopf(PVbus_idx)],...
      [VAopf VMopf(nPVbus_idx) VMopf(PVbus_idx)]);