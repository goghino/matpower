function [results, success, raw] = ipoptscopf_solver(om, model, mpopt)
%IPOPTOPF_SOLVER  Solves AC optimal power flow with security constraints using IPOPT.
%
%   [RESULTS, SUCCESS, RAW] = IPOPTSCOPF_SOLVER(OM, MODEL, MPOPT)
%
%   Inputs are an OPF model object, SCOPF model and a MATPOWER options struct.
%
%   Model is a struct with following fields:
%       .cont Containts a list of contingencies
%       .index Contains functions to handle proper indexing of SCOPF variables
%           .getGlobalIndices
%           .getLocalIndicesOPF
%           .getLocalIndicesSCOPF
%
%   Outputs are a RESULTS struct, SUCCESS flag and RAW output struct.
%
%   The internal x that ipopt works with has structure
%   [Va1 Vm1 Qg1 Pg_ref1... VaN VmN QgN Pg_refN Pg] for all contingency scenarios 1..N
%   with corresponding bounds xmin < x < xmax
%
%   We impose nonlinear equality and inequality constraints g(x) and h(x)
%   with corresponding bounds cl < [g(x); h(x)] < cu
%   and linear constraints l < Ax < u.
%
%   See also OPF, IPOPT.

%   MATPOWER
%   Copyright (c) 2000-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.


%% TODO
% need to work more efficiently with sparse indexing during construction
% of global hessian/jacobian

% how to account for the sparse() leaving out zeros from the sparse
% structure? We want to have exactly same structure across scenarios

%%----- initialization -----
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

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
[vv, ll, nn] = get_idx(om);

cont = model.cont;


%% problem dimensions
nb = size(bus, 1);          %% number of buses
ng = size(gen, 1);          %% number of gens
nl = size(branch, 1);       %% number of branches
ns = size(cont, 1);         %% number of scenarios (nominal + ncont)

% reference bus and generator
ref_bus = find(mpc.bus(:,BUS_TYPE) == REF);
ref_gen = find(mpc.gen(:,BUS_I) == ref_bus); %index of gen connected to ref_bus

% indices of local parts of solution vector x = [VA VM PG QG]
[VAscopf, VMscopf, PGscopf, QGscopf] = model.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = model.index.getLocalIndicesOPF(mpc);

%% build admittance matrices for nominal case
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% bounds on optimization vars xmin <= x <= xmax 
[x0, xmin, xmax] = getv(om); %returns standard OPF form [Va Vm Pg Qg]
xmax = xmax + 1e-10;

% replicate bounds for all scenarios and append global limits
xl = xmin([VAopf VMopf QGopf PGopf(ref_gen)]); %local variables
xg = xmin(PGopf([1:ref_gen-1 ref_gen+1:ng])); %global variables
xmin = [repmat(xl, [ns, 1]); xg];

xl = xmax([VAopf VMopf QGopf PGopf(ref_gen)]); %local variables
xg = xmax(PGopf([1:ref_gen-1 ref_gen+1:ng])); %global variables
xmax = [repmat(xl, [ns, 1]); xg];

%% try to select an interior initial point based on bounds
if mpopt.opf.init_from_mpc ~= 1
    ll = xmin; uu = xmax;
    ll(xmin == -Inf) = -1e10;               %% replace Inf with numerical proxies
    uu(xmax ==  Inf) =  1e10;
    x0 = (ll + uu) / 2;                     %% set x0 mid-way between bounds
    k = find(xmin == -Inf & xmax < Inf);    %% if only bounded above
    x0(k) = xmax(k) - 1;                    %% set just below upper bound
    k = find(xmin > -Inf & xmax == Inf);    %% if only bounded below
    x0(k) = xmin(k) + 1;                    %% set just above lower bound
    
    % adjust voltage angles to match reference bus
    Varefs = bus(ref_bus, VA) * (pi/180);
    for i = 0:ns-1
        idx = model.index.getGlobalIndices(mpc, ns, i);
        x0(idx(VAscopf)) = Varefs(1);      %% set 1 if
    end
end

%% find branches with flow limits
il_ = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
il = [1:nl]';               %% we assume every branch has implicit bounds
                             % TODO insert default limits to branches that
                             % do not satisfy condition above
nl2 = length(il);           %% number of constrained lines

if size(il_, 1) ~= nl2
   error('Not all branches have specified RATE_A field.'); 
end


%% build linear constraints l <= A*x <= u
% for all PV buses: 0 <= Vm0 - Vmi <= 0 for all contingencies i

busPV = find(mpc.bus(:,BUS_TYPE) == 2); %find PV buses
nPV = size(busPV,1); %number of PV buses
nPart = 2*nb+ng+1; %number of local variables for each scenario [Va Vm Qg Pg_ref]

%build "identity" matrix for Vm0 (nominal case)
A1 = sparse(1:nPV*(ns-1),...
            kron(ones(ns-1,1),busPV) + nb,... replicate column indices and add offest nb
            1,...
            nPV*(ns-1),...
            nPart,...
            nPV*(ns-1));
%build "-identity" matrix for -Vmi (all contingency scenarios)
An = sparse(1:nPV*(ns-1),...
            kron(ones(ns-1,1),busPV) + kron([0:ns-2]',nPart*ones(nPV,1)) + nb,... skip previous local variables (excl. nominal case)
            -1,...
            nPV*(ns-1),...
            nPart*(ns-1),...
            nPV*(ns-1));

A = [A1 An sparse(size(A1,1), ng-1)]; %put together constraints for local and global variables
l = zeros(size(A,1),1);
u = zeros(size(A,1),1);


%% build local connectivity matrices
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
Cl = Cf + Ct;                                   %% for each line - from & to 
Cb = Cl' * Cl + speye(nb);                      %% for each bus - contains adjacent buses
Cl2 = Cl(il, :);                                %% branches with active flow limit
Cg = sparse(gen(:, GEN_BUS), (1:ng)', 1, nb, ng); %%locations where each gen. resides

% Jacobian for the power flow and branch flow constraints:
%     dVa  dVm  dQg dPg_ref   <- local variables for each scenario
%    | Cb   Cb  0   1 | ('one' at appropriate row, otherwise zeros) 
%    |                |
%    | Cb   Cb  Cg  0 |
%    |                |
%    | Cl   Cl  0   0 | 
%    |                |
%    | Cl   Cl  0   0 |
Js_local = [
    Cb      Cb    sparse(nb, ng)   sparse(ref_bus, 1, 1, nb, 1, 1);
    Cb      Cb     Cg              sparse(nb, 1);
    Cl2     Cl2   sparse(nl2, ng+1);
    Cl2     Cl2   sparse(nl2, ng+1);
];

CG = [
 Cg(:,[1:ref_gen-1 ref_gen+1:ng]); %skip column of ref_gen, it is local var.
 sparse(nb, ng-1);
 sparse(2*nl2, ng-1);
];

Js = kron(eye(ns), Js_local); %replicate jac. w.r.t local variables
Js = [Js kron(ones(ns,1), CG)]; % replicate and append jac w.r.t global variables
Js = [Js; A]; %append linear constraints

% Hessian of lagrangian
% Hs = f(x)_dxx + c(x)_dxx + h(x)_dxx

%  dVa  dVm  dQg dPg_ref
% | Cb   Cb   0    0| (nominal case has dPg_ref 1)
% |                 |
% | Cb   Cb   0    0|
%
% h(x)_dxx = zeros(nl2,:)

Hs_local =[
    Cb  Cb sparse(nb, ng+1);
    Cb  Cb sparse(nb, ng+1);
    sparse(ng+1, 2*nb+ng+1);
];
Hs = kron(eye(ns), Hs_local);
Hs(2*nb+ng+1, 2*nb+ng+1) = 1; %set d2Pg_ref to 1 in nominal case

%set structure of objective hessian (cost = f(Pg))
d2f = eye(ng-1); %d2 wrt Pg except Pg_ref

%append hessian w.r.t global variables to lower right corner
Hs = [Hs              zeros(size(Hs,1), ng-1);
      zeros(ng-1, size(Hs,2)) d2f];
Hs = tril(Hs);

%% set options struct for IPOPT
options.ipopt = ipopt_options([], mpopt);

%% extra data to pass to functions
options.auxdata = struct( ...
    'om',       om, ...
    'cont',     cont, ...
    'index',    model.index, ...
    'mpopt',    mpopt, ...
    'il',       il, ...
    'A',        A, ...
    'Js',       Js, ...
    'Hs',       Hs    );

%% define variable and constraint bounds
options.lb = xmin;
options.ub = xmax;
options.cl = [repmat([zeros(2*nb, 1);  -Inf(2*nl2, 1)], [ns, 1]); l];
options.cu = [repmat([zeros(2*nb, 1); zeros(2*nl2, 1)], [ns, 1]); u+1e10]; %add 1e10 so that ipopt doesn't remove l==u case

%% assign function handles
funcs.objective         = @objective;
funcs.gradient          = @gradient;
funcs.constraints       = @constraints;
funcs.jacobian          = @jacobian;
funcs.hessian           = @hessian;
funcs.jacobianstructure = @(d) Js;
funcs.hessianstructure  = @(d) Hs;

%% run the optimization, call ipopt
if 1 %have_fcn('ipopt_auxdata')
    [x, info] = ipopt_auxdata(x0,funcs,options);
else
    [x, info] = ipopt(x0,funcs,options);
end

if info.status == 0 || info.status == 1
    success = 1;
else
    success = 0;
    display(['Ipopt finished with error: ', num2str(info.status)]);
end

if isfield(info, 'iter')
    meta.iterations = info.iter;
else
    meta.iterations = [];
end

idx_nom = model.index.getGlobalIndices(mpc, ns, 0); %evaluate cost of nominal case (only Pg/Qg are relevant) 
f = opf_costfcn(x(idx_nom([VAscopf VMscopf PGscopf QGscopf])), om);

% %% update solution data for nominal senario and global vars
% Va = x(vv.i1.Va:vv.iN.Va);
% Vm = x(vv.i1.Vm:vv.iN.Vm);
% Pg = x(ns*2*nb + (1:ng));
% Qg = x(ns*2*nb + ng + (1:ng));
% V = Vm .* exp(1j*Va);
% 
% %%-----  calculate return values  -----
% %% update voltages & generator outputs
% bus(:, VA) = Va * 180/pi;
% bus(:, VM) = Vm;
% gen(:, PG) = Pg * baseMVA;
% gen(:, QG) = Qg * baseMVA;
% gen(:, VG) = Vm(gen(:, GEN_BUS));
% 
% %% compute branch flows
% [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
% Sf = V(branch(:, F_BUS)) .* conj(Yf * V);  %% cplx pwr at "from" bus, p.u.
% St = V(branch(:, T_BUS)) .* conj(Yt * V);  %% cplx pwr at "to" bus, p.u.
% branch(:, PF) = real(Sf) * baseMVA;
% branch(:, QF) = imag(Sf) * baseMVA;
% branch(:, PT) = real(St) * baseMVA;
% branch(:, QT) = imag(St) * baseMVA;
    
%pack some additional info to output so that we can verify the solution
meta.Ybus = Ybus;
meta.Yf = Yf;
meta.Yt = Yt;
meta.lb = options.lb;
meta.ub = options.ub;
meta.A = A;
    
raw = struct('info', info.status, 'meta', meta);
results = struct('f', f, 'x', x);

%% -----  callback functions  -----
function f = objective(x, d)
mpc = get_mpc(d.om);
ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)

% use nominal case to evaluate cost fcn (only pg/qg are relevant)
idx_nom = d.index.getGlobalIndices(mpc, ns, 0);
[VAscopf, VMscopf, PGscopf, QGscopf] = d.index.getLocalIndicesSCOPF(mpc);

f = opf_costfcn(x(idx_nom([VAscopf VMscopf PGscopf QGscopf])), d.om); %assuming only pg/qg are relevant

function grad = gradient(x, d)
mpc = get_mpc(d.om);
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)

%evaluate grad of nominal case
idx_nom = d.index.getGlobalIndices(mpc, ns, 0);
[VAscopf, VMscopf, PGscopf, QGscopf] = d.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = d.index.getLocalIndicesOPF(mpc);

[f, df] = opf_costfcn(x(idx_nom([VAscopf VMscopf PGscopf QGscopf])), d.om);

grad = zeros(size(x,1),1);
grad(idx_nom(PGscopf)) = df(PGopf); %nonzero only nominal case Pg


function constr = constraints(x, d)
mpc = get_mpc(d.om);
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches
ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)
NCONSTR = 2*nb + 2*nl;

constr = zeros(ns*(NCONSTR), 1);

[VAscopf, VMscopf, PGscopf, QGscopf] = d.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = d.index.getLocalIndicesOPF(mpc);

for i = 0:ns-1
    cont = d.cont(i+1);
    idx = d.index.getGlobalIndices(mpc, ns, i);
    [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
    [hn_local, gn_local] = opf_consfcn(x(idx([VAscopf VMscopf PGscopf QGscopf])), d.om, Ybus, Yf, Yt, d.mpopt, d.il);
    constr(i*(NCONSTR) + (1:NCONSTR)) = [gn_local; hn_local];
end

if ~isempty(d.A)
    constr = [constr; d.A*x]; %append linear constraints
end


function J = jacobian(x, d)
mpc = get_mpc(d.om);
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches
ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)
NCONSTR = 2*nb + 2*nl;

J = sparse(ns*(NCONSTR), size(x,1));

BUS_TYPE = 2;
REF = 3;
GEN_BUS = 1;
ref_bus = find(mpc.bus(:,BUS_TYPE) == REF);
ref_gen = find(mpc.gen(:,GEN_BUS) == ref_bus); %index of gen connected to ref_bus

[VAscopf, VMscopf, PGscopf, QGscopf] = d.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = d.index.getLocalIndicesOPF(mpc);

for i = 0:ns-1
    %compute local indices and its parts
    idx = d.index.getGlobalIndices(mpc, ns, i);
    
    cont = d.cont(i+1);
    [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
    [hn, gn, dhn, dgn] = opf_consfcn(x(idx([VAscopf VMscopf PGscopf QGscopf])), d.om, Ybus, Yf, Yt, d.mpopt, d.il);
    dgn = dgn';
    dhn = dhn';
    
    %global variables are ordered to be last in x, need to insert local
    %jacobian accordingly into the global structure
    J(i*NCONSTR + (1:NCONSTR), idx([VAscopf VMscopf QGscopf PGscopf(ref_gen)])) = [dgn(:,[VAopf VMopf QGopf PGopf(ref_gen)]); dhn(:,[VAopf VMopf QGopf PGopf(ref_gen)])]; %dVai dVmi dQg dPg_ref
    J(i*NCONSTR + (1:size(dgn,1)), idx(PGscopf([1:ref_gen-1 ref_gen+1:ng]))) = dgn(:, PGopf([1:ref_gen-1 ref_gen+1:ng]));%dPg 
end
J = [J; d.A]; %append Jacobian of linear constraints


function H = hessian(x, sigma, lambda, d)
mpc = get_mpc(d.om);
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches
ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)
NCONSTR = 2*nb + 2*nl;

H = sparse(size(x,1), size(x,1));

BUS_TYPE = 2;
REF = 3;
GEN_BUS = 1;
ref_bus = find(mpc.bus(:,BUS_TYPE) == REF);
ref_gen = find(mpc.gen(:,GEN_BUS) == ref_bus); %index of gen connected to ref_bus

[VAscopf, VMscopf, PGscopf, QGscopf] = d.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = d.index.getLocalIndicesOPF(mpc);

for i = 0:ns-1
    %compute local indices and its parts
    idx = d.index.getGlobalIndices(mpc, ns, i);
    
    cont = d.cont(i+1);
    [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
    
    lam.eqnonlin   = lambda(i*NCONSTR + (1:2*nb));
    lam.ineqnonlin = lambda(i*NCONSTR + 2*nb + (1:2*nl));
    H_local = opf_hessfcn(x(idx([VAscopf VMscopf PGscopf QGscopf])), lam, sigma, d.om, Ybus, Yf, Yt, d.mpopt, d.il);
    
    %insert hessian w.r.t local variables in into the full hessian
    if i == 0
        %include derivative w.r.t Pg_ref in nominal case, because only this is considered in the cost function
        H(idx([VAscopf VMscopf QGscopf PGscopf(ref_gen)]), idx([VAscopf VMscopf QGscopf PGscopf(ref_gen)])) = H_local([VAopf VMopf QGopf PGopf(ref_gen)], [VAopf VMopf QGopf PGopf(ref_gen)]);
        %hessian w.r.t global variables goes to lower right corner
        H(idx(PGscopf([1:ref_gen-1 ref_gen+1:ng])), idx(PGscopf([1:ref_gen-1 ref_gen+1:ng]))) = H_local(PGopf([1:ref_gen-1 ref_gen+1:ng]), PGopf([1:ref_gen-1 ref_gen+1:ng]));
    else
        %skip Pg_ref in contingency scenarios, it is not used in objective function
        H(idx([VAscopf VMscopf QGscopf]), idx([VAscopf VMscopf QGscopf])) = H_local([VAopf VMopf QGopf], [VAopf VMopf QGopf]);
    end
end

H = tril(H);
