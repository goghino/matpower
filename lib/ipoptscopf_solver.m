function [results, success, raw] = ipoptscopf_solver(om, cont, mpopt)
%IPOPTOPF_SOLVER  Solves AC optimal power flow with security constraints using IPOPT.
%
%   [RESULTS, SUCCESS, RAW] = IPOPTSCOPF_SOLVER(OM, CONT, MPOPT)
%
%   Inputs are an OPF model object, list of branch contingencies CONT
%   and a MATPOWER options struct.
%
%   Outputs are a RESULTS struct, SUCCESS flag and RAW output struct.
%
%   RESULTS is a MATPOWER case struct (mpc) with the usual baseMVA, bus
%   branch, gen, gencost fields, along with the following additional
%   fields:
%       .order      see 'help ext2int' for details of this field
%       .x          final value of optimization variables (internal order)
%       .f          final objective function value
%       .mu         shadow prices on ...
%           .var
%               .l  lower bounds on variables
%               .u  upper bounds on variables
%           .nln
%               .l  lower bounds on nonlinear constraints
%               .u  upper bounds on nonlinear constraints
%           .lin
%               .l  lower bounds on linear constraints
%               .u  upper bounds on linear constraints
%
%   SUCCESS     1 if solver converged successfully, 0 otherwise
%
%   RAW         raw output in form returned by MINOS
%       .xr     final value of optimization variables
%       .info   solver specific termination code
%       .output solver specific output information
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


%% architecture of SCOPF implementation
% Add list of contingencies to options.auxdata and construct
% admittance matrices on fly and construct hessian/jacobian accordingly
% by assembling it from the individual scenarios
% --> need to be more efficient while re-constructing Y matrices

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

%% problem dimensions
nb = size(bus, 1);          %% number of buses
ng = size(gen, 1);          %% number of gens
nl = size(branch, 1);       %% number of branches
ns = size(cont, 1);         %% number of scenarios (nominal + ncont)

%% bounds on optimization vars
[x0, xmin, xmax] = getv(om); %[Va Vm Pg Qg]
xmax = xmax + 1e-10;

% copy bounds for all scenarios and append global PG/QG limits
xl = xmin(1:2*nb);
xg = xmin(2*nb+(1:2*ng));
xmin = [repmat(xl, [ns, 1]); xg];
xl = xmax(1:2*nb);
xg = xmax(2*nb+(1:2*ng));
xmax = [repmat(xl, [ns, 1]); xg];

%% try to select an interior initial point
if mpopt.opf.init_from_mpc ~= 1
    ll = xmin; uu = xmax;
    ll(xmin == -Inf) = -1e10;               %% replace Inf with numerical proxies
    uu(xmax ==  Inf) =  1e10;
    x0 = (ll + uu) / 2;                     %% set x0 mid-way between bounds
    k = find(xmin == -Inf & xmax < Inf);    %% if only bounded above
    x0(k) = xmax(k) - 1;                    %% set just below upper bound
    k = find(xmin > -Inf & xmax == Inf);    %% if only bounded below
    x0(k) = xmin(k) + 1;                    %% set just above lower bound
    Varefs = bus(bus(:, BUS_TYPE) == REF, VA) * (pi/180);
    
    for i = 0:ns-1
        x0(i*2*nb + (1:nb)) = Varefs(1);      %% angles set to first reference angle
    end
end

%% find branches with flow limits
il = [1:nl]';               %% every branch has implicit bounds
                             % find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
                             % TODO insert default limits to branches that
                             % do not satisfy condition above
nl2 = length(il);           %% number of constrained lines

%%-----  run opf  -----
%% build Jacobian and Hessian structure
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
Cl = Cf + Ct;                                   %% for each line - from & to 
Cb = Cl' * Cl + speye(nb);                      %% for each bus - contains adjacent buses
Cl2 = Cl(il, :);                                %% branches with active flow limit
Cg = sparse(gen(:, GEN_BUS), (1:ng)', 1, nb, ng); %%locations where each gen. resides

% Jacobian for the power flow constraints:
%    | Cb   Cb  Cg 0| 
%    |              |
%    | Cb   Cb  0 Cg|
%    |              |
%    | Cl   Cl  0  0| 
%    |              |
%    | Cl   Cl  0  0|
Js_local = [
    Cb      Cb;
    Cb      Cb;
    Cl2     Cl2;
    Cl2     Cl2;
];

CG = [
 Cg sparse(nb, ng);
 sparse(nb, ng) Cg;
 sparse(2*nl2, 2*ng);
];

Js = kron(eye(ns), Js_local);
Js = [Js kron(ones(ns,1), CG)];

% Hessian of lagrangian
% Hs = f(x)_dxx + c(x)_dxx + h(x)_dxx
% | Cb   Cb  0  0| 
% |              |
% | Cb   Cb  0  0|
%
% h(x)_dxx = zeros(nl2,:)

Hs_local =[
    Cb  Cb
    Cb  Cb
];
Hs = kron(eye(ns), Hs_local);
%append hessian w.r.t global variables to lower right corner
xCost = [zeros(2*nb,1); x0(ns*2*nb+1:end)]; %pass global x0, only PG/QG are important, rest fill with zeros
[f, df, d2f] = opf_costfcn(xCost, om);
Hs = [Hs              zeros(size(Hs,1), 2*ng);
      zeros(2*ng, size(Hs,2)) d2f(2*nb+1:end,2*nb+1:end)]; %assuming that d2f has nonzeros only wrt Pg/Qg
Hs = tril(Hs);

%% set options struct for IPOPT
options.ipopt = ipopt_options([], mpopt);

%% extra data to pass to functions
options.auxdata = struct( ...
    'om',       om, ...
    'cont',     cont, ...
    'mpopt',    mpopt, ...
    'il',       il, ...
    'Js',       Js, ...
    'Hs',       Hs    );

%% define variable and constraint bounds
options.lb = xmin;
options.ub = xmax;
options.cl = repmat([zeros(2*nb, 1);  -Inf(2*nl2, 1)], [ns, 1]);
options.cu = repmat([zeros(2*nb, 1); zeros(2*nl2, 1)], [ns, 1]);

%% assign function handles
funcs.objective         = @objective;
funcs.gradient          = @gradient;
funcs.constraints       = @constraints;
funcs.jacobian          = @jacobian;
funcs.hessian           = @hessian;
funcs.jacobianstructure = @(d) Js;
funcs.hessianstructure  = @(d) Hs;

%% run the optimization %TODO call of the ipopt
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
    output.iterations = info.iter;
else
    output.iterations = [];
end

f = opf_costfcn([zeros(2*nb,1); x(ns*2*nb+1:end)], om); %assuming only pg/qg are relevant

%% update solution data for nominal senario and global vars
Va = x(vv.i1.Va:vv.iN.Va);
Vm = x(vv.i1.Vm:vv.iN.Vm);
Pg = x(ns*2*nb + (1:ng));
Qg = x(ns*2*nb + ng + (1:ng));
V = Vm .* exp(1j*Va);

%%-----  calculate return values  -----
%% update voltages & generator outputs
bus(:, VA) = Va * 180/pi;
bus(:, VM) = Vm;
gen(:, PG) = Pg * baseMVA;
gen(:, QG) = Qg * baseMVA;
gen(:, VG) = Vm(gen(:, GEN_BUS));

%% compute branch flows
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
Sf = V(branch(:, F_BUS)) .* conj(Yf * V);  %% cplx pwr at "from" bus, p.u.
St = V(branch(:, T_BUS)) .* conj(Yt * V);  %% cplx pwr at "to" bus, p.u.
branch(:, PF) = real(Sf) * baseMVA;
branch(:, QF) = imag(Sf) * baseMVA;
branch(:, PT) = real(St) * baseMVA;
branch(:, QT) = imag(St) * baseMVA;

results = mpc;
[results.bus, results.branch, results.gen, ...
    results.om, results.x, results.f] = ...
        deal(bus, branch, gen, om, [x(1:2*nb); x(ns*2*nb+(1:2*ng))], f);
    
raw = struct('xr', x, 'info', info.status, 'output', output);

%-----  callback functions  -----
function f = objective(x, d)
mpc = get_mpc(d.om);
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches
ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)

f = opf_costfcn([zeros(2*nb,1); x(ns*2*nb+1:end)], d.om); %assuming only pg/qg are relevant

function df = gradient(x, d)
mpc = get_mpc(d.om);
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches
ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)

[f, df] = opf_costfcn([zeros(2*nb,1); x(ns*2*nb+(1:2*ng))], d.om); %assuming only pg/qg are relevant
df = [zeros((ns-1)*2*nb, 1); df];

function constr = constraints(x, d)
mpc = get_mpc(d.om);
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches
ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)
NCONSTR = 2*nb + 2*nl;

constr = zeros(ns*(NCONSTR), 1);

for i = 0:ns-1
    c = d.cont(i+1);
    xl = [x(i*2*nb + (1:2*nb)); x(ns*2*nb + (1:2*ng))]; % extract local [Vai Vmi] and append global [Pg Qg]
    [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, c);
    [hn_local, gn_local] = opf_consfcn(xl, d.om, Ybus, Yf, Yt, d.mpopt, d.il);
    constr(i*(NCONSTR) + (1:NCONSTR)) = [gn_local; hn_local];
end


function J = jacobian(x, d)
mpc = get_mpc(d.om);
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches
ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)
NCONSTR = 2*nb + 2*nl;

J = sparse(ns*(NCONSTR), size(x,1));


for i = 0:ns-1
    c = d.cont(i+1);
    xl = [x(i*2*nb + (1:2*nb)); x(ns*2*nb + (1:2*ng))]; % extract logal [Vai Vmi] and append global [Pg Qg]
    [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, c);
    [hn, gn, dhn, dgn] = opf_consfcn(xl, d.om, Ybus, Yf, Yt, d.mpopt, d.il);
    dgn = dgn';
    dhn = dhn';
    %global variables are ordered to be last in x, need to insert local
    %jacobian accordingly into the global structure
    J(i*NCONSTR + (1:NCONSTR), i*2*nb + (1:2*nb)) = [dgn(:,1:2*nb); dhn(:,1:2*nb)];
    J(i*NCONSTR + (1:2*nb), ns*2*nb + (1:2*ng)) = dgn(:,2*nb+(1:2*ng));
end


function H = hessian(x, sigma, lambda, d)
mpc = get_mpc(d.om);
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches
ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)
NCONSTR = 2*nb + 2*nl;

H = sparse(size(x,1), size(x,1));

for i = 0:ns-1
    c = d.cont(i+1);
    xl = [x(i*2*nb + (1:2*nb)); x(ns*2*nb + (1:2*ng))]; % extract logal [Vai Vmi] and append global [Pg Qg]
    lam.eqnonlin   = lambda(i*NCONSTR + (1:2*nb));
    lam.ineqnonlin = lambda(i*NCONSTR + 2*nb + (1:2*nl));
    [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, c);
    %compute hessian w.r.t local variables and paste in into the global one
    H_local = opf_hessfcn(xl, lam, sigma, d.om, Ybus, Yf, Yt, d.mpopt, d.il);
    H(i*2*nb+(1:2*nb), i*2*nb+(1:2*nb)) = H_local(1:2*nb, 1:2*nb);
end
%hessian w.r.t global variables goes to lower right corner
H(ns*2*nb+1:end, ns*2*nb+1:end) = H_local(2*nb+1:end, 2*nb+1:end);

H = tril(H);

