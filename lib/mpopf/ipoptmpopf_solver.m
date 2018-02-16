function [results, success, raw] = ipoptmpopf_solver(om, model, mpopt)
%IPOPTOPF_SOLVER  Solves AC optimal power flow using IPOPT.
%
%   [RESULTS, SUCCESS, RAW] = IPOPTMPOPF_SOLVER(OM, MPOPT)
%
%   Inputs are an MPOPF model object and a MATPOWER options struct.
%
%   Model is a struct with following fields:
%       .profile Containts a list of load profile over the time horizon
%       .index Contains functions to handle proper indexing of MPOPF variables
%           .getGlobalIndices
%           .getLocalIndicesOPF
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
%           .nln    (deprecated) 2*nb+2*nl - Pmis, Qmis, Sf, St
%               .l  lower bounds on nonlinear constraints
%               .u  upper bounds on nonlinear constraints
%           .nle    nonlinear equality constraints
%           .nli    nonlinear inequality constraints
%           .lin
%               .l  lower bounds on linear constraints
%               .u  upper bounds on linear constraints
%
%   SUCCESS     1 if solver converged successfully, 0 otherwise
%
%   RAW         raw output in form returned by MINOS
%       .xr     final value of optimization variables
%       .pimul  constraint multipliers
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
mpc = om.get_mpc();
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
[vv, ll, nne, nni] = om.get_idx();

%% problem dimensions
Nt = length(model.profile);
nb = Nt * size(bus, 1);          %% number of buses
ng = Nt * size(gen, 1);          %% number of gens
nl = Nt * size(branch, 1);       %% number of branches
ny = om.getN('var', 'y');   %% number of piece-wise linear costs

%number of EQ and INEQ constriants
neq = 2*nb;
niq = 2*nl;

% indices of local OPF solution vector x = [VA VM PG QG]
[VAopf, VMopf, PGopf, QGopf] = model.index.getLocalIndicesOPF(mpc);

%% linear constraints
[A, l, u] = om.params_lin_constraint();

%% bounds on optimization vars
[x0, xmin, xmax] = om.params_var();

% Note that variables with equal upper and lower bounds are removed by IPOPT
% so we add small perturbation to x_u[], we don't want them removed
% because the Schur solver assumes particular structure that would
% be changed by removing variables.
idx = find(xmin == xmax);
xmax(idx) = xmax(idx) + 1e-10;
%exept for the Va at the refernece bus which we want to remove
ref_bus = find(mpc.bus(:,BUS_TYPE)==3);
for i = 0:Nt-1
    idx = model.index.getGlobalIndices(mpc, Nt, i);
    xmax(idx(VAopf(ref_bus))) = xmin(ref_bus);
end

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% try to select an interior initial point, unless requested not to
if mpopt.opf.start < 2
    s = 1e-3;                   %% set init point inside bounds by s
    lb = xmin; ub = xmax;
    lb(xmin == -Inf) = -1e10;   %% replace Inf with numerical proxies
    ub(xmax ==  Inf) =  1e10;
    x0 = (lb + ub) / 2;         %% set x0 mid-way between bounds
    k = find(xmin == -Inf & xmax < Inf);    %% if only bounded above
    x0(k) = xmax(k) - s;                    %% set just below upper bound
    k = find(xmin > -Inf & xmax == Inf);    %% if only bounded below
    x0(k) = xmin(k) + s;                    %% set just above lower bound
    
    % adjust voltage angles to match reference bus
    Varefs = bus(ref_bus, VA) * (pi/180);
    for i = 0:Nt-1
        idx = model.index.getGlobalIndices(mpc, Nt, i);
        x0(idx(VAopf)) = Varefs(1);
    end
    
    if ny > 0
        error('ny > 0 not considered in MPOPF');
        ipwl = find(gencost(:, MODEL) == PW_LINEAR);
        c = gencost(sub2ind(size(gencost), ipwl, NCOST+2*gencost(ipwl, NCOST)));    %% largest y-value in CCV data
        x0(vv.i1.y:vv.iN.y) = max(c) + 0.1 * abs(max(c));
    end
end

%% find branches with flow limits
il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
nl2 = length(il);           %% number of constrained lines
nx = length(x0);

%% replace equality variable bounds with an equality constraint
%% (since IPOPT does not return shadow prices on variables that it eliminates)
kk = find(xmin(nb+1:end) == xmax(nb+1:end));    %% all bounds except ref angles
nk = length(kk);
if nk
    error('nk > 0 not considered in MPOPF');
    kk = kk + nb;               %% adjust index for missing ref angles
    A = [ A; sparse((1:nk)', kk, 1, nk, nx) ];
    l = [ l; xmin(kk) ];
    u = [ u; xmax(kk) ];
    xmin(kk) = -Inf;
    xmax(kk) = Inf;
end

%%-----  run mpopf  -----

%% set options struct for IPOPT
options.ipopt = ipopt_options([], mpopt);

%% extra data to pass to functions
options.auxdata = struct( ...
    'om',       om, ...
    'Ybus',     Ybus, ...
    'Yf',       Yf(il,:), ...
    'Yt',       Yt(il,:), ...
    'profile',  model.profile, ...
    'index',    model.index, ...
    'mpopt',    mpopt, ...
    'il',       il, ...
    'A',        A, ...
    'nA',       size(A, 1), ...
    'neqnln',   neq, ...
    'niqnln',   niq );

%% define variable and constraint bounds
options.lb = xmin;
options.ub = xmax;
options.cl = [zeros(neq, 1);  -Inf(niq, 1); l];
options.cu = [zeros(neq, 1); zeros(niq, 1); u+1e-10];

%% build Jacobian and Hessian structure
randx = rand(size(x0));
[h, g, dh, dg] = opf_consfcn(randx, om, Ybus, Yf(il,:), Yt(il,:), mpopt, il);
Js = [dg'; dh'; A];
lam = struct('eqnonlin', ones(size(dg,2),1), 'ineqnonlin', ones(size(dh,2),1) );
Hs = tril(mpopf_hessfcn(randx, lam, 1, options.auxdata));

%% assign function handles
funcs.objective         = @objective;
funcs.gradient          = @gradient;
funcs.constraints       = @constraints;
funcs.jacobian          = @jacobian;
funcs.hessian           = @hessian;
funcs.jacobianstructure = @(d) Js;
funcs.hessianstructure  = @(d) Hs;
%funcs.jacobianstructure = @jacobianstructure;
%funcs.hessianstructure  = @hessianstructure;

%% run the optimization %TODO call of the ipopt
if 1 %have_fcn('ipopt_auxdata')
    [x, info] = ipopt_auxdata(x0,funcs,options);
else
    [x, info] = ipopt(x0,funcs,options);
end

%% handle outputs
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
f = opf_costfcn(x, om);

%% print value of the constraints
%disp('Final value of the constraints:');
%[hn, gn] = opf_consfcn(x, om, Ybus, Yf, Yt, mpopt, il)

results = mpc;
[results.om, results.x, results.f] = ...
        deal(om, x, f);

raw = struct('xr', x, 'info', info.status, 'output', output);

%-----  callback functions  -----
% opf_costfcn
% opf_consfcn()
% mpopf_hessfcn()
%  are calling the callbacks defined in mpopf_setup.m
%  mpopf_power_flow_fcn/hess()
%  mpopf_branch_flow_fcn/hess()
%  or evaluates costs defined by add_quad_cost() coefficients

function f = objective(x, d)
f = opf_costfcn(x, d.om);

function df = gradient(x, d)
[f, df] = opf_costfcn(x, d.om);

function c = constraints(x, d)
[hn, gn] = opf_consfcn(x, d.om, d.Ybus, d.Yf, d.Yt, d.mpopt, d.il);
if isempty(d.A)
    c = [gn; hn];
else
    c = [gn; hn; d.A*x];
end

function J = jacobian(x, d)
[hn, gn, dhn, dgn] = opf_consfcn(x, d.om, d.Ybus, d.Yf, d.Yt, d.mpopt, d.il);
J = [dgn'; dhn'; d.A];


function H = hessian(x, sigma, lambda, d)
lam.eqnonlin   = lambda(1:d.neqnln);
lam.ineqnonlin = lambda(d.neqnln+(1:d.niqnln));

H = tril(mpopf_hessfcn(x, lam, sigma, d));

% function Js = jacobianstructure(d)
% Js = d.Js;
% 
% function Hs = hessianstructure(d)
% Hs = d.Hs;