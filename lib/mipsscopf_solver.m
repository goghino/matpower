function [success, raw] = mipsscopf_solver(om, cont, mpopt)
%MIPSOPF_SOLVER  Solves AC optimal power flow using MIPS.
%
%   [RESULTS, SUCCESS, RAW] = MIPSOPF_SOLVER(OM, MPOPT)
%
%   Inputs are an OPF model object and a MATPOWER options struct.
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
%       .pimul  constraint multipliers
%       .info   solver specific termination code
%       .output solver specific output information
%
%   See also OPF, MIPS.

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

%% options
opt = mpopt.mips;
opt.verbose = mpopt.verbose;
if opt.feastol == 0
    opt.feastol = mpopt.opf.violation;  %% = MPOPT.opf.violation by default
end
if ~isfield(opt, 'cost_mult') || isempty(opt.cost_mult)
    opt.cost_mult = 1e-4;
end

%% unpack data
mpc = get_mpc(om);
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
[vv, ll, nn] = get_idx(om);

%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ng = size(gen, 1);          %% number of generators
ns = size(cont, 1);         %% number of scenarios (nominal + ncont)

%% linear constraints
[A, l, u] = linear_constraints(om);

%% bounds on optimization vars
[x0, xmin, xmax] = getv(om);

xl = xmin(1:2*nb);
xg = xmin(2*nb+(1:2*ng));
xmin = [repmat(xl, [ns, 1]); xg];

xl = xmax(1:2*nb);
xg = xmax(2*nb+(1:2*ng));
xmax = [repmat(xl, [ns, 1]); xg];

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% try to select an interior initial point
if mpopt.opf.init_from_mpc ~= 1
    ll = xmin; uu = xmax;
    ll(xmin == -Inf) = -1e10;   %% replace Inf with numerical proxies
    uu(xmax ==  Inf) =  1e10;
    x0 = (ll + uu) / 2;         %% set x0 mid-way between bounds
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

%% -----  run opf  -----
auxdata = struct( ...
    'om',       om, ...
    'cont',     cont, ...
    'mpopt',    mpopt );

f_fcn = @(x) objective_fcn(x, auxdata);
gh_fcn = @(x) constraints_fcn(x, auxdata);
hess_fcn = @(x, lambda, cost_mult) hessian_fcn(x, lambda, cost_mult, auxdata);

[x, f, info, output, Lambda] = ...
  mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt, mpc);
  %pmips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt, mpc);
success = (info > 0);

% %% update solution data
% Va = x(vv.i1.Va:vv.iN.Va);
% Vm = x(vv.i1.Vm:vv.iN.Vm);
% Pg = x(vv.i1.Pg:vv.iN.Pg);
% Qg = x(vv.i1.Qg:vv.iN.Qg);
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
% Sf = V(branch(:, F_BUS)) .* conj(Yf * V);  %% cplx pwr at "from" bus, p.u.
% St = V(branch(:, T_BUS)) .* conj(Yt * V);  %% cplx pwr at "to" bus, p.u.
% branch(:, PF) = real(Sf) * baseMVA;
% branch(:, QF) = imag(Sf) * baseMVA;
% branch(:, PT) = real(St) * baseMVA;
% branch(:, QT) = imag(St) * baseMVA;
    
%pack some additional info to output so that we can verify the solution
output.Ybus = Ybus;
output.Yf = Yf;
output.Yt = Yt;
output.lb = xmin;
output.ub = xmax;

raw = struct('x', x, 'info', info, 'output', output);



%% callback routines
%evaluate objective, its gradient and hessian
function [f, df] = objective_fcn(x, d) 
    mpc = get_mpc(d.om);
    nb = size(mpc.bus, 1);          %% number of buses
    ng = size(mpc.gen, 1);          %% number of gens
    nl = size(mpc.branch, 1);       %% number of branches
    ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)
    
    [f, df] = opf_costfcn([zeros(2*nb,1); x(ns*2*nb+(1:2*ng))], d.om); %assuming only pg/qg are relevant
    df = [zeros((ns-1)*2*nb, 1); df];

%evaluate constraints and jacobian
function [hn, gn, dhn, dgn] = constraints_fcn(x, d)
    
    mpc = get_mpc(d.om);
    nb = size(mpc.bus, 1);          %% number of buses
    ng = size(mpc.gen, 1);          %% number of gens
    nl = size(mpc.branch, 1);       %% number of branches
    ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)
    NCONSTR_dgn = 2*nb;
    NCONSTR_dhn = 2*nl;

    dgn = sparse(size(x,1), ns*(NCONSTR_dgn));
    dhn = sparse(size(x,1), ns*(NCONSTR_dhn));

    hn = [];
    gn = [];

    for i = 0:ns-1
        cont = d.cont(i+1);
        xl = [x(i*2*nb + (1:2*nb)); x(ns*2*nb + (1:2*ng))]; % extract local [Vai Vmi] and append global [Pg Qg]
        [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
        [hn_local, gn_local, dhn_local, dgn_local] = opf_consfcn(xl, d.om, Ybus, Yf, Yt, d.mpopt);
        
        gn = [gn; gn_local];
        hn = [hn; hn_local];
        
        %insert into global struct, dgn & dhn are transposed
        dgn(i*2*nb + (1:2*nb), i*NCONSTR_dgn + (1:NCONSTR_dgn)) = dgn_local(1:2*nb, :); %diagonal block
        dgn(ns*2*nb + (1:2*ng), i*NCONSTR_dgn + (1:2*nb)) = dgn_local(2*nb+(1:2*ng), :); %generator coupling
   
        dhn(i*2*nb + (1:2*nb), i*NCONSTR_dhn + (1:NCONSTR_dhn)) = dhn_local(1:2*nb, :); %diagonal block
    end
                
%evaluate hessian of Jacobian    
function H = hessian_fcn(x, lambda, sigma, d)
    mpc = get_mpc(d.om);
    nb = size(mpc.bus, 1);          %% number of buses
    ng = size(mpc.gen, 1);          %% number of gens
    nl = size(mpc.branch, 1);       %% number of branches
    ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)

    H = sparse(size(x,1), size(x,1));

    for i = 0:ns-1
        cont = d.cont(i+1);
        xl = [x(i*2*nb + (1:2*nb)); x(ns*2*nb + (1:2*ng))]; % extract logal [Vai Vmi] and append global [Pg Qg]
        lam.eqnonlin   = lambda.eqnonlin(i*2*nb + (1:2*nb));
        lam.ineqnonlin = lambda.ineqnonlin(i*2*nl + (1:2*nl));
        [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
        %compute hessian w.r.t local variables and insert in into the global one
        H_local = opf_hessfcn(xl, lam, sigma, d.om, Ybus, Yf, Yt, d.mpopt, 1:nl);
        H(i*2*nb+(1:2*nb), i*2*nb+(1:2*nb)) = H_local(1:2*nb, 1:2*nb);
    end
    %hessian w.r.t global variables goes to lower right corner
    H(ns*2*nb+1:end, ns*2*nb+1:end) = H_local(2*nb+1:end, 2*nb+1:end);

    H = tril(H);

