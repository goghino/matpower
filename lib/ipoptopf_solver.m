function [results, success, raw] = ipoptopf_solver(om, mpopt)
%IPOPTOPF_SOLVER  Solves AC optimal power flow using IPOPT.
%
%   [RESULTS, SUCCESS, RAW] = IPOPTOPF_SOLVER(OM, MPOPT)
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
mpc = get_mpc(om);
[baseMVA, bus, gen, branch, gencost] = ...
    deal(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, mpc.gencost);
[vv, ll, nn] = get_idx(om);

%% problem dimensions
nb = size(bus, 1);          %% number of buses
ng = size(gen, 1);          %% number of gens
nl = size(branch, 1);       %% number of branches
ny = getN(om, 'var', 'y');  %% number of piece-wise linear costs

%% linear constraints
[A, l, u] = linear_constraints(om);

%% bounds on optimization vars
[x0, xmin, xmax] = getv(om);

% Note that variables with equal upper and lower bounds are removed by IPOPT
% so we add small perturbation to x_u[], we don't want them removed
% because the Schur solver assumes particular structure that would
% be changed by removing variables.
idx = find(xmin == xmax);
xmax(idx) = xmax(idx) + 1e-10;
%exept for the Va at the refernece bus which we want to remove
ref_bus = find(mpc.bus(:,BUS_TYPE)==3);
xmax(ref_bus) = xmin(ref_bus);

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
    x0(vv.i1.Va:vv.iN.Va) = Varefs(1);  %% angles set to first reference angle
    if ny > 0
        ipwl = find(gencost(:, MODEL) == PW_LINEAR);
    %     PQ = [gen(:, PMAX); gen(:, QMAX)];
    %     c = totcost(gencost(ipwl, :), PQ(ipwl));
        c = gencost(sub2ind(size(gencost), ipwl, NCOST+2*gencost(ipwl, NCOST)));    %% largest y-value in CCV data
        x0(vv.i1.y:vv.iN.y) = max(c) + 0.1 * abs(max(c));
    %     x0(vv.i1.y:vv.iN.y) = c + 0.1 * abs(c);
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
    kk = kk + nb;               %% adjust index for missing ref angles
    A = [ A; sparse((1:nk)', kk, 1, nk, nx) ];
    l = [ l; xmin(kk) ];
    u = [ u; xmax(kk) ];
    xmin(kk) = -Inf;
    xmax(kk) = Inf;
end

%%-----  run opf  -----
%% build Jacobian and Hessian structure
nA = size(A, 1);                %% number of original linear constraints
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
Cf = sparse(1:nl, f, ones(nl, 1), nl, nb);      %% connection matrix for line & from buses
Ct = sparse(1:nl, t, ones(nl, 1), nl, nb);      %% connection matrix for line & to buses
Cl = Cf + Ct;                                   %% for each line - from & to 
Cb = Cl' * Cl + speye(nb);                      %% for each bus - contains adjacent buses
Cl2 = Cl(il, :);                                %% branches with active flow limit
Cg = sparse(gen(:, GEN_BUS), (1:ng)', 1, nb, ng); %%locations where each gen. resides
nz = nx - 2*(nb+ng);
nxtra = nx - 2*nb;

% Jacobian for the power flow constraints:
%      | dPf_da  dPf_dV dPf_dP  dPf_dQ |   | Cb   Cb  Cg 0| 
% Js = |                               | = |              |
%      | dQf_da  dQf_dV dPf_dP  dPf_dQ |   | Cb   Cb  0 Cg|
%       and line active/reactive power flow constraints
Js = [
    Cb      Cb      Cg              sparse(nb,ng)   sparse(nb,nz); %nz - user variables
    Cb      Cb      sparse(nb,ng)   Cg              sparse(nb,nz);
    Cl2     Cl2     sparse(nl2, 2*ng)               sparse(nl2,nz);
    Cl2     Cl2     sparse(nl2, 2*ng)               sparse(nl2,nz);
    A;
];

% Hessian of lagrangian
% Hs = f(x)_dxx + c(x)_dxx + h(x)_dxx
%            | dPf_daa  dPf_daV dPf_daP  dPf_daQ |   | Cb   Cb  0  0| 
% c(x)_dxx = |                                   | = |              |
%            | dQf_dVa  dQf_dVV dPf_dVP  dPf_dVQ |   | Cb   Cb  0  0|
%
% h(x)_dxx = zeros(nl2,:)
[f, df, d2f] = opf_costfcn(x0, om);
Hs = tril(d2f + [
    Cb  Cb  sparse(nb,nxtra);
    Cb  Cb  sparse(nb,nxtra);
    sparse(nxtra,nx);
]);

%% set options struct for IPOPT
options.ipopt = ipopt_options([], mpopt);

%% extra data to pass to functions
options.auxdata = struct( ...
    'om',       om, ...
    'Ybus',     Ybus, ...
    'Yf',       Yf(il,:), ...
    'Yt',       Yt(il,:), ...
    'mpopt',    mpopt, ...
    'il',       il, ...
    'A',        A, ...
    'nA',       nA, ...
    'neqnln',   2*nb, ...
    'niqnln',   2*nl2, ...
    'Js',       Js, ...
    'Hs',       Hs    );

% %% check Jacobian and Hessian structure
% xr                  = rand(size(x0));
% lambda              = rand(2*nb+2*nl2, 1);
% options.auxdata.Js  = jacobian(xr, options.auxdata);
% options.auxdata.Hs  = tril(hessian(xr, 1, lambda, options.auxdata));
% Js1 = options.auxdata.Js;
% options.auxdata.Js = Js;
% Hs1 = options.auxdata.Hs;
% [i1, j1, s] = find(Js);
% [i2, j2, s] = find(Js1);
% if length(i1) ~= length(i2) || norm(i1-i2) ~= 0 || norm(j1-j2) ~= 0
%     error('something''s wrong with the Jacobian structure');
% end
% [i1, j1, s] = find(Hs);
% [i2, j2, s] = find(Hs1);
% if length(i1) ~= length(i2) || norm(i1-i2) ~= 0 || norm(j1-j2) ~= 0
%     error('something''s wrong with the Hessian structure');
% end

%% define variable and constraint bounds
options.lb = xmin;
options.ub = xmax;
options.cl = [zeros(2*nb, 1);  -Inf(2*nl2, 1); l];
options.cu = [zeros(2*nb, 1); zeros(2*nl2, 1); u+1e-10]; %TODO drosos added eps
%% TODO simulating contingency for line 9, case 9. Delete afterwards!!!!
% options.cu(2*nb+9) = Inf;
% options.cu(2*nb+18) = Inf;

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

%disp('Objective value:');
%f

%% update solution data
Va = x(vv.i1.Va:vv.iN.Va);
Vm = x(vv.i1.Vm:vv.iN.Vm);
Pg = x(vv.i1.Pg:vv.iN.Pg);
Qg = x(vv.i1.Qg:vv.iN.Qg);
V = Vm .* exp(1j*Va);

%%-----  calculate return values  -----
%% update voltages & generator outputs
bus(:, VA) = Va * 180/pi;
bus(:, VM) = Vm;
gen(:, PG) = Pg * baseMVA;
gen(:, QG) = Qg * baseMVA;
gen(:, VG) = Vm(gen(:, GEN_BUS));

%% compute branch flows
Sf = V(branch(:, F_BUS)) .* conj(Yf * V);  %% cplx pwr at "from" bus, p.u.
St = V(branch(:, T_BUS)) .* conj(Yt * V);  %% cplx pwr at "to" bus, p.u.
branch(:, PF) = real(Sf) * baseMVA;
branch(:, QF) = imag(Sf) * baseMVA;
branch(:, PT) = real(St) * baseMVA;
branch(:, QT) = imag(St) * baseMVA;

%% line constraint is typically on square of limit
%% so we must fix multipliers
muSf = zeros(nl, 1);
muSt = zeros(nl, 1);
if ~isempty(il)
    if upper(mpopt.opf.flow_lim(1)) == 'P'
        muSf(il) = info.lambda(2*nb+    (1:nl2));
        muSt(il) = info.lambda(2*nb+nl2+(1:nl2));
    else
        muSf(il) = 2 * info.lambda(2*nb+    (1:nl2)) .* branch(il, RATE_A) / baseMVA;
        muSt(il) = 2 * info.lambda(2*nb+nl2+(1:nl2)) .* branch(il, RATE_A) / baseMVA;
    end
end

%% extract shadow prices for equality var bounds converted to eq constraints
%% (since IPOPT does not return shadow prices on variables that it eliminates)
if nk
    lam_tmp = info.lambda(2*nb+2*nl2+nA-nk+(1:nk));
    kl = find(lam_tmp < 0);             %% lower bound binding
    ku = find(lam_tmp > 0);             %% upper bound binding
    info.zl(kk(kl)) = -lam_tmp(kl);
    info.zu(kk(ku)) =  lam_tmp(ku);

    info.lambda(2*nb+2*nl2+nA-nk+(1:nk)) = [];  %% remove these shadow prices
    nA = nA - nk;                               %% reduce dimension accordingly
end


%% update Lagrange multipliers
bus(:, MU_VMAX)  = info.zu(vv.i1.Vm:vv.iN.Vm);
bus(:, MU_VMIN)  = info.zl(vv.i1.Vm:vv.iN.Vm);
gen(:, MU_PMAX)  = info.zu(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_PMIN)  = info.zl(vv.i1.Pg:vv.iN.Pg) / baseMVA;
gen(:, MU_QMAX)  = info.zu(vv.i1.Qg:vv.iN.Qg) / baseMVA;
gen(:, MU_QMIN)  = info.zl(vv.i1.Qg:vv.iN.Qg) / baseMVA;
bus(:, LAM_P)    = info.lambda(nn.i1.Pmis:nn.iN.Pmis) / baseMVA;
bus(:, LAM_Q)    = info.lambda(nn.i1.Qmis:nn.iN.Qmis) / baseMVA;
branch(:, MU_SF) = muSf / baseMVA;
branch(:, MU_ST) = muSt / baseMVA;

%% package up results
nlnN = getN(om, 'nln');

%% extract multipliers for nonlinear constraints
kl = find(info.lambda(1:2*nb) < 0);
ku = find(info.lambda(1:2*nb) > 0);
nl_mu_l = zeros(nlnN, 1);
nl_mu_u = [zeros(2*nb, 1); muSf; muSt];
nl_mu_l(kl) = -info.lambda(kl);
nl_mu_u(ku) =  info.lambda(ku);

%% extract multipliers for linear constraints
lam_lin = info.lambda(2*nb+2*nl2+(1:nA));   %% lambda for linear constraints
kl = find(lam_lin < 0);                     %% lower bound binding
ku = find(lam_lin > 0);                     %% upper bound binding
mu_l = zeros(nA, 1);
mu_l(kl) = -lam_lin(kl);
mu_u = zeros(nA, 1);
mu_u(ku) = lam_lin(ku);

mu = struct( ...
  'var', struct('l', info.zl, 'u', info.zu), ...
  'nln', struct('l', nl_mu_l, 'u', nl_mu_u), ...
  'lin', struct('l', mu_l, 'u', mu_u) );

results = mpc;
[results.bus, results.branch, results.gen, ...
    results.om, results.x, results.mu, results.f] = ...
        deal(bus, branch, gen, om, x, mu, f);

pimul = [ ...
  results.mu.nln.l - results.mu.nln.u;
  results.mu.lin.l - results.mu.lin.u;
  -ones(ny>0, 1);
  results.mu.var.l - results.mu.var.u;
];
raw = struct('xr', x, 'pimul', pimul, 'info', info.status, 'output', output);


%-----  callback functions  -----
function f = objective(x, d)
xnew = x;
% mpc = get_mpc(d.om);
% xnew = check_ramps(x,mpc);
f = opf_costfcn(xnew, d.om);

function df = gradient(x, d)
xnew = x;
% mpc = get_mpc(d.om);
% xnew = check_ramps(x,mpc);
[f, df] = opf_costfcn(xnew, d.om);

function c = constraints(x, d)
xnew = x; 
% mpc = get_mpc(d.om);
% xnew = check_ramps(x,mpc);

% xnew = 1:size(x,1); %for case9 SCOPF
% xnew = xnew';

%xnew = [1:24]';
[hn, gn] = opf_consfcn(xnew, d.om, d.Ybus, d.Yf, d.Yt, d.mpopt, d.il);
if isempty(d.A)
    c = [gn; hn];
else
    c = [gn; hn; d.A*xnew];
end

function J = jacobian(x, d)
xnew = x;

%case9
%xnew = [[1:9]' ; [1:9]'; [1.300000000050000,1.550000000050000,1.400000000050000,0.5e-10,0.5e-10,0.5e-10]']

%xnew = [0:343]'; %case118
%xnew = [0:3227]'; %case_pegase1354
%xnew = [0:35501]'; %case_pegase13659

% mpc = get_mpc(d.om);
% xnew = check_ramps(x,mpc);

%xnew = [1:24]'; %SCOPF, [Jbase(:,1:18) zeros(36,18) Jbase(:,19:end); zeros(34,18) Jcont8]

[hn, gn, dhn, dgn] = opf_consfcn(xnew, d.om, d.Ybus, d.Yf, d.Yt, d.mpopt, d.il);
J = [dgn'; dhn'; d.A];


function H = hessian(x, sigma, lambda, d)
xnew = x;
% mpc = get_mpc(d.om);
% xnew = check_ramps(x,mpc);
lam.eqnonlin   = lambda(1:d.neqnln);
lam.ineqnonlin = lambda(d.neqnln+(1:d.niqnln));

% reference values for C++ case118
% xnew = [0:343]'; %ones(344,1);
% lam.eqnonlin = [0:235]'; % ones(236,1);
% lam.ineqnonlin = [236:607]'; %ones(372,1);
% sigma = 1;
 
% reference values for C++ case_pegase1354
% xnew = [0:3227]'; %ones(344,1);
% lam.eqnonlin = [0:2707]'; % ones(236,1);
% lam.ineqnonlin = [2708:5571]'; %ones(372,1);
% sigma = 1;

% reference values for C++ case_pegase13659
% xnew = [0:35501]'; %ones(344,1);
% lam.eqnonlin = [0:27317]'; % ones(236,1);
% lam.ineqnonlin = []'; %ones(372,1);
% sigma = 1;

% reference values for C++ case9
%xnew = [[1:9]' ; [1:9]'; [1.300000000050000,1.550000000050000,1.400000000050000,0.5e-10,0.5e-10,0.5e-10]']
%xnew = [1:24]';

% xnew = [1:9 1:9 1:6]';
% 
% lam.eqnonlin = ones(d.neqnln,1);
% 
% lam.ineqnonlin = ones(d.niqnln,1);

H = tril(opf_hessfcn(xnew, lam, sigma, d.om, d.Ybus, d.Yf, d.Yt, d.mpopt, d.il));
%writecsr('/Users/Juraj/Desktop/mtlb_hess.csr',H,1);

% function Js = jacobianstructure(d)
% Js = d.Js;
% 
% function Hs = hessianstructure(d)
% Hs = d.Hs;