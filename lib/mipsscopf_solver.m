function [results, success, raw] = mipsscopf_solver(om, model, mpopt)
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

cont = model.cont;

%% problem dimensions
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of branches
ng = size(gen, 1);          %% number of generators
ns = size(cont, 1);         %% number of scenarios (nominal + ncont)

% get indices of REF gen and of REF/PV buses
[REFgen_idx, nREFgen_idx] = model.index.getREFgens(mpc);
[REFbus_idx,nREFbus_idx] = model.index.getXbuses(mpc,3);%3==REF
[PVbus_idx, nPVbus_idx] = model.index.getXbuses(mpc,2);%2==PV

% indices of local OPF solution vector x = [VA VM PG QG]
[VAscopf, VMscopf, PGscopf, QGscopf] = model.index.getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = model.index.getLocalIndicesOPF(mpc);

%% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% linear constraints
[A, l, u] = linear_constraints(om);

%% bounds on optimization vars
[x0, xmin, xmax] = getv(om); %returns standard OPF form [Va Vm Pg Qg]

% Note that variables with equal upper and lower bounds are removed by IPOPT
% so we add small perturbation to x_u[], we don't want them removed
% because the Schur solver assumes particular structure that would
% be changed by removing variables.
% Exept for the Va at the refernece bus which we want to remove.
idx = find(xmin == xmax);
xmax(idx) = xmax(idx) + 1e-10;
xmax(REFbus_idx) = xmin(REFbus_idx);

% replicate bounds for all scenarios and append global limits
xl = xmin([VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)]); %local variables
xg = xmin([VMopf(PVbus_idx) PGopf(nREFgen_idx)]); %global variables
xmin = [repmat(xl, [ns, 1]); xg];

xl = xmax([VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)]); %local variables
xg = xmax([VMopf(PVbus_idx) PGopf(nREFgen_idx)]); %global variables
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
    Varefs = bus(REFbus_idx, VA) * (pi/180);
    for i = 0:ns-1
        idx = model.index.getGlobalIndices(mpc, ns, i);
        x0(idx(VAscopf)) = Varefs(1);
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

%% -----  run opf  -----
auxdata = struct( ...
    'om',       om, ...
    'cont',     cont, ...
    'index',    model.index, ...
    'il',       il, ...
    'mpopt',    mpopt );

f_fcn = @(x) objective_fcn(x, auxdata);
gh_fcn = @(x) constraints_fcn(x, auxdata);
hess_fcn = @(x, lambda, cost_mult) hessian_fcn(x, lambda, cost_mult, auxdata);

[x, f, info, output, Lambda] = ...
  mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn, opt, mpc);
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
meta.Ybus = Ybus;
meta.Yf = Yf;
meta.Yt = Yt;
meta.lb = xmin;
meta.ub = xmax;
meta.A = A;
    
raw = struct('info', info, 'meta', meta);
results = struct('f', f, 'x', x);

%% callback routines
%evaluate objective and its gradient
function [f, df] = objective_fcn(x, d) 
    mpc = get_mpc(d.om);
    ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)

    % use nominal case to evaluate cost fcn (only pg/qg are relevant)
    idx_nom = d.index.getGlobalIndices(mpc, ns, 0);
    [VAscopf, VMscopf, PGscopf, QGscopf] = d.index.getLocalIndicesSCOPF(mpc);
    [VAopf, VMopf, PGopf, QGopf] = d.index.getLocalIndicesOPF(mpc);

    [f, dfOPF] = opf_costfcn(x(idx_nom([VAscopf VMscopf PGscopf QGscopf])), d.om);

    df = zeros(size(x,1),1);
    df(idx_nom(PGscopf)) = dfOPF(PGopf); %nonzero only nominal case Pg

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

    gn = zeros(ns*NCONSTR_dgn,1);
    hn = zeros(ns*NCONSTR_dhn,1);
    
    % get indices of REF gen and PV bus
    [REFgen_idx, nREFgen_idx] = d.index.getREFgens(mpc);
    [PVbus_idx, nPVbus_idx] = d.index.getXbuses(mpc,2);%2==PV
    
    [VAscopf, VMscopf, PGscopf, QGscopf] = d.index.getLocalIndicesSCOPF(mpc);
    [VAopf, VMopf, PGopf, QGopf] = d.index.getLocalIndicesOPF(mpc);

    for i = 0:ns-1
        cont = d.cont(i+1);
        idx = d.index.getGlobalIndices(mpc, ns, i);
        [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);
        [hn_local, gn_local, dhn_local, dgn_local] = opf_consfcn(x(idx([VAscopf VMscopf PGscopf QGscopf])), d.om, Ybus, Yf, Yt, d.mpopt, d.il);
        
        gn(i*NCONSTR_dgn + (1:NCONSTR_dgn),1) = gn_local;
        hn(i*NCONSTR_dhn + (1:NCONSTR_dhn)) = hn_local;
        
         %jacobian wrt local variables
        dgn(idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx)]), i*NCONSTR_dgn + (1:NCONSTR_dgn)) = dgn_local([VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)],:);
        %jacobian wrt global variables
        dgn(idx([VMscopf(PVbus_idx) PGscopf(nREFgen_idx)]), i*NCONSTR_dgn + (1:NCONSTR_dgn)) = dgn_local([VMopf(PVbus_idx) PGopf(nREFgen_idx)],:);
        
        %jacobian wrt local variables
        dhn(idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx)]), i*NCONSTR_dhn + (1:NCONSTR_dhn)) = dhn_local([VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)],:);
        %jacobian wrt global variables
        dhn(idx([VMscopf(PVbus_idx) PGscopf(nREFgen_idx)]), i*NCONSTR_dhn + (1:NCONSTR_dhn)) = dhn_local([VMopf(PVbus_idx) PGopf(nREFgen_idx)], :);
    end
                
%evaluate hessian of Lagrangian    
function H = hessian_fcn(x, lambda, sigma, d)
    mpc = get_mpc(d.om);
    nb = size(mpc.bus, 1);          %% number of buses
    ng = size(mpc.gen, 1);          %% number of gens
    nl = size(mpc.branch, 1);       %% number of branches
    ns = size(d.cont, 1);           %% number of scenarios (nominal + ncont)
    NCONSTR = 2*nb + 2*nl;

    H = sparse(size(x,1), size(x,1));
    % get indices of REF gen and PV bus
    [REFgen_idx, nREFgen_idx] = d.index.getREFgens(mpc);
    [PVbus_idx, nPVbus_idx] = d.index.getXbuses(mpc,2);%2==PV

    [VAscopf, VMscopf, PGscopf, QGscopf] = d.index.getLocalIndicesSCOPF(mpc);
    [VAopf, VMopf, PGopf, QGopf] = d.index.getLocalIndicesOPF(mpc);

    for i = 0:ns-1
        %compute local indices and its parts
        idx = d.index.getGlobalIndices(mpc, ns, i);

        cont = d.cont(i+1);
        [Ybus, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch, cont);

        lam.eqnonlin   = lambda.eqnonlin(i*2*nb + (1:2*nb));
        lam.ineqnonlin = lambda.ineqnonlin(i*2*nl + (1:2*nl));
        H_local = opf_hessfcn(x(idx([VAscopf VMscopf PGscopf QGscopf])), lam, sigma, d.om, Ybus, Yf, Yt, d.mpopt, d.il);

        % H_ll (PG_ref relevant only in nominal case, added to global part)
         H(idx([VAscopf VMscopf(nPVbus_idx) QGscopf]), idx([VAscopf VMscopf(nPVbus_idx) QGscopf])) =...
                H_local([VAopf VMopf(nPVbus_idx) QGopf], [VAopf VMopf(nPVbus_idx) QGopf]);

        % H_lg and H_gl (PG parts are implicitly zero, could leave them out)
        H(idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx)]), idx([VMscopf(PVbus_idx) PGscopf(nREFgen_idx)])) = ...
                H_local([VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)], [VMopf(PVbus_idx) PGopf(nREFgen_idx)]);
        H(idx([VMscopf(PVbus_idx) PGscopf(nREFgen_idx)]), idx([VAscopf VMscopf(nPVbus_idx) QGscopf PGscopf(REFgen_idx)])) = ...
                H_local([VMopf(PVbus_idx) PGopf(nREFgen_idx)], [VAopf VMopf(nPVbus_idx) QGopf PGopf(REFgen_idx)]);

        % H_gg hessian w.r.t global variables (and PG_ref_0) 
        if i == 0
            % H_pg at non-reference gens, these are global variables
            H(idx([PGscopf(nREFgen_idx)]), idx([PGscopf(nREFgen_idx)])) = ...
                H_local([PGopf(nREFgen_idx)], [PGopf(nREFgen_idx)]);

            % H_pgref is local variable for nominal scenario, but used in f()
            H(idx([PGscopf(REFgen_idx)]), idx([PGscopf(REFgen_idx)])) = ...
                H_local([PGopf(REFgen_idx)], [PGopf(REFgen_idx)]);
        end

        %each scenario contributes to hessian w.r.t global VM variables at PV buses
        H(idx([VMscopf(PVbus_idx)]), idx([VMscopf(PVbus_idx)])) = ...
            H(idx([VMscopf(PVbus_idx)]), idx([VMscopf(PVbus_idx)])) + ...
            H_local([VMopf(PVbus_idx)], [VMopf(PVbus_idx)]);
    end

    H = tril(H);