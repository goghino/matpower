function [results, success, info] = ...
    scopf(mpc, cont, mpopt)
%SCOPF  Solves an optimal power flow with security constraints.
%   [RESULTS, SUCCESS] = SCOPF(MPC, CONT, MPOPT)
%
%   Returns either a RESULTS struct and an optional SUCCESS flag, or individual
%   data matrices, the objective function value and a SUCCESS flag. In the
%   latter case, there are additional optional return values. See Examples
%   below for the possible calling syntax options.
%
%   Examples:
%       Output argument options:
%
%       [results, success, info] = scopf(mpc, cont, mpopt)
%
%   See also RUNOPF, DCOPF, UOPF, CASEFORMAT.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
t0 = clock;         %% start timer

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

%% add zero columns to bus, gen, branch for multipliers, etc if needed
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections
if size(mpc.bus,2) < MU_VMIN
  mpc.bus = [mpc.bus zeros(nb, MU_VMIN-size(mpc.bus,2)) ];
end
if size(mpc.gen,2) < MU_QMIN
  mpc.gen = [ mpc.gen zeros(ng, MU_QMIN-size(mpc.gen,2)) ];
end
if size(mpc.branch,2) < MU_ANGMAX
  mpc.branch = [ mpc.branch zeros(nl, MU_ANGMAX-size(mpc.branch,2)) ];
end

%%-----  convert to internal numbering, remove out-of-service stuff  -----
mpc = ext2int(mpc);

%%-----  construct OPF model object  -----
om = scopf_setup(mpc, mpopt);

%%----- build scopf model -----
%pass index functions to solvers in order to properly construct x and evaluate callbacks
index = struct('getGlobalIndices', @getGlobalIndices, ...
               'getLocalIndicesSCOPF', @getLocalIndicesSCOPF, ...
               'getLocalIndicesOPF', @getLocalIndicesOPF);
           
scopf_m = struct('cont', cont, 'index', index);

%%-----  execute the OPF  -----
[results, success, raw] = scopf_execute(om, scopf_m, mpopt);
info = raw.info;

%% verify feasibility of the results 
ns = size(cont, 1);         %% number of scenarios (nominal + ncont)

[VAscopf, VMscopf, PGscopf, QGscopf] = getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = getLocalIndicesOPF(mpc);

BUS_TYPE = 2;
REF = 3;
GEN_BUS = 1;
ref_bus = find(mpc.bus(:,BUS_TYPE) == REF);
ref_gen = find(mpc.gen(:,GEN_BUS) == ref_bus); %index of gen connected to ref_bus

%compute bus voltages and check bounds
for i = 1:ns
    %get contingency
    c = cont(i);
    
    %compute local indices and its parts
    idx = getGlobalIndices(mpc, ns, i-1);
        
    %extract local solution vector [Va Vm Pg Qg]
    xl = results.x(idx([VAscopf VMscopf PGscopf QGscopf]));
    xg = xl([PGopf QGopf]); %[Pg Qg]
    
    %update admittance matrices
    [Ybus, Yf, Yt] = updateYbus(mpc.branch, raw.meta.Ybus, raw.meta.Yf, raw.meta.Yt, c);
    
    %check power generation bounds
    err = find(xg < raw.meta.lb(idx([PGscopf QGscopf])));
    if (not(isempty(err)))
        error('violated lower Pg/Qg limit');
    end
    err = find(xg > raw.meta.ub(idx([PGscopf QGscopf])));
    if (not(isempty(err)))
        error('violated upper Pg/Qg limit');
    end
    
    %bus voltages magnitude bounds p.u.
    err = find(xl(VMopf) < raw.meta.lb(idx(VMscopf)));
    if (not(isempty(err)))
        error('violated lower Vm limit %e', max(mpc.bus(err,VMIN) - xl(VMopf(err))));
    end
    
    err = find(xl(VMopf) > raw.meta.ub(idx(VMscopf)));
    if (not(isempty(err)))
        error('violated upper Vm limit %e', max(xl(VMopf(err)) - mpc.bus(err,VMAX)));
    end
    
    %bus voltage angles limits only reference bus Val = Va = Vau
    err_lb = abs(xl(VAopf(ref_bus)) * (pi/180) - raw.meta.lb(idx(VAscopf(ref_bus))));
    err_ub = abs(xl(VAopf(ref_bus)) * (pi/180) - raw.meta.ub(idx(VAscopf(ref_bus))));
    if (err_lb > 1e-5 || err_ub > 1e-5)
        error('violated Va limits on reference bus %e %e', err_lb, err_ub);
    end
    
    %power flows equations and branch power flows
    [hn_local, gn_local] = opf_consfcn(xl, om, Ybus, Yf, Yt, mpopt);
    
    %g(x) = 0, g(x) = V .* conj(Ybus * V) - Sbus;
    err = find(abs(gn_local) > 1e-4);
    if (not(isempty(err)))
        error('violated PF equations %e', max(abs(gn_local(err))));
    end
    
    %h(x) <= 0, h(x) = Sf .* conj(Sf) - flow_max.^2
    %h(x) for lines with contingency is - flow_max.^2 which satisfy limits implicitly
    err = find(hn_local > 0);
    if(not(isempty(err)))
        error('violated branch power flow limits %e', max(abs(hn_local(err))));
    end
    
    %linear constraints
    if ~isempty(raw.meta.A)
        lin_constr = raw.meta.A * results.x;
        err = find(abs(lin_constr) > 2e-1);
        if(not(isempty(err)))
            error('violated linear constraints %e', max(abs(lin_constr(err))));
        end
    end
end

%% -----  DO NOT revert to original ordering, we are returnting SCOPF solution, not OPF  -----
%results = int2ext(results);

%% -----  helper functions  ----- 
function idx = getGlobalIndices(mpc, ns, i)
% returns indices of local OPF variables of sceanrio i in vector x_ipopt
% OPF variables are ordered local first, global variables then: [Va Vm Qg Pg_ref] [Pg]
% scenarios i are indexed 0..NS-1
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches

nPart = 2*nb+ng+1; %number of local variables for each scenario

li1 = i*nPart + (1:2*nb); %indices of local variables [Va Vm] of scenario i
li2 = i*nPart + 2*nb + (1:ng); %indices of local variables [Qg] of scenario i
li3 = i*nPart + 2*nb + ng + 1; %index of local variable Pg_ref at BUS_ref
gi = ns*nPart + (1:ng-1); %indices of global variables [Pg], Pg not at BUS_ref

idx = [li1 li2 li3 gi]; %return in order [Va Vm Qg Pg_ref] [Pg]

function [VAi, VMi, PGi, QGi] = getLocalIndicesSCOPF(mpc)
%extracts OPF variables from local SCOPF variables vector x
%usage: x_local = getGlobalIndices(..., ..., scenario_i)
%       x(x_local([VAi VMi PGi QGi]))
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches

BUS_TYPE = 2;
REF = 3;
GEN_BUS = 1;
ref_bus = find(mpc.bus(:,BUS_TYPE) == REF);
ref_gen = find(mpc.gen(:,GEN_BUS) == ref_bus); %index of gen connected to ref_bus

VAi = 1:nb;
VMi = nb + (1:nb);
PGi = 2*nb + ng + 1 + [(1:ref_gen-1), 0 ,(ref_gen:ng-1)];
QGi = 2*nb + (1:ng);

function [VAi, VMi, PGi, QGi] = getLocalIndicesOPF(mpc)
%extracts variables from OPF variables vector x
%usage: x([VAi VMi PGi QGi])
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches

VAi = 1:nb;
VMi = nb + (1:nb);
PGi = 2*nb + (1:ng);
QGi = 2*nb + ng + (1:ng);
