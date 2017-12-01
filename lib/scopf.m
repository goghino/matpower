function [results, success] = ...
    scopf(mpc, cont, mpopt, tol)
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
ns = size(cont, 1);         %% number of scenarios (nominal + ncont)

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

%%----- analyze given list of contingencies or create them automatically
if ns >= 1 && cont(1) > 0
    cont_filt = cont;
    
    %%-- remove contingencies outside bounds of the mpc.branch
    cont_filt(cont_filt > size(mpc.branch,1)) = [];
    N = length(cont_filt);
    
    %TODO - add warning when removing contingency from given list!!
else
    %cont specifies number of contingencies to create
    N = abs(cont);
    
    % generate ns-1 contingencies
    nbranch = size(mpc.branch,1);
    cont_filt = [1:nbranch]';
end

%%-- remove branches causing islands or isolated buses
lines = findIslandBranches(mpc);
cont_filt = setdiff(cont_filt, lines);

%%-- remove duplicate lines
duplicates = findDuplicateBranches(mpc);
cont_filt = setdiff(cont_filt, duplicates);

%%-- remove lines causing Qg violations
critical = findQgCritical(mpc, 150); %max 150% violation wrt Qmax, Qmin
cont_filt = setdiff(cont_filt, critical);

cont = [-1; cont_filt(1:min(N,length(cont_filt)))]; %add nominal case and N contingencies
ns = length(cont);

%%-----  construct OPF model object  -----
om = scopf_setup(mpc, mpopt);

%%----- build scopf model -----
%pass index functions to solvers in order to properly construct x and evaluate callbacks
index = struct('getGlobalIndices', @getGlobalIndices, ...
               'getLocalIndicesSCOPF', @getLocalIndicesSCOPF, ...
               'getLocalIndicesOPF', @getLocalIndicesOPF, ...
               'getXbuses', @getXbuses, ...
               'getREFgens', @getREFgens);
           
scopf_aux = struct('cont', cont, 'index', index);

%%-----  execute the SCOPF  -----
t0 = clock;         %% start timer

[results, success, raw] = scopf_execute(om, scopf_aux, mpopt);
info = raw.info;

et = etime(clock, t0);      %% compute elapsed time

%% verify feasibility of the results 

[VAscopf, VMscopf, PGscopf, QGscopf] = getLocalIndicesSCOPF(mpc);
[VAopf, VMopf, PGopf, QGopf] = getLocalIndicesOPF(mpc);

[REFbus_idx,nREFbus_idx] = getXbuses(mpc,3);%3==REF

TOL_EQ = 1e-5;
TOL_LIN = 1e-6;

if nargin==4
    TOL_EQ = tol;
    TOL_LIN = tol;
end

%verigy feasibility and check bounds
if (mpopt.verbose >= 1)
    fprintf('\n--------------------------------------------------------------\n');
    fprintf('--------------------------------------------------------------\n');
    fprintf('Verifing feasibility of the SCOPF solution with tolerance %e.\n', TOL_EQ);
    errors = 0;
    for i = 1:ns
        fprintf('\tscenario %d ...\n', i-1);

        %get contingency
        c = cont(i);

        %compute local indices and its parts
        idx = getGlobalIndices(mpc, ns, i-1);

        %extract local solution vector [Va Vm Pg Qg]
        xl = results.x(idx([VAscopf VMscopf PGscopf QGscopf]));

        %update admittance matrices
        [Ybus, Yf, Yt] = updateYbus(mpc.branch, raw.meta.Ybus, raw.meta.Yf, raw.meta.Yt, c);

        %check power generation bounds l < [Pg Qg] < u
        err = find(xl([PGopf QGopf]) < raw.meta.lb(idx([PGscopf QGscopf])));
        if (~isempty(err))
            fprintf('violated %d lower Pg/Qg limits %e\n', length(err), max(raw.meta.lb(idx([PGscopf QGscopf])) - xl([PGopf QGopf])));
            errors = errors + 1;
        end
        err = find(xl([PGopf QGopf]) > raw.meta.ub(idx([PGscopf QGscopf])));
        if (~isempty(err))
            fprintf('violated %d upper Pg/Qg limits %e\n', length(err), max(xl([PGopf QGopf]) - raw.meta.ub(idx([PGscopf QGscopf]))));
            errors = errors + 1;
        end

        %bus voltages magnitude bounds p.u. l < Vm < u
        err = find(xl(VMopf) < raw.meta.lb(idx(VMscopf)));
        if (~isempty(err))
            fprintf('violated %d lower Vm limits %e\n', length(err), max(raw.meta.lb(idx(VMscopf)) - xl(VMopf)));
            errors = errors + 1;
        end

        err = find(xl(VMopf) > raw.meta.ub(idx(VMscopf)));
        if (~isempty(err))
            fprintf('violated %d upper Vm limits %e\n', length(err), max(xl(VMopf) - raw.meta.ub(idx(VMscopf))));
            errors = errors + 1;
        end

        %bus voltage angles limits only reference bus l = Va = u
        err_lb = abs(xl(VAopf(REFbus_idx)) - raw.meta.lb(idx(VAscopf(REFbus_idx))));
        err_ub = abs(xl(VAopf(REFbus_idx)) - raw.meta.ub(idx(VAscopf(REFbus_idx))));
        if (err_lb > TOL_EQ)
            fprintf('violated lower Va limit on reference bus %e\n', err_lb);
            errors = errors + 1;
        end
        if (err_ub > TOL_EQ)
            fprintf('violated upper Va limit on reference bus %e\n', err_ub);
            errors = errors + 1;
        end

        %evaluate power flows equations and branch power flows g(x), h(x)
        [hn_local, gn_local] = opf_consfcn(xl, om, Ybus, Yf, Yt, mpopt);

        %g(x) = 0, g(x) = V .* conj(Ybus * V) - Sbus;
        err = find(abs(gn_local) > TOL_EQ);
        if (~isempty(err))
            fprintf('violated %d PF equations with max %e\n', length(err), max(abs(gn_local(err))));
            errors = errors + 1;
        end

        %h(x) <= 0, h(x) = Sf .* conj(Sf) - flow_max.^2
        %h(x) for lines with contingency is (- flow_max.^2) which satisfy limits implicitly
        err = find(hn_local > 0);
        if(not(isempty(err)))
            fprintf('violated %d branch power flow limits %e\n', length(err), max(abs(hn_local(err))));
            errors = errors + 1;
        end

        %linear constraints
        if (~isempty(raw.meta.A))
            lin_constr = raw.meta.A * results.x;
            err = find(abs(lin_constr) > TOL_LIN);
            if(not(isempty(err)))
                fprintf('violated %d linear constraints %e\n', length(err), max(abs(lin_constr(err))));
                errors = errors + 1;
            end
        end

        if (errors == 0)
            fprintf('\tPASSED\n');
        else
            fprintf('\tFAILED with %d errors\n', errors);
            errors = 0;
        end
    end
    fprintf('\n--------------------------------------------------------------\n');
    fprintf('--------------------------------------------------------------\n');
end

%% -----  DO NOT revert to original ordering, we are returnting SCOPF solution, not OPF  -----
%results = int2ext(results);

%% prepare output

if nargout <= 2
    results.et = et;
    results.success = success;
    results.raw = raw;
else
    error('scopf.m: Incorrect number of output arguments\n')
end

%% -----  helper functions  ----- 
function idx = getGlobalIndices(mpc, ns, i)
% returns indices of local OPF variables of sceanrio i in vector x_ipopt
% OPF variables are ordered local first, global variables then: [Va Vm Qg Pg_ref] [Pg]
% scenarios i are indexed 0..NS-1
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches

[PVbus_idx, nPVbus_idx] = getXbuses(mpc,2);%2==PV
nPV = size(PVbus_idx,1); %number of PV buses
nnPV = nb - nPV; %number of non PV buses

nPart = nb+nnPV+ng+1; %number of local variables for each scenario [Va Vm_npv Qg Pg_ref]

li1 = i*nPart + (1:nb); %indices of local variables [Va] of scenario i
li2 = i*nPart + nb + (1:nnPV); %indices of local variables [Vm not at BUS_PV] of scenario i
li3 = i*nPart + nb + nnPV + (1:ng); %indices of local variables [Qg] of scenario i
li4 = i*nPart + nb + nnPV + ng + 1; %index of local variable [Pg] at BUS_ref
gi1 = ns*nPart + (1:nPV); %indices of global variables [Vm] Vm at BUS_PV
gi2 = ns*nPart + nPV + (1:ng-1); %indices of global variables [Pg], Pg not at BUS_ref

idx = [li1 li2 li3 li4 gi1 gi2]; %return in order [Va Vm_npv Qg Pg_ref] [Vm_pv Pg_nref]

function [Xbus_idx, nXbus_idx] = getXbuses(mpc, type)
%returns indices of buses with specified type and its complement to the full bus set
BUS_TYPE = 2;
Xbus_idx = find(mpc.bus(:,BUS_TYPE) == type);
nXbus_idx = find(mpc.bus(:,BUS_TYPE) ~= type);

function [REFgen_idx, nREFgen_idx] = getREFgens(mpc)
%returns indices of generators connected to reference bus and its
%complement to the full generator set
BUS_TYPE = 2;
REF = 3;
GEN_BUS = 1;
REFbus_idx = find(mpc.bus(:,BUS_TYPE) == REF);
REFgen_idx = find(mpc.gen(:,GEN_BUS) == REFbus_idx); %index of gen connected to ref_bus
nREFgen_idx = find(mpc.gen(:,GEN_BUS) ~= REFbus_idx); %index of gens not connected to ref_bus


function [VAi, VMi, PGi, QGi] = getLocalIndicesSCOPF(mpc)
%extracts OPF variables from local SCOPF variables vector x
%usage: x_local = getGlobalIndices(..., ..., scenario_i)
%       x(x_local([VAi VMi PGi QGi]))
nb = size(mpc.bus, 1);          %% number of buses
ng = size(mpc.gen, 1);          %% number of gens
nl = size(mpc.branch, 1);       %% number of branches

[REFgen_idx, nREFgen_idx] = getREFgens(mpc);
[PVbus_idx, nPVbus_idx] = getXbuses(mpc, 2);%2==PV
[REFbus_idx, nREFbus_idx] = getXbuses(mpc,3);%3==ref

nPV = size(PVbus_idx,1); %number of PV buses
nnPV = size(nPVbus_idx,1); %number of non PV buses

nPart = nb+nnPV+ng+1; %number of local variables for each scenario [Va Vm_npv Qg Pg_ref]

%do some tests
if (nPV + nnPV ~= nb)
   error('PV and nonPV buses are not equal to the whole bus set'); 
end
if (length(REFbus_idx) > 1)
   error('multiple REF buses'); 
end

if(length(REFgen_idx) > 1)
    error('reference bus has multiple generators')
end
if(isempty(REFgen_idx))
    error('reference bus has no connected generators')
end

tmp = [find(mpc.bus(mpc.gen(:,1), 2) == 1); %1==PQ
       find(mpc.bus(mpc.gen(:,1), 2) == 4)];%4==isolated
if(~isempty(tmp))
    error('generator conected to non-PV bus (not considering REF buses)')
end

VAi = 1:nb;
VMi = zeros(1,nb); VMi(nPVbus_idx) = nb + (1:nnPV); VMi(PVbus_idx) = nPart + (1:nPV);
PGi = zeros(1,ng); PGi(nREFgen_idx) = nPart+nPV + (1:ng-1); PGi(REFgen_idx) = nb+nnPV+ng + (1);
QGi = nb + nnPV + (1:ng);

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
