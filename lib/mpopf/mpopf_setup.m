function om = mpopf_setup(mpc, mpopf_aux, mpopt)
%OPF  Constructs an MPOPF model object from a MATPOWER case struct.
%   OM = OPF_SETUP(MPC, MPOPF_AUX, MPOPT)
%
%   Assumes that MPC is a MATPOWER case struct with internal indexing,
%   all equipment in-service, etc.
%
%   See also OPF, EXT2INT, OPF_EXECUTE.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% options
dc  = strcmp(upper(mpopt.model), 'DC');
alg = upper(mpopt.opf.ac.solver);
use_vg = mpopt.opf.use_vg;

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

%% define flag to indicate whether we are tied to legacy formulation
%% implemented by legacy MINOS and PDIPM solvers (e.g. with hard-coded
%% costs and constrain
if strcmp(alg, 'MINOPF') || strcmp(alg, 'PDIPM') || ...
        strcmp(alg, 'TRALM') || strcmp(alg, 'SDPOPF')
    error('legacy_formulation not supported in MPOPF');
    legacy_formulation = 1;
else
    legacy_formulation = 0;
end

%% get loading profile
profile = mpopf_aux.profile;

%% get OPF indices
[VAopf, VMopf, PGopf, QGopf] = mpopf_aux.index.getLocalIndicesOPF(mpc);

%% data dimensions
Nt    = length(profile);    %% number of time periods
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections
nnle = 0;                   %% number of nonlinear user-defined equality cons
nnli = 0;                   %% number of nonlinear user-defined inequality cons
if isfield(mpc, 'A')
  nlin = size(mpc.A, 1);    %% number of linear user constraints
else
  nlin = 0;
end
if isfield(mpc, 'N')
  nw = size(mpc.N, 1);      %% number of general cost vars, w
else
  nw = 0;
end

if dc
  error('DC not supported in MPOPF');
else    %% AC
  if use_vg     %% adjust bus voltage limits based on generator Vg setpoint
    %% gen connection matrix, element i, j is 1 if, generator j at bus i is ON
    Cg = sparse(mpc.gen(:, GEN_BUS), (1:ng)', mpc.gen(:, GEN_STATUS) > 0, nb, ng);
    Vbg = Cg * sparse(1:ng, 1:ng, mpc.gen(:, VG), ng, ng);
    Vmax = max(Vbg, [], 2); %% zero for non-gen buses, else max Vg of gens @ bus
    ib = find(Vmax);                %% buses with online gens
    Vmin = max(2*Cg - Vbg, [], 2);  %% same as Vmax, except min Vg of gens @ bus
    Vmin(ib) = 2 - Vmin(ib);

    if use_vg == 1      %% use Vg setpoint directly
        mpc.bus(ib, VMAX) = Vmax(ib);   %% max set by max Vg @ bus
        mpc.bus(ib, VMIN) = Vmin(ib);   %% min set by min Vg @ bus
    elseif use_vg > 0 && use_vg < 1     %% fractional value
        %% use weighted avg between original Vmin/Vmax limits and Vg
        mpc.bus(ib, VMAX) = (1-use_vg) * mpc.bus(ib, VMAX) + use_vg * Vmax(ib);
        mpc.bus(ib, VMIN) = (1-use_vg) * mpc.bus(ib, VMIN) + use_vg * Vmin(ib);
    else
        error('opf_setup: option ''opf.use_vg'' (= %g) cannot be negative or greater than 1', use_vg);
    end
  end
  if isfield(mpc, 'user_constraints')
    if isfield(mpc.user_constraints, 'nle')
      for k = 1:length(mpc.user_constraints.nle)
        nnle = nnle + mpc.user_constraints.nle{k}{2};
      end
    end
    if isfield(mpc.user_constraints, 'nli')
      for k = 1:length(mpc.user_constraints.nli)
        nnli = nnli + mpc.user_constraints.nli{k}{2};
      end
    end
  end
end

%% convert single-block piecewise-linear costs into linear polynomial cost
pwl1 = find(mpc.gencost(:, MODEL) == PW_LINEAR & mpc.gencost(:, NCOST) == 2);
% p1 = [];
if ~isempty(pwl1)
  x0 = mpc.gencost(pwl1, COST);
  y0 = mpc.gencost(pwl1, COST+1);
  x1 = mpc.gencost(pwl1, COST+2);
  y1 = mpc.gencost(pwl1, COST+3);
  m = (y1 - y0) ./ (x1 - x0);
  b = y0 - m .* x0;
  mpc.gencost(pwl1, MODEL) = POLYNOMIAL;
  mpc.gencost(pwl1, NCOST) = 2;
  mpc.gencost(pwl1, COST:COST+1) = [m b];
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
Pg   = gen(:, PG) / baseMVA;
Pmin = gen(:, PMIN) / baseMVA;
Pmax = gen(:, PMAX) / baseMVA;
if ~dc
  Vm   = bus(:, VM);
  Qg   = gen(:, QG) / baseMVA;
  Qmin = gen(:, QMIN) / baseMVA;
  Qmax = gen(:, QMAX) / baseMVA;
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

if dc               %% DC model
  error('DC not supported in MPOPF');
else                %% AC model
  %% more problem dimensions
  nv    = nb;           %% number of voltage magnitude vars
  nq    = ng;           %% number of Qg vars
  q1    = 1+ng;         %% index of 1st Qg column in Ay

  %% find branches with flow limits
  il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
  nl2 = length(il);         %% number of constrained lines

  %% build admittance matrices
  [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

  %% dispatchable load, constant power factor constraints
  [Avl, lvl, uvl]  = makeAvl(baseMVA, gen);
  if(size(Avl,1) > 0)
     error('Not supported in MPOPF: dispatchable load, constant power factor constraints'); 
  end
  
  %% generator PQ capability curve constraints
  [Apqh, ubpqh, Apql, ubpql, Apqdata] = makeApq(baseMVA, gen);
  if(size(Apqh,1) > 0 || size(Apql,1) > 0)
     error('Not supported in MPOPF: generator PQ capability curve constraints'); 
  end

  user_vars = {'Va', 'Vm', 'Pg', 'Qg'};
  ycon_vars = {'Pg', 'Qg', 'y'};
  
  %% nonlinear cost functions
  if ~isempty(ip3)
    cost_Pg = @(x)opf_gen_cost_fcn(x, baseMVA, pcost, ip3, mpopt);
  end
  if ~isempty(qcost) && ~isempty(iq3)
    cost_Qg = @(x)opf_gen_cost_fcn(x, baseMVA, qcost, iq3, mpopt);
  end
end

%% voltage angle reference constraints
Vau = Inf(nb, 1);
Val = -Vau;
Vau(refs) = Va(refs);
Val(refs) = Va(refs);

%% branch voltage angle difference limits
[Aang, lang, uang, iang]  = makeAang(baseMVA, branch, nb, mpopt);
if(size(Aang,1) > 0)
 error('Not supported in MPOPF: branch voltage angle difference limits'); 
end

%% basin constraints for piece-wise linear gen cost variables
if (strcmp(alg, 'PDIPM') && mpopt.pdipm.step_control) || strcmp(alg, 'TRALM')
  %% SC-PDIPM or TRALM, no CCV cost vars
  ny = 0;
  Ay = sparse(0, ng+nq);
  by =[];
else
  ipwl = find(gencost(:, MODEL) == PW_LINEAR);  %% piece-wise linear costs
  ny = size(ipwl, 1);   %% number of piece-wise linear cost vars
  [Ay, by] = makeAy(baseMVA, ng, gencost, 1, q1, 1+ng+nq);
end
if any(gencost(:, MODEL) ~= POLYNOMIAL & gencost(:, MODEL) ~= PW_LINEAR)
    error('opf_setup: some generator cost rows have invalid MODEL value');
end

%% more problem dimensions
nxOPF = (nb+nv + ng+nq);
nx    = Nt * nxOPF;  %% number of standard MPOPF control variables
if nlin
  nz = size(mpc.A, 2) - nx; %% number of user z variables
  if (nz < 0 && size(mpc.A, 2) ~= nxOPF) %% Storate A is supplied in OPF format
    error('opf_setup: user supplied A matrix must have at least %d columns.', nx);
  end
else
  nz = 0;               %% number of user z variables
  if nw                 %% still need to check number of columns of N
    if size(mpc.N, 2) ~= nx
      error('opf_setup: user supplied N matrix must have %d columns.', nx);
    end
  end
end

%% construct MPOPF model object
om = opf_model(mpc);
if ~isempty(pwl1)
  om.userdata.pwl1 = pwl1;
end
if dc
  error('DC not supported in MPOPF');
else

  %% user data
  om.userdata.Apqdata = Apqdata;
  om.userdata.iang = iang;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%---            MPOPF ordering             ---%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%---         Variables ordering            ---%%
  % Va1 Va2 ... Vm1 Vm2 ... Pg1 Pg2 ... Qg1 Qg2   %
  %%---        Equality cons ordering         ---%%
  % PFre1 PFre2 ... PFim1 PFim2 ...               %  
  %%---       Inequality cons ordering        ---%%
  % BFf1 BFf2 ... BFt1 BFt2 ...                   %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% optimization variables
  for i = 1:Nt
    om.add_var(strcat('Va',int2str(i)), nb, Va, Val, Vau);
  end
  for i = 1:Nt
    om.add_var(strcat('Vm',int2str(i)), nb, Vm, bus(:, VMIN), bus(:, VMAX));
  end
  for i = 1:Nt
    om.add_var(strcat('Pg',int2str(i)), ng, Pg, Pmin, Pmax);
  end
  for i = 1:Nt
    om.add_var(strcat('Qg',int2str(i)), ng, Qg, Qmin, Qmax);
  end

%% nonlinear constraints - mismatch

%Eval the mismatch eqs for the whole time horizon with scaled PD, QD
%The function iterates for each timeperiod and puts the resulting matrices
%in the global Jacobian/Hessian with proper offsets
fcn_mis = @(x, time_period)mpopf_power_balance_fcn(x, time_period, mpc, mpopf_aux, Ybus, mpopt);
hess_mis = @(x, lam, time_period)mpopf_power_balance_hess(x, lam, time_period, mpc, mpopf_aux, Ybus, mpopt);

isEq = 1;
om.add_nln_constraint({'Pmis','Qmis'}, [Nt*nb; Nt*nb], isEq, fcn_mis, hess_mis);

%% nonlinear constraints - branch flow limits 

%eval the flow ineqs for the whole time horizon
%The function iterates for each timeperiod and puts the resulting matrices
%in the global Jacobian/Hessian with proper offsets
fcn_flow = @(x, time_period)mpopf_branch_flow_fcn(x, time_period, mpc, mpopf_aux, Yf(il, :), Yt(il, :), il, mpopt);
hess_flow = @(x, lam, time_period)mpopf_branch_flow_hess(x, lam, time_period, mpc, mpopf_aux, Yf(il, :), Yt(il, :), il, mpopt);  

isEq = 0;
om.add_nln_constraint({'Sf','St'}, [Nt*nl;Nt*nl], isEq, fcn_flow, hess_flow);

%% linear constraints
%     om.add_lin_constraint('PQh', Apqh, [], ubpqh, {'Pg', 'Qg'});      %% npqh
%     om.add_lin_constraint('PQl', Apql, [], ubpql, {'Pg', 'Qg'});      %% npql
%     om.add_lin_constraint('vl',  Avl, lvl, uvl,   {'Pg', 'Qg'});      %% nvl
%     om.add_lin_constraint('ang', Aang, lang, uang, {'Va'});           %% nang

%% polynomial generator costs
  for i = 1:Nt
    if ~legacy_formulation
        %% quadratic/linear generator costs
        if ~isempty(cpg)
          om.add_quad_cost(strcat('polPg',int2str(i)), Qpg, cpg, kpg, {strcat('Pg',int2str(i))});
        end
        if ~isempty(cqg)
          om.add_quad_cost(strcat('polQg',int2str(i)), Qqg, cqg, kqg, {strcat('Qg',int2str(i))});
        end

        %% higher order polynomial generator costs
        if ~isempty(ip3)
          om.add_nln_cost(strcat('polPg',int2str(i)), 1, cost_Pg, {strcat('Pg',int2str(i))});
        end
        if ~isempty(qcost) && ~isempty(iq3)
          om.add_nln_cost(strcat('polQg',int2str(i)), 1, cost_Qg, {strcat('Qg',int2str(i))});
        end
    end
  end
end

%% y vars, constraints for piece-wise linear gen costs
if ny > 0
  error('y vars not supported in MPOPF');
  om.add_var('y', ny);
  om.add_lin_constraint('ycon', Ay, [], by, ycon_vars);             %% ncony
  if dc || ~legacy_formulation
    om.add_quad_cost('pwl', [], ones(ny, 1), 0, {'y'});
  end
end

%% add user vars, constraints and costs (as specified via A, ..., N, ...)
if nz > 0
  error('User vars not supported in MPOPF');
  om.add_var('z', nz, z0, zl, zu);
  user_vars{end+1} = 'z';
end
if nlin
  %linear constraints from mpc, storage limits
  % -Einit <  T*P1          < Emax - Einit
  % -Einit <  T*(P1+P2)     < Emax - Einit
  % -Einit <  T*(P1+P2+...) < Emax - Einit
  % where P1 is Pg of [disch charg]
  om.add_lin_constraint('usr', mpc.A, lbu, ubu);         %% nlin
end
if nnle
  for k = 1:length(mpc.user_constraints.nle)
    nlc = mpc.user_constraints.nle{k};
    fcn  = eval(['@(x)' nlc{3} '(x, nlc{6}{:})']);
    hess = eval(['@(x, lam)' nlc{4} '(x, lam, nlc{6}{:})']);
    om.add_nln_constraint(nlc{1:2}, 1, fcn, hess, nlc{5});
  end
end
if nnli
  for k = 1:length(mpc.user_constraints.nli)
    nlc = mpc.user_constraints.nli{k};
    fcn  = eval(['@(x)' nlc{3} '(x, nlc{6}{:})']);
    hess = eval(['@(x, lam)' nlc{4} '(x, lam, nlc{6}{:})']);
    om.add_nln_constraint(nlc{1:2}, 0, fcn, hess, nlc{5});
  end
end
if nw
  user_cost.N = mpc.N;
  user_cost.Cw = Cw;
  if ~isempty(fparm)
    user_cost.dd = fparm(:, 1);
    user_cost.rh = fparm(:, 2);
    user_cost.kk = fparm(:, 3);
    user_cost.mm = fparm(:, 4);
  end
  if ~isempty(H)
    user_cost.H = H;
  end
  om.add_legacy_cost('usr', user_cost, user_vars);
end

%% execute userfcn callbacks for 'formulation' stage
om = run_userfcn(userfcn, 'formulation', om, mpopt);

%% implement legacy user costs using quadratic or general non-linear costs
cp = om.params_legacy_cost();   %% construct/fetch the parameters
[N, H, Cw, rh, mm] = deal(cp.N, cp.H, cp.Cw, cp.rh, cp.mm);
[nw, nx] = size(N);
if nw
    if any(cp.dd ~= 1) || any(cp.kk)    %% not simple quadratic form
        if dc                           %% (includes "dead zone" or
            if any(cp.dd ~= 1)          %%  quadratic "penalty")
                error('opf_setup: DC OPF can only handle legacy user-defined costs with d = 1');
            end
            if any(cp.kk)
                error('opf_setup: DC OPF can only handle legacy user-defined costs with no "dead zone", i.e. k = 0');
            end
        elseif ~legacy_formulation
            %% use general nonlinear cost to implement legacy user cost
            legacy_cost_fcn = @(x)opf_legacy_user_cost_fcn(x, cp);
            om.add_nln_cost('usr', 1, legacy_cost_fcn);
        end
    else                                %% simple quadratic form
        %% use a quadratic cost to implement legacy user cost
        if dc || ~legacy_formulation
            %% f = 1/2 * w'*H*w + Cw'*w, where w = diag(mm)*(N*x - rh)
            %% Let: MN = diag(mm)*N
            %%      MR = M * rh
            %%      HMR  = H  * MR;
            %%      HtMR = H' * MR;
            %%  =>   w = MN*x - MR
            %% f = 1/2 * (MN*x - MR)'*H*(MN*x - MR) + Cw'*(MN*x - MR)
            %%   = 1/2 * x'*MN'*H*MN*x +
            %%          (Cw'*MN - 1/2 * MR'*(H+H')*MN)*x +
            %%          1/2 * MR'*H*MR - Cw'*MR
            %%   = 1/2 * x'*Q*w + c'*x + k
    
            [N, H, Cw, rh, mm] = deal(cp.N, cp.H, cp.Cw, cp.rh, cp.mm);
            nw = size(N, 1);            %% number of general cost vars, w
            M    = sparse(1:nw, 1:nw, mm, nw, nw);
            MN   = M * N;
            MR   = M * rh;
            HMR  = H  * MR;
            HtMR = H' * MR;
            Q = MN' * H * MN;
            c = full(MN' * (Cw - 1/2*(HMR+HtMR)));
            k = (1/2 * HtMR - Cw)' * MR;
            om.add_quad_cost('usr', Q, c, k);
        end
    end
end