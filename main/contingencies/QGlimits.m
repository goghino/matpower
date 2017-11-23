clear; close all;
define_constants;

setenv('OMP_NUM_THREADS', '1')
%setenv('IPOPT_WRITE_MAT','1');
addpath('/Users/Juraj/Documents/Code/PowerGrid/matrices/'); %readcsr

[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%% Select and configure the solver
IPOPT = 1;
if (IPOPT == 1)
    %for further options see ipopt.opt
    mpopt = mpoption('opf.ac.solver', 'IPOPT', 'verbose', 2);
else
    mpopt = mpoption('opf.ac.solver', 'MIPS', 'verbose', 3);
    mpopt.mips.max_it  = 100;
    mpopt.mips.max_stepsize = 1e12;
    mpopt.mips.feastol = 1e-4;
    mpopt.mips.gradtol = 1e-2;
    mpopt.mips.comptol = 1e-2;
    mpopt.mips.costtol = 1e-2;
    %mpopt.mips.linsolver='PARDISO';
end

%% load MATPOWER case struct, see help caseformat
%mpc = loadcase('case9');
mpc = loadcase('case118');
%mpc = loadcase('case1354pegase');
%mpc = loadcase('case9241pegase');
%mpc = case231swiss;

%insert missing branch limits
no_limit = find(mpc.branch(:,RATE_A) < 1e-10);
mpc.branch(no_limit,RATE_A) = 1e5;

%% run OPF and extract x*
[MVAbase, bus, gen, gencost, branch, f, success, et] = runopf(mpc, mpopt);
xOPF = [bus(:,VA); bus(:,VM); gen(:,PG); gen(:,QG)];
mpcOPF = mpc;
mpcOPF.bus = bus;
mpcOPF.gen = gen;
mpcOPF.branch = branch;


%% run PF with part of the x* for each cont and analyze x_c
cont = [1:6 8 10:65 67:112 114:132 135:175 178:182 185:186]; %case 118
%cont = [13 51   114   116   118   137   141   155   160   185]; %problematic
critical = 100; %percentage of the QG output w.r.t QG upper bound
problematic_cont = []; %contingencies that violate critical QG limit

% list of corrections, c: bus1 bus2 ... busnc
corrections = sparse(length(cont), size(mpc.bus,1));

for c = cont
    %copy the OPF mpc so that we can make changes to it
    mpc = mpcOPF;
    
    %remove branch i
    if(c > 0)
        mpc.branch(c,BR_STATUS) = 0;
    end
    
    %run PF with OPF solution and modified grid with a contingency c
    mpopt.pf.enforce_q_lims = 2;
    [MVAbase_c, bus_c, gen_c, branch_c, success_c, et_c] = runpf(mpc, mpopt);
    
    %identify QG violation and critical violations
    idx = find(gen_c(:, QG) > gen_c(:, QMAX));
    violations_g = [idx gen_c(idx, QG) gen_c(idx, QMAX) gen_c(idx, QG)./gen_c(idx, QMAX).*100];%gen_idx gen_QG gen_QMAX percentage
    critical_g = find(violations_g(:,4) > critical);
    critical_b = gen_c(idx(critical_g), GEN_BUS);
    
    if ~isempty(critical_g)
       problematic_cont = [problematic_cont c]; 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform iterative change of PV -> PQ buses until there are no violations
    while ~isempty(critical_g)
        
        % 1. set Qg to binding limit
        
        % 2. [Pd, Qd] -= [Pg, Qg] (one at a time) and turn off the generator
        
        % 3. change the type of problematic buses from PV to PQ
        mpc.bus(critical_b, BUS_TYPE) = PQ;
        corrections(c, critical_b) = 1;
        
        % X. What if we overwrite REF bus? [ref, pv, pq] = bustypes(bus, gen);

        %run PF with OPF solution and modified grid with a contingency c and
        %changed bus type in problematic PV buses
        [MVAbase_c, bus_c, gen_c, branch_c, success_c, et_c] = runpf(mpc, mpopt);

        %identify QG violation and critical violations
        idx = find(gen_c(:, QG) > gen_c(:, QMAX));
        violations_g = [idx gen_c(idx, QG) gen_c(idx, QMAX) gen_c(idx, QG)./gen_c(idx, QMAX).*100];%gen_idx gen_QG gen_QMAX percentage
        critical_g = find(violations_g(:,4) > critical);
        critical_b = gen_c(idx(critical_g), GEN_BUS);
    end
    
    % 4. restore injections from limited gens (those at Q limits)
    %    [Pd, Qd] += [Pg, Qg] and turn on the generators
    %    adjust voltage angles to make original ref bus correct
end

spy(corrections)

%%
%problematic_cont200 =
%    13    51   114   116   118   137   141   155   160   185

%%TODO
%1. x-axis contingency (excluding island and representing multiple lines as one)
%2. y-asix violation of Qmin, Qmax, Vm_min/max when using OPF solution with
%   contingent network