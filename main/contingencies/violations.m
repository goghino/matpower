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
mpc = loadcase('case118');
%mpc = loadcase('case1354pegase');
%mpc = loadcase('case9241pegase');

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

%% identify set of feasible contingencies
nbrch = size(mpc.branch,1);
cont = 1:nbrch;
cont_filt = [];

% 1. that do not create islands in the network
fprintf('\n\nFiltering contingencies that create islands or isolated buses...\n');
for i = 1:nbrch

    %get contingency
    c = cont(i);
    
    %remove branch for non-nominal case
    mpc_test = mpc;
    if(c > 0)
        mpc_test.branch(c,:) = [];
    end
    
    %test for islands or isolated buses
    [islands, isolated] = find_islands(mpc_test);
    
    if (~isempty(isolated))
        fprintf('\tRemoving branch %d leaves bus %d isolated\n', c, isolated);
        continue;
    end
    
    if (size(islands,2) > 1)
       fprintf('\tRemoving branch %d creates %d separate islands in the network\n', c, size(islands,2)); 
       continue;
    end
    
    cont_filt = [cont_filt c];
end

% 2. exclude multiple lines, or remove them all
fprintf('\n\nIdentifying duplicate lines...\n');
duplicates = findDuplicateBranches(mpc)
cont_filt = setdiff(cont_filt, duplicates);

fprintf('Considering %f%% of the branches as contingencies.\n',length(cont_filt)/nbrch*100);


%% Analyze the PF violations after the contingency using the X_OPF
QG_max = [];
QG_min = [];
Vm_max = [];
Vm_min = [];
for c = cont_filt
    %copy the OPF mpc so that we can make changes to it
    mpc = mpcOPF;
    
    %remove branch i
    if(c > 0)
        mpc.branch(c,BR_STATUS) = 0;
    end
    
    %run PF with OPF solution and modified grid with a contingency c
    %mpopt.pf.enforce_q_lims = 1;
    [MVAbase_c, bus_c, gen_c, branch_c, success_c, et_c] = runpf(mpc, mpopt);
    
    %test if XXX changes from OPF solution
    % Does not change: Pg
    % Changes: Va, Vm
    % Violates: Qg
    diff = bus_c(:, VM) - mpc.bus(:, VM);
    if(any(abs(diff) > 0))
       fprintf('Changed values: %d out of %d\n', sum(abs(diff) > 0), size(diff,1)); 
    end
    
    %identify QG MAX and MIN violation
    idx = find(gen_c(:, QG) > gen_c(:, QMAX));
    n = length(idx);
    if(n > 0)
        QG_max = [QG_max; repmat(c, [n,1]) idx ((gen_c(idx, QG) - gen_c(idx, QMAX))./gen_c(idx, QMAX)) .* 100];%cont gen_idx percentage
    end
    
    idx = find(gen_c(:, QG) < gen_c(:, QMIN));
    n = length(idx);
    if(n > 0)
        QG_min = [QG_min; repmat(c, [n,1]) idx ((gen_c(idx, QMIN) - gen_c(idx, QG))./gen_c(idx, QMIN)) .* 100];%cont gen_idx percentage
    end
    
%     %identify Vm MAX and MIN violations
%     idx = find(bus_c(:, VM) > bus_c(:, VMAX));
%     n = length(idx);
%     if(n > 0)
%         Vm_max = [Vm_max; repmat(c, [n,1]) idx ((bus_c(idx, VM) - bus_c(idx, VMAX))./bus_c(idx, VMAX)) .* 100];%cont gen_idx percentage
%     end
%     
%     idx = find(bus_c(:, VM) < bus_c(:, VMIN));
%     n = length(idx);
%     if(n > 0)
%         Vm_min = [Vm_min; repmat(c, [n,1]) idx ((bus_c(idx, VMIN) - bus_c(idx, VM))./bus_c(idx, VMIN)) .* 100];%cont gen_idx percentage
%     end
end

sz = 20;

ax1 = subplot(2,1,1);
scatter(ax1, QG_max(:,1),QG_max(:,3),sz,'filled');
title('QG Upper bound violations');
xlabel('Contingency');
ylabel('Violation (% max)');

ax2 = subplot(2,1,2);
scatter(ax2, QG_min(:,1),abs(QG_min(:,3)),sz,'filled');
title('QG Lower bound violations');
xlabel('Contingency');
ylabel('Violation (% max)');


% ax3 = subplot(2,2,1);
% scatter(ax3, Vm_max(:,1),Vm_max(:,3),sz,'filled');
% title('Vm Upper bound violations');
% xlabel('Contingency');
% ylabel('Violation (% max)');
% 
% ax4 = subplot(2,2,2);
% scatter(ax4, Vm_min(:,1),abs(Vm_min(:,3)),sz,'filled');
% title('Vm Lower bound violations');
% xlabel('Contingency');
% ylabel('Violation (% max)');
%% Try to rerun OPF to see if the violations can be removed and the network is feasible 
% it might happen that the network will be feasible for individual contingencies,
% but not in the case they co-exist in the contingency set