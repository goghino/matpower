function [lines, violations] = findQgCritical(mpc, branches, limit)
%findQgCritical  Identifies lines that result into critical Qg violation afer removal.
%   [Branches] = findQgCritical(MPC, BRANCHES, LIMIT)
%
%   Returns list of branches that create Qg violatios higher than limit after
%   their removal from the network. The Qg analysis is performed with
%   respect to the OPF solution.
%   
%   Examples:
%       Branches = findQgCritical(mpc, branches, limit);
%
%   See also findIslandBranches, findDuplicateBranches.

define_constants;

mpopt0 = mpoption('verbose', 0, 'out.all', 0);
if have_fcn('ipopt')
    mpopt = mpoption(mpopt0, 'opf.ac.solver', 'IPOPT');
else
    mpopt = mpoption(mpopt0, 'opf.ac.solver', 'MIPS');
end

lines = [];
violations = [];

%% run OPF and extract x*
[MVAbase, bus, gen, gencost, branch, f, success, et] = runopf(mpc, mpopt);
if(success ~= 1)
   error('OPF finished with error %d\n', success); 
end

xOPF = [bus(:,VA); bus(:,VM); gen(:,PG); gen(:,QG)];
mpcOPF = mpc;
mpcOPF.bus = bus;
mpcOPF.gen = gen;
mpcOPF.branch = branch;
mpcOPF.order.state = 'e';

%% identify QG violations
for ci = 1:length(branches)
    c = branches(ci);
    
    % copy the OPF mpc so that we can make changes to it
    mpc_test = mpcOPF;
    
    % remove branch i
    mpc_test.branch(c,BR_STATUS) = 0;
    
    % run PF with OPF solution and modified grid with a contingency c
    
    %mpopt.pf.enforce_q_lims = 1;
    warning('off', 'MATLAB:singularMatrix');
    warning('off', 'MATLAB:nearlySingularMatrix');
    [MVAbase_c, bus_c, gen_c, branch_c, success_c, et_c] = runpf(mpc_test, mpopt);
    
    %identify QG MAX and MIN violation
    idx = find(gen_c(:, QMAX) == 0);
    gen_c(idx, QMAX) = 1e-5;
    idx = find(gen_c(:, QG) > gen_c(:, QMAX));
    viol = ((gen_c(idx, QG) - gen_c(idx, QMAX))./gen_c(idx, QMAX)) .* 100;
    if (any(viol > limit))
        lines = [lines; c];
        violations = [violations; max(viol)];
    end
    
    idx = find(gen_c(:, QMIN) == 0);
    gen_c(idx, QMIN) = 1e-5;
    idx = find(gen_c(:, QG) < gen_c(:, QMIN));
    viol = ((gen_c(idx, QMIN) - gen_c(idx, QG))./gen_c(idx, QMIN)) .* 100;
    if (any(viol > limit))
        lines = [lines; c];
        violations = [violations; max(viol)];
    end
end

end