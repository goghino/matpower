function [lines, violations] = findQgCritical(mpc, branches, limit)
%findQgCritical  Identifies lines that result into critical Qg violation afer removal.
%   [Branches] = findQgCritical(MPC, BRANCHES, LIMIT)
%
%   Returns list of branches that create Qg violatios higher than limit after
%   their removal from the network. The Qg analysis is performed with
%   respect to the OPF solution.
%
%   Limit represents relative violation to QMAX or QMIN.
%   violation = (QG - QMAX) / QMAX >= limit
%   violation = (QMIN - QG) / QMIN >= limit
%   If the Qmax or Qmin is 0, the violation is considered
%   for situations if abs(QG - Qmin) > 1e-2
%   
%   Examples:
%       Branches = findQgCritical(mpc, branches, limit);
%
%   See also findIslandBranches, findDuplicateBranches.

define_constants;

mpopt0 = mpoption('verbose', 0, 'out.all', 0, 'opf.start', 3);
if have_fcn('ipopt')
    mpopt = mpoption(mpopt0, 'opf.ac.solver', 'IPOPT');
else
    mpopt = mpoption(mpopt0, 'opf.ac.solver', 'MIPS');
end

lines = [];
violations = [];

%% run OPF and extract x*
[mpcOPF, success] = runopf(mpc, mpopt);
if(success ~= 1)
   error('OPF finished with error %d in findQgCritical\n', success); 
end

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
    idxMAX = find(gen_c(:, QG) > gen_c(:, QMAX) & gen_c(:, QMAX) ~= 0);
    violMAX = (gen_c(idxMAX, QG) - gen_c(idxMAX, QMAX)) ./ gen_c(idxMAX, QMAX);
    idx = find(gen_c(:, QMAX) == 0); %handles situation of QMAX = 0 (prevents division by zero)
    if (any(abs(gen_c(idx, QMAX) - gen_c(idx, QG)) > 1e-2))
        violMAX = [violMAX; limit];
    end
    
    idxMIN = find(gen_c(:, QG) < gen_c(:, QMIN) & gen_c(:, QMIN) ~= 0);
    violMIN = (gen_c(idxMIN, QMIN) - gen_c(idxMIN, QG)) ./ gen_c(idxMIN, QMIN);
    idx = find(gen_c(:, QMIN) == 0); %handles situation of QMIN = 0 (prevents division by zero)
    if (any(abs(gen_c(idx, QMIN) - gen_c(idx, QG)) > 1e-2))
        violMIN = [violMIN; limit];
    end
    
    if (any(violMAX >= limit))
        lines = [lines; c];
        violations = [violations; max(violMAX)];
    elseif (any(violMIN >= limit))
        lines = [lines; c];
        violations = [violations; max(violMIN)];
    end
end

end