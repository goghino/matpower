function [feasible] = findOPFfeasible(mpc, branches)
%findOPFfeasible  Identifies lines that result into infeasible OPF problem afer removal.
%   [Lines] = findOPFfeasible(MPC, BRANCHES)
%
%   
%   Examples:
%       Branches = findOPFfeasible(mpc, branches);
%
%   See also findIslandBranches, findDuplicateBranches.

%Script runs OPFs for all the specified contingencies and tries
%to identify which contingencies make OPF problem infeasible

define_constants;

printInfo = 0;

cont = branches;
ns = length(cont);

IPOPTresults = zeros(ns,1);
IPOPTiters = zeros(ns,1);

if(printInfo == 1)
    fprintf('Contingency Status Iterations\n');
end
for i = 1:ns
        %update mpc first by removing a line
        c = cont(i);
        mpc_tmp = mpc;
        if(c > 0)
            mpc_tmp.branch(c,BR_STATUS) = 0;
        end
    
        %% initialization mode
        %%  1 = default starting point
        %%  2 = starting point taken directly from mpc
        %%  3 = AC power flow solution used as starting point
        mpopt0 = mpoption('verbose', 0, 'out.all', 0, 'opf.start', 3);
        mpopt_tmp = mpoption(mpopt0, 'opf.ac.solver', 'IPOPT');
        
        %run opf
        [results, success] = runopf(mpc_tmp, mpopt_tmp);
        
        IPOPTresults(i) = success;
        IPOPTiters(i) = results.raw.output.iterations;

        if(printInfo == 1)
        fprintf('%d %d %d\n', c, success, results.raw.output.iterations);	
        end
end

idx = find(IPOPTresults == 1);
feasible = cont(idx);

end
