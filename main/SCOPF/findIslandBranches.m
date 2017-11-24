function lines = findIslandBranches(mpc)
%findIslandBranches  Identifies lines that create islands or isolated buses.
%   [Branches] = findIslandBranches(MPC)
%
%   Returns list of branches that would create islands or isolated
%   buses when removed from the network.
%   
%   Output argument options:
%       Branches = findIslandBranches(mpc);
%
%   See also findDupliateBranches.

    nbranch = size(mpc.branch,1);
    lines = [];

    for l = 1:nbranch
        
        mpc_test = mpc;
        mpc_test.branch(l,:) = [];
        
        %%-- remove branches causing islands or isolated buses
        [islands, isolated] = find_islands(mpc_test);
        
        if (~isempty(isolated) || size(islands,2) > 1)
            lines = [lines; l];
        end
    end 
end