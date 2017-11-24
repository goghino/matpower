function lines = findDuplicateBranches(mpc)
%findDuplicateBranches  Identifies duplicate lines that connect the same two buses.
%   [Duplicates] = findDuplicateBranches(MPC)
%
%   Returns list of duplicates, where each row lists branches that
%   connect the same two buses. Zeros are used to align the resulting
%   matrix to be of size M x N, where M is size of the set of multi-line
%   connections and N is the maximum number of the lines connecting two
%   buses.
%
%   Examples:
%       Consider 4 bus network with following branches:
%        1: [1 2]
%        2: [1 2]
%        3: [2 3]
%        4: [3 4]
%        5: [3 4]
%        6: [3 4]
%        7: [4 1]
%       The result will be:
%        Duplicates = [1 2 0; 4 5 6]
%   
%   Output argument options:
%       duplicates = findDupliateBranches(mpc);
%
%   See also findIslandBranches.

    [F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
        TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
        ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

    lines = zeros(5,2);

    nbrch = size(mpc.branch,1);
    nbus = size(mpc.bus,1);
    i = 1;

    for b = 1:nbrch
        
        %skip branch if it is already identified
        if(ismember(b, lines))
            continue;
        end

        f = mpc.branch(b, F_BUS);
        t = mpc.branch(b, T_BUS);

        %find branches with the same t->f ends (or reversed f->t)
        idx = find((mpc.branch(:,F_BUS) == f) & (mpc.branch(:,T_BUS) == t));
        idx_rev = find((mpc.branch(:,F_BUS) == t) & (mpc.branch(:,T_BUS) == f));

        duplicates = [idx idx_rev];

        if (length(duplicates) > 1)
            lines(i,1:length(duplicates)) = duplicates;
            i = i + 1;
        end
    end
end

% for b = 1:nbus
%    idx1 = find(mpc.branch(:, F_BUS) == b);
%    idx2 = find(mpc.branch(:, T_BUS) == b);
%
%    brchs = [idx1; idx2];
%    buses = [mpc.branch(idx1, T_BUS); mpc.branch(idx2, F_BUS)];
%
%    uniBuses = unique(buses);
%
%    for i = uniBuses
%        if sum(buses == i) > 1
%
%        end
%    end
% end