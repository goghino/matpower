clear; close all;

addpath( ...
    '/home/juraj/matpower/lib', ...
    '/home/juraj/matpower/lib/t', ...
    '/home/juraj/matpower/data', ...
    '/home/juraj/matpower/mips/lib', ...
    '/home/juraj/matpower/mips/lib/t', ...
    '/home/juraj/matpower/most/lib', ...
    '/home/juraj/matpower/most/lib/t', ...
    '/home/juraj/matpower/mptest/lib', ...
    '/home/juraj/matpower/mptest/lib/t', ...
    '/home/juraj/matpower-3.12.8', ...
    '-end' );

setenv('OMP_NUM_THREADS', '1')

define_constants;
%%
cases = {
%     'case9Q',
%     'case30pwl',
%     'case24_ieee_rts',
%     'case_RTS_GMLC',
%     'case30Q',
%     'case1888rte',
%     'case1951rte',
%     'case2736sp',
%     'case2848rte',
%     'case2737sop',
%     'case2746wop',
%     'case2746wp',
%     'case21k',
%     'case42k'
%     'case2868rte',
%     'case3012wp',
%     'case3120sp',
%     'case3375wp',
%     'case6468rte',
%     'case6470rte',
%     'case6495rte',
%     'case6515rte',
%     'case5',
%     'case6ww',
%     'case9',
%     'case14',
%     'case89pegase', %EXIT: Converged to a point of local infeasibility. Problem may be infeasible.
%     'case9241pegase',
%     'case13659pegase',
%     'case_ACTIVSg10k',
%-------------------
     'case14',
     'case30',
     'case39',
     'case57',
     'case89pegase',
     'case118',
     'case300',
     'case1354pegase',
     'case2383wp',
     'case2869pegase',
%     'case9241pegase',
%     'case13659pegase', %cannot find any feasible contingency
%     'case_ACTIVSg200', %did not converge
%     'case_ACTIVSg500', %did not converge
%     'case_ACTIVSg2000', %did not converge
%     'case1888rte', %no connected gen to ref bus
%     'case1951rte', %no connected gen to ref bus
%     'case2736sp', %multiple generators at ref bus
%     'case2737sop', %multiple generators at ref bus
%     'case2746wop',%multiple generators at ref bus
%     'case2746wp' %multiple generators at ref bus
};


mpopt0 = mpoption('verbose', 0, 'out.all', 0);
mpopt = {
      mpoption(mpopt0, 'opf.ac.solver', 'IPOPT', 'opf.init_from_mpc', 0), % flat start
      mpoption(mpopt0, 'opf.ac.solver', 'IPOPT', 'opf.init_from_mpc', 1), % local PF
      mpoption(mpopt0, 'opf.ac.solver', 'IPOPT', 'opf.init_from_mpc', 2), % nominal OPF
      mpoption(mpopt0, 'opf.ac.solver', 'IPOPT', 'opf.init_from_mpc', 3)  % local OPF
      mpoption(mpopt0, 'opf.ac.solver', 'IPOPT', 'opf.init_from_mpc', 4)  % nominal OPF + 1 cont
};
%     mpoption(mpopt0, 'opf.ac.solver', 'MIPS', 'mips.step_control', 0),
%     mpoption(mpopt0, 'opf.ac.solver', 'MIPS', 'mips.step_control', 1),
%     mpoption(mpopt0, 'opf.ac.solver', 'FMINCON'),
%%
nc = length(cases);
maxlen = 0;
for c = 1:nc
    if length(cases{c}) > maxlen
        maxlen = length(cases{c});
    end
end

%% run benchmarks
sts = zeros(nc, 16); % [nb, ng, nl, nlc, ndc, ni, nib, p,q, P,L,Q] and [X, G, H, A]

fprintf('\nResults%s       Ipopt-flat        Ipopt-pf-local     Ipopt-opf-nom     Ipopt-opf-local    total secs\n',  repmat(' ', 1, maxlen+5-length('Results')));
fprintf(  '       %s    iter   C   secs     iter   C   secs    iter   C   secs    iter   C   secs      elapsed \n',   repmat(' ', 1, maxlen+5-length('Results')));
fprintf(  '-------%s   ------- - --------  ------- - -------- ------- - -------- ------- - -------- -----------\n', repmat(' ', 1, maxlen+5-length('Results')));
na = length(mpopt);
res.success = zeros(nc, na);
res.it = zeros(nc, na);
res.et = zeros(nc, na);
res.nc = zeros(nc,na);

t0 = tic;
for c = 1:nc
    spacers = repmat(' ', 1, maxlen+3-length(cases{c}));
    fprintf('%s %s', cases{c}, spacers);
    mpc = loadcase(cases{c});
    if isfield(mpc, 'dcline')
        mpc = toggle_dcline(mpc, 'on');
    end
    %insert missing branch limits
    no_limit = find(mpc.branch(:,RATE_A) < 1e-10);
    mpc.branch(no_limit,RATE_A) = 1e5;
    if isfield(mpc, 'gencost')
        for a = 1:length(mpopt)
            %------------------------------------
            r = runscopf(mpc, -20, mpopt{a}, 1e-4);
            %------------------------------------
            sts(c, 13:16) = [r.raw.meta.lenX r.raw.meta.lenG r.raw.meta.lenH r.raw.meta.lenA];
            res.success(c,a) = r.success;
            res.it(c,a) = r.raw.output.iterations;
            res.et(c,a) = r.et;
	    res.nc(c,a) = length(r.raw.meta.cont) - 1;
            fprintf('  %7d', res.it(c, a));
            if res.success(c,a)
                fprintf('  ');
            else
                fprintf(' x');
            end
            if res.et(c, a) < 0.1
                fprintf(' %8.3f', res.et(c, a));
            elseif res.et(c, a) < 10
                fprintf(' %8.2f', res.et(c, a));
            elseif res.et(c, a) < 100
                fprintf(' %8.1f', res.et(c, a));
            else
                fprintf(' %8.0f', res.et(c, a));
            end
        end
    else
        fprintf('no gencost');
    end
    fprintf('%10.1f\n', toc(t0));
end

%% case statistics
% for c = 1:nc
%     mpc = loadcase(cases{c});
%     nb = size(mpc.bus, 1);
%     ng = size(mpc.gen, 1);
%     nl = size(mpc.branch, 1);
%     nlc = length(find(mpc.branch(:, RATE_A)));
%     ndc = 0;
%     if isfield(mpc, 'dcline')
%         ndc = size(mpc.dcline, 1);
%         mpc = toggle_dcline(mpc, 'on');
%     end
%     pqcost = 'p.';
%     if size(mpc.gencost, 1) == 2*ng
%         pqcost = 'pq';
%     end
%     costtype = '...';
%     if any(mpc.gencost(:, MODEL) == PW_LINEAR)
%         costtype(1) = 'P';
%     end
%     if any(mpc.gencost(:, MODEL) == POLYNOMIAL)
%         k = find(mpc.gencost(:, MODEL) == POLYNOMIAL);
%         if any(mpc.gencost(k, NCOST) == 2)
%             costtype(2) = 'L';
%         end
%         if any(mpc.gencost(k, NCOST) == 3 & mpc.gencost(k, COST) ~= 0)
%             costtype(3) = 'Q';
%         elseif any(mpc.gencost(k, NCOST) == 3 & mpc.gencost(k, COST) == 0)
%             costtype(2) = 'L';
%         end
%         if any(mpc.gencost(k, NCOST) > 3)
%             warning('cost order higher than quadratic');
%         end
%     end
%     [g, i] = find_islands(mpc);
%     sts(c, 1:7) = [nb ng nl nlc ndc length(g) length(i)];
%     sts(c, 8:9) = pqcost;
%     sts(c, 10:12) = costtype;
% end

% fprintf('\n\nCase Statistics%s  nb    ng    nl   nlc    X     G     H    ndc  ni  nib pq PLQ\n', repmat(' ', 1, maxlen-10));
% fprintf(    '---------------%s----- ----- ----- ----- ----- ----- ----- ----- --- --- -- ---\n', repmat(' ', 1, maxlen-10));
% for c = 1:nc
%     spacers = repmat('.', 1, maxlen+3-length(cases{c}));
%     fprintf('%s %s ', cases{c}, spacers);
%     fprintf('%5d %5d %5d %5d %5d %5d %5d %3d %3d %3d %1s%1s %1s%1s%1s\n', sts(c, 1:4), sts(c, 13:15), sts(c, 5:7), sts(c, 8), sts(c, 9), sts(c, 10), sts(c, 11), sts(c, 12));
% end


%% visualize convergence
fprintf('===================== Convergence results ========================\n');
fprintf('x0 \t sucess \t iter \t NC \t time\n');
for c = 1:nc
fprintf('%s\n', cases{c});
    for i = 1:na
        fprintf('%d \t %d \t\t %d \t %d \t %.2f\n', mpopt{i}.opf.init_from_mpc, res.success(c, i), res.it(c, i), res.nc(c,i), res.et(c, i));
    end
end

fprintf('===================== Convergence results ========================\n');
fprintf('case flat localPF nominalOPF localOPF nominalOPF1cont\n');
for c = 1:nc
fprintf('%s ', cases{c});
    for i = 1:na
        fprintf('%d ', res.it(c, i));
    end
    fprintf('\n');
end
