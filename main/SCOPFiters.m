clear; close all;

addpath( ...
    '/Users/Juraj/Documents/Optimization/matpower/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower/data', ...
    '/Users/Juraj/Documents/Optimization/matpower/mips/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower/mips/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower/most/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower/most/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower/mptest/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower/mptest/lib/t', ...
    '-end' );

setenv('OMP_NUM_THREADS', '1')
%setenv('IPOPT_WRITE_MAT','1');
addpath('/Users/Juraj/Documents/Code/PowerGrid/matrices/'); %readcsr

define_constants;
%%
cases = {
%     'case_ieee30',
%     'case_RTS_GMLC',
%     'case5',
%     'case6ww',
    'case9',
%     'case9Q',
%     'case14',
%     'case24_ieee_rts',
%     'case30',
%     'case30pwl',
%     'case30Q',
%     'case39',
%     'case57',
     'case89pegase',
    'case118',
%    'case300',
%     'case1354pegase',
%     'case1888rte',
%     'case1951rte',
%     'case2383wp',
%     'case2736sp',
%     'case2737sop',
%     'case2746wop',
%     'case2746wp',
%     'case2848rte',
%     'case2868rte',
%     'case2869pegase',
%     'case3012wp',
%     'case3120sp',
%     'case3375wp',
%     'case6468rte',
%     'case6470rte',
%     'case6495rte',
%     'case6515rte',
%     'case_ACTIVSg10k',
%     'case_ACTIVSg200',
%     'case_ACTIVSg500',
%     'case_ACTIVSg2000',
%     'case9241pegase',
%     'case13659pegase',
%     'case21k',
%     'case42k'
};

mpopt0 = mpoption('verbose', 0, 'out.all', 0);
mpopt = {
%     mpoption(mpopt0, 'opf.ac.solver', 'MIPS', 'mips.step_control', 0),
%     mpoption(mpopt0, 'opf.ac.solver', 'MIPS', 'mips.step_control', 1),
%     mpoption(mpopt0, 'opf.ac.solver', 'FMINCON'),
      mpoption(mpopt0, 'opf.ac.solver', 'IPOPT', 'opf.init_from_mpc', 0),
      mpoption(mpopt0, 'opf.ac.solver', 'IPOPT', 'opf.init_from_mpc', 1)
};

%%
nc = length(cases);
maxlen = 0;
for c = 1:nc
    if length(cases{c}) > maxlen
        maxlen = length(cases{c});
    end
end

%% case statistics
sts = zeros(nc, 12);   %% nb, ng, nl, nlc, ndc, ni, nib, p,q, P,L,Q
for c = 1:nc
    fprintf('.');
    mpc = loadcase(cases{c});
    nb = size(mpc.bus, 1);
    ng = size(mpc.gen, 1);
    nl = size(mpc.branch, 1);
    nlc = length(find(mpc.branch(:, RATE_A)));
    ndc = 0;
    if isfield(mpc, 'dcline')
        ndc = size(mpc.dcline, 1);
        mpc = toggle_dcline(mpc, 'on');
    end
    pqcost = 'p.';
    if size(mpc.gencost, 1) == 2*ng
        pqcost = 'pq';
    end
    costtype = '...';
    if any(mpc.gencost(:, MODEL) == PW_LINEAR)
        costtype(1) = 'P';
    end
    if any(mpc.gencost(:, MODEL) == POLYNOMIAL)
        k = find(mpc.gencost(:, MODEL) == POLYNOMIAL);
        if any(mpc.gencost(k, NCOST) == 2)
            costtype(2) = 'L';
        end
        if any(mpc.gencost(k, NCOST) == 3 & mpc.gencost(k, COST) ~= 0)
            costtype(3) = 'Q';
        elseif any(mpc.gencost(k, NCOST) == 3 & mpc.gencost(k, COST) == 0)
            costtype(2) = 'L';
        end
        if any(mpc.gencost(k, NCOST) > 3)
            warning('cost order higher than quadratic');
        end
    end
    [g, i] = find_islands(mpc);
    sts(c, 1:7) = [nb ng nl nlc ndc length(g) length(i)];
    sts(c, 8:9) = pqcost;
    sts(c, 10:12) = costtype;
end

%% sort cases
[junk, k] = sort(sts(:,1) + sts(:,2) + sts(:,4) + sts(:,5));    %% nb+ng+nlc+ndc
sorted_cases = cases(k);
sts = sts(k,:);

fprintf('\n\nCase Statistics%s  nb    ng    nl   nlc  ndc  ni nib pq PLQ\n', repmat(' ', 1, maxlen-10));
fprintf(    '---------------%s----- ----- ----- ----- --- --- --- -- ---\n', repmat(' ', 1, maxlen-10));
for c = 1:nc
    spacers = repmat('.', 1, maxlen+3-length(sorted_cases{c}));
    fprintf('%s %s ', sorted_cases{c}, spacers);
    fprintf('%5d %5d %5d %5d %3d %3d %3d %1s%1s %1s%1s%1s\n', sts(c, 1:7), sts(c, 8), sts(c, 9), sts(c, 10), sts(c, 11), sts(c, 12));
end

%%
fprintf('\nResults%s       Ipopt-flat         Ipopt-opf       total secs\n',  repmat(' ', 1, maxlen+5-length('Results')));
fprintf(  '       %s    iter   C   secs     iter   C   secs     elapsed\n',   repmat(' ', 1, maxlen+5-length('Results')));
fprintf(  '-------%s   ------- - --------  ------- - -------- -----------\n', repmat(' ', 1, maxlen+5-length('Results')));
sts = [sts zeros(nc, 4)];   %% nv, nnle, nnli, nlin
na = length(mpopt);
res.success = zeros(nc, na);
res.it = zeros(nc, na);
res.et = zeros(nc, na);

t0 = tic;
for c = 1:nc
    spacers = repmat(' ', 1, maxlen+3-length(sorted_cases{c}));
    fprintf('%s %s', sorted_cases{c}, spacers);
    mpc = loadcase(sorted_cases{c});
    if isfield(mpc, 'dcline')
        mpc = toggle_dcline(mpc, 'on');
    end
    if isfield(mpc, 'gencost')
        for a = 1:length(mpopt)
            %------------------------------------
            r = runscopf(mpc, [], mpopt{a}, 1e-4);
            %------------------------------------
            sts(c, 13:16) = [r.raw.meta.lenX r.raw.meta.lenG r.raw.meta.lenH r.raw.meta.lenA];
            res.success(c,a) = r.success;
            res.it(c,a) = r.raw.output.iterations;
            res.et(c,a) = r.et;
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