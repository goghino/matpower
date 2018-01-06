clc
close all

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

warning('off','all');

define_constants;

%% initialization mode
% 0 - flat start
% 1 - local PF
% 2 - nominal OPF
% 3 - local OPF
mpopt0 = mpoption('opf.init_from_mpc', 2);
mpopt0 = mpoption(mpopt0, 'verbose', 0, 'out.all', 0);

cases = {  
    'case_ieee30',
    'case30',
    'case39',
    'case57',
    'case89pegase',
    'case118',
    'case300',
    'case1354pegase',
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

solvers = {
    %'MIPS (\)',
    %'MIPS (P)',
    %'MIPS-Sc (\)',
    %'MIPS-Sc (P)',
    %'fmincon',
    'Ipopt',
    %'Knitro'
};
mpopt = {
    %mpoption(mpopt0, 'opf.ac.solver', 'MIPS', 'mips.step_control', 0, 'mips.linsolver', ''),
    %mpoption(mpopt0, 'opf.ac.solver', 'MIPS', 'mips.step_control', 0, 'mips.linsolver', 'PARDISO'),
    %mpoption(mpopt0, 'opf.ac.solver', 'MIPS', 'mips.step_control', 1, 'mips.linsolver', ''),
    %mpoption(mpopt0, 'opf.ac.solver', 'MIPS', 'mips.step_control', 1, 'mips.linsolver', 'PARDISO'),
    %mpoption(mpopt0, 'opf.ac.solver', 'FMINCON'), 
    mpoption(mpopt0, 'opf.ac.solver', 'IPOPT'),
    %mpoption(mpopt0, 'opf.ac.solver', 'KNITRO')
};
na = length(solvers);
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

fprintf('\nResults%s', repmat(' ', 1, maxlen+4-length('Results')));
for a = 1:na
    spacers = repmat(' ', 1, round((20-length(solvers{a}))/2));
    fprintf('%-20s', sprintf('%s%s', spacers, solvers{a}) );
end
fprintf(' total secs\n');
fprintf(  '       %s',   repmat(' ', 1, maxlen+4-length('Results')));
for a = 1:na
    fprintf('   iter   *   secs  ');
end
fprintf('   elapsed\n');
fprintf(  '-------%s', repmat(' ', 1, maxlen+4-length('Results')));
for a = 1:na
    fprintf('  ------- - --------');
end
fprintf(  ' -----------\n');
sts = [sts zeros(nc, 4)];   %% nv, nnle, nnli, nlin
res.success = zeros(nc, na);
res.it = zeros(nc, na);
res.et = zeros(nc, na);
failures = zeros(na, 1);
total_et = zeros(na, 1);
total_it = zeros(na, 1);
t0 = tic;
% nc = 15;
s = warning('query', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');
for c = 1:nc
    spacers = repmat('.', 1, maxlen+3-length(sorted_cases{c}));
    fprintf('%s %s', sorted_cases{c}, spacers);
    mpc = loadcase(sorted_cases{c});
    if isfield(mpc, 'dcline')
        mpc = toggle_dcline(mpc, 'on');
    end
    if isfield(mpc, 'gencost')
        rpf = mpc;
        %insert missing branch limits
        no_limit = find(rpf.branch(:,RATE_A) < 1e-10);
        rpf.branch(no_limit,RATE_A) = 1e5;
        for a = 1:na
            %% -- run OPF
            %r = runopf(rpf, mpopt{a});
            
	    %% -- run SCOPF
	    cont = 5;
            r = runscopf(rpf, cont, mpopt{a}, 1e-4);

            %% -- process results
            sts(c, 13:16) = [r.om.var.N r.om.nle.N r.om.nli.N r.om.lin.N];
            res.success(c,a) = r.success;
            res.it(c,a) = r.raw.output.iterations;
            res.et(c,a) = r.et;
            fprintf('  %7d', res.it(c, a));
            if res.success(c,a)
                fprintf('  ');
            else
                fprintf(' x');
                failures(a) = failures(a) + 1;
            end
            total_it(a) = total_it(a) + res.it(c, a);
            total_et(a) = total_et(a) + res.et(c, a);
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
warning(s.state, 'MATLAB:nearlySingularMatrix');

fprintf('\n\nCase Statistics%s  nb    ng    nl   nlc  ndc  nvar   nnle   nnli   nlin   ni nib pq PLQ\n', repmat(' ', 1, maxlen-10));
fprintf(    '---------------%s----- ----- ----- ----- --- ------ ------ ------ ------ --- --- -- ---\n', repmat(' ', 1, maxlen-10));
for c = 1:nc
    spacers = repmat('.', 1, maxlen+3-length(sorted_cases{c}));
    fprintf('%s %s ', sorted_cases{c}, spacers);
    fprintf('%5d %5d %5d %5d %3d%7d%7d%7d%7d %3d %3d %1s%1s %1s%1s%1s\n', ...
        sts(c, [1:5 13:16 6:7]), sts(c, 8), sts(c, 9), sts(c, 10), sts(c, 11), sts(c, 12));
end
fprintf('\n');


fprintf('\nResults%s', repmat(' ', 1, maxlen+4-length('Results')));
for a = 1:na
    spacers = repmat(' ', 1, round((20-length(solvers{a}))/2));
    fprintf('%-20s', sprintf('%s%s', spacers, solvers{a}) );
end
fprintf('\n');
fprintf(  '       %s',   repmat(' ', 1, maxlen+4-length('Results')));
for a = 1:na
    fprintf('   iter   *   secs  ');
end
fprintf('\n');
fprintf(  '-------%s', repmat(' ', 1, maxlen+4-length('Results')));
for a = 1:na
    fprintf('  ------- - --------');
end
fprintf(  '\n');
for c = 1:nc
    spacers = repmat('.', 1, maxlen+3-length(sorted_cases{c}));
    fprintf('%s %s', sorted_cases{c}, spacers);
    for a = 1:na
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
    fprintf('\n');
end

fprintf('\nSummary\n');
fprintf('Solver           Failures  Total ET  Total It\n');
fprintf('---------------  --------  --------  --------\n');
for a = 1:na
    fprintf('%-15s%8d%12g%9d\n', solvers{a}, failures(a), total_et(a), total_it(a));
end

exit
