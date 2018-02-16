load('ipopt_benchmark_new.mat'); 

N = length(bench_params);

mu = zeros(N,1);
sigma = zeros(N,1);
failed = zeros(N,1);

for i = 1:N
    p = bench_params{i};
    mu(i) = mean(p.iterations);
    sigma(i) = std(p.iterations); %normalized by length(iterations)-1
    failed(i) = sum(p.iterations==250);
end

[sorted, idx] = sort(mu,'ascend');
fprintf('Avg_it\tStddev\tfailed\n');
fprintf('%.2f\t%.2f\t%d\n',[sorted, sigma(idx), failed(idx)].');

%% print the best set of params
bench_params{idx(1)}
bench_params{idx(2)}
bench_params{idx(3)}
