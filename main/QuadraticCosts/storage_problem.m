clc
close all

addpath( ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/data', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/mips/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/mips/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/most/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/most/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/mptest/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower-origin/mptest/lib/t', ...
    '-end' );

setenv('OMP_NUM_THREADS', '1')

define_constants;

%% create base case file
mpc = case118;
%mpc = case30; % works with top 2%, not with static = 10
%mpc = case89pegase;
%mpc = case300;

%number of time periods
N = 1;

%% run the benchmark
%Emax = 0.01; %max capacity of the storage relative to PD at given bus
%Rcount = 0.02; % fraction of PD buses that will have storage
Rfirst = 0.0; % relative position of the storage placement when sorted by PD from largest first

results = struct('N', [], 'PFviol', [], 'BranchViol', [], 'EmaxViol', [], 'EminViol', [] );

Emax_list = [100 10 1 0.1 0.01];
Emax_list = [0];
Rcount_list = [0.2 0.1 0.05 0.02 0.01];
Rcount_list = [0.04];

for Emax = Emax_list
    for Rcount = Rcount_list
        %% execute MPOPF
        opt_solution = MPOPF(mpc, N, Emax, Rcount, Rfirst);
        
        %% verify constraints
        % 0 < Einit + T*P1              < Emax 
        % 0 < Einit + T*P1 + T*P2       < Emax 
        % 0 < Einit + T*P1 + T*P2 + ... < Emax 

        nstorage = opt_solution.nstorage;
        N = opt_solution.horizon;
        Pgen_discharge = reshape(opt_solution.gen(opt_solution.id_gen_storages_discharge,PG),[nstorage, N ]);
        Pgen_charge    = reshape(opt_solution.gen(opt_solution.id_gen_storages_charge   ,PG),[nstorage, N ]);
        Pgen_storage = Pgen_discharge+Pgen_charge;
        Pgen_storage_max = repmat(opt_solution.P_storage_max_MW,[1,N]);
        Pgen_storage_min = repmat(opt_solution.P_storage_min_MW,[1,N]);

        %Einit - (Ec + Ed)_1 - (Ec + Ed)_2 - ... - (Ec + Ed)_N
        % where Ec_i = Pci * eps_c and Pc < 0
        %       Ed_i = Pdi / eps_d and Pd > 0
        E_storage = [opt_solution.E_storage_init_MWh,...
                     repmat(opt_solution.E_storage_init_MWh,[1,N]) - ...
                     cumsum(Pgen_discharge./repmat(opt_solution.c_discharge,[1,N]) +...
                            Pgen_charge.*repmat(opt_solution.c_charge,[1,N]),2)];
        E_storage_max = repmat(opt_solution.E_storage_max_MWh,[1,N+1]);

        E = E_storage;
        E_u = opt_solution.E_storage_max_MWh;
        E_l = 0*opt_solution.E_storage_max_MWh;

        %PF constraints and branch PF
        [h, g] = opf_consfcn(opt_solution.x, opt_solution.om);
        gn = 2*size(mpc.bus,1);
        hn = 2*size(mpc.branch,1);

        PFviol = 0;
        BranchViol = 0;
        EmaxViol = 0;
        EminViol = 0;

        for i = 1:N
            E_i = E((1:nstorage), i);
            idx = find(E_i > E_u );
            nviol = length(idx);
            viol = max(E_i(idx) - E_u(idx));
            if (nviol > 0)
               fprintf('[N=%2d] Vioated MAX energy limit for %d storages with error %e\n', i, nviol, viol); 
               EmaxViol = max(EmaxViol, viol);
            end

            idx = find(E_i < E_l );
            nviol = length(idx);
            viol = max(E_l(idx) - E_i(idx));
            if (nviol > 0)
               fprintf('[N=%2d] Vioated MIN energy limit for %d storages with error %e\n', i, nviol, viol); 
               EminViol = max(EminViol, viol);
            end

            g_local = g((1:gn) + (i-1)*gn);
            idx = find(g_local > 1e-5 );
            nviol = length(idx);
            viol = max(g_local(idx));
            if (nviol > 0)
               fprintf('[N=%2d] Vioated %d PF equations with error %e\n', i, nviol, viol); 
               PFviol = max(PFviol, viol);
            end

            h_local = h((1:hn) + (i-1)*hn);
            idx = find(h_local > 0 );
            nviol = length(idx);
            viol = max(h_local(idx));
            if (nviol > 0)
               fprintf('[N=%2d] Vioated %d branch flow limits with error %e\n', i, nviol, viol); 
               BranchViol = max(BranchViol, viol);
            end
        end

        %% save results
        results.N = [results.N opt_solution.raw.output.iterations];
        results.PFviol = [results.PFviol PFviol];
        results.BranchViol = [results.BranchViol BranchViol];
        results.EmaxViol = [results.EmaxViol EmaxViol];
        results.EminViol = [results.EminViol EminViol];
    end
end


%% plot the results
figure;
nEmax = length(Emax_list);
nRcount = length(Rcount_list);
Legend=cell(nEmax,1);
for i = 1:nEmax
    semilogy(1:nRcount, results.PFviol((1:nRcount) + (i-1)*nRcount));
    xlabel('Number of storages [% of load buses]');
    ylabel('PF violation');
    xticks(1:nRcount);
    xticklabels(Rcount_list*100);
    Legend{i} = strcat('E [MWh] = ', num2str(Emax_list(i)), '*PD');
    hold on;
end
title('PF violations')
legend(Legend);


figure;
nEmax = length(Emax_list);
nRcount = length(Rcount_list);
Legend=cell(nEmax,1);
for i = 1:nEmax
    semilogy(1:nRcount, results.EmaxViol((1:nRcount) + (i-1)*nRcount));
    xlabel('Number of storages [% of load buses]');
    ylabel('Emax violation');
    xticks(1:nRcount);
    xticklabels(Rcount_list*100);
    Legend{i} = strcat('E [MWh] = ', num2str(Emax_list(i)), '*PD');
    hold on;
end
title('Emax Violations')
legend(Legend);

%%
%figure; plot(E(1:nstorage_applied:end)); title('Energy level of 1st storage over time'); xlabel('T [h]'); ylabel('E [MWh]')
%figure; plot(E(2:nstorage_applied:end)); title('Energy level of 2st storage over time'); xlabel('T [h]'); ylabel('E [MWh]')

%% verify constraints graphically (Alex's code)

figure();
stairs(0:N, [Pgen_storage Pgen_storage(:,end)  ]', 'LineWidth',2);
hold on;
stairs(0:N, [Pgen_storage_max Pgen_storage_max(:,end)  ]', '--');
stairs(0:N, [Pgen_storage_min Pgen_storage_min(:,end)  ]', '--');
grid on
title('storage network injection power : trajectory (solid) and maximum (dashed)')
xlabel('t [hours]')
ylabel('storage power [MW]')

figure();
plot(0:N, E_storage', 'LineWidth',2);
hold on;
plot(0:N, E_storage_max','--');
grid on
title('storage energy level: trajectory (solid) and maximum (dashed)')
xlabel('t [hours]')
ylabel('storage energy [MWh]')


%% parse x vector
% nb = size(mpc.bus,1);
% ng = size(mpc.gen,1);
% 
% offset = nb*N + nb*N; %skip Va Vm
% 
% Pg_discharge = opt_solution.x(offset + ng*N + (1:nstorage_applied:nstorage_applied*N)) * opt_solution.baseMVA;
% offset = offset + nstorage_applied*N;
% Pg_charge = opt_solution.x(offset + (1:nstorage_applied:nstorage_applied*N)) * opt_solution.baseMVA;
% offset = offset + nstorage_applied*N + ng*N;
% Qg_discharge = opt_solution.x(offset + (1:nstorage_applied*N)) * opt_solution.baseMVA;
% offset = offset + nstorage_applied*N;
% Qg_charge = opt_solution.x(offset + (1:nstorage_applied*N)) * opt_solution.baseMVA;
% 
% figure; plot(Pg_discharge(1:nstorage_applied:end)); title('Pg discharge over time'); xlabel('T [h]'); ylabel('P [MW]')
% figure; plot(Pg_charge(1:nstorage_applied:end));title('Pg charge over time'); xlabel('T [h]'); ylabel('P [MW]')
