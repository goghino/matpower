function plot_storage_results(mpc_solution)

    define_constants

    nstorage = mpc_solution.nstorage;
    N = mpc_solution.horizon;
    Pgen_discharge = reshape(mpc_solution.gen(mpc_solution.id_gen_storages_discharge,PG),[nstorage, N ]);
    Pgen_charge    = reshape(mpc_solution.gen(mpc_solution.id_gen_storages_charge   ,PG),[nstorage, N ]);
    Pgen_storage = Pgen_discharge+Pgen_charge;
    Pgen_storage_max = repmat(mpc_solution.P_storage_max_MW,[1,N]);
    Pgen_storage_min = repmat(mpc_solution.P_storage_min_MW,[1,N]);

    E_storage = [mpc_solution.E_storage_init_MWh, repmat(mpc_solution.E_storage_init_MWh,[1,N]) - cumsum(Pgen_discharge./repmat(mpc_solution.c_discharge,[1,N]) + Pgen_charge.*repmat(mpc_solution.c_charge,[1,N]),2)];
    E_storage_max = repmat(mpc_solution.E_storage_max_MWh,[1,N+1]);


    max_product_charge_discharge_MWMW = max(max(abs(Pgen_discharge.*Pgen_charge)))


    figure;
    subplot(2,1,1)
    plot(Pgen_discharge');
    grid on
    title('storage network injection - discharging')
    xlabel('t [hours]')
    ylabel('storage power [MW]')
    subplot(2,1,2)
    plot(Pgen_charge');
    grid on
    title('storage network injection - charging')
    xlabel('t [hours]')
    ylabel('storage power [MW]')


    figure;
    subplot(2,1,1)
    stairs(0:N, [Pgen_storage Pgen_storage(:,end)  ]', 'LineWidth',2);
    hold on;
    stairs(0:N, [Pgen_storage_max Pgen_storage_max(:,end)  ]', '--');
    stairs(0:N, [Pgen_storage_min Pgen_storage_min(:,end)  ]', '--');
    grid on
    title('storage network injection power : trajectory (solid) and maximum (dashed)')
    xlabel('t [hours]')
    ylabel('storage power [MW]')
    subplot(2,1,2)
    plot(0:N, E_storage', 'LineWidth',2);
    hold on;
    plot(0:N, E_storage_max','--');
    grid on
    title('storage energy level: trajectory (solid) and maximum (dashed)')
    xlabel('t [hours]')
    ylabel('storage energy [MWh]')

    figure;
    subplot(2,1,1)
    stairs(0:N, sum([Pgen_storage Pgen_storage(:,end)])', 'b-', 'LineWidth',2);
    hold on;
    stairs(0:N, sum([Pgen_storage_max Pgen_storage_max(:,end)  ])', 'b--');
    stairs(0:N, sum([Pgen_storage_min Pgen_storage_min(:,end)  ])', 'b--');
    grid on
    title('total storage network injection power : trajectory (solid) and maximum (dashed)')
    xlabel('t [hours]')
    ylabel('storage power [MW]')
    subplot(2,1,2)
    plot(0:N, sum(E_storage)','b-', 'LineWidth',2);
    hold on;
    plot(0:N, sum(E_storage_max)','b--');
    grid on
    title('total storage energy level: trajectory (solid) and maximum (dashed)')
    xlabel('t [hours]')
    ylabel('storage energy [MWh]')

end
