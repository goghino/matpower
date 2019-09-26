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

    %set matching colors for discharge and the corresponding maximum
    map = get(gca,'ColorOrder');
    ns = size(Pgen_discharge, 1); %number of storage devices
    map_n = size(map,1); %size of the color pallete
    map = repmat(map, [ceil(ns/map_n),1]); %make color pallete match the ns
    map(ns+1:end, :) = [];
    set(gca, 'ColorOrder',map, 'NextPlot','ReplaceChildren')

    fig=figure;
    set(fig,'defaultAxesColorOrder',map);
    subplot(2,1,1)
    stairs(Pgen_discharge', 'LineWidth',2);
    hold on;
    stairs(0:N, [Pgen_storage_max Pgen_storage_max(:,end)  ]', '--');
    grid on
    title('Storage network injection - discharging')
    xlabel('t [hours]')
    ylabel('Storage power [MW]')
    set(gca,'FontSize',20);
    subplot(2,1,2)
    stairs(Pgen_charge', 'LineWidth',2);
    hold on;
    stairs(0:N, [Pgen_storage_min Pgen_storage_min(:,end)  ]', '--');
    grid on
    title('Storage network injection - charging')
    xlabel('t [hours]')
    ylabel('Storage power [MW]')
    set(gca,'FontSize',20);


    fig=figure;
    set(fig,'defaultAxesColorOrder',map);
    subplot(2,1,1)
    stairs( [Pgen_storage ]', 'LineWidth',2);
    hold on;
    stairs(0:N, [Pgen_storage_max Pgen_storage_max(:,end)  ]', '--');
    stairs(0:N, [Pgen_storage_min Pgen_storage_min(:,end)  ]', '--');
    grid on
    title('Storage network injection power : trajectory (solid) and maximum (dashed)')
    xlabel('t [hours]')
    ylabel('Storage power [MW]')
    set(gca,'FontSize',20);
    subplot(2,1,2)
    plot(0:N, E_storage', 'LineWidth',2);
    hold on;
    plot(0:N, E_storage_max','--');
    grid on
    title('Storage energy level: trajectory (solid) and maximum (dashed)')
    xlabel('t [hours]')
    ylabel('Storage energy [MWh]')
    set(gca,'FontSize',20);

    fig=figure;
    set(fig,'defaultAxesColorOrder',map);
    subplot(2,1,1)
    stairs(sum([Pgen_storage])', 'b-', 'LineWidth',2);
    hold on;
    stairs(0:N, sum([Pgen_storage_max Pgen_storage_max(:,end)  ])', 'b--');
    stairs(0:N, sum([Pgen_storage_min Pgen_storage_min(:,end)  ])', 'b--');
    grid on
    title('Total storage network injection power : trajectory (solid) and maximum (dashed)')
    xlabel('t [hours]')
    ylabel('Storage power [MW]')
    set(gca,'FontSize',20);
    subplot(2,1,2)
    plot(0:N, sum(E_storage)','b-', 'LineWidth',2);
    hold on;
    plot(0:N, sum(E_storage_max)','b--');
    grid on
    title('Total storage energy level: trajectory (solid) and maximum (dashed)')
    xlabel('t [hours]')
    ylabel('Storage energy [MWh]')
    set(gca,'FontSize',20);

end
