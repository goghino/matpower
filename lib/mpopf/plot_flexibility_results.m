function plot_flexibility_results(mpc_solution)

    define_constants
    
    if(~isfield(mpc_solution, 'storageFlexibility'))
       return; 
    end
    if(mpc_solution.storageFlexibility ~= 1)
        return;
    end

    nstorage = mpc_solution.nstorage;
    N = mpc_solution.horizon;
    
    %% get storage profile
    %skip generators in the mpc.gen structure and extract storage charge/discharge
    Pgen_discharge = reshape(mpc_solution.gen(mpc_solution.id_gen_storages_discharge,PG),[nstorage, N ]);
    Pgen_discharge = [Pgen_discharge Pgen_discharge(:,end)];
    
    Pgen_charge    = reshape(mpc_solution.gen(mpc_solution.id_gen_storages_charge   ,PG),[nstorage, N ]);
    Pgen_charge    = [Pgen_charge Pgen_charge(:,end)];
    
    max_product_charge_discharge_MWMW = max(max(abs(Pgen_discharge.*Pgen_charge)))
    
    Pgen_storage_max = repmat(mpc_solution.P_storage_max_MW,[1,N+1]);
    Pgen_storage_min = repmat(mpc_solution.P_storage_min_MW,[1,N+1]);


    %% get flexibility profile
    flex_locations = mpc_solution.id_gen_storages_discharge + 2*nstorage*N; %skip generators and charge and discharge in mpc.gen
    Pgen_d_up    = reshape(mpc_solution.gen(flex_locations,PG),[nstorage, N ]); %discharge up flexibility
    Pgen_d_up    = [Pgen_d_up Pgen_d_up(:,end)];
    
    flex_locations = flex_locations + nstorage*N; %additionaly skip also up-dis flexibility
    Pgen_d_down  = reshape(mpc_solution.gen(flex_locations,PG),[nstorage, N ]); %discharge down flexibility
    Pgen_d_down  = [Pgen_d_down Pgen_d_down(:,end)];
    
    flex_locations = flex_locations + nstorage*N; %additionaly skip also down-dis flexibility
    Pgen_c_up    = reshape(mpc_solution.gen(flex_locations   ,PG),[nstorage, N ]); %charge up flexibility
    Pgen_c_up    = [Pgen_c_up Pgen_c_up(:,end)];
    
    flex_locations = flex_locations + nstorage*N; %additionaly skip also up-charg flexibility
    Pgen_c_down  = reshape(mpc_solution.gen(flex_locations   ,PG),[nstorage, N ]); %charge down flexibility
    Pgen_c_down    = [Pgen_c_down Pgen_c_down(:,end)];
    
    max_product_dischargeUP_dischargeDOWN = max(max(abs(Pgen_d_up.*Pgen_d_down)))
    max_product_chargeUP_chargeDOWN = max(max(abs(Pgen_c_up.*Pgen_c_down)))
    max_product_dischargeUP_chargeDOWN = max(max(abs(Pgen_d_up.*Pgen_c_down)))
    max_product_chargeUP_dischargeDOWN = max(max(abs(Pgen_c_up.*Pgen_d_down)))
    
    max_product_discharge_chargeDOWN = max(max(abs(Pgen_discharge.*Pgen_c_down)))
    max_product_discharge_chargeUP = max(max(abs(Pgen_discharge.*Pgen_c_up)))
    max_product_charge_dischargeDOWN = max(max(abs(Pgen_charge.*Pgen_d_down)))
    max_product_charge_dischargeUP = max(max(abs(Pgen_charge.*Pgen_d_up)))


    %% Plot output of each storage individually with flexibility and the limits
    figure;
    for s = 1:nstorage
        subplot(nstorage,1,s)
        stairs(Pgen_discharge(s,:)'+Pgen_d_up(s,:)'+Pgen_d_down(s,:)', 'b', 'LineWidth',3); hold on;
        stairs(Pgen_discharge(s,:)', 'b--', 'LineWidth',2); hold on;
        stairs(Pgen_d_up(s,:)','g:', 'LineWidth',2); hold on;
        stairs(Pgen_d_down(s,:)','g--', 'LineWidth',2); hold on;
        stairs([Pgen_storage_max(s,:)]', 'r:', 'LineWidth',2);
        legend('P_{ed}+u_{ed}+d_{ed}', 'P_{ed}', 'u_{ed}', 'd_{ed}', 'P_{max}');
        grid on
        title("Storage network injection - discharging flexibility. Device " + num2str(s))
        xlabel('t [hours]')
        ylabel('Storage power [MW]')
        %set(gca,'FontSize',20);
    end
    
    figure;
    for s = 1:nstorage
        subplot(nstorage,1,s)
        stairs(Pgen_charge(s,:)'+Pgen_c_up(s,:)'+Pgen_c_down(s,:)', 'b', 'LineWidth',3); hold on;
        stairs(Pgen_charge(s,:)', 'b--', 'LineWidth',2); hold on;
        stairs(Pgen_c_up(s,:)','g:', 'LineWidth',2); hold on;
        stairs(Pgen_c_down(s,:)','g--', 'LineWidth',2); hold on;
        stairs([Pgen_storage_min(s,:)]', 'r', 'LineWidth',2);
        legend('P_{ec}+u_{ec}+d_{ec}', 'P_c', 'u_{ec}', 'd_{ec}', 'P_{max}');
        grid on
        title("Storage network injection - charging flexibility. Device " + num2str(s))
        xlabel('t [hours]')
        ylabel('Storage power [MW]')
        %set(gca,'FontSize',20);
    end

%% Plot the flexibility
    figure;
    requirement = [mpc_solution.FlexibilityReq.up' mpc_solution.FlexibilityReq.up(end)];
    subplot(2,1,1)
    stairs(sum(Pgen_d_up,1),'b--', 'LineWidth',2); hold on;
    stairs(sum(Pgen_c_up,1),'b:', 'LineWidth',2); hold on;
    stairs(requirement, 'r:', 'LineWidth',2);
    legend('u_{ed}', 'u_{ec}', 'u_{min}');
    grid on
    title('Cumulative Storage up flexibility (all storage devices)')
    xlabel('t [hours]')
    ylabel('Storage power [MW]')
    %set(gca,'FontSize',20);
    
    requirement = [mpc_solution.FlexibilityReq.down' mpc_solution.FlexibilityReq.down(end)];
    subplot(2,1,2)
    stairs(sum(Pgen_d_down,1),'b--', 'LineWidth',2); hold on;
    stairs(sum(Pgen_c_down,1),'b:', 'LineWidth',2); hold on;
    stairs(requirement, 'r:', 'LineWidth',2);
    legend('d_{ed}', 'd_{ec}', 'd_{max}');
    grid on
    title('Cumulative Storage down flexibility (all storage devices)')
    xlabel('t [hours]')
    ylabel('Storage power [MW]')
    %set(gca,'FontSize',20);

%%  Plot overall storage profile
    %set matching colors for discharge and the corresponding maximum
    map = get(gca,'ColorOrder');
    ns = size(Pgen_discharge, 1); %number of storage devices
    map_n = size(map,1); %size of the color pallete
    map = repmat(map, [ceil(ns/map_n),1]); %make color pallete match the ns
    map(ns+1:end, :) = [];
    set(gca, 'ColorOrder',map, 'NextPlot','ReplaceChildren')
    
    Pgen_storage = Pgen_discharge+Pgen_charge;
    
    fig=figure;
    set(fig,'defaultAxesColorOrder',map);
    subplot(2,1,1)
    stairs(Pgen_storage' , 'LineWidth',3);
    hold on;
    stairs(Pgen_storage_max', '--', 'LineWidth',3);
    stairs(Pgen_storage_min', ':', 'LineWidth',3);
    grid on
    title('Individual storage network injection power : trajectory (solid) and maximum(dashed)/minimum(dotted)')
    xlabel('t [hours]')
    ylabel('Storage power [MW]')
    %set(gca,'FontSize',20);

    %remove the extra added period at the end (added for the plotting purposes)
    Pgen_discharge = Pgen_discharge(:,1:N);
    Pgen_charge = Pgen_charge(:,1:N);
    Pgen_d_up = Pgen_d_up(:,1:N);
    Pgen_d_down = Pgen_d_down(:,1:N);
    Pgen_c_up = Pgen_c_up(:,1:N);
    Pgen_c_down = Pgen_c_down(:,1:N);
    
    E_storage = [mpc_solution.E_storage_init_MWh, repmat(mpc_solution.E_storage_init_MWh,[1,N]) - cumsum((Pgen_discharge)./repmat(mpc_solution.c_discharge,[1,N]) + (Pgen_charge).*repmat(mpc_solution.c_charge,[1,N]),2)];
    E_storage_max = repmat(mpc_solution.E_storage_max_MWh,[1,N+1]);


    subplot(2,1,2)
    plot(E_storage', 'LineWidth',3);
    hold on;
    plot(E_storage_max','--', 'LineWidth',3);
    grid on
    title('Individual storage energy level: trajectory (solid) and maximum (dashed)')
    xlabel('t [hours]')
    ylabel('Storage energy [MWh]')
    %set(gca,'FontSize',20);

end
