function plot_flexibility_results(mpc_solution)

    define_constants
    
    if(~isfield(mpc_solution, 'enableDemandShift'))
       return; 
    end
    if(mpc_solution.enableDemandShift ~= 1)
        return;
    end

    nstorage = mpc_solution.nstorage;
    N = mpc_solution.horizon;
    
    nstorageN = nstorage * N;    
    ngen_storageN =  nstorageN*2;
    ngen_storageFlexibilityN = 0;
    if mpc_solution.storageFlexibility
        ngen_storageFlexibilityN = nstorageN*4; %up/down charg, up/down discharg.
    end
    
    DS_bus_idx = mpc_solution.demandShift.busesID;
    nbus_demandShift = length(mpc_solution.demandShift.busesID); %buses which provide DS
    nbus_demandShiftN = nbus_demandShift*N; %buses which provide DS, overall count
    ngen_demandShiftN = nbus_demandShift*4*N; %response+upFlex and rebound+downFlex
    
    ngen_nostorageN = size(mpc_solution.gen,1) - (ngen_storageN + ngen_storageFlexibilityN + ngen_demandShiftN);
    
    offset = ngen_nostorageN + ngen_storageN + ngen_storageFlexibilityN;

    %% get DS variables Pdu, ud, Pdd, dd
    id_gen_Pdu = offset +                       (1:nbus_demandShiftN);
    id_gen_ud =  offset + nbus_demandShiftN   + (1:nbus_demandShiftN);
    id_gen_Pdd = offset + 2*nbus_demandShiftN + (1:nbus_demandShiftN);
    id_gen_dd =  offset + 3*nbus_demandShiftN + (1:nbus_demandShiftN);

    %skip generators in the mpc.gen structure and extract DS
    Pdu = reshape(mpc_solution.gen(id_gen_Pdu,PG),[nbus_demandShift, N ]);
    Pdu = [Pdu Pdu(:,end)];
    
    ud = reshape(mpc_solution.gen(id_gen_ud,PG),[nbus_demandShift, N ]);
    ud = [ud ud(:,end)];
    
    Pdd = reshape(mpc_solution.gen(id_gen_Pdd,PG),[nbus_demandShift, N ]);
    Pdd = [Pdd Pdd(:,end)];
    
    dd = reshape(mpc_solution.gen(id_gen_dd,PG),[nbus_demandShift, N ]);
    dd = [dd dd(:,end)];
    
    b = reshape(mpc_solution.x(end-nbus_demandShiftN+1:end),[nbus_demandShift, N ]);
    b = [b, b(:,end)];
    
    max_product_Pdu_Pdd_MWMW = max(max(abs(Pdu.*Pdd)))
    max_product_ud_dd_MWMW = max(max(abs(ud.*dd)))
    max_product_Pdu_dd_MWMW = max(max(abs(Pdu.*dd)))
    max_product_ud_Pdd_MWMW = max(max(abs(ud.*Pdd)))
    
    Prsp = repmat(mpc_solution.demandShift.responsePowerMW,[1,N+1]);
    Prb = repmat(mpc_solution.demandShift.reboundPowerMW,[1,N+1]);


    %% Plot output of each DS bus individually with flexibility and the limits
    figure;
    for s = 1:nbus_demandShift
        subplot(nbus_demandShift,1,s)
        stairs(-Pdu(s,:)'+ud(s,:)', 'b', 'LineWidth',3); hold on;
        stairs(-Pdu(s,:)', 'b--', 'LineWidth',2); hold on;
        stairs(ud(s,:)','g:', 'LineWidth',2); hold on;
        stairs([Prsp(s,:)]', 'r:', 'LineWidth',2);
        legend('-P_{du}+u_{d}', '-P_{du}', 'u_{d}', 'P_{rsp}');
        grid on; xticks(1:N);
        title("Response power. Bus " + num2str(DS_bus_idx(s)))
        xlabel('t [hours]')
        ylabel('Power [MW]')
        %set(gca,'FontSize',20);
    end
    
    figure;
    for s = 1:nbus_demandShift
        subplot(nbus_demandShift,1,s)
        stairs(Pdd(s,:)'-dd(s,:)', 'b', 'LineWidth',3); hold on;
        stairs(Pdd(s,:)', 'b--', 'LineWidth',2); hold on;
        stairs(-dd(s,:)','g:', 'LineWidth',2); hold on;
        stairs([Prb(s,:)]', 'r:', 'LineWidth',2);
        legend('P_{dd}-d_{d}', 'P_{dd}', '-d_{d}', 'P_{rb}');
        grid on; xticks(1:N);
        title("Rebound power. Bus " + num2str(DS_bus_idx(s)))
        xlabel('t [hours]')
        ylabel('Power [MW]')
        %set(gca,'FontSize',20);
    end
    
    max_product_binary = max(max(abs(b.*(b-1))))
    
    figure;
    for s = 1:nbus_demandShift
        subplot(nbus_demandShift,1,s)
        stairs(b(s,:)','b', 'LineWidth',2); hold on;
        legend('b_{d}');
        grid on; xticks(1:N);
        title("Binary variable corresponding to the activation of DS at bus" + num2str(DS_bus_idx(s)))
        xlabel('t [hours]')
        %set(gca,'FontSize',20);
    end

    %% Plot the flexibility
    figure;
    requirement = [mpc_solution.FlexibilityReq.up' mpc_solution.FlexibilityReq.up(end)];
    subplot(2,1,1)
    stairs(sum(ud,1),'b--', 'LineWidth',2); hold on;
    stairs(requirement, 'r:', 'LineWidth',2);
    legend('u_{d}', 'u_{min}');
    grid on; xticks(1:N);
    title('Cumulative DS up flexibility (all demands)')
    xlabel('t [hours]')
    ylabel('Power [MW]')
    %set(gca,'FontSize',20);
    
    requirement = [mpc_solution.FlexibilityReq.down' mpc_solution.FlexibilityReq.down(end)];
    subplot(2,1,2)
    stairs(sum(dd,1),'b--', 'LineWidth',2); hold on;
    stairs(requirement, 'r:', 'LineWidth',2);
    legend('d_{d}', 'd_{max}');
    grid on; xticks(1:N);
    title('Cumulative DS down flexibility (all demands)')
    xlabel('t [hours]')
    ylabel('Power [MW]')
    %set(gca,'FontSize',20);
end
