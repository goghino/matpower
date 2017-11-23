function mpcN_opf_storage = create_storage_case_file(mpc,load_scaling_profile, p_storage)

    define_constants
    mpcN_opf = create_multi_time_step_case_file4(mpc,load_scaling_profile);


    nnodes =size(mpc.bus,1);

    id_storage_location = p_storage.id_storage_location;
    E_storage_max_MWh   = p_storage.E_storage_max_MWh;
    E_storage_init_MWh  = p_storage.E_storage_init_MWh ;
    nstorage = length(id_storage_location);

    ngen0 = size(mpc.gen,1);
    ngen_total = ngen0+nstorage*2;

    P_storage_max_MW   =  p_storage.rPmaxEmax_MW_per_MWh*E_storage_max_MWh;
    P_storage_min_MW   =  p_storage.rPminEmax_MW_per_MWh*E_storage_max_MWh;
    if length(p_storage.c_discharge)> 1
        c_discharge        = p_storage.c_discharge;
    else
        c_discharge        = p_storage.c_discharge*ones(nstorage,1);
    end
    if length(p_storage.c_charge)> 1
        c_charge        = p_storage.c_charge;
    else
        c_charge        = p_storage.c_charge*ones(nstorage,1);
    end


    mpcN_opf_storage = add_storage2bNoPeriodic(mpcN_opf,nnodes,id_storage_location,P_storage_max_MW,P_storage_min_MW,E_storage_max_MWh,E_storage_init_MWh,c_discharge,c_charge, 1)
    % max(max(abs(mpcN_opf_storage.bus - mpcN_opf_storage2.bus)))
    % max(max(abs(mpcN_opf_storage.branch - mpcN_opf_storage2.branch)))
    % max(max(abs(mpcN_opf_storage.gen - mpcN_opf_storage2.gen)))
    % max(max(abs(mpcN_opf_storage.gencost - mpcN_opf_storage2.gencost)))
    % max(max(abs(mpcN_opf_storage.A - mpcN_opf_storage2.A)))
    % max(max(abs(mpcN_opf_storage.l - mpcN_opf_storage2.l)))
    % max(max(abs(mpcN_opf_storage.u - mpcN_opf_storage2.u)))


    clear mpcN_opf,nnodes,id_storage_location,P_storage_max_MW,P_storage_min_MW,E_storage_max_MWh,E_storage_init_MWh,c_discharge,c_charge;
    nb = size(mpcN_opf_storage.bus,1);
    ng = size(mpcN_opf_storage.gen,1);
    nl = size(mpcN_opf_storage.branch,1);
    nA = size(mpcN_opf_storage.A,1);
    N = length(load_scaling_profile);

    %[idT1,idT2,T1,T2,idHnz] = compute_sparsity_structure(nb,ng,nl,nA,N);

    % J_ordered = T1*Js*T2;
    % H_ordered = T2'*Hs*T2;


    fprintf('nb, %5d %5d\n', nb, nb/N);
    fprintf('ng, %5d %5d\n', ng, ng/N);
    fprintf('nl, %5d %5d\n', nl, nl/N);
    fprintf('nA, %5d %5d\n', nA, nA/N);
    fprintf('N , %5d\n', N);
    fprintf('nstorage, %5d \n', nstorage);

    adjust_ipopt_opt(N,nb/N,ng/N,nl/N,nstorage)

end
