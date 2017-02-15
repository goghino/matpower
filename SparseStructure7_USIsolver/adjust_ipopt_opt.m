function adjust_ipopt_opt(N,nbus,ngen,nbranch,nstorage)
    return;

    A = regexp( fileread('ipopt.opt'), '\n', 'split');
    A{9} = ['number_of_timesteps    '        num2str(N)];
    A{10} = ['number_of_buses   '        num2str(nbus)];
    A{11} = ['number_of_generators    '        num2str(ngen)];
    A{12} = ['number_of_lines    '        num2str(nbranch)];
    A{13} = ['number_of_storages    '        num2str(nstorage)];
%     A{14} = ['nlp_scaling_method            "none"   '];
%     A{14} = ['line_search_method            "penalty"   '];

    a = A(1:end-1);

    fid = fopen('ipopt.opt', 'w');
    fprintf(fid, '%s\n', a{:});
    fclose(fid);
end
