%Script runs OPFs for all the specified contingencies and tries
%to identify which contingencies make OPF problem infeasible

define_constants;

mpc = case1354pegase;
load cont1354_noQGviol.mat;
cont = cont_filt;
ns = length(cont);

IPOPTresults = zeros(ns,1);
IPOPTiters = zeros(ns,1);

fprintf('Contingency Status Iterations\n');
for i = 1:ns
        %update mpc first by removing a line
        c = cont(i);
        mpc_tmp = mpc;
        if(c > 0)
            mpc_tmp.branch(c,BR_STATUS) = 0;
        end
    
        %run opf
        mpopt0 = mpoption('verbose', 0, 'out.all', 0);
        mpopt_tmp = mpoption(mpopt0, 'opf.ac.solver', 'IPOPT');
        [results, success] = runopf(mpc_tmp, mpopt_tmp);
        
        IPOPTresults(i) = success;
        IPOPTiters(i) = results.raw.output.iterations;

	fprintf('%d %d %d\n', c, success, results.raw.output.iterations);	
end

idx = find(IPOPTresults ~= 1);
cont(idx) = [];
dlmwrite('cont1354.txt',cont-1);
