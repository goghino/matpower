clear; close all;
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

define_constants;
%% Select and configure the solver
SOLVER = 1;
if (SOLVER == 1)
    %for further options see ipopt.opt
    mpopt = mpoption('opf.ac.solver', 'IPOPT', 'verbose', 2);
elseif (SOLVER == 2)
    mpopt = mpoption('opf.ac.solver', 'MIPS', 'verbose', 3);
    mpopt.mips.max_it  = 100;
    mpopt.mips.max_stepsize = 1e12;
    mpopt.mips.feastol = 1e-4;
    mpopt.mips.gradtol = 1e-2;
    mpopt.mips.comptol = 1e-2;
    mpopt.mips.costtol = 1e-2;
    %mpopt.mips.linsolver='PARDISO';
else
    mpopt = mpoption('opf.ac.solver', 'OPTIZELLE');
end

mpopt = mpoption(mpopt, 'opf.init_from_mpc', 0);

%% load MATPOWER case struct, see help caseformat
%mpc = loadcase('case89pegase');
%mpc = loadcase('case118');
%mpc = loadcase('case89pegase');
mpc = loadcase('case1354pegase');
%mpc = loadcase('case9241pegase');
%mpc = case231swiss;

%insert missing branch limits
no_limit = find(mpc.branch(:,RATE_A) < 1e-10);
mpc.branch(no_limit,RATE_A) = 1e5;

%% Specify contingencies case 118
cont = -9999999;

%% Set required tolerance for satysfying PF equations (i.e. equality constraints |gx| < tol) and run SCOPF
tol = 1e-4;
runscopf(mpc, cont, mpopt, tol);

%%
% read-in IPOPT hessian and permute it
% name = ls('mat-ipopt_004-*');
% starts = strfind(name,'mat-ipopt'); %check if there is more mat-ipopt_004s
% if size(starts,2) > 1
%    name = name(1:starts(2)-2); %pick first file 004_
% else
%     name = name(1:end-1); % delete trailing ' '
% end
% 
% A = readcsr(name, 0, 1);
% 
% [P Pinv, npart] = KKTpermute(mpc, length(cont)+1);
% AP = A(P,P');
