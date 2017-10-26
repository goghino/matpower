clear; close all;
define_constants;

setenv('OMP_NUM_THREADS', '1')
%setenv('IPOPT_WRITE_MAT','1');
addpath('/Users/Juraj/Documents/Code/PowerGrid/matrices/'); %readcsr

%% Select and configure the solver
mpopt = mpoption('opf.ac.solver', 'MIPS', 'verbose', 3);
mpopt.mips.max_it=50;
mpopt.mips.feastol = 1e-4;
mpopt.mips.gradtol = 1e-2;
mpopt.mips.comptol = 1e-2;
mpopt.mips.costtol = 1e-2;
%mpopt.mips.linsolver='PARDISO'; %TODOOoooooooooooooooooooooooooooooooooooooooooooo

%mpopt = mpoption('opf.ac.solver', 'IPOPT', 'verbose', 2);
%%for further options see ipopt.opt

%% load MATPOWER case struct, see help caseformat
%mpc = loadcase('case9');
%mpc = loadcase('case118');
%mpc = loadcase('case9241pegase');
mpc = loadcase('case1354pegase');
%mpc = case231swiss;


%insert missing branch limits
no_limit = find(mpc.branch(:,RATE_A) < 1e-10);
mpc.branch(no_limit,RATE_A) = 1e5;

%% Specify contingencies

%cont = [1:6 8 10:112 114:132 135:175 178:182 185:186]; %case118
%cont = [ 4 5 14 15 18 19 20 22 26 32 33 36 37  40 41 42 43 44 46 47 48 49 50 51 54 55 64 72 74 75 79 80 87 99 100]; %pegase 1354
%cont = [1:34 37:92 94:121 124:173 176:203 206:219 222:225 227:320 323:366 368:433 436:16000]; %case pegase 9k
cont = [ 46 47 ]; %error in case1k
%cont = [8];

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