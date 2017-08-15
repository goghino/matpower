clear; close all;
define_constants;

setenv('OMP_NUM_THREADS', '1')
setenv('IPOPT_WRITE_MAT','1');
addpath('/Users/Juraj/Documents/Code/PowerGrid/matrices/'); %readcsr

%load MATPOWER case struct, see help caseformat
mpc = loadcase('case9');

nbus = size(mpc.bus, 1);
nbrch = size(mpc.branch, 1);
ngen = size(mpc.gen, 1);

npart = 2*nbus + 2*nbus + 2*nbrch + 2*nbrch;


%% This script tries to compare solution of IPOPT and MIPS solvers
% and check feasibility of the solution internaly in scopf.m

%list of case9 branches that do not leave isolated bus if cut
%contingencies = [2 3 5 6 8 9];

% specify list of branch contingencies, equals to OPF if empty
cont = [2]'; 
ns = size(cont,1) + 1; %number of scenarios (ncont + nominal)

%% run solvers
tol = 1e-2;
mpopt = mpoption('opf.ac.solver', 'MIPS', 'verbose', 2, ...
                 'mips.max_it', 150, ...
                 'mips.feastol', tol, ...
                 'mips.gradtol', tol, ...
                 'mips.comptol', tol, ...
                 'mips.costtol', tol);
%[results1, success1, info1] = runscopf(mpc, cont, mpopt);

%for solver options see ipopt.opt
mpopt = mpoption('opf.ac.solver', 'IPOPT', 'verbose', 2);
[results2, success2, info2] = runscopf(mpc, cont, mpopt);
%delete *.iajaa; %remove IPOPT matrices

% compare solutions x
[results1.x, results2.x]

results1.f - results2.f