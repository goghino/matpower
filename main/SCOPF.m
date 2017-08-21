clear; close all;
define_constants;

setenv('OMP_NUM_THREADS', '1')
setenv('IPOPT_WRITE_MAT','1');
addpath('/Users/Juraj/Documents/Code/PowerGrid/matrices/'); %readcsr

%mpopt = mpoption('opf.ac.solver', 'MIPS', 'verbose', 2, 'mips.max_it', 100);
mpopt = mpoption('opf.ac.solver', 'IPOPT', 'verbose', 2);

%load MATPOWER case struct, see help caseformat
mpc = loadcase('case118');

%%

cont = [1:10];

runscopf(mpc, cont, mpopt);