clc
close all

addpath( ...
    '/Users/Juraj/Documents/Optimization/matpower/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower/data', ...
    '/Users/Juraj/Documents/Optimization/matpower/mips/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower/mips/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower/most/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower/most/lib/t', ...
    '/Users/Juraj/Documents/Optimization/matpower/mptest/lib', ...
    '/Users/Juraj/Documents/Optimization/matpower/mptest/lib/t', ...
    '-end' );

addpath( ...
    '/Users/Juraj/Documents/Optimization/matpower/lib/mpopf', ...
    '/Users/Juraj/Documents/Optimization/matpower/lib/mpopf/ticinoData', ...
    '-end' );

setenv('OMP_NUM_THREADS', '1')

warning('off','all');

define_constants;
%% create case
mpc = case1354pegase;

%% initialization mode
%  1 = default starting point
%  2 = starting point taken directly from mpc
%  3 = AC power flow solution used as starting point
init_mode = 1;

mpopt = mpoption('verbose', 2, 'out.all', 1);
mpopt = mpoption(mpopt, 'opf.start', init_mode);
mpopt = mpoption(mpopt, 'opf.ac.solver', 'IPOPT');
                            
%%
Rfirst = 0.00; %% relative position of the storage placement when sorted by PD from largest first
Rcount = 0.02; %% fraction of the PD buses that will have storage
Emax = 2; %% max capacity of the storage relative to a PD at given bus
N = 1; %% Number of time periods

[r, SUCCESS] = runmpopf_dummy(mpc, mpopt, N, Emax, Rcount, Rfirst);