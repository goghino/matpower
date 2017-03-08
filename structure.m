close all;

ORIGINAL = 0;
RAMPS = 1;

if ORIGINAL
    nb=118; ng=72; nl=186; ns=9; N=24;
    file = '~/Documents/Code/MatPowerMatrices/case118factor1N24/mat-ipopt_023-01.iajaa';
else
    %mpc = case118
    %-------------
    %number_of_timesteps    24
    %number_of_buses       118
    %number_of_generators   54
    %number_of_lines       186
    %number_of_storages      3
    nb=118; ng=60; nl=186; ns=3; N=24;
    if RAMPS
        file = '~/Documents/Optimization/matpower/SparseStructure7_USIsolver/ramp.iajaa';
    else
        file = '~/Documents/Optimization/matpower/SparseStructure7_USIsolver/noramp.iajaa';
    end
end


pg = PowerGridHessian(file, nb, ng, nl, ns, N, 1);
pg.CheckStructure('nsAInDiagBlocks', RAMPS);