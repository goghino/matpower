close all;

ORIGINAL = 1;

if ORIGINAL
    nb=118; ng=72; nl=186; ns=9; N=24;
    file = '~/Documents/Code/MatPowerMatrices/case118factor1N24/mat-ipopt_023-01.iajaa';
else
    nb=118; ng=72; nl=186; ns=9; N=31;
    file = '~/Documents/Optimization/matpower/SparseStructure7_USIsolver/mat-ipopt_023-01.iajaa';
end


pg = PowerGridHessian(file, nb, ng, nl, ns, N, 1);
pg.CheckStructure('nsAInDiagBlocks');