close all;
clear all;

%spys at the beginning differ in nnz and the very first also in structure block_11 (why)?
files = {...
    '../mat-pardiso_001-phase-22.mtx', ...
    '../mat-pardiso_003-phase-22.mtx', ...
    '../mat-pardiso_005-phase-22.mtx', ...
    '../mat-pardiso_009-phase-22.mtx', ...
    '../mat-pardiso_013-phase-22.mtx'};
%%
for f = files
    A = load(char(f), '-ascii');

    i = A(:,1);
    j = A(:,2);
    a = A(:,3);
    S = sparse(i,j,a);
    figure;
    spy(S);
    set(gca,'fontsize',12)
end