function spyKKT(name, ns, mpc)

    %name = 'mat-ipopt_000-02.iajaa';
    %ns = 9;
    %mpc = case118;

    addpath('/Users/Juraj/Documents/Code/PowerGrid/matrices/'); %readcsr

    A = readcsr(name, 0, 1);

    [P Pinv, npart] = KKTpermute(mpc, ns);
    AP = A(P,P'); 
    if (A - AP(Pinv, Pinv') ~= sparse(size(A,1),size(A,2)))
       error('Inverse permutation does not result in original matrix.') 
    end
    
    spy(A);
    figure; spy(AP);

end