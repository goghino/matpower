function h = binary_complementarity_hess(x, lam)

x = x{1};

h = sum(lam) * sparse(diag(2*ones(size(x,1),1)));

end