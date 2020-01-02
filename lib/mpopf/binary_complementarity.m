function [f, df] = binary_complementarity(x)

x = x{1};

f = x.*(x-1);

if nargout > 1
    df=sparse(diag(2*x - 1));
end

end
