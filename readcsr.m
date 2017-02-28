function [A, prows, pcols, pdata] = readcsr(filename, has_cols, is_sym)

  fin = fopen(filename, 'r');

  n        = fscanf(fin, '%d', 1);
  m = 0;
  if (has_cols)
    m        = fscanf(fin, '%d', 1);
  else
    m = n;
  end
  nonzeros = fscanf(fin, '%d', 1);

  fprintf('number of rows: %d\n', n); 
  fprintf('number of cols: %d\n', m); 
  fprintf('number of nnzr: %d\n', nonzeros); 

  %iA       = zeros(n+1);
  %jA       = zeros(nonzeros);
  %data     = zeros(nonzeros);

  iA       = fscanf(fin, '%d', n+1);
  jA       = fscanf(fin, '%d', nonzeros);
  data     = fscanf(fin, '%f', nonzeros);
  I        = zeros(nonzeros,1);

  fclose(fin); % finished reading the file

  fprintf('file read. Constructing sparse matrix\n');

  first    = 1-iA(1);
  iA       = iA + first;
  jA       = jA + first;
data = data; % + 1.0;

for i=1:n
    for index=iA(i):iA(i+1) - 1
        I(index,1) = i;
    end
end

A = sparse(I, jA, data);
        
if (is_sym)
  A = A+A.'-spdiags(diag(A,0), 0, size(A,1), size(A,2)); 
end
