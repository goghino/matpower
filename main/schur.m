% set up sizes of the grid and the system
mpc = case9;

nb = 9;
nb_PV = 2;
ng = 3;
nl = 9;
ns = 3;

nH = (nb - 1) + (nb-nb_PV) + ng + 1;
nG = 2*nb;
nK = 4*nl;

npart = nH + nG + nK;
N2 = nb_PV + (ng - 1);

%% set up blocks and build a system
H = triu(magic(npart)); H = H + H' - diag(diag(H));
s = sum(H,2); J_T = repmat(s, [1,N2]); J = J_T';

KKT_ = kron(eye(ns), H);
jt = kron([1:ns]', J_T);
KKT = sparse([KKT_ jt; jt' eye(N2)]);

% use inverse permutation to create original KKT
[P Pinv nPart] = KKTpermute(mpc, 3);
KKT_orig = KKT(Pinv,Pinv);
writecsr('/Users/Juraj/Desktop/KKT.csr', triu(KKT_orig), 0);

%reference solution
%x_ref = ones(length(P),1);
x_ref = [1:length(P)]';

% create RHS
%rhs_orig = KKT_orig * x_ref;
rhs_orig = KKT_orig * x_ref;
rhs = rhs_orig(P);

%% 

% Schur complement S = Sc - sum(JH\J_T)
S = eye(N2);
for i = 1:ns
   S = S - (i*J)*(H\(i.*J_T)); 
end

% rhs for Schur system rhs = bs - sum(JH\bi)
rhs_S = rhs(ns*npart+1:end);
for i = 1:ns
   range = ((i-1)*npart+1):i*npart;
   rhs_S = rhs_S - (i*J)*(H\rhs(range,1)); 
end

xs = S \ rhs_S;

% rhs for partitions rhs_i = bi - J_T xs
rhs_i = zeros(npart,ns);
for i = 1:ns
   range = ((i-1)*npart+1):i*npart;
   rhs_i(:,i) = rhs(range,1) - (i.*J_T) * xs;
end

% H * xi = bi - J_T xs
xi = zeros(npart,ns);
for i = 1:ns
   xi(:,i) = H \ rhs_i(:,i);
end

x = [ xi(:); xs];
