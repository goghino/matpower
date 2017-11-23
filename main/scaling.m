clear; close all;
define_constants;

setenv('OMP_NUM_THREADS', '1')
setenv('IPOPT_WRITE_MAT','1');
addpath('/Users/Juraj/Documents/Code/PowerGrid/matrices/'); %readcsr

%load MATPOWER case struct, see help caseformat
mpc = loadcase('case9');

nbus = size(mpc.bus, 1);
nbrch = size(mpc.branch, 1);
ngen = size(mpc.gen, 1);
ns = 2;

%% permute system and form SC

npart = 2*nbus + 2*nbus + 2*nbrch + 2*nbrch;
nJ = 2*ngen;
N = ns*npart + nJ;
part1 = 1:npart;
part2 = npart + (1:npart);
partC = 2*npart + (1:nJ);

name = ls('mat-ipopt_009-*');
starts = strfind(name,'mat-ipopt'); %check if there is more mat-ipopt_004s
if size(starts,2) > 1
   name = name(1:starts(2)-2); %pick first file 004_
else
    name = name(1:end-1); % delete trailing ' '
end

%name = '/Users/Juraj/Desktop/KKT.csr'

A = readcsr(name, 0, 1); 
[P Pinv] = KKTpermute(nbus, nbrch, ngen, ns);
AP = A(P,P');

x = ones(N,1);
b = AP*x;

K1 = AP(part1, part1);
J1 = AP(ns*npart + (1:2*ngen), part1);

K2 = AP(part2, part2);
J2 = AP(ns*npart + (1:2*ngen), part2);

C = AP(partC, partC);

%%
% npart = 5;
% nJ = 2;
% N = 2*npart + nJ;
% part1 = 1:npart;
% part2 = npart + (1:npart);
% partC = 2*npart + (1:nJ);
% 
% K1 = magic(npart);
% J1 = rand(nJ, npart);
% 
% K2 = magic(5);
% J2 = rand(nJ, npart);
% 
% C = magic(nJ);
% 
% A = [K1       zeros(npart) J1';
%      zeros(npart) K2       J2';
%      J1           J2       C  ];
% x = ones(N,1);
% b = A*x;

%%
%% compute schur complement S Without scaling
% S = C - sum(Ji inv(Ki) JiT)
% Si = Ji inv(Ki) JiT

Si_noscale1 = J1 * (K1 \ J1');
Si_noscale2 = J2 * (K2 \ J2');
S_noscale = C - Si_noscale1 - Si_noscale2;

Si_rhs1 = J1 * (K1 \ b(part1));
Si_rhs2 = J2 * (K2 \ b(part2));
xc = S_noscale \ (b(partC) - Si_rhs1 - Si_rhs2);

x1 = K1 \ (b(part1) - J1'*xc);
x2 = K2 \ (b(part2) - J2'*xc);

x_noscale = [x1; x2; xc];
cond(K1)

%%

% With max row scaling
% maxK = max(abs(K1),[],2);
% idx = find(abs(maxK) < 1e-5);
% maxK(idx) = 1; %do not scale zeros
% M = diag(1./sqrt(maxK));

% With diagonal/maxrow scaling
%diagK =  abs(diag(K1));
diagK = max(abs(K1),[],2);
idx = find(diagK < 1e-5);
diagK(idx) = 1; %do not scale zeros
M1 = diag(1./sqrt(diagK));

%diagK = abs(diag(K2)); 
diagK = max(abs(K2),[],2);
idx = find(diagK < 1e-5);
diagK(idx) = 1; %do not scale zeros
M2 = diag(1./sqrt(diagK));

Si_scale1 =  J1 * M1 *((M1*K1*M1) \ (M1*J1'));
Si_scale2 =  J2 * M2 *((M2*K2*M2) \ (M2*J2'));
S_scale = C - Si_scale1 - Si_scale2;

Si_rhs1 = J1 * M1 * (M1*K1*M1 \ M1*b(part1));
Si_rhs2 = J2 * M2 * (M2*K2*M2 \ M2*b(part2));
xc = S_scale \ (b(partC) - Si_rhs1 - Si_rhs2);

x1 = M1 * (M1*K1*M1 \ M1*(b(part1) - J1'*xc));
x2 = M2 * (M2*K2*M2 \ M2*(b(part2) - J2'*xc));

x_scale = [x1; x2; xc];
cond(M1*K1*M1)

%% Direct solve

x_direct = AP\b;

%scale the whole matrix AP

%diagK =  abs(diag(K1));
diagK = max(abs(AP),[],2);
idx = find(diagK < 1e-5);
diagK(idx) = 1; %do not scale zeros
M = diag(1./sqrt(diagK));

x_direct_scaled = M * (M*AP*M \ M*b);


%% 
[norm(x-x_direct) norm(x-x_direct_scaled) norm(x-x_noscale) norm(x-x_scale)]
[norm(b-AP*x_direct)/norm(b) norm(b-AP*x_direct_scaled)/norm(b) norm(b-AP*x_noscale)/norm(b) norm(b-AP*x_scale)/norm(b)]

%% scaling is essentialy this diagonal scale
% A = magic(5);
% b = sum(A,2);
% A\b;
% 
% S = diag(1./sqrt(diag(A)));
% y = (S*A*S)\(S*b);
% x = S*y;