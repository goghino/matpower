clear; close all;
define_constants;

setenv('OMP_NUM_THREADS', '1')
setenv('IPOPT_WRITE_MAT','1');
addpath('/Users/Juraj/Documents/Code/PowerGrid/matrices/'); %readcsr

%mpopt = mpoption('opf.ac.solver', 'MIPS', 'verbose', 2, 'mips.max_it', 100);
mpopt = mpoption('opf.ac.solver', 'IPOPT', 'verbose', 2);

%load MATPOWER case struct, see help caseformat
mpc = loadcase('case118');

nbus = size(mpc.bus, 1);
nbrch = size(mpc.branch, 1);
ngen = size(mpc.gen, 1);

delete *.iajaa;


%%

%list of case9 branches that do not leave isolated bus if cut
%contingencies = [2 3 5 6 8 9];
contingencies = [1:6 8 10:50];

%condition number for the full IPOPT system
cond_KKT = zeros(length(contingencies),1);
%condition number of the diagonal blocks after reordering
cond_block = zeros(length(contingencies),2);
cond_pblock = zeros(length(contingencies),2);

for i = 1:20
    
    % specify list of branch contingencies, equals to OPF if empty
    cont = contingencies(1:2*i);
    ns = length(cont) + 1; %number of scenarios (ncont + nominal)
 
    %runopf(mpc, mpopt);
    runscopf(mpc, cont, mpopt);

    % read-in IPOPT hessian and permute it
    name = ls('mat-ipopt_004-*');
    starts = strfind(name,'mat-ipopt'); %check if there is more mat-ipopt_004s
    if size(starts,2) > 1
       name = name(1:starts(2)-2); %pick first file 004_
    else
        name = name(1:end-1); % delete trailing ' '
    end
    
    A = readcsr(name, 0, 1);
    
    %condition num of full IPOPT KKT system
    cond_KKT(i,1) = cond((A));
    
    [P Pinv, npart] = KKTpermute(mpc, ns);
    AP = A(P,P'); 
    if (A - AP(Pinv, Pinv') ~= sparse(size(A,1),size(A,2)))
       error('Inverse permutation does not result in original matrix.') 
    end
    
    %diagonal blocks after permutation
    for j = 1:ns
        block = AP((j-1)*npart+(1:npart),(j-1)*npart+(1:npart));
        cond_block(i,j) = cond(full(block));
        
        %preprocess block
%         d = abs(max(block,[],2));
%         idx = find(d < 1e-12);
%         d(idx) = 1;
%         S = diag(1./sqrt(d));
%         cond_pblock(i,j) = cond(full(S*block*S));
     
    end
    
    %remove IPOPT matrices
    delete *.iajaa;
end

cond_KKT
cond_block
%cond_pblock

%isolating bus by contingency creates singular Cb, zero on diagonal
% -> singular local hessian and fucked up condition numbers!!!
% -> verify SC for reasonable matrices without isolated bus

%% fix diagonal (for case 9, cont on branch 1)
% idx = find(abs(A) > 1e10);
% fix = [floor(idx / size(A,1))+1 mod(idx, size(A,1))]
% for i = 1:size(fix,1)
% A(fix(i,1), fix(i,2)) = -A(fix(i,1), fix(i,2)+3); %A[1,1] = -A[1,4] => branch1  (1->4)
% end
% cond(full(A))


%% check inertia (n,m,0)
% E = eig(A)
% size(find(E == 0))
% size(find(E > 0))
% size(find(E < 0))

%desired inertia is obtained if dw is sufficiently large
%and the constraint Jacobian dc(xk)T has full rank

%In practice, however, the iteration matrix can become
%so ill-conditioned that the factorization cannot be performed
%successfully, even with very large values of dw and some dc > 0.
%In this case, we give up on the current step computation and switch
%directly to the feasibility restoration phase, hoping that the matrix
%has better properties close to feasible points.

%if dw exceeds dw_max than abort directly to feasibility restoration phase,
%i.e. it was not possible to make the matrix inertia correct!!!

%% TODO:
%1. verify identity with opf and SCOPF with no contingencies
%2. create 1 contingency and compare structure of jac/hes of nominal case
%   vs. contingency, investigate how sparse structure behaves, it it stays
%   constant of if there are missing entries corresponding to a contingency
%3. generate tests for matpower_cpp, keep in mind that the entries are kept
%   there but they have zero value.