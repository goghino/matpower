clear; close all;
define_constants;

setenv('OMP_NUM_THREADS', '1')
setenv('IPOPT_WRITE_MAT','1');
addpath('/Users/Juraj/Documents/Code/PowerGrid/matrices/'); %readcsr

%load MATPOWER case struct, see help caseformat
mpopt = mpoption('opf.ac.solver', 'IPOPT', 'verbose', 2);

%load MATPOWER case struct, see help caseformat
mpc = loadcase('case9');

nbus = size(mpc.bus, 1);
nbrch = size(mpc.branch, 1);
ngen = size(mpc.gen, 1);


%%

%list of case9 branches that do not leave isolated bus if cut
contingencies = [2 3 5 6 8 9];

%condition number for the full IPOPT system
cond_KKT = zeros(nbrch,1);
%condition number of the small diagonal blocks before reordering
cond_hess = zeros(nbrch,2);
%condition number of the big diagonal blocks after reordering
cond_block = zeros(nbrch,2);
cond_pblock = zeros(nbrch,2);
%condition of part of the diagonal block from cond_block
cond_bblock1 = zeros(nbrch,2); %SC w.r.t smaller block
cond_bblock2 = zeros(nbrch,2); %SC w.r.t bigger block
cond_bblock3 = zeros(nbrch,2); %SC w.r.t bigger block

for i = contingencies
    
    % specify list of branch contingencies, equals to OPF if empty
    cont = [i]'; 
    ns = size(cont,1) + 1; %number of scenarios (ncont + nominal)

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
    [P Pinv] = KKTpermute(nbus, nbrch, ngen, ns);
    AP = A(P,P'); 
    
%     figure; spy(AP)
%     figure; spy(A)

    %condition num of full IPOPT KKT system
    cond_KKT(i,1) = cond(full(A));

    %local hessians wrt primal variables for nominal and contingency
    nH = 2*nbus;
    for j = 1:ns
        cond_hess(i,j) = cond(full(A((j-1)*nH+(1:nH),(j-1)*nH+(1:nH))));
    end
    
    %diagonal blocks after permutation
    npart = 2*nbus + 2*nbus + 2*nbrch + 2*nbrch;
    nppart = 2*nbus + 2*nbus + 2*nbrch;
    for j = 1:ns
        block = AP((j-1)*npart+(1:npart),(j-1)*npart+(1:npart));
        cond_block(i,j) = cond(full(block));
        
        %preprocess block
        d = abs(diag(block));
        idx = find(d < 1e-12);
        d(idx) = 1;
        S = diag(1./sqrt(d));
        cond_pblock(i,j) = cond(full(S*block*S));
        
        %break down the block into more blocks
        block11 = block(1:nppart,1:nppart);
        block22 = block(nppart+(1:2*nbrch),nppart+(1:2*nbrch));
        B = block(nppart+(1:2*nbrch),1:nppart);
        BT = block(1:nppart,nppart+(1:2*nbrch));
        
        SC = block11 - BT*(block22\B); %SC w.r.t block 1,1
        
        cond_bblock1(i,j) = cond(full(block11));
        cond_bblock2(i,j) = cond(full(SC));
        
        %do 3-rd level of SC
        block11 = SC(1:18,1:18);
        block22 = SC(19:end, 19:end);
        BT = SC(1:18, 19:end);
        B = SC(19:end, 1:18);
        
        SC = block22 - B*(block11\BT);
        cond_bblock3(i,j) = cond(full(SC));
    end
    
    %remove IPOPT matrices
    delete *.iajaa;
end

cond_KKT
cond_hess
cond_block
cond_pblock
cond_bblock1
cond_bblock2

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