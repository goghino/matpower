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


%% specify list of branch contingencies, equals to OPF if empty
cond_KKT = zeros(nbrch,1);
cond_block = zeros(nbrch,2);

for i = 1:nbrch
    cont = [i]'; 
    ns = size(cont,1) + 1; %number of scenarios (ncont + nominal)

    %results = runopf(mpc, mpopt);
    results = runscopf(mpc, cont, mpopt);

    % inspect IPOPT hessian
    A = readcsr('mat-ipopt_005-01.iajaa', 0, 1); 
    [P Pinv] = KKTpermute(nbus, nbrch, ngen, ns);
    AP = A(P,P'); 
    
%     figure; spy(AP)
%     figure; spy(A)

    cond_KKT(i,1) = cond(full(A));

    %local hessians wrt primal variables for nominal and contingency
    nH = 2*nbus;
    cond_block(i,1) = cond(full(A(1:nH,1:nH)));
    cond_block(i,2) = cond(full(A(nH+(1:nH),nH+(1:nH))));
    
    %remove IPOPT matrices
    delete *.iajaa;
end

cond_KKT
cond_block

%isolating bus by contingency creates singular Cb, zero on diagonal
% -> singular local hessian

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
%4. could it be that Jac of constraints is singular? reason could be same
%   branch flow limits and very similar voltages on the attached branches
%   so that rows becomes almost the same