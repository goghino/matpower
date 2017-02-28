classdef PowerGridHessian < handle
    
    properties
        
        KKT              % the original matrix
        H                % the reordered matrix
        b                % the rhs
        x                % the solution vector
        makeSpyPlots     % for debugging (not recommended for large problems)
        iperm            % inverse permutation
        
        nrows            % rows of KKT Hessian
        npart            % rows of each diagonal block
        nb               % number of buses
        ng               % number of generators
        nl               % number of lines
        ns               % number of storages
        N                % number of timesteps
        
        %%% ------------------------------------------- %%%
        %%% H              stored for Schur application %%%
        %%% ------------------------------------------- %%%
        A                % H diagonal blocks
        B                % H Bs
        D                % H D
        S                % Schur-complement
        
        %%% ------------------------------------------- %%%
        %%% Preconditioner stored for Schur application %%%
        %%% ------------------------------------------- %%%
        PA               % preconditioners diagonal blocks
        PB               % preconditioners Bs
        PD               % preconditioners D
        
    end
    
    
    
    methods
        
        
        function this = PowerGridHessian(filename, nb, ng, nl, ns, N, show_spy_plots)
            
            close all;
            clc;
            
            this.KKT          =[];
            this.H            =[];
            this.x            =[];
            this.b            =[];
            this.iperm        =[];
            
            this.nb           = nb;
            this.ng           = ng;
            this.nl           = nl;
            this.ns           = ns;
            this.N            = N;
            this.makeSpyPlots = show_spy_plots;
            
            % read the Hessian stored in the file filename, in symmtric CSR format
            
            tstart = tic;
            this.KKT          = readcsr(filename, 0, 1);
            
            %x1               = IPOPTHessian\b;
            t1                = toc(tstart);
            
            fprintf('time needed for reading Hessian is %f10.2\n',t1);
            this.nrows        = size(this.KKT,1);
            this.x            = rand(this.nrows, 1);
            this.b            = this.KKT*this.x;
            
            
            
            if (this.makeSpyPlots == true)
                spy(this.KKT, 'k');
            end
            
            %  [  H     0     J^T     P^T   ]
            %  [                            ]
            %  [  0     K     0      -S     ]
            %  [                            ]
            %  [  J     0     0       0     ]
            %  [                            ]
            %  [  P    -S     0       0     ]
            
            
            
            %                                                      -1
            %          [ H    J^T ]       [  0   P^T ] [  K   -S  ]   [  0    0  ]
            %  Schur = [          ]   -   [          ] [          ]   [          ]
            %          [ J    0   ]       [  0   0   ] [ -S    0  ]   [  P    0  ]
            %
            %
            
            %                                                        -1
            %          [ H    J^T ]       [  0   P^T ] [  0        -S   ]   [  0    0  ]
            %  Schur = [          ]   -   [          ] [   -1   -1   -1 ]   [          ]
            %          [ J    0   ]       [  0   0   ] [ -S   -S  K S   ]   [  P    0  ]
            %
            %
            
            %                                                        -1
            %          [ H    J^T ]       [  0   P^T ] [  0        -S   ]   [  0    0  ]
            %  Schur = [          ]   -   [          ] [   -1   -1   -1 ]   [          ]
            %          [ J    0   ]       [  0   0   ] [ -S   -S  K S   ]   [  P    0  ]
            
            
            %                                               -1
            %          [ H    J^T ]       [  0   P^T ] [  -S  P      0   ]
            %  Schur = [          ]   -   [          ] [   -1   -1       ]
            %          [ J    0   ]       [  0   0   ] [ -S  K S  P  0   ]
            
            
            %                                    -1   -1
            %          [ H    J^T ]       [ -P^T S  K S  P   0   ]
            %  Schur = [          ]   -   [                      ]
            %          [ J    0   ]       [  0               0   ]
            
            
            
            %                     -1   -1            T
            %          [  H + P^T S  K S  P         J   ]
            %  Schur = [                                ]
            %          [  J                         0   ]
            
            
            
            %  [  H     J^T     P^T    0     ]
            %  [                             ]
            %  [  J     0       0      0     ]
            %  [                             ]
            %  [  P     0       0     -S     ]
            %  [                             ]
            %  [  0     0      -S      K     ]
            
            
            
            fprintf('call now this.SolveSchur(reordering), where reordering is one of the following...\n');
            fprintf('reordering options\n');
            fprintf('------------------\n');
            fprintf('noAInDiagBlocks\n');
            fprintf('nsAInDiagBlocks\n');
            fprintf('2nsAInDiagBlocks\n');
            fprintf('------------------\n');
            
            
            %fprintf('The difference of the 2 solutions: || x1 - x2 ||_2 = %25.16e\n', norm(x1-x2)/norm(x));
            %fprintf('The difference of the 2 solutions: || x  - x1 ||_2 = %25.16e\n', norm(x -x1)/norm(x));
            %fprintf('The difference of the 2 solutions: || x  - x2 ||_2 = %25.16e\n', norm(x -x2)/norm(x));
            %fprintf('The difference of the 2 solutions: ||r|| / ||b||   = %25.16e\n', norm(IPOPTHessian*(x -x2))/norm(b));
            %fprintf('The difference of the 2 solutions:         ||b||   = %25.16e\n', norm(b));
            
        end
        
        
        function y = applyMatrix(this, x)
            
            %  [  A     B  ]
            %  [           ]
            %  [  B'    D  ]
            %
            % y  := ( D - B' * A^{-1} * B )*x
            
            N     = this.N;
            npart = this.npart;
            
            % S := D
            nrows  = size(this.H,1);
            I      = N*npart+1:nrows;
            y      = this.H(I, I)*x;
            
            
            % form the Schur block by block
            for n=0:N-1
                I    = n*npart+1:(n+1)*npart;
                H_nn = this.H(I, I);
                BnT  = this.H(I, N*npart+1:end);
                Bn   = this.H(N*npart+1:end, I);
                y = y - Bn * ( H_nn \ (BnT*x) );
            end
            
            %A  = this.H(1:N*npart, 1:N*npart);
            %B  = this.H(1:N*npart, N*npart+1:end);
            %S  = this.H(N*npart+1:end, N*npart+1:end) - B'*(A\B);
            %y = S*x;
        end
        
        
        
        
        
        function x = applyPrec(this, b)
            
            
            N     = this.N;
            npart = this.npart;
            
            
            H_nn  = zeros(npart, npart);
            I     = this.index(npart, N);
            % I     = index(npart, 1);
            H_p   = this.H(I, I);
            for n=0:N-1
                I    = n*npart+1:(n+1)*npart;
                H_nn = H_nn + this.H(I, I)/N;
            end
            
            H_nn = H_p;
            % S  =  D - \sum_n Bn' * Ann^{-1} * Bn
            S       = this.H(N*npart+1:end, N*npart+1:end);
            for n=0:N-1
                I    = this.index(npart, n+1);
                %H_nn = this.H(I, I);
                BnT  = this.H(I, N*npart+1:end);
                Bn   = this.H(N*npart+1:end, I);
                S  = S - Bn * (H_nn \ BnT);
            end
            
            x = S \ b;
        end
        
        
        
        
        
        function x = applyGlobalPrec(this, b)
            
            
            N     = this.N;
            npart = this.npart;
            
            
            H_nn  = zeros(npart, npart);
            I     = this.index(npart, N);
            % I     = index(npart, 1);
            H_p   = this.H(I, I);
            for n=0:N-1
                I    = n*npart+1:(n+1)*npart;
                H_nn = H_nn + this.H(I, I)/N;
            end
            
            H_nn = H_p;
            
            % S  =  D - \sum_n Bn' * Ann^{-1} * Bn
            S       = this.H(N*npart+1:end, N*npart+1:end);
            for n=0:N-1
                I    = this.index(npart, n+1);
                %H_nn = this.H(I, I);
                BnT  = this.H(I, N*npart+1:end);
                Bn   = this.H(N*npart+1:end, I);
                S  = S - Bn * (H_nn \ BnT);
            end
            
            I    = this.index(npart, N);
            x(I) = S \ b(I);
            
            for n=0:N-1
                I    = this.index(npart, n+1);
                BnT  = this.H(I, N*npart+1:end);
                y    = b(I) - BnT*x2;
                x(I) = this.H(I, I) \ y;
            end
            
            
        end
        
        
        
        
        function J = index(this, npart, i)
            
            J = (i-1)*npart + 1:i*npart;
            
        end
        
        function Scale = ScaleSystem(this)
            
            S = speye(this.nrows);
            for i=1:1
                fprintf('computing scaling  ...\n');
                s   = (max(abs(this.KKT)));
                s(s < 1e-12) = 1.0;
                s   = 1.0./s;
                
                fprintf('before the scaling kappa(KKT) %e\n', condest(this.KKT));
                %s                 = load('MC64mat-ipopt_025-01.scales')';
                %s                 = load('MC64330mat-ipopt_002-01.iajaa')';
                %Scale             = speye(this.nrows, this.nrows);
                Scale             = spdiags(s', 0, this.nrows, this.nrows);
                fprintf('scaling the matrix ...\n');
                this.KKT          = Scale*this.KKT;
                fprintf('scaling the rhs b  ...\n');
                this.b            = Scale*this.b;
                fprintf('after  the scaling kappa(KKT) %e\n', condest(this.KKT));
                S = Scale*S;
            end
            Scale = S;
        end
        
        
        
        
        function [H, x] = SolveGlobal(this, reordering)
            
            %Scale = ScaleSystem(this);
            %sc=full(diag(Scale));
            %fprintf('max and min scales = %25.16f %25.16f\n',max(sc),min(sc));
            p = computePermutation(this, reordering);
            %p  = amd(IPOPTHessian);
            fprintf('size of permutation array is %d and size of rhs b is %d\n',length(p), length(this.b));
            
            b        = this.b;
            b(p)  = this.b;
            this.H  = this.reorder(p, this.KKT);
            
            if (this.makeSpyPlots == true)
                figure; spy(this.H);
            end
            
            % reorder it back using the inverse permutation
            H0 = this.reorderBack(p, this.H);
            
            
            % and test if we get the same matrix
            y = rand(this.nrows, 1);
            fprintf('norm of H0 and KKT = ||H0*y - KKT*y|| = %25.16e\n', norm(this.KKT*y - H0*y));
            
            fprintf('H0 is %d by %d\n', this.nrows, this.nrows);
            fprintf('H0 is %d by %d\n', size(this.H,1), size(this.H,2));
            
            N     = this.N;
            npart = this.npart;
            
            P     = sparse(this.nrows, this.nrows);
            S     = this.H(N*npart+1:end, N*npart+1:end);
            nS    = size(S,1);
            
            H_nn = sparse(npart,npart);
            I    = 1:npart; %(N-1)*npart+1:N*npart; %1:npart;
            H_00 = this.H(I, I);
            B    = this.H(I, N*npart+1:end);
            B0   = this.H(I, N*npart+1+this.ns:N*npart+4*this.ns);
            alpha = sparse(npart,npart);
            P = this.H;
             if (this.makeSpyPlots == true)
                figure; spy([H_00  B; B' , sparse(nS,nS)]);
            end
            %H_00 = H_nn / N;
            for n=0:N-1
                I             = n*npart+1:(n+1)*npart;
                H_n           = this.H(I, I);
                for k=1:npart
                  alpha(k,k)  = sum(abs(H_n(k,:))/ abs(H_00(k,:)));
                end

                P(I, I)       = alpha*H_00;
                fprintf('||H_%2d - alpha*H_00||_1 = %e  relative  %e\n', n, norm(P(I,I) - H_n, 'fro'), norm(P(I,I) - H_n, 'fro')/norm(H_n, 'fro'));
            end
            
            
            fprintf('condition numbers kappa(KKT), kappa(H), kappa(P) %e %e %e\n', condest(this.KKT), condest(this.H), condest(P));
            x0 = this.x;
            x  = zeros(this.nrows, 1);
            [x,flag,relres,iter,resvec] = gmres(this.H,b,135,1e-4,1,P); %MATLAB native
            %[x,flag,relres,iter,resvec] = cgs(this.H,b,1e-16,150,P); % MATLAB native
            %[x,flag,relres,iter,resvec] = bicgstab(this.H,b,1e-16,150,P); % MATLAB native
            %fprintf('iterative method finished with flag %4d\n', flag);
            fprintf('iter %4d  %25.16e\n', [1:length(resvec) ; resvec' ]);
            y = x(p);
            %x = Scale*y;
            %y=x;
            fprintf('|x - x0|      = %26.15e\n', norm(y - x0));
            fprintf('|x - x0|/|x0| = %26.15e\n', norm(y - x0)/norm(x0));
            fprintf('|Hx - b|      = %26.15e\n', norm(this.H*x - b));
            fprintf('|Hx - b|/|b|  = %26.15e\n', norm(this.H*x - b)/norm(b));
            fprintf('|Ax - b|      = %26.15e\n', norm(this.KKT*y - this.b));
            fprintf('|Ax - b|/|b|  = %26.15e\n', norm(this.KKT*y - this.b)/norm(this.b));
            fprintf('|Ax0- b|      = %26.15e\n', norm(this.KKT*x0 - this.b));
            fprintf('|Ax0- b|/|b|  = %26.15e\n', norm(this.KKT*x0 - this.b)/norm(this.b));
            
            %    Computing Schur after reordering KKT-component-wise, and then
            %    reordering again per timestep
            
        end
        
        
        
        
        
        function [H, x] = CheckStructure(this, reordering)
            
            %Scale = ScaleSystem(this);
            %sc=full(diag(Scale));
            %fprintf('max and min scales = %25.16f %25.16f\n',max(sc),min(sc));
            p = computePermutation(this, reordering);
            %p  = amd(IPOPTHessian);
            fprintf('size of permutation array is %d and size of rhs b is %d\n',length(p), length(this.b));
            
            ns      = this.ns;
            b       = this.b;
            b(p)    = this.b;
            this.H  = this.reorder(p, this.KKT);
            
            if (this.makeSpyPlots == true)
                figure; spy(this.H,'k');
            end
            
            % reorder it back using the inverse permutation
            H0 = this.reorderBack(p, this.H);
            
            
            % and test if we get the same matrix
            y = rand(this.nrows, 1);
            fprintf('norm of H0 and KKT = ||H0*y - KKT*y|| = %25.16e\n', norm(this.KKT*y - H0*y));
            
            fprintf('H0 is %d by %d\n', this.nrows, this.nrows);
            fprintf('H0 is %d by %d\n', size(this.H,1), size(this.H,2));
            
            N     = this.N;
            npart = this.npart;
            
            P     = sparse(this.nrows, this.nrows);
            S     = this.H(N*npart+1:end, N*npart+1:end);
            nS    = size(S,1);
            
            H_nn = sparse(npart,npart);
            I    = 1:npart; %(N-1)*npart+1:N*npart; %1:npart;
            H_00 = this.H(I, I);
            B    = this.H(I, N*npart+1:end);
            B0   = this.H(I, N*npart+1+this.ns:N*npart+4*this.ns);
            alpha = sparse(npart,npart);
            P = this.H;
             if (this.makeSpyPlots == true)
                figure; spy([H_00  B; B' , sparse(nS,nS)], 'k');
            end
            %H_00 = H_nn / N;
            for n=0:N-1
                I             = n*npart+1:(n+1)*npart;
                H_n           = this.H(I, I);
                B             = this.H(I, N*npart+1:end);
                %spy(B);
                %pause;
                for i=1:N+1
                  B_i         = B(:,     (i-1)*ns+1:i*ns);
                  fprintf('|B_%2d| = %26.15e\n', i, norm(full(B_i)));
                end
                fprintf('============== n = %d\n', n);
                %for k=1:N
                %  B           = this.H(I, N*npart+1:end);
                %end
                %fprintf('|B - B0|      = %26.15e\n', norm(y - x0));
            end           
            %    Computing Schur after reordering KKT-component-wise, and then
            %    reordering again per timestep
            
        end
        
        
        
        
        % nb=118; ng=72; nl=186; N=24; ns=9; pg=PowerGridHessian(file, % nb,ng,nl,ns,N, false)
        % file='matrices/case118factor1N24/mat-ipopt_033-01.iajaa'
        function out = SolveSchur(this, reordering)
            clc;
            
            % Scale = ScaleSystem(this);
            % sc=full(diag(Scale));
            % fprintf('max and min scales = %25.16f %25.16f\n',max(sc),min(sc));
            p = computePermutation(this, reordering);
            b       = this.b;
            b(p)    = this.b;
            this.H  = this.reorder(p, this.KKT);
            
            spy(this.H);
            if (this.makeSpyPlots == true)
                figure; spy(this.H, 'k');
            end
            
            return;
            % reorder it back using the inverse permutation
            H0 = this.reorderBack(p, this.H);
            
            % and test if we get the same matrix
            y = rand(this.nrows, 1);
            fprintf('norm of H0 and KKT = ||H0*y - KKT*y|| = %25.16e\n', norm(this.KKT*y - H0*y));
            
            fprintf('H0 is %d by %d\n', this.nrows, this.nrows);
            fprintf('H0 is %d by %d\n', size(this.H,1), size(this.H,2));
            
            N     = this.N;
            npart = this.npart;
            
            
            b1    = b(1:N*npart,1);
            b2    = b(N*npart+1:end,1);
            
            % r2 := b2 - B^T A^{-1} b1
            % ------------------------
            r2    = b2;
            % form the Schur block by block
            S   = this.H(N*npart+1:end, N*npart+1:end);
            Sp  = this.H(N*npart+1:end, N*npart+1:end);
            I   = (N-1)*npart+1:N*npart;
            H_p = this.H(I, I);
            for n=0:N-1
                I    = n*npart+1:(n+1)*npart;
                H_nn = this.H(I, I);
                BnT  = this.H(I, N*npart+1:end);
                Bn   = this.H(N*npart+1:end, I);
                
                if (n == 0)
                    D = [H_nn BnT; Bn 0*S];
                    spy(D);pause
                end
                
                r2 = r2 - Bn * ( H_nn \ b1(I));
                
                S  = S  - Bn * (H_nn \ BnT);
                Sp = Sp - Bn * (H_p \ BnT);
            end
            
            x2 = zeros(length(r2),1);
            
            [S1, Sp1] = this.formSchur();
            %[x2,error_norm,iter,flag] = gmres(@(x)applyMatrix(this,x),x2,r2,S,100,1,1e-15);
            %[x2,error_norm,iter,flag] = gmres(@(x)applyMatrix(this,x),x2,r2,@(x)applyPrec(this,x),50,1,1e-15);
            %[x2,error_norm,iter,flag] = cgs(@(x)applyMatrix(this,x),x2,r2,@(x)applyPrec(this,x),30,1e-16);
            
            
            [x2,flag,relres,iter,resvec] = gmres(S,r2,12,1e-14,1,Sp);
            %[x2,flag,relres,iter,resvec] = gmres(this.S,r2,100,1e-14,1,@(x)applyPrec(this,x));
            %[x2,flag,relres,iter,resvec] = cgs(@(x)applyMatrix(this,x),r2,1e-14,100,@(x)applyPrec(this, x));
            %[obj.x,flag,relres,iter,resvec] = pcg(obj.A,obj.b,1e-8,100,@(x)precondition(obj, x));
            fprintf('iterative method finished with flag %4d\n', flag);
            fprintf('iter %4d  %25.16e\n', [1:length(resvec) ; resvec' ]);
            
            
            % x1 := A^{-1} ( b1 - B^T x2 )
            % ----------------------------
            x1    = b1;
            for n=0:N-1
                I     = n*npart+1:(n+1)*npart;
                BnT   = this.H(I, N*npart+1:end);
                H_nn  = this.H(I, I);
                x1(I) = H_nn \ ( b1(I) - BnT * x2 );
            end
            x     = [x1; x2];
            x0    = this.x;
            y     = x(p);
            %x     = Scale*y;
            x     = y;
            
            fprintf('|x - x0|      = %26.15e\n', norm(x - x0));
            fprintf('|x - x0|/|x0| = %26.15e\n', norm(x - x0)/norm(x0));
            fprintf('|Ax - b|      = %26.15e\n', norm(this.KKT*y - this.b));
            fprintf('|Ax - b|/|b|  = %26.15e\n', norm(this.KKT*y - this.b)/norm(this.b));
            %    Computing Schur after reordering KKT-component-wise, and then
            %    reordering again per timestep
            
        end
        
        
        
        function [S, Sp] = formSchur(this)
            
            N     = this.N;
            npart = this.npart;
            ns    = this.ns;
            nS    = size(this.H,1) - N*npart;
            fprintf('The size of the Schur complement is %d x %d\n', nS, nS);
            S     = this.H(N*npart+1:end, N*npart+1:end);
            Sp    = this.H(N*npart+1:end, N*npart+1:end);
            
            I     = (N-1)*npart+1:N*npart;
            H_p   = this.H(I, I);
            
            for n=0:N-1
                I    = n*npart+1:(n+1)*npart;
                H_nn = this.H(I, I);
                BnT  = this.H(I, N*npart+1:end);
                Bn   = this.H(N*npart+1:end, I);
                S  = S  - Bn * (H_nn \ BnT);
                Sp = Sp - Bn * (H_p \ BnT);
            end
            
            S0     = this.H(N*npart+1:end, N*npart+1:end);
            S0p    = this.H(N*npart+1:end, N*npart+1:end);
            
            m = size(S0,1);
            S0(ns+1:m-ns, ns+1:m-ns) = zeros(m-2*ns, m-2*ns);
            S0(ns+1:m-ns, ns+1:m-ns) = zeros(m-2*ns, m-2*ns);
            
            
            for n=0:N-1
                C    = SchurComplement(this, n);
                S0   = S0  -  C;
            end
            
            
            b = rand(m,1);
            [x,flag,relres,iter,resvec] = gmres(S,b,100,1e-14,1,S0);
            fprintf('iter %4d  %25.16e\n', [1:length(resvec) ; resvec' ]);
            fprintf('|| b - S * x ||_2 = %25.16e\n', norm(b-S*x)/norm(b));
            fprintf('k(S0\\S)        = %25.16e\n', cond(S0\S));
            fprintf('||S ||_inf     = %25.16e\n', norm(S,  inf));
            fprintf('||S0||_inf     = %25.16e\n', norm(S0, inf));
            fprintf('||S - S0||_inf = %25.16e\n', norm(S-S0, inf));
            
            [L, D] = this.SchurBlockLDL(S(2*ns+1:end, 2*ns+1:end));
            %[L, D] = this.SchurLDL(S(2*ns+1:end, 2*ns+1:end));
        end
        
        
        function x = SolveSchurComplement(this, S, b)
            
            ns      = this.ns;
            nblocks = size(S, 1) / ns;
            
            % [  A11    A12 ]
            % |             |
            % [  A21    H   ]
            
            A11     = S(1:2*ns, 1:2*ns);
            A12     = S(1:2*ns, 2*ns+1:end);
            A21     = S(2*ns+1:end, 1:2*ns);
            
            H = zeros((nblocks-2)*ns, 3*ns);
            
            %     [ H11 H12 H12 H12 H12] [ X1 ]   [ R1 ]
            %     | H12 H21 H22 H22 H22| | X2 |   | R2 |
            % H = | H12 H22 H31 H32 H32| | X3 | = | R3 |
            %     | H12 H22 H32 H41 H42| | X4 |   | R4 |
            %     [ H12 H22 H32 H42 H51] [ X4 ]   [ R4 ]
            
            %     [ H11 H12 ]
            %     | H21 H22 |
            % H = | H31 H32 |
            %     | H42 H42 |
            %     [ H51     ]
            
            I1 = 1:ns;
            I2 = I1 + ns;
            I3 = I2 + ns;
            for i=3:nblocks
                I                    = (i-1)*ns+1:i*ns;
                S_ii                 = S(I, I);
                S_iip1               = S(I, I + ns);
                
                % Hi1
                H(I-2*ns, I1)        = Sii;
                % Hi2
                H(I-2*ns, I2)        = Siip1;
            end
            
            for i=1:nblocks-2
                I        = (i-1)*ns+1:i*ns;
                
                % Hi3    = Hi1 \ Hi2
                H(I, I3) = H(I, I1) \ H(I, I2);
                % A21i   = Hi1 \ A21i
                A21(I,:) = H(I, I1) \ A21(I, :);
                
                for j=i+1:nblocks-2
                    J         = (j-1)*ns+1:j*ns;
                    
                    % Hj -= Hi2*Hi3
                    H(J, :)   = H(J, :)   - H(I, I2)*H(I, I3);
                    A21(J, :) = A21(J, :) - H(I, I2)*A21(I, :);
                end
            end
            
            %     [  I  H1  H1  H1  H1 ] [ X1 ]   [ R1 ]
            % ~   |     I   H2  H2  H2 | | X2 |   | R2 |
            % H = |         I   H3  H3 | | X3 | = | R3 |
            %     |             I   H4 | | X4 |   | R4 |
            %     [                 I  ] [ X4 ]   [ R4 ]
            
            %           _   _   _   _
            %     [ I   H12 H12 H12 H12] [ X1 ]   [ H11^-1 R1 ]
            %     |     H21 H22 H22 H22| | X2 |   | R2 - H12 H11^-1 R1 |
            % H = |     H22 H31 H32 H32| | X3 | = | R3 - H12 H11^-1 R1 |
            %     |     H22 H32 H41 H42| | X4 |   | R4 - H12 H11^-1 R1 |
            %     [     H22 H32 H42 H51] [ X4 ]   [ R4 - H12 H11^-1 R1 ]
            
            %     [ H11 H12  H11^-1 H12 ]
            %     | K21 K22  K21^-1 K22 |
            % D = | L31 L32  L31^-1 L32 |
            %     | M41 M42  M41^-1 M42 |
            %     [ P51      P51        ]
            
            
            for i=nblocks-3:-1:3
                I            = (i-1)*ns+1:i*ns;
                D            = zeros(ns, ns);
                for j=i+1:nblocks-2
                    J = (j-1)*ns+1:j*ns;
                    D = D + A21(J, :);
                end
                A21(I, :) = A21(I, :) - H(I, I3)*D;
            end
        end
        
        
        
        function [B] = SchurBlockLU(this, A)
            
            ns = this.ns;
            n  = size(A,1) / ns;
            
            for k=1:n-1
                J        = k+1:n
                A(k+1,k) = A(k+1,k) / A(k,k);
                for j=k+1:n
                    %A(j,k) = A(j,k) / A(k,k);
                    A(j,k) = A(k+1,k);
                end
                
                a      = A(k+1,k)*A(k,k+1);
                D      = repmat(a, n-k);
                A(J,J) = A(J,J) - D; %A(J,k) * A(k,J);
            end
            
            L = A;
            
        end
        
        
        
        
        
        function C = SchurComplement(this, n)
            N     = this.N;
            npart = this.npart;
            ns    = this.ns;
            nS    = size(this.H,1) - N*npart;
            
            I    = n*npart+1:(n+1)*npart;
            H_nn = this.H(I, I);
            BnT  = this.H(I, N*npart+1:end);
            Bn   = this.H(N*npart+1:end, I);
            
            C = zeros(nS, nS);
            b3=[BnT(:, ns+1:2*ns)  BnT(:, (n+1)*ns+1:(n+3)*ns)];
            b = b3' * (H_nn \ b3);
            
            
            b11 = b(1:ns,        1:ns);
            b12 = b(1:ns,        ns+1:2*ns);
            b13 = b(1:ns,        2*ns+1:3*ns);
            b31 = b(2*ns+1:3*ns, 1:ns);
            b22 = b(ns+1:2*ns,   ns+1:2*ns);
            b23 = b(ns+1:2*ns,   2*ns+1:3*ns);
            b32 = b(2*ns+1:3*ns, ns+1:2*ns);
            b33 = b(2*ns+1:3*ns, 2*ns+1:3*ns);
            
            C(ns+1:2*ns, ns+1:2*ns)                     = b11;
            C((n+1)*ns+1:(n+2)*ns, (n+1)*ns+1:(n+2)*ns) = b22;
            
            C(ns+1:2*ns, (n+1)*ns+1:(n+2)*ns) = b12;
            C((n+1)*ns+1:(n+2)*ns, ns+1:2*ns) = b12';
            
            nblocks = nS / ns;
            for k=n:nblocks-3
                J1          = (k+2)*ns+1:(k+3)*ns;
                C(J1, ns+1:2*ns) = b31;
                C(ns+1:2*ns, J1) = b13;
                
                C(J1, (n+1)*ns+1:(n+2)*ns) = b32;
                C((n+1)*ns+1:(n+2)*ns, J1) = b23;
                for l=n:nblocks-3
                    J2        = (l+2)*ns+1:(l+3)*ns;
                    C(J1, J2) = b33;
                end
            end
        end
        
        
        function p = computePermutation(this, reordering)
            
            nb     = this.nb; % number of buses
            ng     = this.ng; % number of generators
            nl     = this.nl; % number of lines
            ns     = this.ns; % number of storages
            N      = this.N;  % number of timesteps
            
            
            NB     = nb * N;
            NG     = ng * N;
            NL     = nl * N;
            
            NH     = 2*NB + 2*NG;
            NJ     = 2*NB;
            NP     = 2*NL + (N+1)*ns;
            NK     = 2*NL + (N+1)*ns;
            
            
            %  [  H     0     J^T     P^T   ]
            %  [                            ]
            %  [  0     K     0      -S     ]
            %  [                            ]
            %  [  J     0     0       0     ]
            %  [                            ]
            %  [  P    -S     0       0     ]
            
            %  [  H     J^T     P^T    0     ]
            %  [                             ]
            %  [  J     0       0      0     ]
            %  [                             ]
            %  [  P     0       0     -S     ]
            %  [                             ]
            %  [  0     0      -S      K     ]
            
            p      = [];
            %%%%%%%%%
            % H  part
            % %%%%%%%
            
            if (strcmp(reordering, 'noAInDiagBlocks') == 1)
                for i = 0:N-1
                    p = [p   i*nb+[1:nb]   NB+i*nb+[1:nb]   2*NB+i*ng+[1:ng]   2*NB+NG+i*ng+[1:ng] ]; % H
                    p = [p   NH+NK+i*nb+[1:nb]       NH+NK+NB+i*nb+[1:nb]    ];                       % J
                    p = [p   NH+NK+NJ+i*nl+[1:nl]    NH+NK+NJ+NL+i*nl+[1:nl] ];                       % P
                    p = [p   NH+i*nl+[1:nl]          NH+NL+i*nl+[1:nl]       ];                       % K
                end
                p = [p   NH+2*NL+1:NH+NK    NH+NK+NJ+2*NL+1:NH+NK+NJ+NP ];
                this.npart = 2*nb + 2*ng + 2*nb + 4*nl;
            end
            
            if (strcmp(reordering, 'nsAInDiagBlocks') == 1)
                for i = 0:N-1
                    p = [p   i*nb+[1:nb]   NB+i*nb+[1:nb]   2*NB+i*ng+[1:ng]   2*NB+NG+i*ng+[1:ng] ]; % H
                    p = [p   NH+NK+i*nb+[1:nb]       NH+NK+NB+i*nb+[1:nb]    ];                       % J
                    p = [p   NH+NK+NJ+i*nl+[1:nl]    NH+NK+NJ+NL+i*nl+[1:nl] ];                       % P
                    p = [p   NH+i*nl+[1:nl]          NH+NL+i*nl+[1:nl]       ];                       % K
                    p = [p   NH+2*NL+1 + i*ns:NH+2*NL+(i+1)*ns   ];
                end
                p = [p   NH+2*NL+N*ns+1:NH+NK    NH+NK+NJ+2*NL+1:NH+NK+NJ+NP ];
                this.npart = 2*nb + 2*ng + 2*nb + 4*nl + ns;
            end
            
            if (strcmp(reordering, '2nsAInDiagBlocks') == 1)
                for i = 0:N-1
                    p = [p   i*nb+[1:nb]   NB+i*nb+[1:nb]   2*NB+i*ng+[1:ng]   2*NB+NG+i*ng+[1:ng] ]; % H
                    p = [p   NH+NK+i*nb+[1:nb]       NH+NK+NB+i*nb+[1:nb]    ];                       % J
                    p = [p   NH+NK+NJ+i*nl+[1:nl]    NH+NK+NJ+NL+i*nl+[1:nl] ];                       % P
                    p = [p   NH+i*nl+[1:nl]          NH+NL+i*nl+[1:nl]       ];                       % K
                    p = [p   NH+2*NL+1 + i*ns:NH+2*NL+(i+1)*ns   ];
                    p = [p   NH+NK+NJ+2*NL+ns + i*ns+1:NH+NK+NJ+2*NL+ns + (i+1)*ns   ];               % Ps
                end
                p=[p   NH+2*NL+N*ns+1:NH+NK    NH+NK+NJ+2*NL+1:NH+NK+NJ+2*NL+ns ];
                this.npart = 2*nb + 2*ng + 2*nb + 4*nl + 2*ns;
            end
            
            %  [  H     0     J^T     P^T   ]
            %  [                            ]
            %  [  0     K     0      -S     ]
            %  [                            ]
            %  [  J     0     0       0     ]
            %  [                            ]
            %  [  P    -S     0       0     ]
            
            if (strcmp(reordering, 'Klaus') == 1)
                p=[2*NB+1:NH+NK   1:2*NB   NH+NK+1:NH+NK+NJ  NH+NK+NJ+1:NH+NK+NJ+NP ];
            end
            
            if (strcmp(reordering, 'Klaus1') == 1)
                p=[2*NB+2*NB+2*NL:-1:2*NB+1 1:2*NB 4*NB+2*NL+1:4*NB+2*NL+2*(N+1)*ns ];
            end
            
            this.iperm=p;
            for i=1:length(p)
                p(this.iperm(i)) = i;
            end
            
            fprintf('Permutation0 array has numel #%d %d\n',length(p));
            
            
        end
        
        
        
        
        
        function p = computePermutation1(nrows, nb, ng, nl, ns, N)
            
            
            NA     = ns * (N+1);
            NB     = nb * N;
            NG     = ng * N;
            NL     = nl * N;
            
            NH     = 2*NB + 2*NG;
            NJ     = 2*NB;
            NP     = 2*NL + (N+1)*ns;
            NK     = 2*NL + (N+1)*ns;
            
            
            %  [  H     0     J^T     P^T   ]
            %  [                            ]
            %  [  0     K     0      -S     ]
            %  [                            ]
            %  [  J     0     0       0     ]
            %  [                            ]
            %  [  P    -S     0       0     ]
            
            %  [  H     J^T     P^T    0     ]
            %  [                             ]
            %  [  J     0       0      0     ]
            %  [                             ]
            %  [  P     0       0     -S     ]
            %  [                             ]
            %  [  0     0      -S      K     ]
            
            p      = [];
            %%%%%%%%%
            % H  part
            % %%%%%%%
            
            %    for i = 0:N-1
            %      p = [p   i*nb+[1:nb]   NB+i*nb+[1:nb]   2*NB+i*ng+[1:ng]   2*NB+NG+i*ng+[1:ng] ]; % H
            %      p = [p   NH+NK+i*nb+[1:nb]       NH+NK+NB+i*nb+[1:nb]    ];                       % J
            %      p = [p   NH+NK+NJ+i*nl+[1:nl]    NH+NK+NJ+NL+i*nl+[1:nl] ];                       % P
            %      p = [p   NH+i*nl+[1:nl]          NH+NL+i*nl+[1:nl]       ];                       % K
            %    end
            %    p=[p   NH+2*NL+1:NH+NK    NH+NK+NJ+2*NL+1:NH+NK+NJ+NP ];
            
            p = [p   1:NH   NH+NK+1:NH+NK+NJ   NH+NK+NJ+1:NH+NK+NJ+NP   NH+1:NH+NK ]; % H
            iperm=p;
            for i=1:length(p)
                p(iperm(i)) = i;
            end
            
            fprintf('Permutation1 array has numel #%d %d\n',length(p));
            
            
        end
        
        
        
        function p = computePermutation2(nrows, nb, ng, nl, ns, N)
            
            
            NA     = ns * (N+1);
            NB     = nb * N;
            NG     = ng * N;
            NL     = nl * N;
            
            NH     = 2*NB + 2*NG;
            NJ     = 2*NB;
            NP     = 2*NL + (N+1)*ns;
            
            
            %  [  H     0     J^T     P^T   ]
            %  [                            ]
            %  [  0     K     0      -S     ]
            %  [                            ]
            %  [  J     0     0       0     ]
            %  [                            ]
            %  [  P    -S     0       0     ]
            
            %  [  H     J^T     P^T    0     ]
            %  [                             ]
            %  [  J     0       0      0     ]
            %  [                             ]
            %  [  P     0       0     -S     ]
            %  [                             ]
            %  [  0     0      -S      K     ]
            
            p      = [];
            %%%%%%%%%
            % H  part
            % %%%%%%%
            
            for i = 0:N-1
                p = [p   i*nb+[1:nb]   NB+i*nb+[1:nb]   2*NB+i*ng+[1:ng]   2*NB+NG+i*ng+[1:ng] ]; % H
                p = [p   NH+i*nb+[1:nb]       NH+NB+i*nb+[1:nb]    ];                             % J
                p = [p   NH+NJ+i*nl+[1:nl]    NH+NJ+NL+i*nl+[1:nl] ];                             % P
                p = [p   NH+NJ+2*NL+i*ns+1:NH+NJ+2*NL+(i+1)*ns ];                               % P
            end
            p = [p   NH+NJ+2*NL+N*ns+1:NH+NJ+NP];                                             % P
            
            iperm=p;
            for i=1:length(p)
                p(iperm(i)) = i;
            end
            
            fprintf('Permutation1 array has numel #%d %d\n',length(p));
            
            
        end
        
        
        
        function B = reorder(this, permutation, A)
            [rows,cols,vals] = find(A);
            
            nonzeros       = length(rows);
            for n          = 1:nonzeros
                i            = rows(n);
                j            = cols(n);
                I            = permutation(i);
                J            = permutation(j);
                perm_rows(n) = I;
                perm_cols(n) = J;
                perm_vals(n) = vals(n);
            end
            
            B = sparse(perm_rows, perm_cols, perm_vals);
            
        end
        
        
        
        function B = reorderBack(this, permutation, A)
            [rows,cols,vals] = find(A);
            
            n = size(A, 1);
            fprintf('size of A is %d by %d\n',n,n);
            for i = 1:n
                inverse_permutation(permutation(i)) = i;
            end
            
            
            
            nonzeros       = length(rows);
            for n          = 1:nonzeros
                i            = rows(n);
                j            = cols(n);
                I            = inverse_permutation(i);
                J            = inverse_permutation(j);
                perm_rows(n) = I;
                perm_cols(n) = J;
                perm_vals(n) = vals(n);
            end
            
            B = sparse(perm_rows, perm_cols, perm_vals);
            
        end
        
        
        
        
        function out = SolveKlaus(this)
            
            figure; spy(this.KKT, 'k');
            p = computePermutation(this, 'Klaus');
            fprintf('finished computing the permutation\n');
            %b       = this.b;
            %b(p)    = this.b;
            H  = this.reorder(p, this.KKT);
            
            figure; spy(H, 'b');
            %if (this.makeSpyPlots == true)
            %  figure; spy(this.H, 'b');
            %end
            
            nb     = this.nb; % number of buses
            ng     = this.ng; % number of generators
            nl     = this.nl; % number of lines
            ns     = this.ns; % number of storages
            N      = this.N;  % number of timesteps
            
            
            NB     = nb * N;
            NG     = ng * N;
            NL     = nl * N;
            
            NH     = 2*NB + 2*NG;
            NJ     = 2*NB;
            NP     = 2*NL + (N+1)*ns;
            NK     = 2*NL + (N+1)*ns;
            
            
            D      = H(1:2*NG+2*NL, 1:2*NG+2*NL);
            figure; spy(D,'c');
            S      = H(1:2*NG+2*NL, 2*NG+2*NL+1:end);
            figure; spy(S,'m');
            D_1S   = D\S;
            
            H(2*NG+2*NL+1:end, 2*NG+2*NL+1:end) = H(2*NG+2*NL+1:end, 2*NG+2*NL+1:end) -S'*D_1S;
            figure; spy(D_1S,'r');
            figure; spy(H,'g');
            H1 = H(2*NG+2*NL+1:end, 2*NG+2*NL+1:end);
            figure; spy(H1,'b');
            
            p = computePermutation(this, 'Klaus1');
            H2  = this.reorder(p, H1);
            figure; spy(H2, 'k');
            
            D = H2(1:2*NL+2*NB-(N+1)*ns, 1:2*NL+2*NB-(N+1)*ns);
            figure; spy(D,'c');
            B = H2(1:2*NL+2*NB-(N+1)*ns, 2*NB+2*NL+1:2*NB+2*NL+2*NB);
            figure; spy(B,'c');
            d=spdiags(abs(D));
            I = find(d < 1e-10);
            fprintf('the number of zeros is %d\n',length(I));
            save -ascii 'd.txt' I;
            fprintf('minmax %f %f\n',min(d),max(d));
            S = -B'*(D\B);
            figure; spy(S,'r');
            A = H2(2*NB+2*NL+1:end, 2*NB+2*NL+1:end);
            figure; spy(A,'b');
            A = A + S;
            figure; spy(A,'k');
            
        end
        
        
    end  % methods
    
end     % PowerGridHessian
