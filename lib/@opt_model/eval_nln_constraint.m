function [g, dg] = eval_nln_constraint(om, x, iseq, time_period)
%EVAL_NLN_CONSTRAINT  Builds and returns full set of nonlinear constraints.
%   [G, DG] = OM.EVAL_NLN_CONSTRAINT(X, ISEQ)
%   Builds a full set of nonlinear equality or inequality constraints and
%   their gradients for a given value of the optimization vector X,
%   based on constraints added by ADD_NLN_CONSTRAINT.
%
%       g(X) = 0
%       h(X) <= 0
%
%   Example:
%       [g, dg] = om.eval_nln_constraint(x, 1)
%       [h, dh] = om.eval_nln_constraint(x, 0)
%
%   See also OPT_MODEL, ADD_NLN_CONSTRAINT, EVAL_NLN_CONSTRAINT_HESS.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% get constraint type
TIME_NOT_SPECIFIED = 0;
if (nargin < 4)
    time_period = TIME_NOT_SPECIFIED; 
elseif (time_period == 0)
    error('Parameter time_period has to be either negative or positive, cannot be zero!');
end

if iseq         %% equality constraints
    om_nlx = om.nle;
else            %% inequality constraints
    om_nlx = om.nli;
end

%% initialize g, dg
%properly set size of the Hessian, depending if we want only single
%period or whole horizont, based on the paramter time_period
if (time_period <= 0)
    %whole time horizon
    nrows = om_nlx.N;
    ncols = om.var.N;
elseif (time_period > 0)
    %only single period
    mpc = om.get_mpc();
    if (iseq)
        nrows = 2*size(mpc.bus,1); %EQ constr
    else
        nrows = 2*size(mpc.branch,1); %INEQ constr
    end
    ncols = 2*size(mpc.bus,1) + 2*size(mpc.gen,1);    
end
g = NaN(nrows, 1);
dg = sparse(0, ncols);   %% build gradient by stacking

%% calls to substruct() are relatively expensive, so we pre-build the
%% structs for addressing cell and numeric array fields, updating only
%% the subscripts before use
sc = struct('type', {'.', '{}'}, 'subs', {'', 1});  %% cell array field
sn = sc; sn(2).type = '()';                         %% num array field

%% fill in each piece
for k = 1:om_nlx.NS
    name = om_nlx.order(k).name;
    idx  = om_nlx.order(k).idx;
    if isempty(idx)
        if ~isfield(om_nlx.data.fcn, name)
            continue;   %% skip, there is no function handle stored here,
                        %% the function value for this named set was included
                        %% in the value computed by a previous named set
        end
        N = om_nlx.idx.N.(name);    %% number of constraint functions
                                    %% evaluated for this named set
        if isfield(om_nlx.data.include, name)
            N = N + sum(om_nlx.data.include.(name).N);
        end
    else
        % (calls to substruct() are relatively expensive ...
        % sn = substruct('.', name, '()', idx);
        % sc = substruct('.', name, '{}', idx);
        % ... so replace them with these more efficient lines)
        sn(1).subs = name;
        sn(2).subs = idx;
        sc(1).subs = name;
        sc(2).subs = idx;
        N = subsref(om_nlx.idx.N, sn);
    end
    if N                                %% non-zero number of rows
        if isempty(idx)
            fcn = om_nlx.data.fcn.(name);       %% fcn for kth constraint set
            i1 = om_nlx.idx.i1.(name);          %% starting row index
            iN = i1 + N - 1;                    %% ending row index
            vs = om_nlx.data.vs.(name);         %% var sets
        else
            fcn = subsref(om_nlx.data.fcn, sc); %% fcn for kth constraint set
            i1 = subsref(om_nlx.idx.i1, sn);    %% starting row index
            iN = subsref(om_nlx.idx.iN, sn);    %% ending row index
            vs = subsref(om_nlx.data.vs, sc);   %% var sets
        end
        xx = om.varsets_x(x, vs);
        if(time_period == TIME_NOT_SPECIFIED)
            [gk, dgk] = fcn(xx);    %% evaluate kth constraint and gradient
        else
            [gk, dgk] = fcn(xx, time_period);    %% evaluate kth constraint and gradient
        end
        
        %in general, we can add multiple equality constraints in the
        %mpopf_setup, this is handled by the indexing i1..iN, which 
        %gives us global indexes for the current set of equality
        %constraints
        if (time_period <= 0)
            g(i1:iN) = gk;          %% assign kth constraint
        elseif (time_period > 0)
            assert(om_nlx.NS == 2); %[Pmis, Qmis]
            %however, in the case of MPOPF evaluaing single period,
            %we assume only a single set of eq. constraints exist
            %(global for whole time horizon) and we process it without
            %using the i1..iN indexes
            g(1:nrows) = gk;
        end
        
        if isempty(vs)          %% all rows of x
            if(time_period <= 0)
                %full multiperiod horizon
                if size(dgk, 2) == om.var.N
                    dg = [dg; dgk];
                else                %% must have added vars since adding
                                    %% this constraint set
                    dg(i1:iN, 1:size(dgk, 2)) = dgk;
                end
            else
                assert(om_nlx.NS == 2); %[Pmis, Qmis]
                dg = [dg; dgk];
            end
        else                    %% selected rows of x
            jj = om.varsets_idx(vs);    %% column indices for var set
            dgi = sparse(N, om.var.N);
            dgi(:, jj) = dgk;
            dg = [dg; dgi];
        end
    end
end
