function [f, df, d2f] = mpopf_costfcn(x, om, time_period)
%OPF_COSTFCN  Evaluates objective function, gradient and Hessian for OPF.
%   [F, DF, D2F] = OPF_COSTFCN(X, OM)
%
%   Objective function evaluation routine for AC optimal power flow,
%   suitable for use with MIPS or FMINCON. Computes objective function value,
%   gradient and Hessian.
%
%   Inputs:
%     X : optimization vector
%     OM : OPF model object
%
%   Outputs:
%     F   : value of objective function
%     DF  : (optional) gradient of objective function (column vector)
%     D2F : (optional) Hessian of objective function (sparse matrix)
%
%   Examples:
%       f = opf_costfcn(x, om);
%       [f, df] = opf_costfcn(x, om);
%       [f, df, d2f] = opf_costfcn(x, om);
%
%   See also OPF_CONSFCN, OPF_HESSFCN.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- evaluate objective function -----
%% general nonlinear costs
%if we want to evaluate only one time period, we need to
%specify the names of the cost function added in mpopf_setup
%For whole time period horizon, the name is not specified and 
%all cost functions are evaluated
if (time_period > 0)
   %use Pg, Qg cost only if present
    names = cell(0);
    
   name = strcat('polPg',int2str(time_period)); 
   if isfield(om.qdc.idx.N, name)
       names{1} = name;
   end
   
   name = strcat('polQg',int2str(time_period));
   if isfield(om.qdc.idx.N, name)
        names{2} = name;
   end
        
    %get problem dimensions needed to build proper offsets
    mpc = om.get_mpc();
    ng = size(mpc.gen,1);
    nb = size(mpc.bus,1);
    N = 2*nb + 2*ng;
    
    if(om.nlc.NS > 0)
       error('MPOPF does not support nonlinear cost function'); 
    end
end

if nargout == 3
    [f, df, d2f]    = om.eval_nln_cost(x);
    if om.qdc.NS
        if (time_period < 0) %evaluate cost for whole time horizont
            [fq, dfq, d2fq] = om.eval_quad_cost(x);
            f = f + sum(fq);
            df = df + dfq;
            d2f = d2f + d2fq;
        else %evaluate only single time step (named set polPgi + polQgi)
            offset = 2*nb; %offset for df/dPg, skip df/dVa and df/dVm
            f = 0;
            df = zeros(N,1);
            d2f = sparse(N,N);
            for name = names
                [fq, dfq, d2fq] = om.eval_quad_cost(x, name{1});%name is cell, reference it to get str
                f = f + sum(fq);
                df(offset + (1:ng),1) = dfq;
                d2f(offset + (1:ng),offset + (1:ng)) = spdiags(d2fq,0,ng,ng); 
                offset = offset + ng; %offset for df/dQg
            end
        end
    end
elseif nargout == 2
    [f, df]   = om.eval_nln_cost(x);
    if om.qdc.NS
        if (time_period < 0) %evaluate cost for whole time horizont
            [fq, dfq] = om.eval_quad_cost(x);
            f = f + sum(fq);
            df = df + dfq;
        else %evaluate only single time step (named set polPgi + polQgi)
            offset = 2*nb; %offset for df/dPg, skip df/dVa and df/dVm
            f = 0;
            df = zeros(N,1);
            for name = names
                [fq, dfq] = om.eval_quad_cost(x, name{1});%name is cell, reference it to get str
                f = f + sum(fq);
                df(offset + (1:ng),1) = dfq;
                offset = offset + ng; %offset for df/dQg
            end
        end
    end
else
    f  = om.eval_nln_cost(x);
    if om.qdc.NS
        if (time_period < 0) %evaluate cost for whole time horizont
            [fq] = om.eval_quad_cost(x);
            f = f + sum(fq);
        else %evaluate only single time step (named set polPgi + polQgi)
            offset = 2*nb; %offset for df/dPg, skip df/dVa and df/dVm
            f = 0;
            for name = names
                fq = om.eval_quad_cost(x, name{1});%name is cell, reference it to get str
                f = f + sum(fq);
                offset = offset + ng; %offset for df/dQg
            end
        end
    end
end
