function [results, success, raw] = mpopf_execute(om, mpopf_aux, mpopt)
%OPF_EXECUTE  Executes the OPF specified by an OPF model object.
%   [RESULTS, SUCCESS, RAW] = OPF_EXECUTE(OM, MPOPT)
%
%   RESULTS are returned with internal indexing, all equipment
%   in-service, etc.
%
%   See also OPF, OPF_SETUP.

%   MATPOWER
%   Copyright (c) 2009-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

%%-----  setup  -----
%% options
dc  = strcmp(upper(mpopt.model), 'DC');
alg = upper(mpopt.opf.ac.solver);
sdp = strcmp(alg, 'SDPOPF');

%% get indexing
[vv, ll, nne, nni] = om.get_idx();

if mpopt.verbose > 0
    v = mpver('all');
    fprintf('\nMATPOWER Version %s, %s', v.Version, v.Date);
end

%%-----  run DC OPF solver  -----
if dc
  if mpopt.verbose > 0
    fprintf(' -- DC Optimal Power Flow\n');
  end
  error('DC not supported in MPOPF');
  [results, success, raw] = dcopf_solver(om, mpopt);
else
  %%-----  run AC OPF solver  -----
  if mpopt.verbose > 0
    fprintf(' -- AC Optimal Power Flow\n');
  end

  %% ZIP loads?
  if (~isempty(mpopt.exp.sys_wide_zip_loads.pw) && ...
          any(mpopt.exp.sys_wide_zip_loads.pw(2:3))) || ...
          (~isempty(mpopt.exp.sys_wide_zip_loads.qw) && ...
          any(mpopt.exp.sys_wide_zip_loads.qw(2:3)))
    switch alg
    case {'PDIPM', 'TRALM', 'MINOPF', 'SDPOPF'}
      warning('opf_execute: ''%s'' solver does not support ZIP load model. Converting to constant power loads.', alg)
      mpopt = mpoption(mpopt, 'exp.sys_wide_zip_loads', ...
                        struct('pw', [], 'qw', []));
    end
  end

  %% run specific AC OPF solver
  switch alg
    case 'MIPS'
      error('MIPS not supported in MPOPF');
      [results, success, raw] = mipsopf_solver(om, mpopt);
    case 'IPOPT'
      if ~have_fcn('ipopt')
        error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires IPOPT (see http://www.coin-or.org/projects/Ipopt.xml)', alg);
      end
      [results, success, raw] = ipoptmpopf_solver(om, mpopf_aux, mpopt);
%     case 'PDIPM'
%       if mpopt.pdipm.step_control
%         if ~have_fcn('scpdipmopf')
%           error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires SCPDIPMOPF (see http://www.pserc.cornell.edu/tspopf/)', alg);
%         end
%       else
%         if ~have_fcn('pdipmopf')
%           error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires PDIPMOPF (see http://www.pserc.cornell.edu/tspopf/)', alg);
%         end
%       end
%       [results, success, raw] = tspopf_solver(om, mpopt);
%     case 'TRALM'
%       if ~have_fcn('tralmopf')
%         error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires TRALM (see http://www.pserc.cornell.edu/tspopf/)', alg);
%       end
%       [results, success, raw] = tspopf_solver(om, mpopt);
%     case 'MINOPF'
%       if ~have_fcn('minopf')
%         error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires MINOPF (see http://www.pserc.cornell.edu/minopf/)', alg);
%       end
%       [results, success, raw] = mopf_solver(om, mpopt);
%     case 'FMINCON'
%       if ~have_fcn('fmincon')
%         error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires FMINCON (Optimization Toolbox 2.x or later)', alg);
%       end
%       [results, success, raw] = fmincopf_solver(om, mpopt);
%     case 'KNITRO'
%       if ~have_fcn('knitro')
%         error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires KNITRO (see http://www.ziena.com/)', alg);
%       end
%       [results, success, raw] = ktropf_solver(om, mpopt);
%     case 'SDPOPF'
%       if ~have_fcn('yalmip')
%         error('opf_execute: MPOPT.opf.ac.solver = ''%s'' requires YALMIP (see http://users.isy.liu.se/johanl/yalmip/)', alg);
%       end
%       [results, success, raw] = sdpopf_solver(om, mpopt);
    otherwise
      error('opf_execute: MPOPT.opf.ac.solver = ''%s'' is not a valid AC OPF solver selection', alg);
  end
end
if ~isfield(raw, 'output') || ~isfield(raw.output, 'alg') || isempty(raw.output.alg)
    raw.output.alg = alg;
end

end