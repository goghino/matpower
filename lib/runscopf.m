function [results, success, info] = ...
                runscopf(casedata, cont, mpopt, tol)
%RUNOPF  Runs a security constrained optimal power flow.
%   [RESULTS, SUCCESS] = RUNOPF(CASEDATA, MPOPT, CONT, FNAME, SOLVEDCASE)
%
%   Runs a security constrained optimal power flow (AC OPF),
%   optionally returning a RESULTS struct and SUCCESS flag.
%
%   Inputs (all are optional):
%       CASEDATA : either a MATPOWER case struct or a string containing
%           the name of the file with the case data (default is 'case9')
%           (see also CASEFORMAT and LOADCASE)
%       CONT  : list of branch contingencies, empty for standard OPF
%               OR number specifying number of contingencies to be selected
%                  automatically iff CONT < 0 and length(CONT)==1
%       MPOPT : MATPOWER options struct to override default options
%           can be used to specify the solution algorithm, output options
%           termination tolerances, and more (see also MPOPTION).
%       FNAME : name of a file to which the pretty-printed output will
%           be appended
%       SOLVEDCASE : name of file to which the solved case will be saved
%           in MATPOWER case format (M-file will be assumed unless the
%           specified name ends with '.mat')
%
%   Outputs (all are optional):
%       RESULTS : results struct, with the following fields:
%           (all fields from the input MATPOWER case, i.e. bus, branch,
%               gen, etc., but with solved voltages, power flows, etc.)
%           order - info used in external <-> internal data conversion
%           et - elapsed time in seconds
%           success - success flag, 1 = succeeded, 0 = failed
%           (additional OPF fields, see OPF for details)
%       SUCCESS : the success flag can additionally be returned as
%           a second output argument
%
%   Calling syntax options:
%       [RESULTS, SUCCESS] = runscopf(casedata, cont, mpopt, tol);
%
%   Example:
%       results = runscopf('case30',[1 2 3], mpopt, 1e-4);
%
%   See also RUNDCOPF, RUNOPF.

%%-----  initialize  -----
%% default arguments
if nargin < 5
    solvedcase = '';                %% don't save solved case
    if nargin < 4
        fname = '';                 %% don't print results to a file
        if nargin < 3
            mpopt = mpoption;       %% use default options
            if nargin < 1
                casedata = 'case9'; %% default data file is 'case9.m'
                cont = [];          %% by default no contingency
            end
        end
    end
end

%%-----  add nominal case to the list of contingencies  -----
if ~isempty(cont) && (size(cont, 1) > 1 && size(cont, 2) > 1)
    error('List of contingencies must be a vector.');
end

if size(cont, 2) > 1
   cont = cont'; %scopf expects a column vector
end

% we were given a vector, add nominal case in this case;
% else assume negative N, specifying no. of contingencies to be chosen automatically
if ~(length(cont) == 1 && cont < 0)
    if (sum(cont > 0) ~= length(cont))
       error('List of contingencies must contain only positive integers\n'); 
    end
end 

%%-----  run the optimal power flow  -----
[results, success] = scopf(casedata, cont, mpopt, tol);
info = results.raw.info;


end
