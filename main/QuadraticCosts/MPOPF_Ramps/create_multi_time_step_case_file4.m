%Function that replicates mpc case for each time step
%and scales P,Q demand by profile provided as input parameter for every bus
function mpcN      = create_multi_time_step_case_file(mpc0,load_scaling_profile)

%number of ???
N                  = length(load_scaling_profile);

%number of buses
nnodes             = size(mpc0.bus,1);

% Nb
nbranches          = size(mpc0.branch,1);

% Ng
ngens              = size(mpc0.gen,1);

if not(min(mpc0.bus(:,1) == (1:nnodes)'))
    error('buses not properly numbered')
end

%copy some variables into new mpc
mpcN.version       = mpc0.version ;
mpcN.baseMVA       = mpc0.baseMVA;

%replicate buses, real and reactive power demand scaled by given profile
mpcN.bus           = repmat(mpc0.bus,[N,1]); %repeat bus structure as many times as there are time steps
mpcN.bus(:,1)      = 1:(nnodes*N); %number buses, idx 1 ==  bus number (positive integer)
mpcN.bus(:,3:4)    = mpcN.bus(:,3:4).*repmat(kron(load_scaling_profile,ones(nnodes,1)),[1 2]) ;

%replicate branches, add offsets for each timestep to from-to columns
%in order to connect buses only within timestep case (buses numbering: first nnodes for t=1, then t=2 and so on...)
mpcN.branch        = repmat(mpc0.branch,[N,1]);
mpcN.branch(:,1:2) = mpcN.branch(:,1:2) + kron( (0:(N-1))' , ones(nbranches,2)*nnodes ); %from-to bus numbers

%replicate generator, add offset nnodes for generators for each timestep so
%that it belongs to proper bus
mpcN.gen           = repmat(mpc0.gen,[N,1]);
mpcN.gen(:,1)      = mpcN.gen(:,1) + kron( (0:(N-1))' , ones(ngens,1)*nnodes ); %GEN_BUS bus number


% mpcN.gencost       = [ones(ngens*N,1)*[2 0 0 3 0], kron(price_trajectory_CHFperMWTimeUnit,ones(ngens,1)), zeros(ngens*N,1)];
mpcN.gencost = kron(ones(N,1),mpc0.gencost); %% same marginal cost for all time steps

mpcN.load_scaling_profile              = load_scaling_profile;
% mpcN.price_trajectory_CHFperMWTimeUnit = price_trajectory_CHFperMWTimeUnit;
mpcN.id_gen_original_generators        = (1:(N*ngens))';

mpcN.horizon       = N;


end
