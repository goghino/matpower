function mpcN      = create_multi_time_step_case_file(mpc0)

N                  = length(mpc0.load_profile(:,1));

nnodes             = size(mpc0.bus,1);

% Nb
nbranches          = size(mpc0.branch,1);

% Ng
ngens              = size(mpc0.gen,1);

if not(min(mpc0.bus(:,1) == (1:nnodes)'))
    error('buses not properly numbered')
end

mpcN.version       = mpc0.version ;
mpcN.baseMVA       = mpc0.baseMVA;

mpcN.bus           = repmat(mpc0.bus,[N,1]);
mpcN.bus(:,1)      = 1:(nnodes*N); %fix bus IDs

%load can be split into different categories - e.g. residential and
%commertial, each is scaled using different load scaling profile
base_load = mpcN.bus(:,3:4);
mpcN.bus(:,3:4) = zeros(size(mpcN.bus(:,3:4)));
for load_cat=1:size(mpc0.load_profile, 2) %number of cols
  load_scaling_profile = mpc0.load_profile(:, load_cat); %needs to be column vec
  load_ratios = mpc0.load_ratios(:, load_cat);%needs to be column vec
  %sum of different load categories (ratio of base load) scaled by its load scaling profile
  mpcN.bus(:,3:4) = mpcN.bus(:,3:4) + base_load.*repmat(kron(ones(N,1),load_ratios),[1 2]).*repmat(kron(load_scaling_profile,ones(nnodes,1)),[1 2]) ; %Scale PD, QD for each N by load_profile
end

figure;
subplot(2,1,1)
for b=1:nnodes
    plot(mpcN.bus(b:15:end,3)*mpc0.baseMVA, 'LineWidth',2); hold on;
end
title('Pd at the buses');
xlabel('T [hrs]');
ylabel('Load demand [MW]');
subplot(2,1,2)
for b=1:nnodes
    plot(mpcN.bus(b:15:end,3)*mpc0.baseMVA, 'LineWidth',2); hold on;
end
title('Qd at the buses');
xlabel('T [hrs]');
ylabel('Load demand [MW]');

mpcN.branch        = repmat(mpc0.branch,[N,1]);
mpcN.branch(:,1:2) = mpcN.branch(:,1:2) + kron( (0:(N-1))' , ones(nbranches,2)*nnodes ); %fix F_BUS, T_BUS


mpcN.gen           = repmat(mpc0.gen,[N,1]); %fix GEN_BUS
mpcN.gen(:,1)      = mpcN.gen(:,1) + kron( (0:(N-1))' , ones(ngens,1)*nnodes );


mpcN.gencost = kron(ones(N,1),mpc0.gencost); %% same marginal cost for all time steps

mpcN.load_scaling_profile              = load_scaling_profile;

mpcN.id_gen_original_generators        = (1:(N*ngens))';

mpcN.horizon       = N;


end
