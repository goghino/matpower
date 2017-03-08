function xnew = check_ramps(x, mpc)
%CHECK_RAMPS Imposes generator ramping constraints to MPOPF problem.
%
% Feasibility of power generation ramps is checked between time steps
% and if ramp limit is violated the generator powers are adjusted accordingly
% and the modified quantities are used for evaluation of IPOPT callback functions.
% 
% x contains [Va Vm P Q] where P is [Pgen Pc Pd] and Pgen is the real
% power of generator, Pc is charging of storages and Pd is discharging.
%
% mpc holds info about structure of the MPOPF problem

%TODO: remove
x_bk = x;

%get MPOPF problem parameters
N = mpc.horizon;
ns = mpc.nstorage;
NP = size(mpc.gen,1); %size of the P variable of OPF (includes gens and 2*storages)
NG = (NP - N*ns*2);
ng = NG/N;
i2e = mpc.order.gen.i2e;
e2i = mpc.order.gen.e2i;

%convert x to external ordering (gens first, then storages)
offset = size(mpc.bus,1)*2; %skip [Va Vm] variables
xe = x(offset+1:offset+NP);
xe = xe(i2e);

%build constraint matrix for one timestep and global time horizon
A = [-eye(ng), eye(ng)];

%TODO: not necessary, removve, only debug purpose
Ag = [kron(eye(N-1), -eye(ng)), zeros((N-1)*ng,ng)] + [zeros((N-1)*ng,ng), kron(eye(N-1), eye(ng))];
ramps_old = Ag * xe(1:NG);
ro = reshape(ramps_old,ng,N-1);
%----------------------------

%% handle ramp max violations by setting power differences to max acceptable rate
% iterate ofer time and check violations of generators, if there is a
% violation fix it by limiting ramp and move to next iteration in time 
% and fix possible cascade effects of the previous adjustment
for t = 0:N-2
    %compute generator ramps of real power generation
    ramps = A * xe(t*ng+1:(t+2)*ng); %ramps at t, t+1
    %figure;plot(ramps);axis([1,ng,-2,2]);

    %find ramp limits violations
    r_max = mpc.ramp_max;
    ii = find(ramps > r_max);

    %ramp = Pi_tj - Pi_t(j-1) [i-th gen, j-th time step]
    %new value is: Pi_tj = Pi_t(j-1) + r_max
    xe(t*ng + ii + ng) = xe(t*ng + ii) + r_max;
end
% ramps = Ag * xe(1:NG);
% r = reshape(ramps,ng,N-1);
% figure;plot(r');axis([0,N-1,-2,2]);

% %% handle min violations by setting power differences to min acceptable rate
for t = 0:N-2
    %compute generator ramps of real power generation
    ramps = A * xe(t*ng+1:(t+2)*ng); %ramps at t, t+1
    %figure;plot(ramps);axis([1,ng,-2,2]);

    %find ramp limits violations
    r_min = mpc.ramp_min;
    ii = find(ramps < r_min);

    %ramp = Pi_tj - Pi_t(j-1) [i-th gen, j-th time step]
    %new value is: Pi_tj = Pi_t(j-1) + r_min [r_min is already negative]
    xe(t*ng + ii + ng) = xe(t*ng + ii) + r_min;
end

%% ----------------------
%TODO: not necessary, is here just to plot the final result
ramps_new = Ag * xe(1:NG);
rn = reshape(ramps_new,ng,N-1);
%figure;plot(ro','b');hold on; plot(rn','r');
%----------------------------

%convert x to internal ordering
x(offset+1:offset+NP) = xe(e2i);

%return new generator powers
xnew = x;
%xnew = x_bk; %TODO: remove