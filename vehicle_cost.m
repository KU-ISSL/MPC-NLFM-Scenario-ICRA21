function cost = vehicle_cost(u,w,ref,scenario_opt)
%VEHICLE_COST Cost function of scenarios used for optimization.
%   Inputs:
%       u - input sequence over horizon (nu,H)
%       w - disturbance
%       ref - reference signal
%       scenario_opt - structure containing cost terms

xk = scenario_opt.xk;
nx = size(xk,1);
nu = size(u,1);

Q = scenario_opt.costFcn.Q;
if all(size(Q)~=[nx,nx])
    error('''Q'' is not the correct size: [%d, %d]',size(Q,1),size(Q,2));
end

Qf = scenario_opt.costFcn.Qf;
if all(size(Qf)~=[nx,nx])
    error('''Qf'' is not the correct size: [%d, %d]',size(Qf,1),size(Qf,2));
end

R = scenario_opt.costFcn.R;
if all(size(R)~=[nu,nu])
    error('''R'' is not the correct size: [%d, %d]',size(R,1),size(R,2));
end

H   = scenario_opt.H;
N   = scenario_opt.N;
sys = scenario_opt.sys;

% ref can be specified in three ways:
% if ref is not specified, set to zeros
% if ref is specified as a ROW vector, fill matrix for entire horizon
% else, ref should be specified as a (nx,H+1) matrix
if isempty(ref)
    ref = zeros(nx,H+1);
elseif size(ref,2) == nx
    ref = repmat(ref',1,H+1);
elseif all(size(ref)~=[nx, H+1])
    error('''ref'' is not the correct size: [%d, %d]',size(ref,1),size(ref,2));
end

x = zeros(nx,H+2);
% u : (nu,H)
% w : (nw,H,K)  [nw=nx]
x(:,1) = xk;
cost = 0;

for k = 1:N %scenarios
    for i = 1:H+1 %prediction horizon
        if i==H+1
            cost = cost + (x(:,i) - ref(:,i))'*Qf*(x(:,i) - ref(:,i));
        else
            cost = cost + (x(:,i) - ref(:,i))'*Q*(x(:,i) - ref(:,i)) + u(:,i)'*R*u(:,i);
            x(:,i+1) = sys(x(:,i),u(:,i),w(:,i,k));
        end
        
    end
end

end

