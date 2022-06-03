function [uopt,fcost] = Scenario_nonlinear(ukm1,w,ref,scenario_opt)
%SCENARIO_NONLINEAR  Computes optimal control via the Scenario Method for nonlinear systems. 
%
%   [u,f] = SCENARIO_NONLINEAR(ukm1,w,ref,scenario_opt)
%
%   Inputs:
%       ukm1          : previous control solution
%       w             : deterministic disturbance scenarios
%       ref           : reference trajectory/ goal point
%       scenario_opt  : structure containing info for the cost & constraints
%
%   Outputs:
%       uopt   : optimal control sequence
%       fcost  : value function at uopt


% UNPACK
% global x0 H K costFcn nx sys
% x0 = xk;
% nx = size(w,1); %state dim
% H = size(w,2); %horizon
% K = size(w,3); %scenarios
% costFcn = scenario_opt.costFcn;
% sys     = scenario_opt.sys;

%cost
cost = @(u) scenario_opt.costFcn.cost(u',w,ref,scenario_opt);

%constraints
A     = scenario_opt.constraints.A;
b     = scenario_opt.constraints.b;
Aeq   = scenario_opt.constraints.Aeq;
beq   = scenario_opt.constraints.beq;
lb    = scenario_opt.constraints.lb.*ones(size(ukm1));
ub    = scenario_opt.constraints.ub.*ones(size(ukm1));
nlcon = @(u) scenario_opt.constraints.nlcon(u',w,scenario_opt);

% OPTIMIZATION
%use last control solution as warm-start 

opt = optimoptions(@fmincon,'Display','off');
[us,fcost,eflag] = fmincon(cost,ukm1,A,b,Aeq,beq,lb,ub,nlcon);

% RETURN ENTIRE CONTROL SEQUENCE
uopt = us;

end

%% Optimization Function
% function [f] = cost(u)
% global x0 H K costFcn nx sys
% x = zeros(nx,H+1);
% x(:,1) = x0;
% w = zeros(nx,H,K);
% costsum = 0;
% 
% for k = 1:K %scenarios
%     for i = 1:H %prediction horizon
%         costsum = costsum + costFcn.cost(x(i),u(i));
%         x(i+1) = sys(x(i),u(i),w(:,i,k));
%     end
% end
% 
% end
