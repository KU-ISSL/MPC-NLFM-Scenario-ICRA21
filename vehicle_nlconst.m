function [c,ceq] = vehicle_nlconst(u,w,scenario_struct)
%VEHICLE_NLCONST State constraints via nonlinear constraints.s
%   Inputs
%       u : control input over horizon
%       w : disturbance

xk = scenario_struct.xk;
H = scenario_struct.H;
K = scenario_struct.N;
sys = scenario_struct.sys;
sl = scenario_struct.constraints.stateL;
su = scenario_struct.constraints.stateU;

nx = size(xk,1);
% u : (nu,H);

% STATE PREDICTION OVER HORIZON
x = zeros(nx,H+1);
x(:,1) = xk;
% w : (nw,H,K);   [nw=nx]
xcon = zeros(2*nx,H+1,K); %constrained state vector
bounds = zeros(2*nx,H+1,K); %constraints


for k = 1:K %scenarios
    for i = 2:H+1 %prediction horizon
        x(:,i) = sys(x(:,i-1),u(:,i-1),w(:,i-1,k));
        xcon( 1:nx , i , k ) = -x(:,i);
        xcon( 1+nx:2*nx , i , k ) = x(:,i);
    end
    bounds(:,2:end,k) = [-sl;su];
end

% COMPUTE INEQUALITY CONSTRAINT
%c <= 0
c = xcon - bounds;
c = c(:);
ceq = [];

end

