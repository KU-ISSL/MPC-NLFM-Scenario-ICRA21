function f = vehicle_model(x,u,w,Aw,qw)
%VEHICLE_MODEL State transition function for a vehcile based on the bicycle
%model.
%   States:
%       x(1) = r : position
%       x(2) = v : velocity
%       x(3:end) : latent input GP SS
%
%   Inputs:
%       u: applied control
%       w: white noise disturbance
%% System Model
% Original system with unknown latent input

% [dx] = [v*cos(h+beta)]*dt + [qx         ]*[dB]
% |dy|   |v*sin(h+beta)|      |   qy      |
% |dv|   |    a + w    |      |      qv   |
% [dh]   [v/l*sin(beta)]      [         qh]
% beta = atan(0.5*tan(s))
% u = [a; s]
% w = latent input

% system parameters 
Nx = 4;     %state dimension 
l = 0.25;   %half length of car
qx = 0.01;
qy = 0.01;
qv = 0.01; 
qh = deg2rad(0.5);

beta = @(s) atan(0.5*tan(s));   %slip angle eq.
sysf = @(x,u,w) [x(3)*cos(x(4)+beta(u(2)));
                 x(3)*sin(x(4)+beta(u(2)));
                 u(1) + w;
                 x(3)/l*sin(beta(u(2)))];

%% Augmented System
% Augmented system with original states and GP prior states. Augmented
% model

% dg = fa(g)*dt + La*dB
% y = h(g) + r    r~N(0,R)
% h = H*g

dt = 0.2;                        %time step
steps = 1;                       %Euler integration steps (50)
ddt = dt/steps;
Q = blkdiag(diag([qx,qy,qv,qh]), qw); %diffusion

for i=1:steps
    %augmented system model
    x = x + [sysf(x(1:Nx),u,x(Nx+1));     
             Aw*x(Nx+1:end)].*ddt + ...
             sqrt(ddt)*Q*w; 
end
f = x;

end


