% This is an example of using MPC for path planning on a nonlinear LFM
% autonomous vehicle

clear; close all
%% Setup
% [dx] = [v*cos(h+beta)]*dt + [qx         ]*[dB]
% |dy|   |v*sin(h+beta)|      |   qy      |
% |dv|   |    a + w    |      |      qv   |
% [dh]   [v/l*sin(beta)]      [         qh]
% beta = atan(0.5*tan(s))
% u = [a; s]

dt = 0.2;       %sampling time
sig2 = 2^2;    %latent model parameters
l = 4;

nw = 3;     %GP state
nx = nw+4;  %augmented state
ny = 4;     %output state
nu = 2;     %input state

scenario_opt = struct('costFcn',struct,...
                      'constraints',struct,...
                      'xk',[],...
                      'H',[],...
                      'N',[],...
                      'sys',[]);

%% Modeling
% Latent GP Prior Model
% Construct linear SDE for the Guassian Prior. 

% Hyperparamters
% theta = [sig2,l];
theta = [sig2, l];  

% Spectral factorization
[h,qw,nw] = specfactor('matern',theta,2);

% Create state space
% dz/dt = Aw*z + Bw*randn(nw,1)
Aw = zeros(nw);
Aw(1:end-1,2:end) = eye(nw-1);
Aw(end,:) = -h;

Bw = zeros(nw,1);
Bw(end) = 1;

qw = Bw*sqrt(qw)*Bw';

% System Functions
% system parameters inside model function definition
scenario_opt.sys = @(x,u,w) vehicle_model(x,u,w,Aw,qw);

R = diag([0.15^2 0.15^2 0.15^2 deg2rad(0.5)^2]); %Noise covariance
H = [eye(nx-nw), zeros(4,nw)]; %Output function
outputFcn = @(x) H*x + mvnrnd(zeros(ny,1),R)';

x0 = [0; 0; 1; deg2rad(45); 0; 0; 0];
scenario_opt.xk = x0;

%% Optimization
scenario_opt.H = 10; %control horizon

%Compute min required scenarios
eps = 0.05;                 %violation probability
beta = 10e-10;              %reliability level
dvars = nu*scenario_opt.H;  %decision variables
scenario_opt.N = 50; %num_scenarios(eps,beta,dvars); %total number of scenarios

% Cost Function 
scenario_opt.costFcn = struct('Q',eye(nx),...
                              'R',eye(nu),...
                              'Qf',eye(nx),...
                              'cost',[]);
scenario_opt.costFcn.Q = diag([0.5 2 12 0 0 0 0]);
scenario_opt.costFcn.Qf = diag([0.5 2 12 0 0.1 0.1 0.1]);
scenario_opt.costFcn.R  = diag([0.5, 0.75]);
scenario_opt.costFcn.cost = @(u,w,ref,scen_struct) vehicle_cost(u, w, ref, scen_struct);

%% Constraints
scenario_opt.constraints = struct('A',[],...
                                  'b',[],...
                                  'Aeq',[],...
                                  'beq',[],...
                                  'ub',[],...
                                  'lb',[],...
                                  'nlcon',[]);
% State constraints
stateLower = [-1 -0.75 0 -pi -1000 -1000 -1000];
stateUpper = [70 0.75 8 pi 1000 1000 1000];
scenario_opt.constraints.stateL = repmat(stateLower',1,scenario_opt.H); %same lower bound for entire horizon
scenario_opt.constraints.stateU = repmat(stateUpper',1,scenario_opt.H); %same upper bound for entire horizon
scenario_opt.constraints.nlcon = @(u,w,scen_struct) vehicle_nlconst(u, w, scen_struct);

% Input constraints
inputLower = [-5,-deg2rad(25)];
inputUpper = [5,deg2rad(25)];
scenario_opt.constraints.ub = inputUpper;
scenario_opt.constraints.lb = inputLower;

%% Estimation Setup
% Setup parameters and initial conditions for estimation. Initial state,
% particle filter parameters

t = 0:dt:30;
T = length(t);
%x0 = [0; 0; 0.1; 45deg];
u0 = [0;0];

% INITIALIZATIONS
x = zeros(nx,T);                    %state
y = zeros(ny,T);                    %measurement
xh_p = zeros(nx, T);                %estimate
xsig2_p = zeros(nx, T);             %uncertainty
u = zeros(nu,T);                    %control
w = zeros(nx,scenario_opt.H,scenario_opt.N); %disturbance scenarios

% PARAMETERS 
obs = @(x) H*x;                                 %observation equation
p_obs_noise = @(v) mvnpdf(v, zeros(ny,1), R);   %likelihood equations
p_yk_given_xk = @(yk, xk) p_obs_noise(yk - obs(xk));
gen_x0 = @(x,P) mvnrnd(x,P);                    %initial particle dist.

% BOOTSTRAP PARTICLE FILTER
pf.k                = 1;                    % initial iteration number
pf.Ns               = 6000;                 % number of particles
pf.wk               = zeros(pf.Ns, 1);      % weights array
pf.wkm1             = zeros(pf.Ns, 1);      % weights at last time
pf.particles        = zeros(nx, pf.Ns);     % particles array
pf.xkm1             = zeros(nx, pf.Ns);     % particles at last time
pf.xh0              = x0;                   % initial state prediction
pf.P0 = blkdiag(0.05*eye(nx-nw),0.05*eye(nw)); % initial state uncertainty
pf.gen_x0           = gen_x0;               % function for sampling from initial pdf p_x0
pf.p_yk_given_xk    = p_yk_given_xk;        % function of the observation likelihood PDF p(y[k] | x[k])
pf.resample         = 'systematic_resampling'; % resampling strategy 

%% Control
scenario_opt.xk = x0;
ukm1 = zeros(scenario_opt.H,nu);

x(:,1) = x0;        %initial state
y(:,1) = outputFcn(x(:,1)); %initial measurement
xh_p(:,1) = pf.xh0; %initial estimate
u(:,1) = u0;        %initial control
plan = zeros(scenario_opt.H,nu,T);
cost = [];
active_constraints = zeros(2,scenario_opt.H,T);
track_width = 1; %half width

%plotting 
% car1 = car(x0(1),x0(2),0.5,0.2,x0(3),x0(4));
% car1.color = 'blue';
% save = 0;
% if save
%     set(gca,'nextplot','replacechildren');
%     video = VideoWriter('live_plot1.avi');
%     open(video);
% end


tic
for i = 2:T 
    % Generate Scenarios
    w(:,:,:) = randn(nx,scenario_opt.H,scenario_opt.N);
    
    % Compute reference
    ref = [62.83, 0, 4, 0, 0, 0, 0];
    
    % Compute constraints
    const_center(1,:) = scenario_opt.xk(1):stateUpper(3)*dt:scenario_opt.xk(1)+stateUpper(3)*dt*scenario_opt.H; %xk : v*dt : xk+v*dt*H
    const_center(2,:) = 3*sin(0.2*const_center(1,:));
    scenario_opt.constraints.stateL(2,:) = const_center(2,2:end) - track_width; 
    scenario_opt.constraints.stateU(2,:) = const_center(2,2:end) + track_width;
    active_constraints(:,:,i) = [const_center(2,2:end) - track_width; const_center(2,2:end) + track_width];
    
    % Compute optimal control
    ctrlstart = tic;
    [uopt,fcost] = Scenario_nonlinear(ukm1,w,ref,scenario_opt); 
    controlTime = toc(ctrlstart);
    ukm1 = uopt;
    u(:,i) = uopt(1,:)'; %applied control
    plan(:,:,i) = uopt;
    cost = [cost, fcost];
    %plot_control(scenario_opt.xk,plan(:,:,i),[const_center(1,2:end);active_constraints(:,:,i)],ref,scenario_opt.sys,car1);
    
    % Apply control (simulate dynamics)
    x(:,i) = scenario_opt.sys(x(:,i-1),u(:,i),randn(nx,1));
    
    % Get measurement
    y(:,i) = outputFcn(x(:,i));
    
    % Perform estimation 
    pf.k = i; 
    eststart = tic;
    [xh_p(:,i), xsig2_p(:,i), pf] = pf_control_scenario(@vehicle_model, y(:,i), pf, u(:,i), Aw, qw);
    estTime = toc(eststart);
    scenario_opt.xk = xh_p(:,i);
    
    % Goal Condition
    if norm([xh_p(1,i)-ref(1),xh_p(2,i)-ref(2)]) <= 0.5
        T = i;
        break
    end

end
tot_time = toc;

%% Review
fprintf('Total Control Time: %1.f min %2.f sec\n', floor(tot_time/60), mod(tot_time,60));
fprintf('Average Control Time: %1.f min %2.1f sec\n',floor(tot_time/T/60), mod(tot_time/T,60));
control_error = x(2,:)-3*sin(0.2*x(1,:));
control_error_sum = sum(abs(x(2,:)-3*sin(0.2*x(1,:))));
constraint_error_upper = 3*sin(0.2*x(1,:))+track_width - x(2,:);
constraint_error_lower = x(2,:) - (3*sin(0.2*x(1,:))-track_width);
constraint_idx = (constraint_error_lower < 0) | (constraint_error_upper < 0);
num_violations = sum(constraint_idx);
violation_prob = num_violations/(2*T);
% if save
%     close(video);
% end

%% Movie
figure(11)
xplan = zeros(nx,scenario_opt.H+1);
car1 = car(x0(1),x0(2),0.5,0.2,x0(3),x0(4));
car1.color = 'blue';
save = 0;
if save
    set(gca,'nextplot','replacechildren');
    video = VideoWriter('constraint_test.avi');
    open(video);
end

for i = 2:T
    %plot goal point
    plot(ref(1),ref(2),'k*','MarkerSize',7,'LineWidth',0.9) 
    hold on
    
    %compute constraints
    xplot_const1 = linspace(stateLower(1),stateUpper(1),141);
    plot(xplot_const1,3*sin(0.2*xplot_const1)-track_width,'k')
    plot(xplot_const1,3*sin(0.2*xplot_const1)+track_width,'k')
    %include active_constraints
    xplot_const2 = xh_p(1,i-1):stateUpper(3)*dt:xh_p(1,i-1)+stateUpper(3)*dt*(scenario_opt.H);
    plot(xplot_const2(2:end),active_constraints(1,:,i),'k','LineWidth',1.5);
    plot(xplot_const2(2:end),active_constraints(2,:,i),'k','LineWidth',1.5);
    

    %compute planned trajectory
    xplan(:,1) = xh_p(:,i-1);
    for j = 2:scenario_opt.H+1
        xplan(:,j) = scenario_opt.sys(xplan(:,j-1),plan(j-1,:,i)',zeros(nx,1));
    end
    plot(xplan(1,:),xplan(2,:),'r-')

    %plot estimated current state
    car1 = car1.update(xh_p(1,i-1), xh_p(2,i-1), xh_p(3,i-1), xh_p(4,i-1));
    plot(car1);
    
    hold off
    axis([-1 71 -5 5])
    if save
        frame = getframe(gcf);
        writeVideo(video,frame);
    end
    pause(0.1)
end
if save
    close(video);
end

%% Plotting
t = 0:dt:(T-1)*dt;
close all
% System States
figure(2)
subplot(4,1,1:2) %Position
plot(ref(1),ref(2),'k*','MarkerSize',7,'LineWidth',0.9), hold on %reference
plot(x(1,1:T),x(2,1:T),'g') %real
plot(xh_p(1,1:T),xh_p(2,1:T),'r','LineWidth',1.5) %estimate mean
% plot(xh_p(1,:)+2*sqrt(xsig2_p(1,:)),xh_p(2,:)+2*sqrt(xsig2_p(2,:)),'r--') %estimate uncertainty
% plot(xh_p(1,:)-2*sqrt(xsig2_p(1,:)),xh_p(2,:)-2*sqrt(xsig2_p(2,:)),'r--')
fill([xh_p(1,1:T)+2*sqrt(xsig2_p(1,1:T)) fliplr(xh_p(1,1:T)-2*sqrt(xsig2_p(1,1:T)))],...
     [xh_p(2,1:T)+2*sqrt(xsig2_p(2,1:T)) fliplr(xh_p(2,1:T)-2*sqrt(xsig2_p(2,1:T)))],[1 0.8 0.8],'edgecolor','r')
plot(y(1,1:T),y(2,1:T),'k.','MarkerSize',3) %measurements
plot(linspace(stateLower(1),stateUpper(1)),...
    3*sin(0.2*linspace(stateLower(1),stateUpper(1)))-track_width,'k') %constraints
plot(linspace(stateLower(1),stateUpper(1)),...
    3*sin(0.2*linspace(stateLower(1),stateUpper(1)))+track_width,'k')
plot([stateLower(1), stateLower(1)],[3*sin(0.2*stateLower(1))-track_width,...
                                        3*sin(0.2*stateLower(1))+track_width],'k')
plot([stateUpper(1), stateUpper(1)],[3*sin(0.2*stateUpper(1))-track_width,...
                                        3*sin(0.2*stateUpper(1))+track_width],'k')
ylabel('y(t)  [m]')
xlabel('x(t)  [m]')
ylim([-5 5])
xlim([-1 71])
legend('Reference','Real','Estimated','Uncertainty','Measurements','Location','best')

subplot(4,1,3) %Velocity
plot(t,stateLower(3)*ones(size(t)),'k'), hold on %constraints
plot(t,stateUpper(3)*ones(size(t)),'k')
plot(t,x(3,1:T),'g') %real
plot(t,xh_p(3,1:T),'r') %estimate
plot(t,xh_p(3,1:T)+2*sqrt(xsig2_p(3,1:T)),'r--') %uncertainty
plot(t,xh_p(3,1:T)-2*sqrt(xsig2_p(3,1:T)),'r--')
plot(t,y(3,1:T),'k.','MarkerSize',3) %measurements
xlabel('t  [s]')
ylabel('v(t)  [m/s]')
ylim([stateLower(3)-1 stateUpper(3)+1])

subplot(4,1,4) %Heading
plot(t,rad2deg(stateLower(4))*ones(size(t)),'k'), hold on %constraints
plot(t,rad2deg(stateUpper(4))*ones(size(t)),'k')
plot(t,rad2deg(x(4,1:T)),'g') %real
plot(t,rad2deg(xh_p(4,1:T)),'r') %estimate
plot(t,rad2deg(xh_p(4,1:T)+2*sqrt(xsig2_p(4,1:T))),'r--') %uncertainty
plot(t,rad2deg(xh_p(4,1:T)-2*sqrt(xsig2_p(4,1:T))),'r--')
plot(t,rad2deg(y(4,1:T)),'k.','MarkerSize',3) %measurements
xlabel('t  [s]')
ylabel('\theta(t)  [deg]')
ylim([rad2deg(stateLower(4))-10 rad2deg(stateUpper(4))+10])

sgtitle('System States')

% Latent States
figure(3)
subplot(3,1,1)
plot(t,x(5,1:T),'g'), hold on
plot(t,xh_p(5,1:T),'r')
% stairs(t,u,'k')
ylabel('z_1= w(t)  [m/s^2]')
legend('Real', 'Estimated','Location','best')

subplot(3,1,2)
plot(t,x(6,1:T),'g'), hold on
plot(t,xh_p(6,1:T),'r')
ylabel('z_2')

subplot(3,1,3)
plot(t,x(7,1:T),'g'), hold on
plot(t,xh_p(7,1:T),'r')
xlabel('t  [s]')
ylabel('z_3')

sgtitle('Latent Input States')

% Control Inputs
figure(4)
subplot(2,1,1)
plot(t,inputLower(1)*ones(size(t)),'k'), hold on
plot(t,inputUpper(1)*ones(size(t)),'k')
stairs(t,u(1,1:T),'b')
xlabel('t  [s]')
ylabel('a(t)  [m/s^2]')
ylim([inputLower(1)-1 inputUpper(1)+1])

subplot(2,1,2)
plot(t,rad2deg(inputLower(2))*ones(size(t)),'k'), hold on
plot(t,rad2deg(inputUpper(2))*ones(size(t)),'k')
stairs(t,rad2deg(u(2,1:T)),'b')
xlabel('t  [s]')
ylabel('\psi(t)  [rad]')
ylim([rad2deg(inputLower(2))-1 rad2deg(inputUpper(2))+1])

sgtitle('Control Input')

% Cost
figure(5)
plot(t(2:end),cost)
sgtitle('Cost')
ylabel('J(x,u)')
xlabel('t  [s]')
% axis([0 23 0 1e6])

% Estimation Residual
figure(6)
for plt=1:nx
    subplot(nx,1,plt)
    stem(t,xh_p(plt,1:T)-x(plt,1:T),'.')
%     ylabel('x_%d',plt)
end
xlabel('t  [s]')
sgtitle('Estimation Residuals')

% Presentation figure
figure(7)
subplot(4,2,1:2) %position
plot(ref(1),ref(2),'k*','MarkerSize',7,'LineWidth',0.9), hold on
plot(x(1,1:T),x(2,1:T),'g')
plot(xh_p(1,1:T),xh_p(2,1:T),'r')
plot(linspace(stateLower(1),stateUpper(1)),...
    3*sin(0.2*linspace(stateLower(1),stateUpper(1)))-track_width,'k') %constraints
plot(linspace(stateLower(1),stateUpper(1)),...
    3*sin(0.2*linspace(stateLower(1),stateUpper(1)))+track_width,'k')
plot([stateLower(1), stateLower(1)],[3*sin(0.2*stateLower(1))-track_width,...
                                        3*sin(0.2*stateLower(1))+track_width],'k')
plot([stateUpper(1), stateUpper(1)],[3*sin(0.2*stateUpper(1))-track_width,...
                                        3*sin(0.2*stateUpper(1))+track_width],'k')
ylabel('y(t)  [m]')
xlabel('x(t)  [m]')
ylim([-5 5])
xlim([-1 71])
legend('Goal','Actual State','Estimated State','Constraints','Location','best')

subplot(4,2,3:4) %velocity
plot(t,stateLower(3)*ones(size(t)),'k'), hold on
plot(t,stateUpper(3)*ones(size(t)),'k')
plot(t,x(3,1:T),'g')
plot(t,xh_p(3,1:T),'r')
ylabel('v(t)  [m/s]')
ylim([stateLower(3)-1 stateUpper(3)+1])

subplot(4,2,5:6) %latent
plot(t,x(5,1:T),'g'), hold on
plot(t,xh_p(5,1:T),'r')
ylabel('w(t)  [m/s^2]')

subplot(4,2,7) %acc
plot(t,inputLower(1)*ones(size(t)),'k'), hold on
plot(t,inputUpper(1)*ones(size(t)),'k')
stairs(t,u(1,1:T),'b')
xlabel('t  [s]')
ylabel('a(t)  [m/s^2]')
ylim([inputLower(1)-1 inputUpper(1)+1])

subplot(4,2,8) %steer
plot(t,rad2deg(inputLower(2))*ones(size(t)),'k'), hold on
plot(t,rad2deg(inputUpper(2))*ones(size(t)),'k')
stairs(t,rad2deg(u(2,1:T)),'b')
xlabel('t  [s]')
ylabel('\psi(t)  [deg]')
ylim([rad2deg(inputLower(2))-5 rad2deg(inputUpper(2))+5])

sgtitle('Estimation and Control of Nonlinear LFM')

%% Paper Plots
figure(8) %Position
plot(linspace(0,100),3*sin(0.2*linspace(0,100)),'k--'), hold on %reference
plot(x(1,1:T),x(2,1:T),'g') %real
plot(xh_p(1,1:T),xh_p(2,1:T),'r') %estimate
% plot(y(1,:),y(2,:),'k.','MarkerSize',3) %measurements
plot(linspace(stateLower(1),stateUpper(1)),...
    3*sin(0.2*linspace(stateLower(1),stateUpper(1)))-track_width,'k') %constraints
plot(linspace(stateLower(1),stateUpper(1)),...
    3*sin(0.2*linspace(stateLower(1),stateUpper(1)))+track_width,'k')
plot([stateLower(1), stateLower(1)],[3*sin(0.2*stateLower(1))-track_width,...
                                        3*sin(0.2*stateLower(1))+track_width],'k')
plot([stateUpper(1), stateUpper(1)],[3*sin(0.2*stateUpper(1))-track_width,...
                                        3*sin(0.2*stateUpper(1))+track_width],'k')
ylabel('y(t)  [m]')
xlabel('x(t)  [m]')
title('Position')
ylim([-5 5])
xlim([-1 71])
legend('Reference','Real','Estimated','Location','best')

figure(9) %Latent Effect
% subplot(2,1,1) %latent
colororder({'b','r'})
yyaxis left
plot(t,x(5,1:T),'b--'), hold on %latent input
plot(t,xh_p(5,1:T),'b-') %estimate
ylabel('w(t)  [m/s^2]')
ylim([inputLower(1)-1 inputUpper(1)+1])

yyaxis right %acceleration
stairs(t,u(1,1:T),'r-')
plot(t,inputLower(1)*ones(size(t)),'k-')
plot(t,inputUpper(1)*ones(size(t)),'k-')
xlabel('t  [s]')
ylabel('a(t)  [m/s^2]')
title('a) Latent Input and Acceleration Input')
ylim([inputLower(1)-1 inputUpper(1)+1])
legend('Latent Disturbance', 'Latent Disturbance Estimate', 'Acceleration Control')

% subplot(2,1,2) %velocity
% plot(t,vref*ones(size(t)),'k--'), hold on
% plot(t,x(3,:),'g')
% plot(t,xh_p(3,:),'r')
% plot(t,stateLower(3)*ones(size(t)),'k')
% plot(t,stateUpper(3)*ones(size(t)),'k')
% ylabel('v(t)  [m/s]')
% xlabel('t  [s]')
% title('b) Velocity')
% ylim([stateLower(3)-1 stateUpper(3)+1])
% legend('Reference', 'True State', 'Estimated State')
