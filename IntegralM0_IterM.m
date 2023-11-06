function [M0] = IntegralM0_IterM(forwarding,B1,ni)
% 版本时间：23.7.24
% 功能: RK45 to calculate emission M0
%% 0.initialization assignment
gamma = 0.267518*1e9;        % Magnetic spin ratio
N     = 6.692*1e28;          % [/m^3]
hq    = 1.054571628*1e-34;   % Planck's constant/2*pi [J.s]
K     = 1.3805*1e-23;        % Boltzmann's constant  [J/K]
T0     = 293;                 % absolute temperature  [K]
M0stre = N * gamma^2 * hq^2 / (4 * K * T0) .* forwarding.MB0.Hsum(ni);  % magnetization vector strength
Ln = length(ni);
rx = zeros(Ln,1);
ry = zeros(Ln,1);                                            
rz = M0stre.*ones(Ln,1);
T = [forwarding.model.WT1(ni),forwarding.model.WT2(ni)];
t1 = forwarding.Pulse.t_90;
%% 
T1 = T(:,1);
T2 = T(:,2);
fT    = forwarding.TLoop.fT;        % RF frequency
T_turns = forwarding.TLoop.turn;    % Number of turns of transmitting/receiving coil
I_trans = forwarding.TLoop.I_trans; % Transmitting/receiving coil energizing current
delt_B0 = forwarding.MB0.Hsum  - 2*pi*fT /gamma; % ΔB0N
BT_P = 0.5*I_trans* T_turns*B1;%B+_T 
BT_P(isnan(BT_P)) = 0;
alpha = atan2(delt_B0(ni),abs(BT_P(ni)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Beff = [abs(BT_P(ni)),zeros(Ln,1),forwarding.B0.Hsum(ni)];
param.M0ini = reshape([rx,ry,rz]',3*Ln,1);
param.T = [T1,T2];
param.Beff1 = Beff;
param.fT = forwarding.TLoop.fT;
tspan = [0 t1];
N = 1;
warning off
options = odeset('MaxStep', 1e-8,'RelTol',1e-14,'AbsTol',1e-14,'Stats','off');
dM0_t = ode45(@(t,M00)Solve_AllBloch_labNorelax(t,M00,param,N),tspan,param.M0ini,options);
M0n1 = dM0_t.y;

M0.rx = M0n1(1:3:end,:);
M0.ry = M0n1(2:3:end,:);
M0.rz = M0n1(3:3:end,:);
M0.rxy = sqrt(M0n1.^2 +  M0n1.^2);
M0.rsum = sqrt(M0.rx.^2 +  M0.ry.^2 + M0.rz.^2);
M0.t = dM0_t.x';



