function [M0] = IntegralM0_GM(forwarding,B1,ni,t)
% Data：23.7.24
% function: Geometric analysis method to calculate emission M0
% Input:
%   B1: RF field 
%   ni: Which voxels are involved in the calculation(usually all,You can also calculate certain voxels individually)
%   t ：sequentially
%   forwarding.Pulse.t_90 :Pulse excitation time
%   forwarding.TLoop.fT :Transmit pulse frequency
% Output:
%   M0.rx, M0.ry,M0.rz : Magnetization Vector Tripartite
%   M0.t :iterative time series
%% 0.initialization assignment
gamma = 0.267518*1e9;        % Magnetic spin ratio
N     = 6.692*1e28;          % [/m^3]
hq    = 1.054571628*1e-34;   % Planck's constant/2*pi [J.s]
K     = 1.3805*1e-23;        % Boltzmann's constant  [J/K]
T0     = 293;                 % absolute temperature  [K]
M0stre = N * gamma^2 * hq^2 / (4 * K * T0) .* forwarding.B0.Hsum(ni);  % magnetization vector strength
Ln = length(ni);
rx = zeros(Ln,1);
ry = zeros(Ln,1);                                            
rz = M0stre.*ones(Ln,1);
%% 
fT    = forwarding.TLoop.fT;        % RF frequency
RW1 = 2*pi*fT;
delt_B0 = forwarding.B0.Hsum  - 2*pi*fT /gamma; % ΔB0N
alpha = atan2(delt_B0(ni),abs(B1(ni)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
rsum(:,1) = sqrt(rx(:, 1).^2 +  ry(:, 1).^2 + rz(:, 1).^2);
rxy(:, 1) = sqrt(rx(:, 1).^2 +  ry(:, 1).^2);  
%%%%%%%%%%% Beff 
Beff_x = abs(B1(ni)); 
Beff_y = zeros(Ln,1); 
Beff_z = delt_B0(ni); 
Beff_xyz = sqrt(Beff_x.^2 + Beff_y.^2 + Beff_z.^2);
%%%%%%%%%%% cosξ
coskexi = (Beff_x.*rx(:, 1) + Beff_y.*ry(:, 1) + Beff_z.*rz(:, 1))...    
    ./ (Beff_xyz.*rsum(:, 1));
%%%%%%%%%%% r'
rT = rsum(:, 1).*sqrt(1-coskexi.^2);
%%%%%%%%%%% 
Beffa_x = rsum(:, 1) .* coskexi .* cos(alpha);
Beffa_z = rsum(:, 1) .* coskexi .* sin(alpha);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% X'OY
%%%%%%%%%%% β
beta = atan2(ry(:, 1).*abs(sin(alpha)),round(rx(:, 1)-Beffa_x,10)).*(alpha~=0)...
    + atan2(ry(:, 1),-rz(:, 1)).*(alpha==0);
tt = t';
Lt = length(tt);
%%%%%%%%%%% θ
Alpha = kron(alpha,tt);
theta = gamma * kron(Beff_xyz,tt).*(Alpha>=0) - gamma * kron(Beff_xyz,tt).*(Alpha<0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Rx = kron(Beffa_x,ones(1,Lt)) + rT .* cos(kron(beta,ones(1,Lt)) - theta) .* abs(sin(alpha));
Ry = rT .* sin(kron(beta,ones(1,Lt)) - theta);
Rz = kron(Beffa_z,ones(1,Lt)) - rT .* cos(kron(beta,ones(1,Lt)) - theta) .* cos(alpha).*(alpha>=0) ...
     + rT .* cos(kron(beta,ones(1,Lt)) - theta) .* cos(alpha).*(alpha<0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
M0.rxy = sqrt(Rx.^2 +  Ry.^2);
M0.rsum = sqrt(Rx.^2 +  Ry.^2 + Rz.^2);
phi1 = atan2(Ry,Rx) - kron(ones(Ln,1),RW1*tt);
M0.rx = M0.rxy .* cos(phi1);
M0.ry = M0.rxy .* sin(phi1);
M0.rz = Rz;
M0.rx(:,1) = rx;
M0.ry(:,1) = ry;
M0.rz(:,1) = rz;
M0.rxy(:,1) = rxy;
M0.sum(:,1) = rsum;
M0.t = tt';
end


