function [M0] = IntegralM0_RM(forwarding,B1,ni,t)
% Data：24.2.23
% function: Rotation Matrix method to calculate emission M0
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
gamma = 0.267518e9;        % Magnetic spin ratio
N     = 6.692e28;          % [/m^3]
hq    = 1.054571628e-34;   % Planck's constant/2*pi [J.s]
K     = 1.3805e-23;        % Boltzmann's constant  [J/K]
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rotation Matrix Method
%%%%%%%%%%% arbitrary rotation axis 
n_x = abs(B1(ni)); 
n_y = zeros(Ln,1); 
n_z = delt_B0(ni); 
n_xyz = sqrt(n_x.^2 + n_y.^2 + n_z.^2);
nx = n_x./n_xyz;
ny = n_y./n_xyz;
nz = n_z./n_xyz;
tt = t';
Lt = length(tt);
%%%%%%%%%%% θ
theta = -gamma .* kron(n_xyz,tt)';
%%%%%%%%%%% rotation matrix 
% R_n = [nx.^2.*(1-cos(theta))+cos(theta)      , nx.*ny*(1-cos(theta))+nz*sin(theta)  , nx.*nz.*(1-cos(theta))-ny*sin(theta);
%        nx.*ny.*(1-cos(theta))-nz.*sin(theta), ny.^2.*(1-cos(theta))+cos(theta)     , ny.*nz.*(1-cos(theta))-nx*sin(theta);
%        nx.*nz.*(1-cos(theta))+ny.*sin(theta), ny.*nz.*(1-cos(theta))-nx*sin(theta) , nz.^2.*(1-cos(theta))+cos(theta)];
R_n = [nx.^2.*(1-cos(theta))+cos(theta)      , nx.*ny.*(1-cos(theta))-nz.*sin(theta)  , nx.*nz.*(1-cos(theta))+ny.*sin(theta);
       nx.*ny.*(1-cos(theta))+nz.*sin(theta) , ny.^2.*(1-cos(theta))+cos(theta)       , ny.*nz.*(1-cos(theta))-nx.*sin(theta);
       nx.*nz.*(1-cos(theta))-ny.*sin(theta) , ny.*nz.*(1-cos(theta))+nx.*sin(theta)  , nz.^2.*(1-cos(theta))+cos(theta)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
M_0 = [rx;ry;rz];
Rx = [R_n(1:Lt,:)*M_0]';
Ry = [R_n(Lt+1:2*Lt,:)*M_0]';
Rz = [R_n(2*Lt+1:3*Lt,:)*M_0]';
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
M0.t = tt';
end


