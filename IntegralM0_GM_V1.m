function [M0] = IntegralM0_GM_V1(forwarding,M0_ini,t0)
% Data：24.2.23
% function: Geometric analysis method to calculate emission M0
% Input:
%   M0_ini: Initial magnetization vector 
%   t0: Pulse excitation time
%   forwarding.TLoop.fT : Pulse emission frequency
%   forwarding.B1Hsum : RF field
%   forwarding.B0.Hsum :Static magnetic field strength
% Output:
%   M0.rx, M0.ry,M0.rz : Magnetization Vector Tripartite
%   M0.t :iterative time series
%% 0.initialization assignment
gamma = 0.267518e9;        % Magnetic spin ratio
Ln = size(M0_ini,1);
rx = M0_ini(:,1);
ry = M0_ini(:,2);                                            
rz = M0_ini(:,3);
%% 
B1 = forwarding.B1Hsum;
B0 = forwarding.B0.Hsum;
fT    = forwarding.TLoop.fT;        % RF frequency
RW1 = 2*pi*fT;
delt_B0 = B0  - 2*pi*fT /gamma; % ΔB0N
alpha = atan2(delt_B0,abs(B1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
rsum(:,1) = sqrt(rx(:, 1).^2 +  ry(:, 1).^2 + rz(:, 1).^2);
rxy(:, 1) = sqrt(rx(:, 1).^2 +  ry(:, 1).^2);  
%%%%%%%%%%% Beff 
Beff_x = abs(B1); 
Beff_y = zeros(Ln,1); 
Beff_z = delt_B0; 
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
beta = atan2(ry(:, 1).*abs(sin(alpha)),rx(:, 1)-Beffa_x).*(alpha~=0)...
    + atan2(ry(:, 1),-rz(:, 1)).*(alpha==0);
tt = t0';
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
phi1 = atan2(Ry,Rx) - kron(RW1,tt);
M0.rx = M0.rxy .* cos(phi1);
M0.ry = M0.rxy .* sin(phi1);
M0.rz = Rz;
M0.t = tt';
end


