function [M0] = IntegralM0_RM_V1(forwarding,M0_ini,t0)
% Data：24.2.23
% function: Rotation Matrix method to calculate emission M0
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rotation Matrix Method
%%%%%%%%%%% arbitrary rotation axis 
n_x = abs(B1); 
n_y = zeros(Ln,1); 
n_z = delt_B0; 
n_xyz = sqrt(n_x.^2 + n_y.^2 + n_z.^2);
nx = n_x./n_xyz;
ny = n_y./n_xyz;
nz = n_z./n_xyz;
tt = t0';
Lt = length(tt);
%%%%%%%%%%% θ
theta = -gamma .* kron(n_xyz,tt);
Nx = repmat(nx,1,length(tt));
Ny = repmat(ny,1,length(tt));
Nz = repmat(nz,1,length(tt));
%%%%%%%%%%% rotation matrix 
R_nx = [reshape(Nx.^2.*(1-cos(theta))+cos(theta),Ln*Lt,1),...
        reshape(Nx.*Ny.*(1-cos(theta))-Nz.*sin(theta),Ln*Lt,1),... 
        reshape(Nx.*Nz.*(1-cos(theta))+Ny.*sin(theta),Ln*Lt,1)];
R_ny = [reshape(Nx.*Ny.*(1-cos(theta))+nz.*sin(theta),Ln*Lt,1),...
        reshape(Ny.^2.*(1-cos(theta))+cos(theta),Ln*Lt,1),...
        reshape(Ny.*Nz.*(1-cos(theta))-Nx.*sin(theta),Ln*Lt,1)];
R_nz = [reshape(Nx.*Nz.*(1-cos(theta))-ny.*sin(theta),Ln*Lt,1),...
        reshape(Ny.*Nz.*(1-cos(theta))+Nx.*sin(theta),Ln*Lt,1),...
        reshape(Nz.^2.*(1-cos(theta))+cos(theta),Ln*Lt,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
M_0 = repmat([rx,ry,rz],Lt,1);
Rx = reshape(sum(R_nx.*M_0,2),Ln,Lt);
Ry = reshape(sum(R_ny.*M_0,2),Ln,Lt);
Rz = reshape(sum(R_nz.*M_0,2),Ln,Lt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
M0.rxy = sqrt(Rx.^2 +  Ry.^2);
M0.rsum = sqrt(Rx.^2 +  Ry.^2 + Rz.^2);
phi1 = atan2(Ry,Rx) - kron(RW1,tt);
M0.rx = M0.rxy .* cos(phi1);
M0.ry = M0.rxy .* sin(phi1);
M0.rz = Rz;
M0.t = tt';
end


