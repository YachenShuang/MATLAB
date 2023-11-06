function dM0_t = Solve_AllBloch_labNorelax(t,M0,param,N)
% 日期：23.7.29
% 功能：Solving the Bloch equation in the lab coordinate system, Beff = B0 + B1 when there is no relaxation excitation, and Beff = B0 when there is relaxation.
% param :  
% param.Beff： Beff 
% param.T   ： relaxation time
gamma = 267518000.000000;
if(N == 1)
    Bx = param.Beff1(:,1)*cos(-2*pi*param.fT*t);
    By = param.Beff1(:,1)*sin(-2*pi*param.fT*t);
    Bz = param.Beff1(:,3);
else
    Bx = param.Beff2(:,1);
    By = param.Beff2(:,2);
    Bz = param.Beff2(:,3);
end
Mx = M0(1:3:end);
My = M0(2:3:end);
Mz = M0(3:3:end);
dMx = gamma*(My.*Bz - Mz.*By);
dMy = gamma*(Mz.*Bx - Mx.*Bz);
dMz = gamma*(Mx.*By - My.*Bx);
dmt = reshape([dMx,dMy,dMz]',length(dMx)*3,1);
dM0_t = dmt;
end